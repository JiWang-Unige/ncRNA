#!/usr/bin/env python3
"""Bio-validation of chr20 lnc predictions outside GENCODE v50 (docs/phase1_results/PREREG_BIOVAL.md, frozen 7c7bc88;
card bioval-data, user 2026-11-03).  Descriptive; nothing here feeds back into any model choice.

  bioval.py [--out outputs/PHASE1/BIOVAL]  -> <out>/RESULT.json

Implementation choices the pre-registration leaves open, fixed here before the analysis is run (committed with this
file):
  - predictions are read from the P1.5 batch-1 chr20 pred.gff3 of the 3 main runs (lnc conf >= tau*, >= 200 nt);
  - a gene's span is [first exon start, last exon end]; v50 gene spans use all transcripts of the gene;
  - "transcript not in v50" (S4): a multi-exon transcript whose intron chain equals that of a same-strand v50
    transcript, or a single-exon transcript overlapping a same-strand v50 single-exon transcript, counts as v50;
  - controls must, like the candidates, have no same-strand exon overlap with any v50 gene; S_inter / S_other
    controls satisfy the stratum's region rule; no N base inside the placed span; no overlap with the span of any
    prediction (lnc or mRNA) of any of the 3 runs; up to 20,000 draws per placement, failures are reported;
  - consensus loci: connected components of candidates linked by same-strand exon overlap that contain candidates
    of >= 2 runs; a locus (and its k-th control) is supported when any member (member's k-th control) is;
  - "known elsewhere": same-strand exon overlap (>= 1 nt) with a LncBook 2.1 transcript exon / an RNAcentral block;
  - Snaptron junctions with strand other than + / - are ignored.
"""
import argparse
import bisect
import glob
import gzip
import json
import os
import sys
import time
from collections import defaultdict

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
from ncrna_annot.eval import ruler  # noqa: E402

RUNS = ("g1b-42-ke", "g1b-123-ke", "g1b-456-ke")
PRED = "outputs/PHASE1/P15/chr20/{run}/pred.gff3"
V50 = "data/human/gencode.v50.annotation.gtf.gz"
RAW = "data/bioval/raw"
CHROM, FLANK, NCTRL, MAXTRY = "chr20", 5000, 100, 20000
NAMES = (CHROM, "20")  # UCSC / Ensembl names


def _open(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


# ---------------------------------------------------------------- evidence
class Points:
    """Same-strand point evidence (TSS / polyA sites), 1-based."""
    def __init__(self, pts):
        self.p = {s: np.array(sorted(v), np.int64) for s, v in pts.items()}

    def near(self, strand, lo, hi):
        a = self.p.get(strand, np.zeros(0, np.int64))
        return bisect.bisect_right(a, hi) - bisect.bisect_left(a, lo) > 0


def read_cage(path):
    pts = defaultdict(list)
    with _open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] in NAMES and f[5] in "+-":
                pts[f[5]].append(int(f[6]) + 1)
    return Points(pts)


def read_polyasite(path):
    pts = defaultdict(list)
    with _open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        ic, ist, ir = head.index("chrom"), head.index("strand"), head.index("rep")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[ic] in NAMES and f[ist] in "+-":
                pts[f[ist]].append(int(f[ir]) + 1)
    return Points(pts)


def read_snaptron(files):
    seen, rows = set(), []
    for fn in files:
        with open(fn) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[0].startswith("DataSource") or len(f) < 12 or f[1] in seen:
                    continue
                seen.add(f[1])
                if f[6] not in "+-" or f[2] != CHROM:
                    continue
                rows.append((f[6], int(f[3]), int(f[4]), int(f[10]), int(f[11])))  # strand, start, end, samples, cov
    return rows


class Junctions:
    def __init__(self, rows, offset, v50_introns):
        by = defaultdict(list)
        for st, s, e, n, cov in rows:
            s, e = s + offset, e + offset
            by[st].append((s, e, n, cov, (st, s, e) in v50_introns))
        self.j = {st: sorted(v) for st, v in by.items()}
        self.starts = {st: [x[0] for x in v] for st, v in self.j.items()}

    def inside(self, strand, a, b):
        js, ss = self.j.get(strand, []), self.starts.get(strand, [])
        return js[bisect.bisect_left(ss, a):bisect.bisect_right(ss, b)]


def gtf_transcripts(path):
    tx = {}
    with _open(path) as fh:
        for line in fh:
            if line.split("\t", 1)[0] not in NAMES:
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] != "exon" or f[6] not in "+-":
                continue
            tid = f[8].split('transcript_id "', 1)[1].split('"', 1)[0]
            tx.setdefault(tid, (f[6], []))[1].append((int(f[3]), int(f[4])))
    return {t: (st, ruler._merge(sorted(ex))) for t, (st, ex) in tx.items()}


def bed12_blocks(path):
    out = []
    with _open(path) as fh:
        for line in fh:
            if line.split("\t", 1)[0] not in NAMES:
                continue
            f = line.rstrip("\n").split("\t")
            s, st = int(f[1]), f[5]
            if st not in "+-":
                continue
            if len(f) >= 12 and f[9]:
                sizes = [int(x) for x in f[10].strip(",").split(",")]
                offs = [int(x) for x in f[11].strip(",").split(",")]
                out.append((st, [(s + o + 1, s + o + z) for o, z in zip(offs, sizes)]))
            else:
                out.append((st, [(s + 1, int(f[2]))]))
    return out


class ExonIx:
    """Same-strand exon lookup over (strand, exons) items: which items have an exon overlapping a query exon."""
    BIN = 100_000

    def __init__(self, items):
        self.items = items
        self.bins = defaultdict(list)
        for i, (st, ex) in enumerate(items):
            for s, e in ex:
                for b in range(s // self.BIN, e // self.BIN + 1):
                    self.bins[(st, b)].append((s, e, i))

    def overlaps(self, strand, exons, frac=None):
        """frac None: any overlap; else overlap >= frac * min(len query exon, len item exon)."""
        for s, e in exons:
            for b in range(s // self.BIN, e // self.BIN + 1):
                for gs, ge, _ in self.bins.get((strand, b), ()):
                    o = min(e, ge) - max(s, gs) + 1
                    if o > 0 and (frac is None or o >= frac * min(e - s + 1, ge - gs + 1)):
                        return True
        return False


# ---------------------------------------------------------------- support
def support(st, exons, ev, ref=False):
    """Dict of S1-S4 (+ sensitivity variants) for one structure.  ref=True drops the two "not in v50" exclusions."""
    a, b = exons[0][0], exons[-1][1]
    tss, tes = (a, b) if st == "+" else (b, a)
    d = 1 if st == "+" else -1
    r = {}
    for w in (100, 500, 1000):
        r[f"S1_{w}"] = ev["cage"].near(st, tss - w, tss + w) if ev["cage"] else False
    r["S3_500"] = ev["pas"].near(st, tes - 500, tes + 500) if ev["pas"] else False
    lo, hi = sorted((tes - 50 * d, tes + 200 * d))
    r["S3_win"] = ev["pas"].near(st, lo, hi) if ev["pas"] else False
    js = [j for j in ev["junc"].inside(st, a, b) if j[1] <= b and (ref or not j[4])] if ev["junc"] else []
    intr = set(ruler.introns(exons))
    r["S2_5"] = any(n >= 5 and cov >= 20 for _, _, n, cov, _ in js)
    r["S2_10"] = any(n >= 10 and cov >= 20 for _, _, n, cov, _ in js)
    r["S2_exact"] = any(n >= 5 and cov >= 20 and (s, e) in intr for s, e, n, cov, _ in js)
    s4 = {}
    for k in ("flair", "chess"):
        ix = ev[k + "_all" if ref else k]
        s4[k] = ix.overlaps(st, exons, 0.5) if ix else False
        chains = ev[k + "_chains_all" if ref else k + "_chains"]
        s4[k + "_chain"] = bool(len(exons) > 1 and chains is not None and (st, tuple(ruler.introns(exons))) in chains)
    r["S4"] = s4["flair"] or s4["chess"]
    r["S4_flair"] = s4["flair"]
    r["S4_chain"] = s4["flair_chain"] or s4["chess_chain"]
    r["primary"] = r["S2_5"] or r["S4"] or (r["S1_500"] and r["S3_500"])
    r["S1_only"] = r["S1_500"] and not r["primary"]
    r["known_lncbook"] = ev["lncbook"].overlaps(st, exons) if ev["lncbook"] else None
    r["known_rnacentral"] = ev["rnacentral"].overlaps(st, exons) if ev["rnacentral"] else None
    return r


KEYS = ("primary", "S1_100", "S1_500", "S1_1000", "S2_5", "S2_10", "S2_exact", "S3_500", "S3_win", "S4", "S4_flair",
        "S4_chain", "S1_only")


def rate(rows, k):
    v = [r[k] for r in rows if r.get(k) is not None]
    return float(np.mean(v)) if v else None


# ---------------------------------------------------------------- placement
def place_controls(cands, L, ok_start, rng):
    """For each candidate, NCTRL placements (st, exons) or None.  ok_start(a, b, strand, exons) checks a placement."""
    out = []
    for c in cands:
        a0, b0 = c["exons"][0][0], c["exons"][-1][1]
        span = b0 - a0
        rel = [(s - a0, e - a0) for s, e in c["exons"]]
        pl = []
        for _ in range(NCTRL):
            for _t in range(MAXTRY):
                a = int(rng.integers(1, L - span))
                st = "+" if rng.random() < 0.5 else "-"
                ex = [(a + s, a + e) for s, e in rel]
                if ok_start(a, a + span, st, ex, c["stratum"]):
                    pl.append((st, ex))
                    break
            else:
                pl.append(None)
        out.append(pl)
    return out


def summarize(cands, ctrl_support, flip_support):
    res = {"n": len(cands)}
    for k in KEYS:
        obs = rate([c["support"] for c in cands], k)
        reps = []
        for i in range(NCTRL):
            rows = [cs[i] for cs in ctrl_support if cs[i] is not None]
            reps.append(rate(rows, k) if rows else None)
        reps = [x for x in reps if x is not None]
        cm = float(np.mean(reps)) if reps else None
        res[k] = {"cand": obs, "ctrl_mean": cm, "ratio": (obs / cm) if obs is not None and cm else None,
                  "p_emp": float(np.mean([x >= obs for x in reps])) if reps and obs is not None else None,
                  "flip": rate(flip_support, k)}
    for k in ("known_lncbook", "known_rnacentral"):
        res[k] = rate([c["support"] for c in cands], k)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="outputs/PHASE1/BIOVAL")
    a = ap.parse_args()
    t0 = time.time()
    out = f"{ROOT}/{a.out}"
    os.makedirs(out, exist_ok=True)
    man = json.load(open(f"{ROOT}/data/bioval/MANIFEST.json"))
    have = lambda n: n in man["files"]  # noqa: E731
    meta = {"started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "manifest_fetched": man["fetched"]}

    # truth
    truth = ruler.read_gtf(f"{ROOT}/{V50}", CHROM)
    tindex = ruler.ExonIndex(truth)
    gspan = [(g["strand"], min(s for ex in g["tx"].values() for s, _ in ex), max(e for ex in g["tx"].values() for _, e in ex))
             for g in truth.values()]
    v50_introns = {(g["strand"],) + i for g in truth.values() for ex in g["tx"].values() for i in ruler.introns(ruler._merge(ex))}
    v50_chains = {(g["strand"], tuple(ruler.introns(ruler._merge(ex)))) for g in truth.values() for ex in g["tx"].values()
                  if len(ex) > 1}
    v50_single = ExonIx([(g["strand"], ruler._merge(ex)) for g in truth.values() for ex in g["tx"].values()
                         if len(ruler._merge(ex)) == 1])
    starts = sorted((s, e) for _, s, e in gspan)
    gs_starts = [s for s, _ in starts]
    max_len = max(e - s for s, e in starts)

    def near_gene(a, b):  # any v50 gene (either strand) within [a - FLANK, b + FLANK]
        lo, hi = a - FLANK, b + FLANK
        i = bisect.bisect_left(gs_starts, lo - max_len)
        for s, e in starts[i:bisect.bisect_right(gs_starts, hi)]:
            if s <= hi and e >= lo:
                return True
        return False

    # evidence
    ev = {}
    ev["cage"] = read_cage(f"{ROOT}/{RAW}/fantom5_cage_peaks.bed.gz") if have("fantom5_cage_peaks.bed.gz") else None
    ev["pas"] = read_polyasite(f"{ROOT}/{RAW}/polyasite2_GRCh38.96.tsv.gz") if have("polyasite2_GRCh38.96.tsv.gz") else None
    sn_files = sorted(glob.glob(f"{ROOT}/{RAW}/snaptron/chr20_*.tsv"))
    rows = read_snaptron(sn_files) if sn_files else []
    meta["snaptron"] = {"files": len(sn_files), "junctions": len(rows)}
    if rows:  # 4b: coordinate convention from v50 introns only
        frac = {}
        for off in (-1, 0, 1):
            js = {(st, s + off, e + off) for st, s, e, _, _ in rows}
            frac[off] = len(v50_introns & js) / len(v50_introns)
        best = max(frac, key=frac.get)
        meta["snaptron_offset_check"] = {"frac_v50_introns_matched": {str(k): v for k, v in frac.items()}, "offset": best}
        if any(frac[o] * 2 > frac[best] for o in frac if o != best):
            json.dump(meta, open(f"{out}/RESULT.json", "w"), indent=1)
            sys.exit(f"Snaptron offset check failed: {frac}")
        ev["junc"] = Junctions(rows, best, v50_introns)
    else:
        ev["junc"] = None
    for k, fn in (("flair", "gtex_v9_flair_filter_transcripts.gtf.gz"), ("chess", "chess3.1.3.GRCh38.gtf.gz")):
        if not have(fn):
            ev[k] = ev[k + "_all"] = ev[k + "_chains"] = ev[k + "_chains_all"] = None
            continue
        tx = list(gtf_transcripts(f"{ROOT}/{RAW}/{fn}").values())
        in50 = [(len(ex) > 1 and (st, tuple(ruler.introns(ex))) in v50_chains) or (len(ex) == 1 and v50_single.overlaps(st, ex))
                for st, ex in tx]
        nov = [t for t, v in zip(tx, in50) if not v]
        ev[k], ev[k + "_all"] = ExonIx(nov), ExonIx(tx)
        ev[k + "_chains"] = {(st, tuple(ruler.introns(ex))) for st, ex in nov if len(ex) > 1}
        ev[k + "_chains_all"] = {(st, tuple(ruler.introns(ex))) for st, ex in tx if len(ex) > 1}
        meta[k] = {"transcripts_chr20": len(tx), "not_v50": len(nov)}
    ev["lncbook"] = ExonIx([(st, ex) for st, ex in gtf_transcripts(f"{ROOT}/{RAW}/lncbook2.1_GRCh38.gtf.gz").values()]) \
        if have("lncbook2.1_GRCh38.gtf.gz") else None
    ev["rnacentral"] = ExonIx(bed12_blocks(f"{ROOT}/{RAW}/rnacentral_homo_sapiens.GRCh38.bed.gz")) \
        if have("rnacentral_homo_sapiens.GRCh38.bed.gz") else None
    meta["available"] = {k: ev[k] is not None for k in ("cage", "pas", "junc", "flair", "chess", "lncbook", "rnacentral")}

    # predictions and candidates
    codes = np.load(f"{ROOT}/data/p1/chr20_codes.npy", mmap_mode="r")
    L = len(codes)
    ncum = np.concatenate([[0], np.cumsum(np.asarray(codes) == 4)])
    preds, cands = {}, {}
    for r in RUNS:
        preds[r] = ruler.read_pred_gff3(f"{ROOT}/{PRED.format(run=r)}", CHROM)
        cands[r] = []
        for p in preds[r]:
            if p["cls"] != "lnc" or tindex.hits(p["strand"], p["exons"]):
                continue
            ex = [tuple(x) for x in p["exons"]]
            cands[r].append({"run": r, "strand": p["strand"], "exons": ex,
                             "stratum": "S_other" if near_gene(ex[0][0], ex[-1][1]) else "S_inter"})
        cands[r].sort(key=lambda c: (c["exons"][0][0], c["strand"]))
    pspans = sorted((p["exons"][0][0], p["exons"][-1][1]) for r in RUNS for p in preds[r])
    ps_starts = [s for s, _ in pspans]
    pmax = max(e - s for s, e in pspans)

    def hits_pred(a, b):
        i = bisect.bisect_left(ps_starts, a - pmax)
        return any(s <= b and e >= a for s, e in pspans[i:bisect.bisect_right(ps_starts, b)])

    def ok_start(a, b, st, ex, stratum):
        if b > L or ncum[b] - ncum[a - 1] > 0 or hits_pred(a, b) or tindex.hits(st, ex):
            return False
        return near_gene(a, b) == (stratum == "S_other")

    rng = np.random.default_rng(0)
    res = {"meta": meta, "runs": {}}
    all_ctrl = {}
    for r in RUNS:
        cs = cands[r]
        for c in cs:
            c["support"] = support(c["strand"], c["exons"], ev)
        pl = place_controls(cs, L, ok_start, rng)
        ctrl = [[support(x[0], x[1], ev) if x else None for x in p] for p in pl]
        flip = [support("-" if c["strand"] == "+" else "+", c["exons"], ev) for c in cs]
        all_ctrl[r] = ctrl
        rr = {"n_candidates": len(cs), "placement_failures": sum(x is None for p in pl for x in p), "strata": {}}
        for sname in ("all", "S_inter", "S_other"):
            idx = [i for i, c in enumerate(cs) if sname == "all" or c["stratum"] == sname]
            rr["strata"][sname] = summarize([cs[i] for i in idx], [ctrl[i] for i in idx], [flip[i] for i in idx])
        res["runs"][r] = rr
        print(r, json.dumps({k: rr["strata"]["all"][k]["cand"] for k in ("primary", "S2_5", "S4")}), flush=True)
    res["mean"] = {s: {k: (float(np.mean([res["runs"][r]["strata"][s][k]["cand"] for r in RUNS
                                          if res["runs"][r]["strata"][s][k]["cand"] is not None]))
                          if any(res["runs"][r]["strata"][s][k]["cand"] is not None for r in RUNS) else None)
                       for k in KEYS} for s in ("all", "S_inter", "S_other")}
    res["mean_ctrl"] = {s: {k: (float(np.mean([res["runs"][r]["strata"][s][k]["ctrl_mean"] for r in RUNS
                                               if res["runs"][r]["strata"][s][k]["ctrl_mean"] is not None]))
                               if any(res["runs"][r]["strata"][s][k]["ctrl_mean"] is not None for r in RUNS) else None)
                            for k in KEYS} for s in ("all", "S_inter", "S_other")}

    # consensus loci (secondary)
    flat = [(r, i, c) for r in RUNS for i, c in enumerate(cands[r])]
    ix = ExonIx([(c["strand"], c["exons"]) for _, _, c in flat])
    parent = list(range(len(flat)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for j, (_, _, c) in enumerate(flat):
        for s, e in c["exons"]:
            for b in range(s // ExonIx.BIN, e // ExonIx.BIN + 1):
                for gs, ge, k in ix.bins.get((c["strand"], b), ()):
                    if gs <= e and s <= ge:
                        parent[find(j)] = find(k)
    comp = defaultdict(list)
    for j in range(len(flat)):
        comp[find(j)].append(j)
    loci = [m for m in comp.values() if len({flat[j][0] for j in m}) >= 2]
    obs = float(np.mean([any(flat[j][2]["support"]["primary"] for j in m) for m in loci])) if loci else None
    reps = []
    for k in range(NCTRL):
        v = []
        for m in loci:
            got = [all_ctrl[flat[j][0]][flat[j][1]][k] for j in m]
            got = [g for g in got if g is not None]
            if got:
                v.append(any(g["primary"] for g in got))
        reps.append(float(np.mean(v)) if v else None)
    reps = [x for x in reps if x is not None]
    res["consensus"] = {"n_loci": len(loci), "primary": obs, "ctrl_mean": float(np.mean(reps)) if reps else None,
                        "p_emp": float(np.mean([x >= obs for x in reps])) if reps and obs is not None else None}

    # positive references (4c)
    lnc_genes = [g for g in truth.values() if g["cls"] == "lnc"]
    ref_g = [support(g["strand"], ruler._merge(sorted(x for ex in g["tx"].values() for x in ex)), ev, ref=True) for g in lnc_genes]
    lnc_ix = ruler.ExonIndex({k: g for k, g in truth.items() if g["cls"] == "lnc"})
    tp = [p for r in RUNS for p in preds[r] if p["cls"] == "lnc" and lnc_ix.hits(p["strand"], p["exons"])]
    ref_tp = [support(p["strand"], [tuple(x) for x in p["exons"]], ev, ref=True) for p in tp]
    res["reference"] = {"v50_lnc_genes": {"n": len(lnc_genes), **{k: rate(ref_g, k) for k in KEYS}},
                        "true_positive_preds_3runs": {"n": len(tp), **{k: rate(ref_tp, k) for k in KEYS}}}
    meta["total_s"] = round(time.time() - t0)
    json.dump(res, open(f"{out}/RESULT.json", "w"), indent=1)
    json.dump([{k: c[k] for k in ("run", "strand", "exons", "stratum", "support")} for r in RUNS for c in cands[r]],
              open(f"{out}/candidates.json", "w"))
    print(json.dumps({"mean_all": res["mean"]["all"], "consensus": res["consensus"], "total_s": meta["total_s"]}), flush=True)


if __name__ == "__main__":
    main()
