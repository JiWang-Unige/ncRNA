#!/usr/bin/env python3
"""S4 (PLAN_BATCH2.md C; card batch2-downloads 2026-11-04): do the zero-shot lnc calls of v1.9 on chicken 5 and
zebrafish 22 (P1.6 reads, final predictions conf >= tau*) land on human lncRNA genes after mapping to hg38 through
UCSC chains, more often than position-shuffled copies of themselves?  An RNA-seq-independent precision estimate.

Mapping (own chain reader, tested in tests/test_s4_synteny.py): every exon base is carried through the aligned
blocks of the chains whose target is the species chromosome; a prediction is lifted when >= 50% of its exon bases
land on one (hg38 chromosome, orientation) and the mapped span is <= 2 Mb.  The lifted strand flips when the chain's
query strand is '-'.  Hit = the lifted span overlaps a same-strand v50 lnc gene body (mRNA bodies reported too).
Control: each prediction shifted 100 times (seed 0) to a uniform random place on its chromosome, exon structure
kept, < 10% N in the shifted exons.  Strata: truth-hit (strict L1 hit to a species lnc gene) / novel (no same-strand
exon overlap with any annotated species gene) / other.

  s4_synteny.py --out outputs/PHASE1/S4/S4.json
"""
import argparse
import gzip
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1"))
from ncrna_annot.eval import ruler  # noqa: E402

SPECIES = {  # target: (chain, chain target name, codes, truth GTF, seqid)
    "chicken_5": ("data/batch2/chains/GCF_016699485.2ToHg38.over.chain.gz", "NC_052536.1", "data/p1/chicken_5_codes.npy",
                  "data/chicken/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.116.gtf.gz", "5"),
    "zebrafish_22": ("data/batch2/chains/danRer11ToHg38.over.chain.gz", "chr22", "data/p1/zebrafish_22_codes.npy",
                     "data/zebrafish/Danio_rerio.GRCz11.116.gtf.gz", "22"),
}
RUNS = ("g1b-42-ke", "g1b-123-ke", "g1b-456-ke")
N_SHUF, MIN_FRAC, MAX_SPAN = 100, 0.5, 2_000_000


def read_chains(fh, tname):
    """Aligned blocks of all chains with target `tname`: arrays sorted by target start (0-based, half open):
    t0, q0 (query + strand coordinate of the block's first base, see qmap), size, qname index, qstrand (+1 / -1),
    qsize, plus the list of query names.  Overlapping target blocks keep the higher-scoring chain."""
    blocks, qnames, cur = [], {}, None
    for line in fh:
        f = line.split()
        if not f:
            cur = None
            continue
        if f[0] == "chain":
            score, tn, tsz, tst, ts, te, qn, qsz, qst, qs, qe = int(f[1]), f[2], int(f[3]), f[4], int(f[5]), int(f[6]), f[7], int(f[8]), f[9], int(f[10]), int(f[11])
            assert tst == "+"
            cur = [score, ts, qs, qnames.setdefault(qn, len(qnames)), 1 if qst == "+" else -1, qsz] if tn == tname else None
            continue
        if cur is None:
            continue
        size = int(f[0])
        score, t, q, qi, qs_, qsz = cur
        blocks.append((t, q, size, qi, qs_, qsz, score))
        if len(f) == 3:
            cur[1], cur[2] = t + size + int(f[1]), q + size + int(f[2])
    blocks.sort(key=lambda b: b[0])
    kept, dropped = [], 0
    for b in blocks:   # target-overlapping blocks (rare in UCSC over.chain files): keep the higher-scoring chain
        if kept and b[0] < kept[-1][0] + kept[-1][2]:
            dropped += 1
            if b[6] > kept[-1][6]:
                kept[-1] = b
            continue
        kept.append(b)
    a = np.array([b[:6] for b in kept], np.int64).reshape(-1, 6)
    return {"t0": a[:, 0], "q0": a[:, 1], "size": a[:, 2], "qi": a[:, 3], "qstrand": a[:, 4], "qsize": a[:, 5],
            "qnames": [n for n, _ in sorted(qnames.items(), key=lambda kv: kv[1])], "dropped_overlaps": dropped}


def map_interval(ch, s, e):
    """0-based half-open target interval -> list of (qname index, qstrand, q_start, q_end) on the query + strand."""
    out = []
    i = max(np.searchsorted(ch["t0"], s, "right") - 1, 0)
    while i < len(ch["t0"]) and ch["t0"][i] < e:
        t0, sz = ch["t0"][i], ch["size"][i]
        a, b = max(s, t0), min(e, t0 + sz)
        if a < b:
            q0 = ch["q0"][i] + (a - t0)
            q1 = q0 + (b - a)
            if ch["qstrand"][i] < 0:   # query coordinates count on the reverse strand
                q0, q1 = ch["qsize"][i] - q1, ch["qsize"][i] - q0
            out.append((int(ch["qi"][i]), int(ch["qstrand"][i]), int(q0), int(q1)))
        i += 1
    return out


def lift(ch, exons, strand):
    """exons 1-based inclusive -> (qname index, hg38 strand, start, end) 0-based half open, or None."""
    pieces, total = [], 0
    for s, e in exons:
        total += e - s + 1
        pieces += map_interval(ch, s - 1, e)
    if not pieces:
        return None
    by = {}
    for qi, qs, a, b in pieces:
        by.setdefault((qi, qs), []).append((a, b))
    (qi, qs), iv = max(by.items(), key=lambda kv: sum(b - a for a, b in kv[1]))
    if sum(b - a for a, b in iv) < MIN_FRAC * total:
        return None
    lo, hi = min(a for a, _ in iv), max(b for _, b in iv)
    if hi - lo > MAX_SPAN:
        return None
    hs = strand if qs > 0 else ("-" if strand == "+" else "+")
    return qi, hs, lo, hi


class BodyIndex:
    """Same-strand gene-body overlap test on hg38 (bodies 1-based inclusive -> 0-based half open)."""
    def __init__(self, truth, cls):
        self.iv = {}
        for g in truth.values():
            if g["cls"] != cls:
                continue
            s = min(a for ex in g["tx"].values() for a, _ in ex) - 1
            e = max(b for ex in g["tx"].values() for _, b in ex)
            self.iv.setdefault((g["chrom"], g["strand"]), []).append((s, e))
        self.ix = {}
        for k, v in self.iv.items():
            v.sort()
            st = np.array([a for a, _ in v])
            en = np.maximum.accumulate(np.array([b for _, b in v]))
            self.ix[k] = (st, en)

    def hit(self, chrom, strand, s, e):
        k = (chrom, strand)
        if k not in self.ix:
            return False
        st, en = self.ix[k]
        i = np.searchsorted(st, e, "left") - 1     # last body starting before e
        return bool(i >= 0 and en[i] > s)


def read_hg38(path):
    """GENCODE v50 genome-wide: {gene_id: {chrom, strand, cls, tx}} (ruler.read_gtf per chromosome is too slow)."""
    genes = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] != "exon":
                continue
            a = ruler._attrs_gtf(f[8])
            g = genes.setdefault(a["gene_id"], {"chrom": f[0], "strand": f[6], "cls": ruler.gene_class(a.get("gene_type", "")), "tx": {}})
            g["tx"].setdefault(a["transcript_id"], []).append((int(f[3]), int(f[4])))
    return genes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/p16_v19.json")
    ap.add_argument("--hg38", default="data/human/gencode.v50.annotation.gtf.gz")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    import p15_score as P
    t0 = time.time()
    cfg = json.load(open(f"{ROOT}/{a.config}"))
    tau = cfg["decode"]["tau"]
    hg = read_hg38(f"{ROOT}/{a.hg38}")
    lnc_ix, mrna_ix = BodyIndex(hg, "lnc"), BodyIndex(hg, "mrna")
    rng = np.random.default_rng(0)
    res = {"tau": tau, "n_shuffles": N_SHUF, "species": {}}
    for sp, (chain, tname, codes_p, gtf, seqid) in SPECIES.items():
        with gzip.open(f"{ROOT}/{chain}", "rt") as fh:
            ch = read_chains(fh, tname)
        codes = np.load(f"{ROOT}/{codes_p}", mmap_mode="r")
        L = len(codes)
        isn = np.asarray(codes) == 4
        cn = np.concatenate([[0], np.cumsum(isn)])
        truth = ruler.read_gtf(f"{ROOT}/{gtf}", seqid)
        lnc_sp = {g: v for g, v in truth.items() if v["cls"] == "lnc"}
        idx_l, idx_all = ruler.ExonIndex(lnc_sp), ruler.ExonIndex(truth)
        out = {"chain_blocks": int(len(ch["t0"])), "aligned_bp": int(ch["size"].sum()), "dropped_overlaps": ch["dropped_overlaps"],
               "chrom_len": L, "runs": {}}

        def score(preds):
            n = lifted = hl = hm = 0
            for strand, ex in preds:
                n += 1
                r = lift(ch, ex, strand)
                if r is None:
                    continue
                lifted += 1
                qn = ch["qnames"][r[0]]
                hl += lnc_ix.hit(qn, r[1], r[2], r[3])
                hm += mrna_ix.hit(qn, r[1], r[2], r[3])
            return {"n": n, "lifted": lifted, "lnc_hit": hl, "mrna_hit": hm,
                    "lift_rate": lifted / n if n else None, "lnc_hit_rate": hl / lifted if lifted else None,
                    "mrna_hit_rate": hm / lifted if lifted else None}

        def shuffled(ex):
            s0, e0 = ex[0][0], ex[-1][1]
            span = e0 - s0 + 1
            for _ in range(50):
                st = int(rng.integers(1, L - span))
                d = st - s0
                sh = [(s + d, e + d) for s, e in ex]
                nb = sum(cn[e] - cn[s - 1] for s, e in sh)
                if nb < 0.1 * sum(e - s + 1 for s, e in sh):
                    return sh
            return None
        for run in RUNS:
            genes = P.load_genes(f"{ROOT}/outputs/PHASE1/P16/{sp}/{run}/genes.json")
            preds = [(p["strand"], sorted(p["exons"])) for p in P.final_preds(genes, tau) if p["cls"] == "lnc"]
            strata = {"truth_hit": [], "novel": [], "other": []}
            for strand, ex in preds:
                if idx_l.hits(strand, ex):
                    strata["truth_hit"].append((strand, ex))
                elif not idx_all.hits(strand, ex):
                    strata["novel"].append((strand, ex))
                else:
                    strata["other"].append((strand, ex))
            rr = {}
            for k, ps in list(strata.items()) + [("all", preds)]:
                obs = score(ps)
                ctrl = []
                for _ in range(N_SHUF):
                    sh = [(s, x) for s, x in ((s, shuffled(ex)) for s, ex in ps) if x is not None]
                    ctrl.append(score(sh)["lnc_hit_rate"] or 0.0)
                ctrl = np.array(ctrl)
                o = obs["lnc_hit_rate"]
                rr[k] = {**obs, "ctrl_lnc_hit_rate_mean": float(ctrl.mean()), "ctrl_p95": float(np.percentile(ctrl, 95)),
                         "enrichment": (o / ctrl.mean()) if (o is not None and ctrl.mean() > 0) else None,
                         "p_emp": float((np.sum(ctrl >= (o or 0)) + 1) / (N_SHUF + 1))}
            out["runs"][run] = rr
            print(json.dumps({sp: run, **{k: (v["n"], v["lifted"], round(v["lnc_hit_rate"] or 0, 3), round(v["ctrl_lnc_hit_rate_mean"], 3),
                                                round(v["enrichment"] or 0, 2), round(v["p_emp"], 3)) for k, v in rr.items()}}), flush=True)
        novel = [out["runs"][r]["novel"] for r in RUNS]
        n_ok = sum(1 for v in novel if v["p_emp"] < 0.01)
        out["novel_mean"] = {k: float(np.mean([v[k] for v in novel if v[k] is not None])) for k in ("lnc_hit_rate", "ctrl_lnc_hit_rate_mean", "lift_rate")}
        out["novel_mean"]["enrichment"] = out["novel_mean"]["lnc_hit_rate"] / out["novel_mean"]["ctrl_lnc_hit_rate_mean"] if out["novel_mean"]["ctrl_lnc_hit_rate_mean"] else None
        out["novel_runs_p_below_0.01"] = n_ok
        out["verdict_novel_synteny_support"] = bool(out["novel_mean"]["enrichment"] is not None and out["novel_mean"]["enrichment"] >= 2 and n_ok == 3)
        res["species"][sp] = out
    res["seconds"] = round(time.time() - t0)
    os.makedirs(os.path.dirname(f"{ROOT}/{a.out}"), exist_ok=True)
    json.dump(res, open(f"{ROOT}/{a.out}", "w"), indent=1)
    print(json.dumps({sp: {k: v[k] for k in ("novel_mean", "novel_runs_p_below_0.01", "verdict_novel_synteny_support")} for sp, v in res["species"].items()}), flush=True)


if __name__ == "__main__":
    main()
