"""M1 weeks 2-4 (docs/m1/PLAN_M1_W24.md): shared helpers -- WTC11 truth, exclusion zone, evidence (ssCAGE / QuantSeq) with
the evidence-amount gradient, S50 / L1 scoring.  Development uses chr18 only; test chromosomes are read once after the freeze.
"""
import gzip
import json
import os
import sys
from collections import defaultdict

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
for d in ("src", "scripts/phase1", "scripts/phase1/b1", "scripts/m1"):
    sys.path.insert(0, os.path.join(ROOT, d))
W = f"{ROOT}/outputs/M1/W24"
GEO = f"{ROOT}/data/wtc11/geo"
ENC = f"{ROOT}/data/wtc11/encode"
CAGE_SAMPLES = ("GSM5625727_WTC11-1_rep1", "GSM5625728_WTC11-1_rep2", "GSM5625729_WTC11-1_rep3",
                "GSM5625730_WTC11-3_rep1", "GSM5625731_WTC11-3_rep2", "GSM5625732_WTC11-3_rep3")
QUANT = ("ENCFF603KSV", "ENCFF937BYX", "ENCFF004EPK")
CAGE_MIN = 2          # pooled (possibly thinned) tags for a cluster to count as evidence
QUANT_MIN_FRAC = 0.67  # QuantSeq site supported by >= 2 of 3 replicates


# ------------------------------------------------------------------ truth
def load_truth(chrom):
    d = json.load(open(f"{W}/truth/{chrom}.json"))
    lnc = [x for x in d["loci"] if x["k"] == "l" and not x["excluded"]]
    ts = sorted([(x["id"], x["strand"], x["start"], x["end"]) for x in lnc], key=lambda t: t[2])
    return {"ts": ts, "lnc": {x["id"]: x for x in lnc}, "all": d["loci"], "zone": [tuple(z) for z in d["exclusion_zone"]]}


def in_zone(zone, s, e):
    """1-based closed [s, e] vs 0-based half-open zone intervals."""
    return any(a < e and s - 1 < b for a, b in zone)


# ------------------------------------------------------------------ evidence
def _cage_clusters(chrom):
    cl = {}
    for k, smp in enumerate(CAGE_SAMPLES):
        with gzip.open(f"{GEO}/{smp}.count.bed.gz", "rt") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[0] != chrom:
                    continue
                key = (int(f[1]), int(f[2]), f[5])
                cl.setdefault(key, np.zeros(len(CAGE_SAMPLES), np.int64))[k] = int(f[6])
    ctss = defaultdict(int)
    for smp in CAGE_SAMPLES:
        with gzip.open(f"{GEO}/{smp}.ctss.bed.gz", "rt") as fh:
            for line in fh:
                if not line.startswith(chrom + "\t"):
                    continue
                f = line.rstrip("\n").split("\t")
                if int(f[4]):
                    ctss[(int(f[2]), f[5])] += int(f[4])   # 1-based position
    out = []
    for (s, e, st), cnt in sorted(cl.items()):
        best, mode = -1, (s + 1 + e) // 2
        for x in range(s + 1, e + 1):
            c = ctss.get((x, st), 0)
            if c > best:
                best, mode = c, x
        out.append({"start": s + 1, "end": e, "strand": st, "tss": mode, "counts": cnt.tolist()})
    return out


def _quant_sites(chrom):
    sites = {}
    for acc in QUANT:
        with gzip.open(f"{ENC}/{acc}.bed.gz", "rt") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[0] != chrom:
                    continue
                pos, st = int(f[3].split(":")[1]), f[5]
                sites[(pos, st)] = max(sites.get((pos, st), 0.0), float(f[6]))
    return [{"pos": p, "strand": st, "frac": fr} for (p, st), fr in sorted(sites.items())]


def evidence_raw(chrom):
    p = f"{W}/evidence/{chrom}.json"
    if not os.path.exists(p):
        os.makedirs(f"{W}/evidence", exist_ok=True)
        json.dump({"cage": _cage_clusters(chrom), "quant": _quant_sites(chrom)}, open(p, "w"))
    return json.load(open(p))


def evidence(chrom, frac=1.0, seed=0):
    """Evidence at one gradient level: CAGE tags binomially thinned per sample (frac), clusters with pooled >= CAGE_MIN;
    QuantSeq sites with >= 2/3 replicate support, each kept with probability frac.  Seed 0..2; frac 1.0 is deterministic."""
    raw = evidence_raw(chrom)
    rng = np.random.default_rng([seed, int(round(frac * 1000))])
    tss = {"+": [], "-": []}
    for c in raw["cage"]:
        n = np.asarray(c["counts"])
        k = int((rng.binomial(n, frac) if frac < 1 else n).sum())
        if k >= CAGE_MIN:
            tss[c["strand"]].append((c["tss"], k))
    pas = {"+": [], "-": []}
    for q in raw["quant"]:
        if q["frac"] >= QUANT_MIN_FRAC and (frac >= 1 or rng.random() < frac):
            pas[q["strand"]].append(q["pos"])
    return {"tss": {s: np.array(sorted(v), np.int64).reshape(-1, 2) for s, v in tss.items()},
            "pas": {s: np.array(sorted(v), np.int64) for s, v in pas.items()}}


# ------------------------------------------------------------------ scoring
def score(truth, preds):
    """preds: dicts with strand, exons [(s, e)] 1-based closed.  Predictions overlapping the exclusion zone are dropped.
    S50: one-to-one reciprocal span overlap >= 0.5 (experiment 2's matcher); L1: span overlaps a truth lnc locus span."""
    import simple_det as D
    ts, zone = truth["ts"], truth["zone"]
    pp = [p for p in preds if not in_zone(zone, p["exons"][0][0], p["exons"][-1][1])]
    m50 = D.s50(ts, pp) if ts and pp else set()
    by = defaultdict(list)
    for gid, st, s, e in ts:
        by[st].append((s, e, gid))
    hit, nhit = set(), 0
    for p in pp:
        s, e = p["exons"][0][0], p["exons"][-1][1]
        h = {g for a, b, g in by[p["strand"]] if a <= e and s <= b}
        hit |= h
        nhit += bool(h)
    nt, n = len(ts), len(pp)
    return {"n_truth": nt, "n_pred": n, "S50_sn": len(m50) / nt if nt else None, "S50_pr": len(m50) / n if n else None,
            "L1_sn": len(hit) / nt if nt else None, "L1_pr": nhit / n if n else None, "S50_ids": sorted(m50)}


def coverage(truth, ev, tol_tss=100, tol_pas=100):
    """Fraction of evaluable truth lnc loci with a CAGE TSS within tol of any truth transcript start, a QuantSeq site within
    tol of any transcript end, and both."""
    n, a, b, ab = 0, 0, 0, 0
    for x in truth["lnc"].values():
        st = x["strand"]
        starts = [t["exons"][0][0] if st == "+" else t["exons"][-1][1] for t in x["tx"]]
        ends = [t["exons"][-1][1] if st == "+" else t["exons"][0][0] for t in x["tx"]]
        T = ev["tss"][st][:, 0] if len(ev["tss"][st]) else np.zeros(0, np.int64)
        P = ev["pas"][st]
        ha = any(len(T) and np.abs(T - s).min() <= tol_tss for s in starts)
        hb = any(len(P) and np.abs(P - e).min() <= tol_pas for e in ends)
        n += 1
        a += ha
        b += hb
        ab += ha and hb
    return {"n": n, "tss": a / n, "pas": b / n, "both": ab / n}


# ------------------------------------------------------------------ "unknown" loci for precision (development option)
def observed_any(chrom):
    """Spans (strand, s, e) of every transcript with >= 1 read in PacBio rep2 or rep3 (the truth replicates; rep1 is the
    long-read baselines' input and is left out so that their own finds are not excused), all novelty classes except
    'Genomic', cached.  Used to drop predictions that sit on WTC11-transcribed loci below the truth bar from the
    precision denominator (they are neither confirmed nor refuted)."""
    p = f"{W}/evidence/{chrom}.observed_rep23.json"
    if os.path.exists(p):
        return [tuple(x) for x in json.load(open(p))]
    import re
    spans = []
    for q, g in (("ENCFF527KWD.tsv", "ENCFF798DWV.gtf.gz"), ("ENCFF486EKE.tsv", "ENCFF422CGK.gtf.gz")):
        want = set()
        with open(f"{ENC}/{q}") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                f = dict(zip(hdr, line.rstrip("\n").split("\t")))
                if f["transcript_novelty"] != "Genomic" and int(f[hdr[-1]]) >= 1:
                    want.add(f["annot_transcript_id"])
        with gzip.open(f"{ENC}/{g}", "rt") as fh:
            for line in fh:
                if not line.startswith(chrom + "\t"):
                    continue
                f = line.split("\t", 8)
                if f[2] == "transcript" and re.search(r'transcript_id "([^"]+)"', f[8]).group(1) in want:
                    spans.append((f[6], int(f[3]), int(f[4])))
    spans = sorted(set(spans))
    json.dump(spans, open(p, "w"))
    return spans


def score_unknown(truth, preds, chrom):
    """score() plus a precision that leaves out predictions not S50-matched to the truth that overlap (same strand) a
    transcript seen in rep2 or rep3 -- 'unknown' -- instead of counting them as false."""
    r = score(truth, preds)
    obs = defaultdict(list)
    for st, s, e in observed_any(chrom):
        obs[st].append((s, e))
    obs = {k: np.array(sorted(v), np.int64).reshape(-1, 2) for k, v in obs.items()}
    pp = [p for p in preds if not in_zone(truth["zone"], p["exons"][0][0], p["exons"][-1][1])]
    n_tp = round(r["S50_pr"] * r["n_pred"]) if r["n_pred"] else 0
    unk = 0
    for p in pp:
        a = obs.get(p["strand"], np.zeros((0, 2), np.int64))
        s, e = p["exons"][0][0], p["exons"][-1][1]
        if len(a) and ((a[:, 0] <= e) & (s <= a[:, 1])).any():
            unk += 1
    # TP predictions also overlap observed transcripts; subtract them from the unknown count
    unk = max(unk - n_tp, 0)
    den = len(pp) - unk
    r["n_unknown"] = unk
    r["S50_pr_known"] = n_tp / den if den > 0 else None
    return r
