#!/usr/bin/env python3
"""PLAN_DEDEND sections 2-4: end localisation (Loc1k, +-20 kb) by Puffin (starts) and PolyaID (ends), CPU.

Truth, evidence and CpG as end_metrics.ends_table (GENCODE v50 Ensembl_canonical, FANTOM5 CAGE +-200 bp, PolyASite +-100 bp,
promoter CpG o/e TSS +-500 bp).  CpG quartiles are computed over the lncRNA starts of each chromosome group
(chr8+chr9 = Puffin held-out; chr18, chr20 = in Puffin training).

  dedend_metrics.py [--puffin outputs/PHASE1/DEDEND/puffin] [--polya outputs/PHASE1/DEDEND/polya] -> <dedend>/DEDEND.json
"""
import argparse
import json
import os
import sys

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
for d in ("src", "scripts/phase1", "scripts/phase1/s3", "scripts/phase1/ends"):
    sys.path.insert(0, os.path.join(ROOT, d))
import end_metrics as EM  # noqa: E402
import s3_features as S3  # noqa: E402

HALF, TOL = 20_000, 1000
GROUPS = {"heldout_chr8_chr9": ["chr8", "chr9"], "chr18": ["chr18"], "chr20": ["chr20"]}


def rows_of(c):
    EM.CHROM = c
    EM.TRUTH[c] = f"{ROOT}/data/human/gencode.v50.annotation.gtf.gz"
    seq = S3.read_fasta(f"{ROOT}/data/human/{c}.fa").upper()
    rows = EM.ends_table(seq)
    for r in rows:
        r["chrom"] = c
    return rows


def puffin_dist(track, r):
    x, L = r["TSS"], len(track)
    a, b = max(0, x - HALF), min(L, x + HALF + 1)
    w = np.asarray(track[a:b, 0 if r["s"] == "+" else 1], np.float32)
    return abs(a + int(np.argmax(w)) - x)


def polya_dist(P, r, primary=True):
    pos, pid, st = P[r["s"]]
    x = r["TES"]
    i, j = np.searchsorted(pos, x - HALF), np.searchsorted(pos, x + HALF, side="right")
    p, q, s = pos[i:j], pid[i:j].astype(np.float32), st[i:j]
    if not len(p):
        return None
    if primary:
        hi = q >= 0.75
        k = np.flatnonzero(hi)[np.argmax(s[hi])] if hi.any() else int(np.argmax(q))
    else:
        k = int(np.lexsort((-np.nan_to_num(s, nan=-1e9), -q))[0])
    return abs(int(p[k]) - x)


def summ(rows, key, sel):
    d = [r[key] for r in rows if sel(r) and r.get(key) is not None]
    return {"n": len(d), "loc1k": float(np.mean(np.array(d) <= TOL)) if d else None,
            "loc200": float(np.mean(np.array(d) <= 200)) if d else None,
            "median_dist": float(np.median(d)) if d else None}


def block(rows, key, end):
    out = {}
    for k, cls in (("l", "lnc"), ("m", "mrna")):
        rr = [r for r in rows if r["k"] == k]
        res = {"all": summ(rr, key, lambda r: True), "with_evidence": summ(rr, key, lambda r: r[f"{end}_ev"]),
               "without_evidence": summ(rr, key, lambda r: not r[f"{end}_ev"])}
        if k == "l" and end == "TSS":
            q = np.quantile([r["cpg"] for r in rr], [0.25, 0.5, 0.75])
            edges = [-np.inf, *q, np.inf]
            for qi in range(4):
                res[f"cpg_q{qi + 1}"] = summ(rr, key, lambda r, lo=edges[qi], hi=edges[qi + 1]: lo < r["cpg"] <= hi)
            res["cpg_q1q2"] = summ(rr, key, lambda r, hi=edges[2]: r["cpg"] <= hi)
        out[f"{cls}_{end}"] = res
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--puffin", default="outputs/PHASE1/DEDEND/puffin")
    ap.add_argument("--polya", default="outputs/PHASE1/DEDEND/polya")
    a = ap.parse_args()
    res = {"plan": "docs/phase1_results/PLAN_DEDEND.md", "puffin": {}, "polya": {}}
    rows = {}
    for c in ("chr8", "chr9", "chr18", "chr20"):
        rows[c] = rows_of(c)
        f = f"{ROOT}/{a.puffin}/{c}.npy"
        if os.path.exists(f + ".done"):
            tr = np.load(f, mmap_mode="r")
            for r in rows[c]:
                r["puffin"] = puffin_dist(tr, r)
        if c in ("chr18", "chr20"):
            fs = [f"{ROOT}/{a.polya}/{c}_{s}.npz" for s in "+-"]
            if all(map(os.path.exists, fs)):
                P = {}
                for s, fn in zip("+-", fs):
                    z = np.load(fn)
                    P[s] = (z["pos"], z["polyaid"], z["strength"])
                for r in rows[c]:
                    r["polya"] = polya_dist(P, r)
                    r["polya_secondary"] = polya_dist(P, r, primary=False)
        print(c, len(rows[c]), flush=True)
    for g, cs in GROUPS.items():
        rr = [r for c in cs for r in rows[c]]
        if any("puffin" in r for r in rr):
            res["puffin"][g] = block(rr, "puffin", "TSS")
        if any("polya" in r for r in rr):
            res["polya"][g] = {"primary": block(rr, "polya", "TES"), "secondary": block(rr, "polya_secondary", "TES")}
    v = {}
    H = res["puffin"].get("heldout_chr8_chr9")
    if H:
        t0 = H["mrna_TSS"]["with_evidence"]["loc1k"]
        t1 = H["lnc_TSS"]["cpg_q1q2"]["loc1k"]
        t2 = H["lnc_TSS"]["with_evidence"]["loc1k"]
        v["T0_valid"] = {"mrna_TSS_cage_loc1k": t0, "valid": bool(t0 >= 0.60)}
        v["T1"] = {"lnc_TSS_lowCpG_loc1k": t1, "outcome": "A" if t1 >= 0.25 else "B" if t1 <= 0.12 else "grey"}
        v["T2"] = {"lnc_TSS_cage_loc1k": t2, "outcome": "A" if t2 >= 0.50 else "B" if t2 <= 0.25 else "grey"}
    if all(c in res["polya"] for c in ("chr18", "chr20")):
        p = {c: res["polya"][c]["primary"]["lnc_TES"]["with_evidence"]["loc1k"] for c in ("chr18", "chr20")}
        v["P0_mrna_TES_polyA_loc1k"] = {c: res["polya"][c]["primary"]["mrna_TES"]["with_evidence"]["loc1k"] for c in ("chr18", "chr20")}
        v["P1"] = {"lnc_TES_polyA_loc1k": p,
                   "outcome": "A" if min(p.values()) >= 0.30 else "B" if max(p.values()) <= 0.15 else "grey"}
    res["verdict"] = v
    print(json.dumps(v, indent=1))
    os.makedirs(f"{ROOT}/outputs/PHASE1/DEDEND", exist_ok=True)
    json.dump(res, open(f"{ROOT}/outputs/PHASE1/DEDEND/DEDEND.json", "w"), indent=1)


if __name__ == "__main__":
    main()
