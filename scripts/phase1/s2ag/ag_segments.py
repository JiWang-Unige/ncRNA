#!/usr/bin/env python3
"""S2-AG analysis B (PLAN_S2_AG.md section 3): lnc locus baseline from AlphaGenome's predicted RNA-seq.

Per strand, bases with signal >= t form segments; gaps <= 50 bp are merged, segments < 200 bp dropped; every segment is
one single-exon lnc prediction (AlphaGenome does not tell coding from non-coding, so segments on mRNAs count as wrong
under strict Pr).  t is chosen on chr18 only: candidates = quantiles 0.50 ... 0.995 of the chr18 signal over all bases
of both strands; the tau rule of protocol section 16 (largest lnc L1 Sn among t with strict Pr >= 0.43; if none, the
largest strict Pr, reported as such).  Then chr20 and chr9 are read once at that t.  Ruler: p15_score.score_preds (P.15),
bootstrap 2,000 for the reported reads.

  ag_segments.py --out outputs/PHASE1/S2AG/SEGMENTS.json
"""
import argparse
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1"))
QS = (0.50, 0.60, 0.70, 0.80, 0.85, 0.90, 0.95, 0.975, 0.99, 0.995)
GAP, MIN_LEN, PR_MIN, N_BOOT = 50, 200, 0.43, 2000


def segments(mask):
    """bool [L] -> list of [s, e) runs, gaps <= GAP merged, runs < MIN_LEN dropped."""
    d = np.diff(np.concatenate([[0], mask.view(np.int8), [0]]))
    s, e = np.flatnonzero(d == 1), np.flatnonzero(d == -1)
    out = []
    for a, b in zip(s, e):
        if out and a - out[-1][1] <= GAP:
            out[-1][1] = b
        else:
            out.append([a, b])
    return [(int(a), int(b)) for a, b in out if b - a >= MIN_LEN]


def preds_at(rna, t):
    out = []
    for col, strand in ((0, "+"), (1, "-")):
        for a, b in segments(rna[:, col] >= t):
            out.append({"id": f"{strand}{a}", "strand": strand, "cls": "lnc", "exons": [(a + 1, b)]})  # 1-based closed
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default="outputs/PHASE1/S2AG/auto/raw")
    ap.add_argument("--v50", default="data/human/gencode.v50.annotation.gtf.gz")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    import p15_score as P
    from ncrna_annot.eval import ruler
    t0 = time.time()
    load = lambda c: np.load(f"{ROOT}/{a.raw}/ag_{c}_rna.npy").astype(np.float32)
    rna18 = load("chr18")
    qv = np.quantile(rna18.ravel(), QS)
    truth18 = ruler.read_gtf(f"{ROOT}/{a.v50}", "chr18")
    grid = []
    for q, t in zip(QS, qv):
        f, _ = P.score_preds(truth18, preds_at(rna18, float(t)), boot=0)
        l1 = f["lnc"]["L1_locus"]
        grid.append({"q": q, "t": float(t), "n_pred": f["n_pred_lnc"], "sn": l1["sn"], "pr_strict": l1["pr_strict"]})
        print(json.dumps(grid[-1]), flush=True)
    ok = [g for g in grid if (g["pr_strict"] or 0) >= PR_MIN]
    pick = max(ok, key=lambda g: g["sn"]) if ok else max(grid, key=lambda g: g["pr_strict"] or 0)
    res = {"plan": "docs/phase1_results/PLAN_S2_AG.md section 3", "grid_chr18": grid, "pick": pick,
           "pick_rule": "max Sn s.t. strict Pr >= 0.43" if ok else "no t reaches strict Pr 0.43: max strict Pr", "reads": {}}
    del rna18
    for chrom in ("chr18", "chr20", "chr9"):
        truth = ruler.read_gtf(f"{ROOT}/{a.v50}", chrom)
        f, _ = P.score_preds(truth, preds_at(load(chrom), pick["t"]), boot=N_BOOT)
        res["reads"][chrom] = P.norm(f)
        l1 = f["lnc"]["L1_locus"]
        print(json.dumps({"chrom": chrom, "n_pred": f["n_pred_lnc"], "sn": l1["sn"], "pr_strict": l1["pr_strict"]}), flush=True)
    res["total_s"] = round(time.time() - t0)
    json.dump(res, open(f"{ROOT}/{a.out}", "w"), indent=1)


if __name__ == "__main__":
    main()
