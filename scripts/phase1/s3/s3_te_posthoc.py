#!/usr/bin/env python3
"""S3-TE post hoc (written after S3TE.json was read; descriptive, not part of PLAN_S3_TE.md): is the lower detection
of lnc genes whose canonical TSS lies in a TE explained by promoter CpG?  Per v50 lnc gene: CpG observed / expected in
canonical TSS +- 500 bp (s3_features.cpg_oe); CpG-high = o/e >= 0.6 (the usual CpG-island ratio).  Reports median
o/e by TSS-in-TE status and, within CpG-high and CpG-low genes, the Sn difference TSS-in-TE minus not (three-run mean
hits as in s3_te.py; 2,000 gene resamples, numpy default_rng(0)).

  s3_te_posthoc.py --out outputs/PHASE1/S3TE/S3TE_posthoc.json
"""
import argparse
import json
import os
import sys

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1/s3"))
import s3_te as T  # noqa: E402
import s3_features as F  # noqa: E402

ruler = T.ruler


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rmsk", default="data/te/rmsk.txt.gz")
    ap.add_argument("--out", default="outputs/PHASE1/S3TE/S3TE_posthoc.json")
    a = ap.parse_args()
    out = {"note": "post hoc, descriptive; CpG-high = canonical TSS +- 500 bp CpG o/e >= 0.6", "chroms": {}}
    for chrom, fn in (("chr20", T.hits_chr20), ("chr9", T.hits_chr9)):
        preds, _ = fn()
        truth = ruler.read_gtf(f"{ROOT}/{T.V50}", chrom)
        L = len(np.load(f"{ROOT}/data/p1/{chrom}_codes.npy", mmap_mode="r"))
        mask, _ = T.te_mask(f"{ROOT}/{a.rmsk}", chrom, L)
        canon = T.canonical(chrom)
        tab = T.gene_table(truth, mask, canon)
        seq = F.read_fasta(f"{ROOT}/data/human/{chrom}.fa")
        gids = sorted(tab)
        lidx = ruler.ExonIndex({k: truth[k] for k in gids})
        H = []
        for run, ps in sorted(preds.items()):
            got = set()
            for p in ps:
                got |= set(lidx.hits(p["strand"], p["exons"]))
            H.append(np.array([g in got for g in gids], float))
        H = np.mean(H, 0)
        oe = []
        for g in gids:
            gg = truth[g]
            tid = canon.get(g) if canon.get(g) in gg["tx"] else max(gg["tx"], key=lambda t: sum(b - a_ + 1 for a_, b in gg["tx"][t]))
            ex = sorted(gg["tx"][tid])
            tss = ex[0][0] - 1 if gg["strand"] == "+" else ex[-1][1] - 1
            oe.append(F.cpg_oe(seq[max(tss - 500, 0):tss + 501]))
        oe = np.array(oe)
        te = np.array([tab[g]["tss_te"] for g in gids])
        hi = oe >= 0.6
        rng = np.random.default_rng(0)
        cells = {"cpg_high": hi, "cpg_low": ~hi}
        boot = {k: [] for k in cells}
        for _ in range(2000):
            i = rng.integers(0, len(gids), len(gids))
            for k, m in cells.items():
                a1, a0 = m[i] & te[i], m[i] & ~te[i]
                boot[k].append(H[i][a1].mean() - H[i][a0].mean() if a1.any() and a0.any() else np.nan)
        res = {"median_cpg_oe": {"tss_te": float(np.median(oe[te])), "tss_not_te": float(np.median(oe[~te]))},
               "frac_cpg_high": {"tss_te": float(hi[te].mean()), "tss_not_te": float(hi[~te].mean())}, "cells": {}}
        for k, m in cells.items():
            res["cells"][k] = {"n_tss_te": int((m & te).sum()), "n_tss_not_te": int((m & ~te).sum()),
                               "sn_tss_te": float(H[m & te].mean()), "sn_tss_not_te": float(H[m & ~te].mean()),
                               "delta": float(H[m & te].mean() - H[m & ~te].mean()), "delta_ci95": np.nanpercentile(boot[k], [2.5, 97.5]).tolist()}
        out["chroms"][chrom] = res
        print(chrom, json.dumps(res), flush=True)
    json.dump(out, open(f"{ROOT}/{a.out}", "w"), indent=1)


if __name__ == "__main__":
    main()
