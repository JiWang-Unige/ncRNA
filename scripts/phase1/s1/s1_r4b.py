#!/usr/bin/env python3
"""S1 R4 second half (PLAN_S1.md section 8 deferral, PLAN_BATCH2.md A): on chr18, which multi-exon v50 lnc genes did
the S1b A2 decodes (truth splice sites, 3 main runs) recover, and how used are their introns in GTEx?
Gene usage = over the gene's multi-exon transcripts, the largest weakest-intron gtexv2 coverage_sum (exact junction
match, 0 if absent; Snaptron samples_count >= 5 as downloaded).  Recovered = a final lnc prediction (conf >= tau*)
on the gene's strand has exactly the intron chain of one of the gene's transcripts.  Descriptive only.

  s1_r4b.py --out outputs/PHASE1/S1/R4B_chr18.json
"""
import argparse
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s1_analyze as A  # noqa: E402
from s1_analyze import ROOT, ruler  # noqa: E402
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1"))
import p15_score as P  # noqa: E402

LAYERS = (0, 1, 10, 100, 1000)
RUNS = ("g1b-42-ke", "g1b-123-ke", "g1b-456-ke")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", default="outputs/PHASE1/DIAG_L2/chr18_v50.gtf")
    ap.add_argument("--snaptron", default="data/batch2/snaptron/chr18_*.tsv")
    ap.add_argument("--s1b", default="outputs/PHASE1/S1B")
    ap.add_argument("--tau", type=float, default=0.6)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    truth = ruler.read_gtf(f"{ROOT}/{a.truth}", "chr18")
    jn = A.read_snaptron(sorted(glob.glob(f"{ROOT}/{a.snaptron}")), "chr18")
    jcov = {(int(s), int(e), st): int(c) for s, e, st, c in zip(jn["s"], jn["e"], jn["strand"], jn["cov"])}
    genes = {}
    for gid, g in truth.items():
        if g["cls"] != "lnc":
            continue
        chains = [ruler.introns(ruler._merge(sorted(ex))) for ex in g["tx"].values()]
        chains = [c for c in chains if c]
        if chains:
            u = max(min(jcov.get((s, e, g["strand"]), 0) for s, e in c) for c in chains)
            genes[gid] = {"strand": g["strand"], "chains": {tuple(c) for c in chains}, "usage": u}
    res = {"n_multi_exon_lnc_genes": len(genes), "layers": ["0", "1-9", "10-99", "100-999", ">=1000"], "runs": {}}
    rec_any = {gid: 0 for gid in genes}
    for run in RUNS:
        pred = P.final_preds(P.load_genes(f"{ROOT}/{a.s1b}/{run}_A2/genes.json"), a.tau)
        have = {(p["strand"], ruler.introns(ruler._merge(sorted(p["exons"])))) for p in pred if p["cls"] == "lnc"}
        rec = {gid for gid, g in genes.items() if any((g["strand"], c) in have for c in g["chains"])}
        for gid in rec:
            rec_any[gid] += 1
        u_rec = np.array([genes[g]["usage"] for g in genes if g in rec])
        u_not = np.array([genes[g]["usage"] for g in genes if g not in rec])
        lay = lambda u: (np.bincount(np.searchsorted(LAYERS, u, "right") - 1, minlength=5) / max(len(u), 1)).tolist()
        res["runs"][run] = {"n_recovered": int(len(u_rec)), "frac_by_layer_recovered": lay(u_rec), "frac_by_layer_not": lay(u_not),
                            "median_usage_recovered": float(np.median(u_rec)) if len(u_rec) else None,
                            "median_usage_not": float(np.median(u_not)) if len(u_not) else None}
    u = np.array([genes[g]["usage"] for g in genes])
    lay = np.searchsorted(LAYERS, u, "right") - 1
    r = np.array([rec_any[g] / len(RUNS) for g in genes])
    res["recovery_rate_by_layer_mean_over_runs"] = [float(r[lay == k].mean()) if (lay == k).any() else None for k in range(5)]
    res["n_genes_by_layer"] = np.bincount(lay, minlength=5).tolist()
    json.dump(res, open(f"{ROOT}/{a.out}", "w"), indent=1)
    print(json.dumps(res), flush=True)


if __name__ == "__main__":
    main()
