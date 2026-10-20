#!/usr/bin/env python3
"""0.5 scoring (PROTOCOL v2). Run on ont after the Modal runs; fetches outputs from the volume.

Primary: L1 Sn on the delta set (v50 lncRNA genes with no same-strand exon overlap with any v41 gene,
n=473 on chr20), seed-averaged per gene; U-ev minus F with a paired bootstrap over delta genes.
Guard: strict L1 Pr under v50 (segments overlapping any v50 lncRNA gene / all segments).
Predictions are unstranded (unstranded labels), so overlaps ignore strand for every arm.
"""
import json
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ruler  # noqa: E402

ROOT = os.path.expanduser("~/ncrna-annot")
OUT = f"{ROOT}/outputs/PHASE0/0.5_unknown"
ARMS, SEEDS = ["F", "Uev", "Upu"], [42, 123, 456]
if not os.path.isdir(OUT):
    subprocess.run([os.path.expanduser("~/envs/modal/bin/modal"), "volume", "get", "ncrna-annot",
                    "/outputs/PHASE0/0.5_unknown", f"{ROOT}/outputs/PHASE0/"],
                   env=dict(os.environ, MODAL_PROFILE="ji-wang-unige"), check=True)
new = ruler.read_gtf(f"{ROOT}/data/human/gencode.v50.annotation.gtf.gz", "chr20")
old = ruler.read_gtf(f"{ROOT}/data/human/gencode.v41.annotation.gtf.gz", "chr20")
delta = sorted(ruler.delta_genes(new, old))
unstr = lambda gs: {k: dict(v, strand=".") for k, v in gs.items()}
lnc_new = unstr({g: v for g, v in new.items() if v["cls"] == "lnc"})
lnc_old = unstr({g: v for g, v in old.items() if v["cls"] == "lnc"})
pc_new = unstr({g: v for g, v in new.items() if v["cls"] == "mrna"})
idx = {"delta": ruler.ExonIndex({g: lnc_new[g] for g in delta}), "v50": ruler.ExonIndex(lnc_new),
       "v41": ruler.ExonIndex(lnc_old), "pc": ruler.ExonIndex(pc_new)}

res, hitvec = {}, {}
for arm in ARMS:
    H = []
    for s in SEEDS:
        d = f"{OUT}/{arm}_s{s}"
        preds = ruler.read_pred_gff3(f"{d}/pred.gff3", "chr20")
        hits = {k: set() for k in idx}
        seg_v50 = seg_pc_only = 0
        for p in preds:
            h = {k: ix.hits(".", p["exons"]) for k, ix in idx.items()}
            for k in idx:
                hits[k] |= h[k]
            seg_v50 += bool(h["v50"])
            seg_pc_only += bool(h["pc"]) and not h["v50"]
        run = json.load(open(f"{d}/run.json"))
        n = max(len(preds), 1)
        res[f"{arm}_s{s}"] = {
            "n_segments": len(preds), "pred_bp": run["pred_bp"], "threshold": run["threshold"],
            "best_epoch": run["best_epoch"], "val_f1_v41": run["val_f1_v41"],
            "delta_L1_sn": len(hits["delta"]) / len(delta),
            "v50_lnc_L1_sn": len(hits["v50"]) / len(lnc_new), "v41_lnc_L1_sn": len(hits["v41"]) / len(lnc_old),
            "strict_pr_v50": seg_v50 / n, "frac_segments_in_mrna_only": seg_pc_only / n}
        H.append([g in hits["delta"] for g in delta])
    hitvec[arm] = np.mean(np.array(H, float), 0)  # seed-averaged per delta gene

summary = {"n_delta": len(delta), "runs": res, "arms": {}}
for arm in ARMS:
    rr = [res[f"{arm}_s{s}"] for s in SEEDS]
    summary["arms"][arm] = {k: [round(float(np.mean([r[k] for r in rr])), 4), round(float(np.std([r[k] for r in rr])), 4)]
                            for k in ("delta_L1_sn", "v50_lnc_L1_sn", "v41_lnc_L1_sn", "strict_pr_v50",
                                      "frac_segments_in_mrna_only", "n_segments", "val_f1_v41")}
rng = np.random.default_rng(0)
for arm in ("Uev", "Upu"):
    diff = hitvec[arm] - hitvec["F"]
    boots = [diff[rng.integers(0, len(diff), len(diff))].mean() for _ in range(1000)]
    pr_drop = summary["arms"]["F"]["strict_pr_v50"][0] - summary["arms"][arm]["strict_pr_v50"][0]
    summary[f"{arm}_minus_F"] = {"delta_L1_sn_diff": round(float(diff.mean()), 4),
                                 "ci95": [round(float(np.percentile(boots, 2.5)), 4), round(float(np.percentile(boots, 97.5)), 4)],
                                 "strict_pr_drop": round(float(pr_drop), 4)}
g = summary["Uev_minus_F"]
summary["GATE_0.5"] = {"criterion": "Uev-F >= 0.05 and CI95 low > 0 and strict Pr drop <= 0.05",
                       "pass": bool(g["delta_L1_sn_diff"] >= 0.05 and g["ci95"][0] > 0 and g["strict_pr_drop"] <= 0.05)}
json.dump(summary, open(f"{OUT}/SUMMARY_0.5.json", "w"), indent=1)
print(json.dumps({k: v for k, v in summary.items() if k != "runs"}, indent=1))
