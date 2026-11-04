#!/usr/bin/env python3
"""v2.0 section 4.2 (protocol section 17 + amendment 1): verdicts of the three pre-registered replications on chr9.
  R-S1  (outputs/PHASE1/V20/rep/S1_chr9.json, judge = own S10k): |matched gap| <= 0.03; AUROC of both classes rises
        by >= 0.05 from the lowest usage layer (1-9 reads) to the highest (>= 1,000), layer AUROC = positive-weighted mean
        over the samples quintiles; R4 median gtexv2 reads of v50 lnc introns < 1/5 of mRNA introns.
  R-S1b (rep/S1B/<run>_A{1,2}/result.json): three-run mean lnc L2 Sn of A2 > A1 and A1 < 0.05.
  R-S3  (rep/S3_chr9.json, v2.0 predictions): promoter CpG and log gene span both "associated" (PLAN_S3 section 4).
-> rep/REPLICATE.json"""
import json
import os

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
REP = f"{ROOT}/outputs/PHASE1/V20/rep"
RUNS = ("v20-g1b-42-ke", "v20-g1b-123-ke", "v20-g1b-456-ke")


def layer_auroc(r1, cls, layer):
    st, n = r1[cls]["auroc_by_stratum"], r1[cls]["n_pos_by_stratum"]
    cells = [(st[layer * 5 + q], n[layer * 5 + q]) for q in range(5) if st[layer * 5 + q] is not None and n[layer * 5 + q] > 0]
    return sum(a * w for a, w in cells) / sum(w for _, w in cells)


def main():
    out = {}
    s1 = json.load(open(f"{REP}/S1_chr9.json"))
    r1 = s1["models"]["own"]["R1"]
    rise = {c: layer_auroc(r1, c, 3) - layer_auroc(r1, c, 0) for c in ("lnc", "mrna")}
    ratio = s1["R4"]["lnc_all"]["median_reads"] / s1["R4"]["mrna_all"]["median_reads"]
    gap = r1["gap"]["matched"]
    out["R-S1"] = {"matched_gap": gap, "matched_gap_ci95": r1["gap"]["matched_ci95"], "layer_rise": rise, "usage_ratio_lnc_over_mrna": ratio,
                   "replicated": bool(abs(gap) <= 0.03 and min(rise.values()) >= 0.05 and ratio < 0.2)}
    l2 = {arm: float(np.mean([json.load(open(f"{REP}/S1B/{r}_{arm}/result.json"))["final"]["lnc"]["L2_intron_chain"]["sn"] for r in RUNS]))
          for arm in ("A1", "A2", "A3")}
    out["R-S1b"] = {"lnc_L2_sn_mean": l2, "replicated": bool(l2["A2"] > l2["A1"] and l2["A1"] < 0.05)}
    s3 = json.load(open(f"{REP}/S3_chr9.json"))["R1"]
    out["R-S3"] = {"cpg_oe": s3["cpg_oe"], "log_span": s3["log_span"],
                   "replicated": bool(s3["cpg_oe"]["associated"] and s3["log_span"]["associated"])}
    json.dump(out, open(f"{REP}/REPLICATE.json", "w"), indent=1)
    print(json.dumps({k: v["replicated"] for k, v in out.items()}), flush=True)


if __name__ == "__main__":
    main()
