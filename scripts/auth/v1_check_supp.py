#!/usr/bin/env python3
"""AUTH v1 step 2, supplement written after reading CHECK.json (post hoc, descriptive only): truth-end coverage by QuantSeq
sites at lower support thresholds, on the training chromosomes.  GM12878 has one biological QuantSeq replicate with two
technical replicates (frac 0.5 / 1.0), so the frozen filter frac >= 0.67 keeps only sites seen in both technical reps.
  v1_check_supp.py -> outputs/AUTH_V1/CHECK_SUPP.json
"""
import json
import os
import sys

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, f"{ROOT}/scripts/auth")
import auth_lib as A                                  # noqa: E402

L = A.L
DIRS = {"wtc11": f"{ROOT}/outputs/AUTH", "h1": f"{ROOT}/outputs/AUTH_V1/h1", "gm12878": f"{ROOT}/outputs/AUTH_V1/gm12878"}
out = {}
for cell, d in DIRS.items():
    A.OUT = d
    res = {}
    for thr in (0.0, 0.5, 0.67):
        n = hit = sites = 0
        for c in A.TRAIN:
            tr = A.load_truth(c)
            raw = A.evidence_raw(c)
            pas = {st: np.array(sorted(q["pos"] for q in raw["quant"] if q["strand"] == st and q["frac"] >= thr - 1e-9), np.int64)
                   for st in "+-"}
            sites += sum(len(v) for v in pas.values())
            cv = L.coverage(tr, {"tss": {st: np.zeros((0, 2), np.int64) for st in "+-"}, "pas": pas})
            n += cv["n"]
            hit += cv["pas"] * cv["n"]
        res[f"frac>={thr}"] = {"sites_train": sites, "truth_end_coverage": hit / n}
    out[cell] = res
A.OUT = DIRS["wtc11"]
json.dump(out, open(f"{ROOT}/outputs/AUTH_V1/CHECK_SUPP.json", "w"), indent=1)
print(json.dumps(out, indent=1))
