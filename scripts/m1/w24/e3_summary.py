#!/usr/bin/env python3
"""chr18 development: E3 vs evidence-only comparators across the evidence gradient (mean over seeds; model rows also over
the three runs).  -> outputs/M1/W24/dev/E3_SUMMARY_chr18.json + table"""
import glob, json, os
from collections import defaultdict
import numpy as np
D = os.path.expanduser("~/ncrna-annot/outputs/M1/W24/dev")
acc = defaultdict(lambda: defaultdict(list))
for f in sorted(glob.glob(f"{D}/E3_chr18_f*_s*.json")):
    r = json.load(open(f))
    for k, v in r["mean"].items():
        acc[k][r["frac"]].append(v)
keys = [k for k in acc if "_d" not in k]
out = {}
for k, byf in acc.items():
    out[k] = {str(fr): {m: float(np.mean([x[m] for x in v])) for m in ("S50_sn", "S50_pr_known", "n_pred")} | {"n_seeds": len(v)}
              for fr, v in sorted(byf.items())}
json.dump(out, open(f"{D}/E3_SUMMARY_chr18.json", "w"), indent=1)
show = [k for k in out if k.startswith("d50000") or k.startswith("d20000_K1_ev") or k.startswith("d100000_K1_ev")]
show = [k for k in show if ("ev-" in k and "_K3_" in k) or ("hy-" in k) or any(t in k for t in ("K3_th-1.0_model", "K3_th0.0_model+tib", "K3_th0.1_model+tib"))]
print(f"{'setting':34s} | " + " | ".join(f"f{fr}: S50  prK  n" for fr in ("0.1", "0.3", "1.0")))
for k in sorted(show):
    cells = []
    for fr in ("0.1", "0.3", "1.0"):
        x = out[k].get(fr)
        cells.append(f"{x['S50_sn']:.3f} {x['S50_pr_known']:.3f} {x['n_pred']:4.0f}" if x else "—")
    print(f"{k:34s} | " + " | ".join(cells))
