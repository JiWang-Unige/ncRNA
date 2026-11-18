#!/usr/bin/env python3
"""Summarise chr18 E2 development outputs: mean over the three v2.0 runs, per (condition, frac, beta, tol, tau).
  e2_summary.py [pattern] -> outputs/M1/W24/dev/E2_chr18.json and a table on stdout"""
import glob, json, os, sys
from collections import defaultdict
import numpy as np
W = os.path.expanduser("~/ncrna-annot/outputs/M1/W24")
pat = sys.argv[1] if len(sys.argv) > 1 else "*"
acc = defaultdict(list)
for f in glob.glob(f"{W}/e2/chr18/*/{pat}.json"):
    if f.endswith(".lnc.json"):
        continue
    r = json.load(open(f))
    key0 = (r["cond"], r["frac"], r["seed"], r["beta"], r["tol"])
    for tau, s in r["by_tau"].items():
        acc[key0 + (float(tau),)].append((r["run"], s))
out = []
for k, v in sorted(acc.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][3], kv[0][4], -kv[0][5], kv[0][2])):
    m = {x: float(np.mean([s[x] for _, s in v if s[x] is not None])) for x in ("S50_sn", "S50_pr", "S50_pr_known", "L1_sn", "n_pred")}
    out.append({"cond": k[0], "frac": k[1], "seed": k[2], "beta": k[3], "tol": k[4], "tau": k[5], "n_runs": len(v), **m,
                "per_run_S50_sn": {r: s["S50_sn"] for r, s in v}})
json.dump(out, open(f"{W}/dev/E2_chr18.json", "w"), indent=1)
print(f"{'cond':6s} {'frac':>4s} {'sd':>2s} {'beta':>4s} {'tol':>3s} {'tau':>5s} {'n':>1s} {'S50':>6s} {'prK':>6s} {'pr':>6s} {'L1':>6s} {'npred':>6s}")
for o in out:
    print(f"{o['cond']:6s} {o['frac']:4.1f} {o['seed']:2d} {o['beta']:4.0f} {o['tol']:3d} {o['tau']:5.2f} {o['n_runs']:1d} {o['S50_sn']:6.3f} {o['S50_pr_known']:6.3f} {o['S50_pr']:6.3f} {o['L1_sn']:6.3f} {o['n_pred']:6.0f}")
