#!/usr/bin/env python3
"""chr18 development: keep only E2 lnc calls whose 5' end lies within w of an ssCAGE TSS ('tss'), whose 3' end lies within
w of a QuantSeq site ('pas'), or either / both; then score (precision with 'unknown' left out).
  e2_filter.py <run> <tag> [frac seed] -> stdout"""
import json, sys
import numpy as np
import w24_lib as L
run, tag = sys.argv[1], sys.argv[2]
frac, seed = (float(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 else (1.0, 0)
C = "chr18"
truth = L.load_truth(C)
ev = L.evidence(C, frac, seed)
lnc = json.load(open(f"{L.W}/e2/{C}/{run}/{tag}.lnc.json"))
def near(arr, x, w):
    return bool(len(arr)) and np.abs(arr - x).min() <= w
out = {}
for tau in (0.36, 0.2, 0.0):
    for w in (200, 1000):
        for rule in ("tss", "pas", "either", "both"):
            keep = []
            for p in lnc:
                if p["conf"] < tau:
                    continue
                st, s, e = p["strand"], p["exons"][0][0], p["exons"][-1][1]
                five, three = (s, e) if st == "+" else (e, s)
                T = ev["tss"][st][:, 0] if len(ev["tss"][st]) else np.zeros(0, np.int64)
                a, b = near(T, five, w), near(ev["pas"][st], three, w)
                if {"tss": a, "pas": b, "either": a or b, "both": a and b}[rule]:
                    keep.append(p)
            r = L.score_unknown(truth, keep, C)
            out[f"tau{tau}_w{w}_{rule}"] = {k: round(v, 3) if isinstance(v, float) else v for k, v in r.items() if k != "S50_ids"}
for k, v in sorted(out.items(), key=lambda kv: -(kv[1]["S50_sn"] or 0))[:14]:
    print(k, v)
