#!/usr/bin/env python3
"""0.5 heads on ont (fallback while the Modal workspace is disabled). Same logic: train_05_core.run.

  train_05_local.py [--arms F,Uev,Upu] [--seeds 42,123,456] [--smoke]
Writes outputs/PHASE0/0.5_unknown/<arm>_s<seed>/ and one JSONL line per finished run to logs/p05_local.jsonl.
"""
import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import train_05_core  # noqa: E402

ROOT = os.path.expanduser("~/ncrna-annot")
ap = argparse.ArgumentParser()
ap.add_argument("--arms", default="F,Uev,Upu")
ap.add_argument("--seeds", default="42,123,456")
ap.add_argument("--smoke", action="store_true")
a = ap.parse_args()
feat = f"{ROOT}/data/p05_feat"
out = f"{ROOT}/outputs/PHASE0/0.5_unknown" + ("_smoke_local" if a.smoke else "")
for arm in a.arms.split(","):
    for s in map(int, a.seeds.split(",")):
        if os.path.exists(f"{out}/{arm}_s{s}/run.json"):
            continue
        t = time.time()
        r = train_05_core.run(arm, s, f"{ROOT}/data/p05", feat, out, a.smoke)
        with open(f"{ROOT}/logs/p05_local.jsonl", "a") as fh:
            fh.write(json.dumps({"arm": arm, "seed": s, "best": r[2], "s": round(time.time() - t)}) + "\n")
