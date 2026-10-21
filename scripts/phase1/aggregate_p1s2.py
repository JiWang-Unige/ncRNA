#!/usr/bin/env python3
"""P1-S2 aggregation and pre-registered decision rules (protocol v1.4 section 11). Runs on ont.

Fetches castwj kernel outputs (ncrna-p1s2-<arm>-<seed>) into outputs/PHASE1/P1S2/<arm>-<seed>/,
tabulates S3 (and S1 from arm R) per seed, means, seed-paired differences and the rule outcomes.
"""
import json
import os
import subprocess
import sys

import numpy as np

ROOT = os.path.expanduser("~/ncrna-annot")
OUT = f"{ROOT}/outputs/PHASE1/P1S2"
ARMS, SEEDS = ("r", "s2", "s4"), (42, 123, 456)


def load(arm, seed, fetch=True):
    d = f"{OUT}/{arm}-{seed}"
    f = f"{d}/p1s2_{arm}_{seed}/result.json"
    if not os.path.exists(f) and fetch:
        os.makedirs(d, exist_ok=True)
        subprocess.run([os.path.expanduser("~/bin/kaggle-as"), "castwj", "kernels", "output",
                        f"castwj/ncrna-p1s2-{arm}-{seed}", "-p", d], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return json.load(open(f)) if os.path.exists(f) else None


def row(x):
    return {"auprc": x["auprc_lnc"], "f1": x["best_locus_f1"], "sn_L1": x["lnc_L1"]["sn"], "pr_L1": x["lnc_L1"]["pr"],
            "L2_sn": x["lnc_L2"]["sn"], "n_pred": x["n_pred_lnc"], "mrna_sn": x["mrna_L1"]["sn"], "mrna_pr": x["mrna_L1"]["pr"]}


def main():
    res = {}
    for arm in ARMS:
        for s in SEEDS:
            r = load(arm, s)
            if r is None:
                continue
            res[(arm.upper(), s)] = row(r["S3"])
            if arm == "r":
                res[("S1", s)] = row(r["S1"])
            res[(arm.upper(), s)]["best_epoch"] = r["best_epoch"]
    table = {f"{a}_{s}": v for (a, s), v in res.items()}
    arms = sorted({a for a, _ in res})
    means = {a: {k: float(np.mean([res[(a, s)][k] for s in SEEDS if (a, s) in res and res[(a, s)][k] is not None]))
                 for k in ("auprc", "f1", "sn_L1", "pr_L1", "L2_sn", "n_pred")} for a in arms}
    complete = all((a, s) in res for a in ("R", "S1", "S2", "S4") for s in SEEDS)

    def paired(a, b):  # b - a per seed
        d = [res[(b, s)]["auprc"] - res[(a, s)]["auprc"] for s in SEEDS if (a, s) in res and (b, s) in res]
        return {"per_seed": d, "mean": float(np.mean(d)) if d else None,
                "same_sign": bool(d) and (all(x > 0 for x in d) or all(x < 0 for x in d))}

    rules = {}
    if complete:
        d1, d2, d3 = paired("R", "S1"), paired("R", "S2"), paired("R", "S4")
        rules = {
            "rule1_R_vs_S1": {**d1, "decision": "re-examine (G2 design)" if d1["mean"] > 0.01 and d1["same_sign"] else "keep G1 + S3"},
            "rule2_R_vs_S2": {**d2, "decision": "lnc-only + coding offset" if d2["mean"] > 0.01 and d2["same_sign"] else "keep union grammar"},
            "rule3_R_vs_S4": {**d3, "decision": "switch main model to F" if d3["mean"] > 0.01 and d3["same_sign"] else "keep U-ev"},
        }
        for k, v in rules.items():
            if not v["same_sign"]:
                v["note"] = "seeds disagree in sign -> indistinguishable (rule 4)"
    out = {"complete": complete, "table": table, "means": means, "rules": rules}
    json.dump(out, open(f"{OUT}/SUMMARY_P1S2.json", "w"), indent=1)
    print("complete:", complete)
    for a in arms:
        m = means[a]
        print(f"{a:3s} AUPRC {m['auprc']:.4f}  bestF1 {m['f1']:.3f}  Sn {m['sn_L1']:.3f}  Pr {m['pr_L1']:.3f}  L2 {m['L2_sn']:.4f}  n_pred {m['n_pred']:.0f}",
              " seeds:", [round(res[(a, s)]["auprc"], 4) for s in SEEDS if (a, s) in res])
    for k, v in rules.items():
        print(k, {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in v.items() if kk != "per_seed"},
              [round(x, 4) for x in v["per_seed"]])


if __name__ == "__main__":
    main()
