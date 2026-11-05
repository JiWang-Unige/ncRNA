#!/usr/bin/env python3
"""v3.0 scoring (protocol section 18, execution detail 10; DRAFT_V30 section 5) of the reads of configs/v30.json.

  v30_score.py --config configs/v30.json   -> <out_dir>/SCORE_V30.json

Per test set and run: the run's arm reports lnc genes whose calibrated probability (the arm's calibration model on
the v20_calib.py features; raw confidence if the arm did not keep its calibration) is >= the arm's tau_v20; mRNA
genes are not filtered.  Ruler: p15_score.score_preds (P1.5).
Hierarchical bootstrap, 2,000 resamples, numpy default_rng(0) per test set: resample the truth lnc genes, then within
each arm draw 3 of its runs with replacement; statistic = mean over the drawn runs of lnc L1 Sn over the drawn genes.
Comparisons HM - H6 (main) and HM - H (secondary).  Strict Pr (reported for every set, a criterion only for the mouse
set): per resample, runs drawn with replacement and each run's reported lnc genes resampled.
Main verdict over the zero-shot sets (group B and the C sets with in_verdict): mean HM - H6 Sn difference >= +0.03,
CI lower bound > 0 on at least half of the sets, and the HM reported lnc genes summed over these sets (three-run means)
<= 1.3 x the H6 sum.  Human guard (chr18 development reads, each arm at its own tau_v20): HM Sn and strict Pr each
>= H6 - 0.02.  Nothing is tuned here.
"""
import argparse
import json
import math
import os
import sys
import time

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1/v20"))
import v20_calib as C  # noqa: E402
import p15_score as P  # noqa: E402
from ncrna_annot.eval import ruler  # noqa: E402

N_BOOT, ARMS = 2000, ("H", "H6", "HM")


def arm_scores(arm, genes, post, codes, tab):
    cal = arm["calibration"]
    rows = C.gene_rows(genes, post, codes, tab)
    X = np.array([[r[f] for f in cal["model"]["features"]] for r in rows], float).reshape(-1, len(cal["model"]["features"]))
    p = C.apply(cal["model"], X) if cal["keep"] else X[:, cal["model"]["features"].index("conf")]
    return {r["k"]: float(s) for r, s in zip(rows, p)}


def ci(x):
    return np.percentile(x, [2.5, 97.5]).tolist()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/v30.json")
    a = ap.parse_args()
    t0 = time.time()
    cfg_path = f"{ROOT}/{a.config}"
    cfg = json.load(open(cfg_path))
    out_dir = f"{ROOT}/{cfg['out_dir']}"
    tabs = {}
    for arm in cfg["arms"].values():
        h = arm["calibration"]["hexamer_table"]
        assert P.sha256(f"{ROOT}/{h}") == arm["calibration"]["hexamer_table_sha256"], "hexamer table changed"
        tabs[h] = json.load(open(f"{ROOT}/{h}"))
    res = {"config": a.config, "config_sha256": P.sha256(cfg_path), "script_sha256": P.sha256(os.path.abspath(__file__)),
           "arms": {k: {"tau_v20": v["tau_v20"], "calibration_kept": v["calibration"]["keep"], "runs": v["runs"]} for k, v in cfg["arms"].items()},
           "sets": {}}
    for t, tg in cfg["targets"].items():
        if t == "chr18":
            continue
        truth = ruler.read_gtf(f"{ROOT}/{tg['truth']}", tg["truth_chrom"])
        lnc = {k: v for k, v in truth.items() if v["cls"] == "lnc"}
        tids = sorted(lnc)
        lidx = ruler.ExonIndex(lnc) if lnc else None
        codes = np.load(f"{ROOT}/{tg['codes']}", mmap_mode="r")
        st = {"group": tg["group"], "in_verdict": bool(tg["in_verdict"]) and len(tids) > 0, "n_truth_lnc": len(tids),
              "n_truth_mrna": sum(v["cls"] == "mrna" for v in truth.values()), "arms": {}}
        hits, prs = {}, {}
        for arm_name in ARMS:
            arm = cfg["arms"][arm_name]
            tab = tabs[arm["calibration"]["hexamer_table"]]
            per, H, R = {}, [], []
            for run in arm["runs"]:
                rd = f"{out_dir}/{t}/{run}"
                meta = json.load(open(f"{rd}/read_meta.json"))
                assert meta["config_sha256"] == res["config_sha256"], f"{t} {run}: read with another config"
                genes = P.load_genes(f"{rd}/genes.json")
                post = np.load(f"{rd}/post.npy", mmap_mode="r")
                sc = arm_scores(arm, genes, post, codes, tab)
                gv = [dict(g, conf=sc.get(k, g["conf"])) for k, g in enumerate(genes)]
                preds = P.final_preds(gv, arm["tau_v20"])
                f, _ = P.score_preds(truth, preds, boot=0)
                rep = [g for g in gv if g["type"] == "l" and g["conf"] >= arm["tau_v20"]]
                rec = {"n_reported_lnc": len(rep), "n_reported_mrna": sum(p["cls"] == "mrna" for p in preds),
                       "mrna_sn": f["mrna"]["L1_locus"]["sn"], "mrna_pr_strict": f["mrna"]["L1_locus"]["pr_strict"]}
                if lidx is not None:
                    got, y_rep = set(), []
                    for g in rep:
                        h = lidx.hits(g["strand"], g["exons"])
                        got |= set(h)
                        y_rep.append(1.0 if h else 0.0)
                    H.append(np.array([x in got for x in tids], float))
                    R.append(np.array(y_rep))
                    lrows = [(k, g) for k, g in enumerate(genes) if g["type"] == "l"]
                    yl = np.array([1.0 if lidx.hits(g["strand"], g["exons"]) else 0.0 for _, g in lrows])
                    rec.update(sn=f["lnc"]["L1_locus"]["sn"], pr_strict=f["lnc"]["L1_locus"]["pr_strict"],
                               pr_lenient=f["lnc"]["L1_locus"]["pr_lenient"], L2_sn=f["lnc"]["L2_intron_chain"]["sn"],
                               ece=C.ece(np.array([sc.get(k, g["conf"]) for k, g in lrows]), yl) if len(lrows) else None)
                per[run] = rec
            mean = lambda key: float(np.mean([r[key] or 0.0 for r in per.values()])) if all(key in r for r in per.values()) else None
            st["arms"][arm_name] = {"runs": per, "mean": {k: mean(k) for k in ("sn", "pr_strict", "pr_lenient", "L2_sn", "ece", "n_reported_lnc",
                                                                                  "n_reported_mrna", "mrna_sn", "mrna_pr_strict")}}
            hits[arm_name], prs[arm_name] = H, R
        if tids:
            rng = np.random.default_rng(0)
            n = len(tids)
            boot = {"HM-H6": [], "HM-H": []}
            boot_pr = {"HM-H6": [], "HM-H": []}
            for _ in range(N_BOOT):
                idx = rng.integers(0, n, n)
                s, p = {}, {}
                for arm_name in ARMS:
                    runs = rng.integers(0, 3, 3)
                    s[arm_name] = float(np.mean([hits[arm_name][r][idx].mean() for r in runs]))
                    vals = []
                    for r in runs:
                        y = prs[arm_name][r]
                        vals.append(y[rng.integers(0, len(y), len(y))].mean() if len(y) else 0.0)
                    p[arm_name] = float(np.mean(vals))
                for k, (x, y) in (("HM-H6", ("HM", "H6")), ("HM-H", ("HM", "H"))):
                    boot[k].append(s[x] - s[y])
                    boot_pr[k].append(p[x] - p[y])
            m = {k: st["arms"][k]["mean"] for k in ARMS}
            st["delta"] = {k: {"sn": m[x]["sn"] - m[y]["sn"], "sn_ci95": ci(boot[k]),
                               "pr_strict": m[x]["pr_strict"] - m[y]["pr_strict"], "pr_strict_ci95": ci(boot_pr[k]),
                               "n_reported_ratio": m[x]["n_reported_lnc"] / max(m[y]["n_reported_lnc"], 1e-9)}
                           for k, (x, y) in (("HM-H6", ("HM", "H6")), ("HM-H", ("HM", "H")))}
        res["sets"][t] = st
        print(json.dumps({"set": t, "sn": {k: st["arms"][k]["mean"]["sn"] for k in ARMS}, "delta": st.get("delta", {}).get("HM-H6")}), flush=True)
    Z = [t for t, s in res["sets"].items() if s["in_verdict"] and s["group"] in ("B", "C")]
    d = [res["sets"][t]["delta"]["HM-H6"] for t in Z]
    n_pos = sum(x["sn_ci95"][0] > 0 for x in d)
    ratio = (sum(res["sets"][t]["arms"]["HM"]["mean"]["n_reported_lnc"] for t in Z)
             / max(sum(res["sets"][t]["arms"]["H6"]["mean"]["n_reported_lnc"] for t in Z), 1e-9))
    g = {k: cfg["arms"][k]["chr18_dev"]["v20"] for k in ("H6", "HM")}
    res["verdict"] = {
        "zero_shot_sets": Z,
        "mean_delta_sn_HM_H6": float(np.mean([x["sn"] for x in d])) if d else None,
        "sets_ci_lower_above_0": n_pos, "sets_needed": math.ceil(len(Z) / 2),
        "reported_ratio_HM_H6": ratio,
        "main": bool(d) and float(np.mean([x["sn"] for x in d])) >= 0.03 and n_pos >= math.ceil(len(Z) / 2) and ratio <= 1.3,
        "human_guard_chr18": {"H6": g["H6"], "HM": g["HM"],
                              "pass": g["HM"]["sn"] >= g["H6"]["sn"] - 0.02 and g["HM"]["pr_strict"] >= g["H6"]["pr_strict"] - 0.02},
        "mouse_chr19": res["sets"].get("mouse_chr19", {}).get("delta")}
    res["seconds"] = round(time.time() - t0)
    json.dump(res, open(f"{out_dir}/SCORE_V30.json", "w"), indent=1)
    print(json.dumps(res["verdict"]), flush=True)


if __name__ == "__main__":
    main()
