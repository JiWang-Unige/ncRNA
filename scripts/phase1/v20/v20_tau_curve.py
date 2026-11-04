#!/usr/bin/env python3
"""v2.0 chr18 development, exploratory (written after CALIB.json was read): Sn / strict Pr of the baseline (raw
confidence) and of v2.0 (calibrated probability, CALIB.json model) on chr18 for tau 0.20-0.90 in steps of 0.01,
three-run means.  chr18 only; nothing here touches chr9.  -> outputs/PHASE1/V20/dev/TAU_CURVE.json"""
import json
import os
import sys

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1/v20"))
import v20_calib as C  # noqa: E402

P, ruler = C.P, C.ruler


def main():
    cal = json.load(open(f"{C.DEV}/CALIB.json"))
    tab = json.load(open(C.HEX))
    truth = ruler.read_gtf(f"{ROOT}/outputs/PHASE1/DIAG_L2/chr18_v50.gtf", "chr18")
    codes = np.load(f"{ROOT}/data/p1/chr18_codes.npy", mmap_mode="r")
    arms = {"baseline": {}, "v20": {}}
    for run in C.RUNS:
        rd = f"{C.DEV}/chr18/{run}"
        genes = P.load_genes(f"{rd}/genes.json")
        rows = C.gene_rows(genes, np.load(f"{rd}/post.npy", mmap_mode="r"), codes, tab)
        X = np.array([[r[f] for f in C.FEATS] for r in rows], float)
        p = C.apply(cal["model"], X)
        sc = {r["k"]: float(s) for r, s in zip(rows, p)}
        arms["baseline"][run] = genes
        arms["v20"][run] = [dict(g, conf=sc.get(k, g["conf"])) for k, g in enumerate(genes)]
    out = {}
    for arm, by in arms.items():
        rows = []
        for t in [round(0.2 + 0.01 * i, 2) for i in range(71)]:
            sn, pr = [], []
            for run, gs in by.items():
                fin, _ = P.score_preds(truth, P.final_preds(gs, t), boot=0)
                sn.append(fin["lnc"]["L1_locus"]["sn"])
                pr.append(fin["lnc"]["L1_locus"]["pr_strict"] or 0.0)
            rows.append({"tau": t, "sn": float(np.mean(sn)), "pr_strict": float(np.mean(pr)), "sn_runs": sn, "pr_runs": pr})
        out[arm] = rows
    json.dump(out, open(f"{C.DEV}/TAU_CURVE.json", "w"), indent=1)
    for arm, rows in out.items():
        print(arm, [(r["tau"], round(r["sn"], 3), round(r["pr_strict"], 3)) for r in rows if round(r["tau"] * 100) % 5 == 0])


if __name__ == "__main__":
    main()
