#!/usr/bin/env python3
"""M1 experiment 3 (docs/m1/PLAN_M1_AG.md section 5): Loc1k of AlphaGenome fold models vs Puffin and a CpG control.

Rows = outputs/M1/AG/tss.json (DEDEND rows with fold assignment); AlphaGenome = outputs/M1/AG/colab_<run>/raw/<model>.npz.
Hit = |argmax - TSS| <= 1 kb inside TSS +- 20 kb.  Primary score: argmax of the mean over the strand's CAGE tracks.
Controls on the same rows: Puffin (DEDEND FANTOM_CAGE track, outputs/PHASE1/DEDEND/puffin/<chrom>.npy, same rule as
dedend_metrics.puffin_dist) and CpG (argmax of the CpG count in a centred 500-bp window).  CpG quartiles over all
lncRNA rows of the chromosome group (as DEDEND, so the Puffin reference values keep their meaning).
Main set: chr8 + chr9 rows that have a fold model; sensitivity: test_only rows; secondary: chr18, chr20.
  ag_fold_metrics.py [--run full] -> outputs/M1/AG/RESULT.json
"""
import argparse
import json
import os
import sys

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
for d in ("src", "scripts/phase1", "scripts/phase1/s3", "scripts/phase1/ends", "scripts/phase1/dedend"):
    sys.path.insert(0, os.path.join(ROOT, d))
import dedend_metrics as DM  # noqa: E402
import s3_features as S3  # noqa: E402

O = f"{ROOT}/outputs/M1/AG"
HALF, TOL, N_BOOT = 20_000, 1000, 2000
GROUPS = {"main_chr8_chr9": ("chr8", "chr9"), "chr18": ("chr18",), "chr20": ("chr20",)}


def cpg_dist(seq, r):
    x = r["TSS"]
    a, b = max(0, x - HALF - 250), min(len(seq), x + HALF + 251)
    s = seq[a:b].upper()
    cg = np.zeros(len(s), np.int32)
    i = s.find("CG")
    while i >= 0:
        cg[i] = 1
        i = s.find("CG", i + 1)
    cs = np.concatenate([[0], np.cumsum(cg)])
    pos = np.arange(max(a, x - HALF), min(b, x + HALF + 1))
    lo = np.clip(pos - 250 - a, 0, len(s))
    hi = np.clip(pos + 251 - a, 0, len(s))
    cnt = cs[hi] - cs[lo]
    return abs(int(pos[int(np.argmax(cnt))]) - x)


def hit(d):
    return None if d is None else float(d <= TOL)


def boot_paired(a, b, seed=0):
    a, b = np.asarray(a, float), np.asarray(b, float)
    rng = np.random.default_rng(seed)
    n = len(a)
    if not n:
        return None
    d = [a[i].mean() - b[i].mean() for i in (rng.integers(0, n, n) for _ in range(N_BOOT))]
    return {"diff": float(a.mean() - b.mean()), "ci95": np.percentile(d, [2.5, 97.5]).tolist()}


def boot_ci(a, seed=0):
    a = np.asarray(a, float)
    rng = np.random.default_rng(seed)
    if not len(a):
        return None
    m = [a[rng.integers(0, len(a), len(a))].mean() for _ in range(N_BOOT)]
    return np.percentile(m, [2.5, 97.5]).tolist()


def summ(rr, key):
    v = [r[key] for r in rr if r.get(key) is not None]
    return {"n": len(v), "loc1k": float(np.mean(v)) if v else None}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default="full")
    a = ap.parse_args()
    rows = json.load(open(f"{O}/tss.json"))
    by_id = {r["id"]: r for r in rows}
    raw = f"{O}/colab_{a.run}/raw"
    for m in ("fold-0", "fold-1", "fold-2", "fold-3"):
        f = f"{raw}/{m}.npz" if a.run != "smoke" else f"{raw}/smoke_{m}.npz"
        if not os.path.exists(f):
            continue
        z = np.load(f)
        for j, rid in enumerate(z["rid"]):
            r = by_id[str(rid)]
            assert r["model"] == m, (rid, r["model"], m)
            r["ag_mean"] = hit(abs(int(z["agg_cage_mean_arg"][j])))
            r["ag_max"] = hit(abs(int(z["agg_cage_max_arg"][j])))
            r["ag_pro"] = hit(abs(int(z["agg_pro_mean_arg"][j])))
            r["ag_best_track"] = hit(int(np.abs(z["cage_arg"][j]).min()))   # oracle over tracks: optimistic, descriptive
    for c in ("chr8", "chr9", "chr18", "chr20"):
        seq = S3.read_fasta(f"{ROOT}/data/human/{c}.fa")
        pf = f"{ROOT}/outputs/PHASE1/DEDEND/puffin/{c}.npy"
        tr = np.load(pf, mmap_mode="r") if os.path.exists(pf + ".done") else None
        for r in rows:
            if r["chrom"] != c or r.get("ag_mean") is None:
                continue
            r["cpg_ctrl"] = hit(cpg_dist(seq, r))
            if tr is not None:
                r["puffin"] = hit(DM.puffin_dist(tr, r))
    res = {"plan": "docs/m1/PLAN_M1_AG.md", "run": a.run, "groups": {}}
    for g, cs in GROUPS.items():
        allg = [r for r in rows if r["chrom"] in cs]
        lnc_all = [r for r in allg if r["k"] == "l"]
        q = np.quantile([r["cpg"] for r in lnc_all], [0.25, 0.5, 0.75])
        rd = [r for r in allg if r.get("ag_mean") is not None]
        sets = {"read": rd, "test_only": [r for r in rd if r["test_only"]]}
        out = {"cpg_quartile_cuts_all_lnc": q.tolist(), "n_rows_group": len(allg), "n_read": len(rd)}
        for sname, ss in sets.items():
            sel = {"G0_mrna_cage": [r for r in ss if r["k"] == "m" and r["TSS_ev"]],
                   "G1_lnc_lowcpg": [r for r in ss if r["k"] == "l" and r["cpg"] <= q[1]],
                   "G2_lnc_cage": [r for r in ss if r["k"] == "l" and r["TSS_ev"]],
                   "lnc_all": [r for r in ss if r["k"] == "l"], "mrna_all": [r for r in ss if r["k"] == "m"]}
            for qi in range(4):
                lo, hi = ([-np.inf, *q, np.inf])[qi], ([-np.inf, *q, np.inf])[qi + 1]
                sel[f"lnc_cpg_q{qi + 1}"] = [r for r in ss if r["k"] == "l" and lo < r["cpg"] <= hi]
            blk = {}
            for k, rr in sel.items():
                blk[k] = {key: summ(rr, key) for key in ("ag_mean", "ag_max", "ag_pro", "ag_best_track", "puffin", "cpg_ctrl")}
                blk[k]["ag_mean"]["ci95"] = boot_ci([r["ag_mean"] for r in rr])
                for ctrl in ("puffin", "cpg_ctrl"):
                    pr = [r for r in rr if r.get(ctrl) is not None]
                    blk[k][f"ag_minus_{ctrl}"] = boot_paired([r["ag_mean"] for r in pr], [r[ctrl] for r in pr])
            out[sname] = blk
        res["groups"][g] = out
    M = res["groups"]["main_chr8_chr9"]["read"]
    g0, g1, g2 = M["G0_mrna_cage"]["ag_mean"], M["G1_lnc_lowcpg"], M["G2_lnc_cage"]["ag_mean"]["loc1k"]
    v = {"G0": {"loc1k": g0["loc1k"], "valid": bool(g0["loc1k"] is not None and g0["loc1k"] >= 0.60)}}
    x = g1["ag_mean"]
    pa, pc = g1["ag_minus_puffin"], g1["ag_minus_cpg_ctrl"]
    if x["loc1k"] >= 0.25 and pa and pc and pa["ci95"][0] > 0 and pc["ci95"][0] > 0:
        o1 = "A"
    elif x["loc1k"] <= 0.12 and x["ci95"][1] < 0.25:
        o1 = "B"
    else:
        o1 = "grey"
    v["G1"] = {"loc1k": x["loc1k"], "ci95": x["ci95"], "minus_puffin": pa, "minus_cpg": pc, "outcome": o1}
    v["G2"] = {"loc1k": g2, "outcome": "A" if g2 >= 0.50 else "B" if g2 <= 0.25 else "grey"}
    res["verdict"] = v
    json.dump(res, open(f"{O}/RESULT{'' if a.run == 'full' else '_' + a.run}.json", "w"), indent=1)
    print(json.dumps(v, indent=1))


if __name__ == "__main__":
    main()
