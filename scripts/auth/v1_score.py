#!/usr/bin/env python3
"""AUTH v1: the one-time H1 transfer test (docs/auth/v1/PLAN_V1.md sections 5.4-7), with the model, standardisation, q
and best rule frozen in outputs/AUTH_V1/dev/MODEL_V1.json (v1_dev.py, WTC11-DPI).

Main (gates): H1 test set = chr2, 6, 10, 16, 22 + chr9, 12, 13, 20 (author 2026-11-24), pooled.  H1 evidence, observed
spans, truth loci and pairs on these chromosomes are generated inside this run only.  Features normalised with H1's own
CAGE total (v1_lib), standardised with the frozen WTC11-DPI mu / sd.
Methods: model / model_rest / P1 / R0-R6 exactly as auth_score.units_all (R2 / R6 with the frozen q).
Gates: T1 model - R0: precision (known) diff >= +0.10, CI low > 0, recall diff >= -0.03;
       T2 model - best rule: precision (known) diff >= +0.05, CI low > 0, recall diff >= -0.03.
Scoring: S50 one-to-one + 'unknown' rule; 1 Mb block bootstrap, 2,000 resamples, seed 0, paired.
Secondary (report only): S1 the frozen model on WTC11-DPI chr2, 6, 10, 16, 22 (AUTH truth) next to AUTH RESULT.json;
  S2 H1 relearning (v1_dev.develop('h1') on H1's 14 training chromosomes, applied to the H1 test set);
  S3 calibration of the frozen model on H1 (in the main table).
  v1_score.py --mode freeze  (only after the author's go-ahead) write outputs/AUTH_V1/FROZEN: sha256 of REQUIRED
  v1_score.py --mode check | --mode test     -> outputs/AUTH_V1/RESULT_V1.json
  v1_score.py --mode dev   (code check, WTC11-DPI training chromosomes in-sample + S2 code path on WTC11) -> dev/SCORE_DEV_V1.json
"""
import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v1_lib as V                                          # noqa: E402
import v1_dev as VD                                         # noqa: E402

A = V.A
import auth_score as AS                                     # noqa: E402  (units_all, verdict gate logic re-used)

FROZEN, RESULT = f"{V.OUTV}/FROZEN", f"{V.OUTV}/RESULT_V1.json"
REQUIRED = ("scripts/auth/v1_lib.py", "scripts/auth/v1_dev.py", "scripts/auth/v1_score.py", "scripts/auth/v1_check.py",
            "scripts/auth/auth_lib.py", "scripts/auth/auth_truth.py", "scripts/auth/auth_dev.py", "scripts/auth/auth_score.py",
            "scripts/m1/w24/w24_lib.py", "scripts/m1/w24/e3_anchor.py", "scripts/m1/w24/wtc11_truth.py", "scripts/m1/simple_det.py",
            "scripts/phase1/b1/b1_strict.py", "scripts/phase1/s3/s3_features.py",
            "outputs/AUTH_V1/dev/MODEL_V1.json", "outputs/AUTH_V1/dev/DEV_V1.json", "outputs/AUTH_V1/CHECK.json")
METHODS = ("model", "model_rest", "P1", "R0", "R1", "R2", "R3", "R4", "R5", "R6")


def check_frozen():
    assert os.path.exists(FROZEN), "outputs/AUTH_V1/FROZEN missing: the author has not frozen v1"
    fz = json.load(open(FROZEN))
    missing = [p for p in REQUIRED if p not in fz]
    assert not missing, f"FROZEN lacks {missing}"
    bad = [p for p, h in fz.items() if A.sha256(f"{A.ROOT}/{p}") != h]
    assert not bad, f"differs from the frozen version: {bad}"
    assert not os.path.exists(RESULT), "the v1 test has already been read once (RESULT_V1.json exists)"
    dev = json.load(open(f"{V.OUTV}/dev/DEV_V1.json"))
    assert dev["D0"]["pass"], "D0 failed: H1 must not be read (plan 6.1)"
    return fz


def evaluate(cell, chroms, M, lib):
    """Main table of one cell on one chromosome set with model M (frozen or relearned)."""
    V.use(cell)
    pc = A.v41_pc(chroms)
    obs = A.observed(chroms)
    lengths = {c: A.chrom_len(c) for c in chroms}
    bid = A.block_ids(chroms, lengths)
    res, starts, ts_by = {m: {} for m in METHODS}, {m: {} for m in METHODS}, {}
    three = {0: 0, 1: 0, 2: 0}
    cal = []
    for c in chroms:
        P = V.load_pairs(cell, c, pc, lib)
        ts = A.load_truth(c)["ts"]
        ts_by[c] = ts
        ocv = A.obs_cover(obs[c])
        acc, s_tw = AS.units_all(P, M)
        for m, a in acc.items():
            res[m][c] = A.score_units(P, a, ts, ocv)
            starts[m][c] = P["s"][a]
        for k, n in zip(*np.unique(P["label"], return_counts=True)):
            three[int(k)] += int(n)
        tw = P["label"] != 1
        cal.append((A.sigmoid(s_tw[tw]), (P["label"][tw] == 2).astype(float)))
        print(cell, c, {m: (res[m][c]["tp"], res[m][c]["n_pred"]) for m in ("model", "R0")}, flush=True)
    T = {m: A.block_table(res[m], ts_by, starts[m], bid) for m in METHODS}
    W = A.boot_weights(len(bid))
    pooled = {m: A.metrics([res[m][c] for c in chroms]) for m in METHODS}
    best = M["best_rule"]
    paired = {"model-R0": A.paired(W, T["model"], T["R0"]), f"model-{best}": A.paired(W, T["model"], T[best]),
              "P1-R0": A.paired(W, T["P1"], T["R0"]), "model_rest-R0": A.paired(W, T["model_rest"], T["R0"])}
    cp = (np.concatenate([x[0] for x in cal]), np.concatenate([x[1] for x in cal]))
    return {"cell": cell, "chroms": list(chroms), "n_truth": sum(len(v) for v in ts_by.values()), "n_blocks": len(bid),
            "pooled": pooled, "per_chrom": {m: {c: A.metrics([res[m][c]]) for c in chroms} for m in METHODS},
            "ci95": {m: {k: np.nanpercentile(v, [2.5, 97.5]).tolist() for k, v in A.boot_metrics(W, T[m]).items()} for m in METHODS},
            "paired": paired, "labels_pairs": {"wrong": three[0], "uncertain": three[1], "true": three[2]},
            "labels_units": {m: {"true": pooled[m]["tp"], "uncertain": pooled[m]["unknown"], "wrong": pooled[m]["wrong"]} for m in METHODS},
            "calibration_pairs_true_vs_wrong": A.calibration(*cp), "_T": T, "_W": W}


def verdict(ev, M):
    best = M["best_rule"]
    t1, t2 = ev["paired"]["model-R0"], ev["paired"][f"model-{best}"]

    def gate(d, thr):
        pk, rc = d["prec_known"], d["recall"]
        return bool(pk["diff"] >= thr and pk["ci95"][0] > 0 and rc["diff"] >= -0.03)
    out = {"T1": gate(t1, 0.10), "T2": gate(t2, 0.05), "best_rule": best, "T1_numbers": t1, "T2_numbers": t2}
    out["conclusion"] = ("T1+T2 pass: the WTC11-DPI model transfers to H1 and beats the simple rules there" if out["T1"] and out["T2"]
                         else "T1 pass, T2 fail: transfer holds, but the rules transfer equally well" if out["T1"]
                         else "T1 fail: what was learned on WTC11 does not carry over to H1")
    return out


def strip(ev):
    return {k: v for k, v in ev.items() if not k.startswith("_")}


def prepare_h1():
    """H1 test-set inputs (inside the one-time run): truth loci, observed spans, DPI evidence."""
    n = V.write_truth("h1", V.H1_TEST)
    cv = V.check_values()["h1"]["truth_lnc_counts"]
    assert sum(n[c] for c in V.TEST) == cv["test"] and sum(n[c] for c in V.SECONDARY) == cv["secondary"], n
    V.write_observed("h1", V.H1_TEST)
    lib = V.write_evidence("h1", V.H1_TEST)
    assert lib == V.check_values()["h1"]["cage_lib_tags"], lib
    return n, lib


def prepare_wtc11_test():
    """S1 inputs: AUTH truth / observed of chr2, 6, 10, 16, 22 (symlinks) and WTC11 DPI evidence."""
    out = V.CELLS["wtc11"]["out"]
    for c in V.TEST:
        for src, dst in ((f"{A.ROOT}/outputs/AUTH/truth/{c}.json", f"{out}/truth/{c}.json"),
                         (f"{A.ROOT}/outputs/AUTH/evidence/{c}.observed_rep23.json", f"{out}/evidence/{c}.observed_rep23.json")):
            assert os.path.exists(src), src
            if not os.path.exists(dst):
                os.symlink(src, dst)
    return V.write_evidence("wtc11", V.TEST)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True, choices=["freeze", "check", "test", "dev"])
    ap.add_argument("--procs", type=int, default=4)
    a = ap.parse_args()
    t0 = time.time()
    M = json.load(open(f"{V.OUTV}/dev/MODEL_V1.json"))
    lib_w = M["normalisation_constants"]["cage_tags_per_sample"]
    if a.mode == "freeze":
        assert not os.path.exists(FROZEN), "already frozen"
        assert json.load(open(f"{V.OUTV}/dev/DEV_V1.json"))["D0"]["pass"], "D0 failed: nothing is frozen (plan 6.1)"
        fz = {p: A.sha256(f"{A.ROOT}/{p}") for p in REQUIRED}
        json.dump(fz, open(FROZEN, "w"), indent=1)
        print(json.dumps(fz, indent=1))
        return
    if a.mode == "dev":
        ev = evaluate("wtc11", V.TRAIN, M, lib_w)
        m2, d2 = VD.develop("wtc11", a.procs)            # S2 code path, run on WTC11 (in-sample code check only)
        same = all(np.allclose(m2["true_vs_wrong"][k], M["true_vs_wrong"][k]) for k in ("mu", "sd", "w")) and m2["q"] == M["q"]
        out = {"mode": "dev", "main_in_sample": strip(ev), "verdict_in_sample": verdict(ev, M), "s2_code_path_reproduces_model": bool(same),
               "seconds": round(time.time() - t0, 1)}
        json.dump(out, open(f"{V.OUTV}/dev/SCORE_DEV_V1.json", "w"), indent=1, default=float)
        print(json.dumps({"s2_code_path_reproduces_model": same, **{k: out["verdict_in_sample"][k] for k in ("T1", "T2", "best_rule")}}))
        return
    fz = check_frozen()
    print("FROZEN ok:", len(fz), "files", flush=True)
    if a.mode == "check":
        return
    with open(f"{V.OUTV}/test_attempts.log", "a") as fh:
        fh.write(time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()) + " --mode test started\n")
    A.open_test()
    V.ensure_seq_link("h1")
    n_h1, lib_h1 = prepare_h1()
    ev = evaluate("h1", V.H1_TEST, M, lib_h1)
    vd = verdict(ev, M)
    # S2: relearning on H1's training chromosomes, applied to the same H1 test set
    m2, d2 = VD.develop("h1", a.procs)
    ev2 = evaluate("h1", V.H1_TEST, m2, lib_h1)
    s2 = {"model_h1": d2, "q": m2["q"], "best_rule_h1": m2["best_rule"], "pooled": ev2["pooled"],
          "relearned-R0": ev2["paired"]["model-R0"],
          "relearned-transferred": A.paired(ev["_W"], ev2["_T"]["model"], ev["_T"]["model"])}
    # S1: the frozen model on WTC11-DPI test chromosomes
    lib_w2 = prepare_wtc11_test()
    assert lib_w2 == lib_w, lib_w2
    ev1 = evaluate("wtc11", V.TEST, M, lib_w)
    auth = json.load(open(f"{A.ROOT}/outputs/AUTH/RESULT.json"))["main"]["pooled"]
    out = {"plan": "docs/auth/v1/PLAN_V1.md", "mode": "test", "model": "outputs/AUTH_V1/dev/MODEL_V1.json", "q": M["q"],
           "best_rule": M["best_rule"], "h1_truth_lnc_per_chrom": n_h1, "h1_cage_tags_per_sample": lib_h1,
           "main_h1": strip(ev), "verdict": vd,
           "coef_standardised": {k: dict(zip(["intercept"] + M["features"], M[k]["w"])) for k in ("true_vs_wrong", "true_vs_rest")},
           "S1_wtc11_dpi_test": {"pooled": ev1["pooled"], "paired": ev1["paired"], "per_chrom": ev1["per_chrom"],
                                 "auth_result_pooled": {m: auth[m] for m in ("model", "R0", "R2") if m in auth}},
           "S2_h1_relearn": s2, "frozen": fz, "seconds": round(time.time() - t0, 1)}
    json.dump(out, open(RESULT, "w"), indent=1, default=float)
    print(json.dumps({k: vd[k] for k in ("T1", "T2", "best_rule", "conclusion")}, indent=1))
    for m in METHODS:
        t = ev["pooled"][m]
        print(f"{m:10s} recall {t['recall']:.3f}  prec_known {t['prec_known']:.3f}  strict {t['prec_strict']:.3f}  n {t['n_pred']}")


if __name__ == "__main__":
    main()
