#!/usr/bin/env python3
"""AUTH: the one-time test scoring of docs/auth/PLAN_AUTH.md (sections 5.4, 6 and 7), pooled over chr2 + chr6 + chr10 +
chr16 + chr22, with the model, q and the best rule fixed in outputs/AUTH/dev/MODEL.json (auth_dev.py).

Methods (same candidates, truth, greedy and GENCODE v41 same-strand coding mask):
  model       per TSS the highest 'true vs wrong' score among the 20 candidate ends, greedy by score, top q per chromosome
  model_rest  the same with the 'true vs rest' model and its q (reported only)
  P1          the model's end, everything else as R0 (greedy by ssCAGE tags, no q): the 'choose the end' contribution
  R0          ev_farthest (W24 frozen rule, v41 mask) -- main comparator
  R1-R6       section 5.4 (R2 / R6 with the q chosen in development); 'best rule' = MODEL.json best_rule
  isoquant / stringtie  PacBio rep1, reference-free (auth_queue.sh), gene spans >= 200 bp, v41 mask; reference only
Scoring: S50 one-to-one (reciprocal >= 0.5) and the W24 'unknown' rule; CI = 1 Mb block bootstrap, 2,000 resamples,
seed 0, the same resamples for every method (paired differences).
Gates (section 6): G1 model - R0: precision (known) diff >= +0.10 with CI low > 0 and recall diff >= -0.03;
  G2 model - best rule: precision (known) diff >= +0.05 with CI low > 0 and recall diff >= -0.03;
  P1 (secondary): P1 - R0 recall diff CI low > 0;  S (report only): precision (known) >= 0.48 and recall >= 0.40.
Also reported: per-chromosome numbers for every method, the three labels before / after filtering, standardised
coefficients, calibration (10 bins), and the secondary set chr9 / chr12 / chr13 / chr20 side by side with W24
(W24 truth with its exclusion zone and half cuts; Tiberius and v41 masks), plus the AUTH truth (no zone) there.

  auth_score.py --mode check  verify outputs/AUTH/FROZEN (every file's sha256) and that RESULT.json does not exist
  auth_score.py --mode prep   check, then extract the test chromosome sequences (IsoQuant reference)
  auth_score.py --mode test   check, then the one-time test read -> outputs/AUTH/RESULT.json
  auth_score.py --mode dev    the same code on the training chromosomes (in-sample, code check only); the W24
                              side-by-side and long-read code paths on chr18 -> outputs/AUTH/dev/SCORE_DEV.json
"""
import argparse
import json
import os
import re
import sys
import time
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import auth_lib as A                                       # noqa: E402
from auth_truth import TRAIN, TEST, SECONDARY              # noqa: E402

FROZEN, RESULT = f"{A.OUT}/FROZEN", f"{A.OUT}/RESULT.json"
REQUIRED = ("scripts/auth/auth_lib.py", "scripts/auth/auth_truth.py", "scripts/auth/auth_dev.py", "scripts/auth/auth_score.py",
            "scripts/auth/auth_queue.sh", "scripts/m1/w24/w24_lib.py", "scripts/m1/w24/e3_anchor.py",
            "scripts/m1/w24/dev_baselines.py", "scripts/m1/w24/w24_score.py", "scripts/m1/w24/wtc11_truth.py",
            "scripts/m1/simple_det.py", "scripts/phase1/b1/b1_strict.py", "scripts/phase1/s3/s3_features.py",
            "outputs/AUTH/dev/MODEL.json")
METHODS = ("model", "model_rest", "P1", "R0", "R1", "R2", "R3", "R4", "R5", "R6")


def check_frozen():
    assert os.path.exists(FROZEN), "outputs/AUTH/FROZEN missing: the author has not frozen the code"
    fz = json.load(open(FROZEN))
    missing = [p for p in REQUIRED if p not in fz]
    assert not missing, f"FROZEN lacks {missing}"
    bad = [p for p, h in fz.items() if A.sha256(f"{A.ROOT}/{p}") != h]
    assert not bad, f"differs from the frozen version: {bad}"
    assert not os.path.exists(RESULT), "the test has already been read once (outputs/AUTH/RESULT.json exists)"
    return fz


def read_longread(gtf, chrom):
    """Gene spans >= 200 bp of a StringTie / IsoQuant GTF (as w24_score.longread)."""
    g = defaultdict(lambda: [None, 10 ** 12, 0])
    for line in open(gtf):
        f = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(f) < 9 or f[0] != chrom or f[2] != "transcript" or f[6] not in "+-":
            continue
        x = g[(re.search(r'gene_id "([^"]+)"', f[8]).group(1), f[6])]
        x[0], x[1], x[2] = f[6], min(x[1], int(f[3])), max(x[2], int(f[4]))
    return [(st, s, e) for st, s, e in g.values() if e - s + 1 >= 200]


def lr_paths(chrom, mode):
    O = f"{A.L.W}/longread/{chrom}" if mode in ("dev", "secondary") else f"{A.OUT}/longread/{chrom}"   # W24 rep1 runs
    return {"stringtie": f"{O}/stringtie.gtf", "isoquant": f"{O}/isoquant/rep1/rep1.transcript_models.gtf"}


def units_all(P, M, masked=None):
    s_tw, s_tr = A.logit(M["true_vs_wrong"], P["X"]), A.logit(M["true_vs_rest"], P["X"])
    q = M["q"]
    acc = {"model": A.method_units(P, "model", score=s_tw, q=q["model"], masked=masked),
           "model_rest": A.method_units(P, "model", score=s_tr, q=q["model_rest"], masked=masked),
           "P1": A.method_units(P, "P1", score=s_tw, masked=masked)}
    for m in ("R0", "R1", "R3", "R4", "R5"):
        acc[m] = A.method_units(P, m, masked=masked)
    acc["R2"] = A.method_units(P, "R2", q=q["R2"], masked=masked)
    acc["R6"] = A.method_units(P, "R6", q=q["R6"], masked=masked)
    return acc, s_tw


def evaluate(chroms, M, pc, obs, mode):
    """Main table on one chromosome set (AUTH truth, v41 mask)."""
    lengths = {c: A.chrom_len(c) for c in chroms}
    bid = A.block_ids(chroms, lengths)
    res = defaultdict(dict)
    starts = defaultdict(dict)
    ts_by = {}
    three = {"pairs": defaultdict(int), "model_before_mask_q": defaultdict(int), "model_after_mask": defaultdict(int)}
    cal_pairs, cal_units = [], []
    lr_ok = True
    for c in chroms:
        truth = A.load_truth(c)
        ts = truth["ts"]
        ts_by[c] = ts
        ocv = A.obs_cover(obs[c])
        P = A.build_pairs(c, pc[c], truth, obs[c])
        acc, s_tw = units_all(P, M)
        for m, a in acc.items():
            res[m][c] = A.score_units(P, a, ts, ocv)
            starts[m][c] = P["s"][a]
        # three labels before / after filtering
        for k, n in zip(*np.unique(P["label"], return_counts=True)):
            three["pairs"][int(k)] += int(n)
        nomask = np.zeros(len(P["s"]), bool)
        for key, a in (("model_before_mask_q", A.method_units(P, "model", score=s_tw, q=1.0, masked=nomask)),
                       ("model_after_mask", A.method_units(P, "model", score=s_tw, q=1.0))):
            r = A.score_units(P, a, ts, ocv)
            for k in (0, 1, 2):
                three[key][k] += int((r["lab"] == k).sum())
        tw = P["label"] != 1
        cal_pairs.append((A.sigmoid(s_tw[tw]), (P["label"][tw] == 2).astype(float)))
        r = res["model"][c]
        cal_units.append((A.sigmoid(s_tw[acc["model"]]), (r["lab"] == 2).astype(float), r["lab"] != 1))
        # long-read references
        paths = lr_paths(c, mode)
        if all(os.path.exists(p) for p in paths.values()):
            pcv = {st: A.cover(pc[c][st]) for st in A.STR}
            for tool, p in paths.items():
                sp = [(st, s, e) for st, s, e in read_longread(p, c) if not A.cov_bp(pcv[st], s, e) > 0]
                sp.sort(key=lambda x: (x[1], x[2]))
                res[tool][c] = A.score_spans(ts, ocv, [x[1] for x in sp], [x[2] for x in sp], [x[0] for x in sp])
                starts[tool][c] = np.array([x[1] for x in sp], np.int64)
        else:
            lr_ok = False
        print(c, {m: (res[m][c]["tp"], res[m][c]["n_pred"]) for m in ("model", "R0")}, flush=True)
    methods = [m for m in res if all(c in res[m] for c in chroms)]
    T = {m: A.block_table(res[m], ts_by, starts[m], bid) for m in methods}
    W = A.boot_weights(len(bid))
    pooled = {m: A.metrics([res[m][c] for c in chroms]) for m in methods}
    per_chrom = {m: {c: A.metrics([res[m][c]]) for c in chroms} for m in methods}
    ci = {m: {k: np.nanpercentile(v, [2.5, 97.5]).tolist() for k, v in A.boot_metrics(W, T[m]).items()} for m in methods}
    best = M["best_rule"]
    pairs_ = {"model-R0": A.paired(W, T["model"], T["R0"]), f"model-{best}": A.paired(W, T["model"], T[best]),
              "P1-R0": A.paired(W, T["P1"], T["R0"]), "model_rest-R0": A.paired(W, T["model_rest"], T["R0"]),
              f"model_rest-{best}": A.paired(W, T["model_rest"], T[best])}
    if "isoquant" in T:
        pairs_["model-isoquant"] = A.paired(W, T["model"], T["isoquant"])
    tp_not_obs = {m: sum(res[m][c]["tp_not_observed"] for c in chroms) for m in methods}
    cp = (np.concatenate([x[0] for x in cal_pairs]), np.concatenate([x[1] for x in cal_pairs]))
    cu_k = [x[2] for x in cal_units]
    cu = (np.concatenate([x[0][k] for x, k in zip(cal_units, cu_k)]), np.concatenate([x[1][k] for x, k in zip(cal_units, cu_k)]))
    names = {0: "wrong", 1: "uncertain", 2: "true"}
    labels = {k: {names[i]: v.get(i, 0) for i in (0, 1, 2)} for k, v in three.items()}
    labels.update({f"final_{m}": {"true": pooled[m]["tp"], "uncertain": pooled[m]["unknown"], "wrong": pooled[m]["wrong"]} for m in methods})
    return {"chroms": list(chroms), "n_truth": sum(len(v) for v in ts_by.values()), "n_blocks": len(bid), "pooled": pooled,
            "ci95": ci, "per_chrom": per_chrom, "paired": pairs_, "labels": labels, "longread_present": lr_ok,
            "tp_not_observed": tp_not_obs,
            "calibration": {"pairs_true_vs_wrong": A.calibration(*cp), "model_units_known": A.calibration(*cu)}}


def verdict(ev, M):
    best = M["best_rule"]
    g1, g2, p1 = ev["paired"]["model-R0"], ev["paired"][f"model-{best}"], ev["paired"]["P1-R0"]
    o = ev["pooled"]["model"]

    def gate(d, thr):
        pk, rc = d["prec_known"], d["recall"]
        return bool(pk["diff"] >= thr and pk["ci95"][0] > 0 and rc["diff"] >= -0.03)
    out = {"G1": gate(g1, 0.10), "G2": gate(g2, 0.05), "P1": bool(p1["recall"]["ci95"][0] > 0),
           "S": bool((o["prec_known"] or 0) >= 0.48 and (o["recall"] or 0) >= 0.40), "best_rule": best,
           "G1_numbers": g1, "G2_numbers": g2, "P1_numbers": p1}
    out["conclusion"] = ("G1+G2 pass: authenticity of cheap end evidence is learnable beyond simple rules" if out["G1"] and out["G2"]
                         else "G1 pass, G2 fail: the rules do the work, not the learning" if out["G1"]
                         else "G1 fail: cheap evidence alone cannot separate real from wrong units")
    return out


def w24_side(chroms, M, pc, mode):
    """Secondary set (test) or chr18 (dev): W24 truth (exclusion zone, half cuts), Tiberius and v41 masks, W24 numbers."""
    import w24_score as WS
    L = A.L
    rows = defaultdict(list)
    for c in chroms:
        truth = L.load_truth(c)
        cuts = WS.cuts_of(c)
        truth["ts"] = [t for t in truth["ts"] if not any(t[2] - 1 < x < t[3] for x in cuts)]
        zone = list(truth["zone"]) + [(x - 1, x + 1) for x in cuts]
        tib = WS.read_tib(f"{L.W}/tiberius/{c}.gtf" if c in ("chr12", "chr13") else f"{A.ROOT}/outputs/PHASE1/BASELINES/tiberius/{c}.gtf", c)
        ocv = A.obs_cover(A.observed([c])[c])
        P = A.build_pairs(c, pc[c])
        for mk, masked in (("tiberius", A.masks(P, tib)), ("v41", A.masks(P))):
            acc, _ = units_all(P, M, masked)
            for m, a in acc.items():
                a = np.array([i for i in a if not L.in_zone(zone, int(P["s"][i]), int(P["e"][i]))], np.int64)
                rows[f"{m}|{mk}"].append(A.score_units(P, a, truth["ts"], ocv))
    out = {k: A.metrics(v) for k, v in rows.items()}
    if mode == "test":
        w = json.load(open(f"{L.W}/RESULT.json"))["pooled"]
        ref = {m: w[f"{m}|f1.0_s0"] for m in ("ev_farthest", "ours", "pair", "isoquant", "stringtie")}
    else:
        g = json.load(open(f"{L.W}/dev/E3_chr18_f1.0_s0.json"))["grid"]
        ref = {"ev_farthest": g["d50000_K3_ev-farthest+tib"]["evidence_only"], "pair": g["d50000_K1_ev-nearest+tib"]["evidence_only"]}
    return {"chroms": list(chroms), "auth": out, "w24_reference": ref,
            "note": "W24 numbers use Tiberius; R0|tiberius must equal W24 ev_farthest, R1|tiberius W24 pair"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True, choices=["check", "prep", "test", "dev"])
    a = ap.parse_args()
    t0 = time.time()
    M = json.load(open(f"{A.OUT}/dev/MODEL.json"))
    if a.mode in ("check", "prep", "test"):
        fz = check_frozen()
        print("FROZEN ok:", len(fz), "files", flush=True)
        if a.mode == "check":
            return
        A.open_test()
        A.extract_seq(TEST)
        if a.mode == "prep":
            return
        with open(f"{A.OUT}/test_attempts.log", "a") as fh:
            fh.write(time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()) + " --mode test started\n")
        import auth_truth as AT
        reps = AT.read_reps(TEST + SECONDARY)
        n = {}
        for c in TEST + SECONDARY:
            loci = AT.build(reps, c)
            AT.write(c, loci)
            n[c] = sum(x["k"] == "l" for x in loci)
        assert sum(n[c] for c in TEST) == AT.EXPECTED["test"], n
        main_set, side = TEST, SECONDARY
        A.extract_seq(SECONDARY)
    else:
        main_set, side = TRAIN, ("chr18",)
    pc = A.v41_pc(tuple(main_set) + tuple(side))
    obs = A.observed(tuple(main_set) + tuple(side))
    ev = evaluate(main_set, M, pc, obs, a.mode)
    if a.mode == "test":
        assert ev["longread_present"], "rep1 StringTie / IsoQuant outputs missing: run auth_queue.sh"
    out = {"plan": "docs/auth/PLAN_AUTH.md", "mode": a.mode, "model": f"{A.OUT}/dev/MODEL.json", "q": M["q"], "best_rule": M["best_rule"],
           "main": ev, "verdict": verdict(ev, M),
           "coef_standardised": {k: dict(zip(["intercept"] + M["features"], M[k]["w"])) for k in ("true_vs_wrong", "true_vs_rest")},
           "w24_side_by_side": w24_side(side, M, pc, a.mode)}
    if a.mode == "test":
        out["secondary_auth_truth_v41"] = evaluate(SECONDARY, M, pc, obs, "secondary")
        out["frozen"] = json.load(open(FROZEN))
    out["seconds"] = round(time.time() - t0, 1)
    path = RESULT if a.mode == "test" else f"{A.OUT}/dev/SCORE_DEV.json"
    json.dump(out, open(path, "w"), indent=1, default=float)
    v = out["verdict"]
    print(json.dumps({k: v[k] for k in ("G1", "G2", "P1", "S", "best_rule", "conclusion")}, indent=1))
    for m in METHODS + ("isoquant", "stringtie"):
        if m in ev["pooled"]:
            t = ev["pooled"][m]
            print(f"{m:10s} recall {t['recall']:.3f}  prec_known {t['prec_known']:.3f}  strict {t['prec_strict']:.3f}  n {t['n_pred']}")


if __name__ == "__main__":
    main()
