#!/usr/bin/env python3
"""AUTH v1 development (docs/auth/v1/PLAN_V1.md sections 5.3, 6.1 and 7 step 3): WTC11-DPI, 14 training chromosomes.

1. Inputs = the caches written by v1_check.py in outputs/AUTH_V1/wtc11_dpi (DPI evidence, AUTH truth and observed spans
   via symlinks, raw pairs).  Checks: the DPI evidence recomputed by v1_lib equals the caches on every training
   chromosome; CAGE tags per sample equal CHECK.json; R0 equals CHECK.json (wtc11_dpi R0_train); the cached raw
   log_cage / log_anti columns equal log1p of the pair tags (asserted in v1_lib.normalise).
2. Features normalised (v1_lib): log_cage and log_anti -> log(1 + CPM), n_samples -> fraction of samples.
3. LOCO exactly as auth_dev.py: L2 logistic 'true vs wrong' ('true vs rest' reported only); q of the model: highest
   precision (known) with pooled LOCO recall >= R0-DPI training recall; R2 / R6 q and best-rule eligibility: >= R0 - 0.02.
4. D0 (plan 6.1): model - R0-DPI precision (known) >= +0.10 on the pooled LOCO units, a q found under the recall bound,
   and model recall >= R0-DPI training recall.  If D0 fails nothing is frozen and H1 is not read.
5. Final models on all 14 chromosomes; coefficients, standardisation (mu / sd of the WTC11-DPI rows), q, best rule and
   the WTC11 normalisation constants -> outputs/AUTH_V1/dev/MODEL_V1.json; details -> DEV_V1.json.
develop(cell) is also used by v1_score.py for the post-freeze secondary analysis S2 (H1 relearning).
  v1_dev.py [--procs 4]
"""
import argparse
import json
import math
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v1_lib as V                          # noqa: E402

A = V.A
import auth_dev as AD                       # noqa: E402  (choose_q, QS, RULES, FIXED, slacks: unchanged AUTH constants)

TRAIN = V.TRAIN
DEVD = f"{V.OUTV}/dev"
Q_CONSTRAINT = "model / model_rest: LOCO recall >= R0-DPI training recall; R2 / R6 and best-rule eligibility: >= R0-DPI - 0.02 (PLAN_V1.md 5.3 / 5.4)"
_P = {}


def fit_one(args):
    h, kind = args
    X, y = [], []
    for c in TRAIN:
        if c == h:
            continue
        lab = _P[c]["label"]
        sel = lab != 1 if kind == "tw" else np.ones(len(lab), bool)
        X.append(_P[c]["X"][sel])
        y.append(lab[sel] == 2)
    return h, kind, A.fit(np.vstack(X), np.concatenate(y))


def develop(cell, procs=4, log=print):
    """LOCO development on the cell's 14 training chromosomes -> (model dict, details dict)."""
    lib = V.check_values()[V.CELLS[cell]["check_key"]]["cage_lib_tags"]
    V.use(cell)
    pc = A.v41_pc(TRAIN)
    P = {c: V.load_pairs(cell, c, pc, lib) for c in TRAIN}
    _P.clear()
    _P.update(P)
    truth = {c: A.load_truth(c)["ts"] for c in TRAIN}
    obs = A.observed(TRAIN)
    ocv = {c: A.obs_cover(obs[c]) for c in TRAIN}
    log("pairs loaded", cell, sum(len(P[c]["s"]) for c in TRAIN))

    with Pool(procs) as pool:
        fits = pool.map(fit_one, [(h, k) for h in TRAIN for k in ("tw", "tr")])
    oof = {k: {h: A.logit(m, P[h]["X"]) for h, kk, m in fits if kk == k} for k in ("tw", "tr")}
    log("LOCO fits done", cell)

    acc = {c: {} for c in TRAIN}
    curves = {m: {} for m in ("R2", "R6", "model", "model_rest")}
    res_fixed = {m: {} for m in AD.FIXED + ("P1",)}
    for c in TRAIN:
        Pc = P[c]
        for m in AD.FIXED:
            acc[c][m] = A.method_units(Pc, m)
            res_fixed[m][c] = A.score_units(Pc, acc[c][m], truth[c], ocv[c])
        acc[c]["P1"] = A.method_units(Pc, "P1", score=oof["tw"][c])
        res_fixed["P1"][c] = A.score_units(Pc, acc[c]["P1"], truth[c], ocv[c])
        for m, kind, sc in (("R2", "R2", None), ("R6", "R6", None), ("model", "model", oof["tw"][c]), ("model_rest", "model", oof["tr"][c])):
            acc[c][m] = A.method_units(Pc, kind, score=sc, q=1.0)
            curves[m][c] = A.q_curve(Pc, acc[c][m], truth[c], ocv[c], AD.QS)
    pooled_fixed = {m: A.metrics(list(r.values())) for m, r in res_fixed.items()}
    r0_recall = pooled_fixed["R0"]["recall"]
    slack = {"R2": AD.RULE_SLACK, "R6": AD.RULE_SLACK, "model": AD.MODEL_SLACK, "model_rest": AD.MODEL_SLACK}
    rec_min, rec_min_model = r0_recall - AD.RULE_SLACK, r0_recall - AD.MODEL_SLACK
    qsel = {m: AD.choose_q(curves[m], TRAIN, r0_recall - slack[m]) for m in ("R2", "R6", "model")}
    qsel["model_rest"] = AD.choose_q(curves["model_rest"], TRAIN, rec_min_model, "prec_strict")
    q_found = {m: v is not None for m, v in qsel.items()}
    qsel = {m: v if v is not None else (1.0, None, len(AD.QS) - 1) for m, v in qsel.items()}
    at_q = {m: {c: curves[m][c][qsel[m][2]] for c in TRAIN} for m in curves}
    table = {m: A.metrics(list(r.values())) for m, r in {**res_fixed, **at_q}.items()}
    per_chrom = {m: {c: A.metrics([r[c]]) for c in TRAIN} for m, r in {**res_fixed, **at_q}.items()}
    elig = {m: table[m]["recall"] >= rec_min - 1e-12 for m in AD.RULES}
    best_rule = max((m for m in AD.RULES if elig[m]), key=lambda m: table[m]["prec_known"])
    curve_tab = {m: {str(q): A.metrics([curves[m][c][qi] for c in TRAIN]) for qi, q in enumerate(AD.QS) if round(q * 100) % 5 == 0}
                 for m in curves}

    lengths = {c: A.chrom_len(c) for c in TRAIN}
    bid = A.block_ids(TRAIN, lengths)

    def k_of(m, c, q):
        n = len(acc[c][m])
        return int(math.ceil(q * n - 1e-9)) if q < 1 else n
    Tb = {m: A.block_table(res_fixed[m], truth, {c: P[c]["s"][acc[c][m]] for c in TRAIN}, bid) for m in ("R0", "P1") + tuple(x for x in AD.FIXED if x != "R0")}
    for m in ("R2", "R6", "model"):
        Tb[m] = A.block_table(at_q[m], truth, {c: P[c]["s"][acc[c][m][:k_of(m, c, qsel[m][0])]] for c in TRAIN}, bid)
    Wb = A.boot_weights(len(bid))
    preview = {"model-R0": A.paired(Wb, Tb["model"], Tb["R0"]), f"model-{best_rule}": A.paired(Wb, Tb["model"], Tb[best_rule]),
               "P1-R0": A.paired(Wb, Tb["P1"], Tb["R0"]), "n_blocks": len(bid)}

    d_prec = table["model"]["prec_known"] - table["R0"]["prec_known"]
    d0 = {"rule": "pooled LOCO: model - R0-DPI precision (known) >= +0.10, a q found under the recall bound, model recall >= R0-DPI training recall",
          "model_prec_known": table["model"]["prec_known"], "R0_prec_known": table["R0"]["prec_known"], "diff_prec_known": d_prec,
          "model_recall": table["model"]["recall"], "R0_recall": r0_recall, "q_found": q_found["model"],
          "pass": bool(q_found["model"] and d_prec >= 0.10 and table["model"]["recall"] >= r0_recall - 1e-12)}

    lab = np.concatenate([P[c]["label"] for c in TRAIN])
    s_tw = np.concatenate([oof["tw"][c] for c in TRAIN])
    cpm = np.concatenate([P[c]["X"][:, A.FI["log_cage"]] for c in TRAIN])
    tw = lab != 1
    diag = {"n_pairs": int(len(lab)), "labels": {"true": int((lab == 2).sum()), "uncertain": int((lab == 1).sum()), "wrong": int((lab == 0).sum())},
            "auc_true_vs_wrong": {"model_oof": A.auc(s_tw[tw], lab[tw] == 2), "cage_cpm_only": A.auc(cpm[tw], lab[tw] == 2)},
            "calibration_oof_true_vs_wrong": A.calibration(A.sigmoid(s_tw[tw]), (lab[tw] == 2).astype(float))}

    fin = {k: m for h, k, m in [fit_one((None, "tw")), fit_one((None, "tr"))]}
    coef = {k: {"intercept": fin[k]["w"][0], **{f: fin[k]["w"][i + 1] for i, f in enumerate(A.FEATURES)}} for k in fin}
    model = {"plan": "docs/auth/v1/PLAN_V1.md section 5.3", "cell": cell, "features": list(A.FEATURES),
             "normalised_features": {"log_cage": "log1p(CPM)", "n_samples": "fraction of samples", "log_anti": "log1p(CPM)"},
             "normalisation_constants": {"cage_tags_per_sample": lib, "n_samples": len(lib)}, "logit_C": A.LOGIT_C,
             "true_vs_wrong": fin["tw"], "true_vs_rest": fin["tr"],
             "q": {"model": qsel["model"][0], "model_rest": qsel["model_rest"][0], "R2": qsel["R2"][0], "R6": qsel["R6"][0]},
             "q_found_under_recall_bound": q_found, "best_rule": best_rule, "R0_train_recall": r0_recall,
             "recall_min_model": rec_min_model, "recall_min_rules": rec_min, "q_constraint": Q_CONSTRAINT,
             "train": list(TRAIN), "created": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    dev = {"cell": cell, "R0_train_recall": r0_recall, "q_selected": {m: v[0] for m, v in qsel.items()}, "q_found_under_recall_bound": q_found,
           "best_rule": best_rule, "rule_eligible": elig, "D0": d0, "loco_pooled": table, "loco_block_bootstrap_preview": preview,
           "loco_per_chrom": per_chrom, "q_curves_every_0.05": curve_tab, "pairs": diag, "coef_standardised": coef,
           "n_truth_train": sum(len(v) for v in truth.values())}
    return model, dev


def checks():
    """Module doc item 1 (WTC11-DPI caches written by v1_check.py)."""
    cv = V.check_values()["wtc11_dpi"]
    ev, lib = V.evidence_dicts("wtc11", TRAIN)
    V.use("wtc11")
    same = {c: json.loads(json.dumps(ev[c])) == json.load(open(f"{A.OUT}/evidence/{c}.json")) for c in TRAIN}
    pc = A.v41_pc(TRAIN)
    res = []
    for c in TRAIN:
        P = V.load_pairs("wtc11", c, pc, lib)
        res.append(A.score_units(P, A.method_units(P, "R0"), A.load_truth(c)["ts"], A.obs_cover(A.observed([c])[c])))
    r0 = A.metrics(res)
    out = {"evidence_equal_check_cache": all(same.values()), "evidence_equal_per_chrom": same,
           "cage_tags_equal_check": lib == cv["cage_lib_tags"],
           "R0_equal_check": all(r0[k] == cv["R0_train"][k] for k in ("recall", "prec_known", "n_pred", "tp", "unknown")),
           "truth_is_auth_truth": all(os.path.realpath(f"{A.OUT}/truth/{c}.json") == f"{A.ROOT}/outputs/AUTH/truth/{c}.json" for c in TRAIN)}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=4)
    a = ap.parse_args()
    t0 = time.time()
    os.makedirs(DEVD, exist_ok=True)
    log = lambda *x: print(f"[{time.time() - t0:7.0f}s]", *x, flush=True)   # noqa: E731
    chk = checks()
    log("checks", json.dumps({k: v for k, v in chk.items() if k != "evidence_equal_per_chrom"}))
    assert all(v for k, v in chk.items() if isinstance(v, bool)), chk
    model, dev = develop("wtc11", a.procs, log)
    dev = {"plan": "docs/auth/v1/PLAN_V1.md section 7 step 3", "checks": chk, "q_constraint": Q_CONSTRAINT, **dev,
           "seconds": round(time.time() - t0, 1)}
    json.dump(model, open(f"{DEVD}/MODEL_V1.json", "w"), indent=1)
    json.dump(dev, open(f"{DEVD}/DEV_V1.json", "w"), indent=1, default=float)
    log("q", dev["q_selected"], "best rule", dev["best_rule"], "D0", json.dumps(dev["D0"]))
    for m in ("R0", "R1", "R2", "R3", "R4", "R5", "R6", "P1", "model", "model_rest"):
        t = dev["loco_pooled"][m]
        log(f"{m:10s} recall {t['recall']:.3f}  prec_known {t['prec_known']:.3f}  strict {t['prec_strict']:.3f}  n {t['n_pred']}")


if __name__ == "__main__":
    main()
