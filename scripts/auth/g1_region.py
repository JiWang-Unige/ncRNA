#!/usr/bin/env python3
"""Post-hoc region-level ablation (ChatGPT round 2, G1 + G2), WTC11 LRGASP evidence, AUTH candidates / truth / scoring.
Each variant = 'true vs wrong' L2 logistic on a feature subset, fitted on all candidates or on legal candidates only
(no same-strand v41 coding overlap); '_legalsel' = the end is chosen among legal candidates only; q chosen exactly as AUTH (pooled LOCO recall >= R0 training recall, best precision
known), final fit on the 14 training chromosomes, then applied to the AUTH test chromosomes (already read once by AUTH:
post-hoc).  'full_all' must reproduce the frozen AUTH model and its test numbers.
  g1_region.py -> outputs/AUTH/posthoc/G1_REGION.json
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import auth_lib as A                       # noqa: E402
import auth_dev as AD                      # noqa: E402

TRAIN, TEST = A.TRAIN, A.TEST
VARIANTS = {"full": tuple(A.FEATURES), "geometry": ("log_len", "end_rank"),
            "signal_density": ("log_quant_signal", "quant_nrep", "log_len", "log_quant_inside", "log_cage_inside"),
            "cage_only": ("log_cage",), "quant_only": ("log_quant_signal",)}


def sub(X, cols):
    return X[:, [A.FI[c] for c in cols]]


def fit(P, chroms, cols, legal):
    X, y = [], []
    for c in chroms:
        lab = P[c]["label"]
        m = lab != 1
        if legal:
            m &= P[c]["X"][:, A.FI["pc_same"]] == 0
        X.append(sub(P[c]["X"][m], cols))
        y.append(lab[m] == 2)
    return A.fit(np.vstack(X), np.concatenate(y))


def score(M, P, cols, legal_sel=False):
    s = A.logit(M, sub(P["X"], cols))
    if legal_sel:            # choose the end among legal (non-coding-overlap) candidates only
        s = np.where(P["X"][:, A.FI["pc_same"]] == 0, s, -1e9)
    return s


def main():
    A.open_test()
    pc = A.v41_pc(TRAIN + TEST)
    P = {c: A.build_pairs(c, pc[c]) for c in TRAIN + TEST}
    truth = {c: A.load_truth(c)["ts"] for c in TRAIN + TEST}
    obs = A.observed(TRAIN + TEST)
    ocv = {c: A.obs_cover(obs[c]) for c in TRAIN + TEST}
    r0_train = A.metrics([A.score_units(P[c], A.method_units(P[c], "R0"), truth[c], ocv[c]) for c in TRAIN])["recall"]
    bid = A.block_ids(TEST, {c: A.chrom_len(c) for c in TEST})
    W = A.boot_weights(len(bid))

    def table(acc_by):
        res = {c: A.score_units(P[c], acc_by[c], truth[c], ocv[c]) for c in TEST}
        T = A.block_table(res, truth, {c: P[c]["s"][acc_by[c]] for c in TEST}, bid)
        return A.metrics(list(res.values())), T
    base, T0 = table({c: A.method_units(P[c], "R0") for c in TEST})
    out = {"R0_test": base, "R0_train_recall": r0_train, "variants": {}}
    for name, cols in VARIANTS.items():
        for legal, lsel in ((False, False), (True, False), (False, True), (True, True)):
            if lsel and name not in ("full", "signal_density"):
                continue
            key = f"{name}_{'legal' if legal else 'all'}" + ("_legalsel" if lsel else "")
            curves = {}
            for h in TRAIN:
                M = fit(P, [c for c in TRAIN if c != h], cols, legal)
                acc = A.method_units(P[h], "model", score=score(M, P[h], cols, lsel), q=1.0)
                curves[h] = A.q_curve(P[h], acc, truth[h], ocv[h], AD.QS)
            q = AD.choose_q(curves, TRAIN, r0_train)
            qv = q[0] if q else 1.0
            loco = A.metrics([curves[c][q[2] if q else len(AD.QS) - 1] for c in TRAIN])
            Mf = fit(P, TRAIN, cols, legal)
            met, T = table({c: A.method_units(P[c], "model", score=score(Mf, P[c], cols, lsel), q=qv) for c in TEST})
            out["variants"][key] = {"cols": list(cols), "q": qv, "q_found": bool(q), "loco": loco, "test": met,
                                    "test_minus_R0": A.paired(W, T, T0)}
            print(key, "q", qv, "test recall %.3f prec %.3f n %d" % (met["recall"], met["prec_known"], met["n_pred"]),
                  "dPrec", round(out["variants"][key]["test_minus_R0"]["prec_known"]["diff"], 3), flush=True)
    os.makedirs(f"{A.OUT}/posthoc", exist_ok=True)
    json.dump(out, open(f"{A.OUT}/posthoc/G1_REGION.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
