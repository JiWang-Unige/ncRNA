#!/usr/bin/env python3
"""AUTH step 3 (docs/auth/PLAN_AUTH.md sections 5 and 7): development on the 14 training chromosomes only.

1. Inputs for the training chromosomes: soft-masked sequence (hg38.fa.gz), evidence, rep2 / rep3 observed spans, v41
   protein_coding genes, candidate pairs with labels and features (auth_lib.build_pairs).
2. Consistency checks against the frozen W24 code: chr18 evidence and observed spans equal the W24 caches; the
   sequence equals data/human/chr18.fa; the pairs equal e3_anchor.candidates(ev, 50 kb, 20); R0 equals
   e3_anchor.select_evidence(candidates(ev, 50 kb, 3), 'farthest', v41) on every training chromosome; R0 with the
   Tiberius mask on chr18 reproduces the W24 development number (dev/E3_chr18_f1.0_s0.json).
3. Leave-one-chromosome-out (LOCO): the L2 logistic model ('true vs wrong'; 'true vs rest' reported only) is fitted on
   13 chromosomes and scores the 14th; q (top fraction of the accepted units, per chromosome) is chosen on the pooled
   out-of-fold units: highest precision (known) with recall >= R0's training recall (ties: larger q) -- plan section
   5.3 as revised by the author on 2026-11-24 before the freeze and before any test number was read (it was R0 - 0.02;
   the rules R2 / R6 and the best-rule eligibility keep R0 - 0.02).  The
   rules R1-R6 are evaluated on the same chromosomes (R2 / R6 choose q the same way); the best rule is the eligible
   one (same recall bound) with the highest precision (known).
4. Final models fitted on all 14 chromosomes; q and the best rule unchanged -> outputs/AUTH/dev/MODEL.json.
  auth_dev.py [--procs 6] -> outputs/AUTH/dev/{DEV.json, MODEL.json}
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
import auth_lib as A                       # noqa: E402
from auth_truth import TRAIN               # noqa: E402

QS = tuple(round(0.05 + 0.01 * i, 2) for i in range(96))
RULES = ("R1", "R2", "R3", "R4", "R5", "R6")
FIXED = ("R0", "R1", "R3", "R4", "R5")
MODEL_SLACK, RULE_SLACK = 0.0, 0.02     # recall bound = R0 training recall - slack (plan 5.3, revised 2026-11-24)
Q_CONSTRAINT = ("model / model_rest: LOCO recall >= R0 training recall; R2 / R6 and best-rule eligibility: >= R0 - 0.02 "
                "(PLAN_AUTH.md 5.3 revision, author 2026-11-24, before the freeze and before any test number was read)")
DEV = f"{A.OUT}/dev"
PC = {}


def prep(c):
    t0 = time.time()
    A.evidence_raw(c)
    P = A.build_pairs(c, PC[c])
    return c, int(len(P["s"])), round(time.time() - t0, 1)


def fit_one(args):
    h, kind = args
    X, y = [], []
    for c in TRAIN:
        if c == h:
            continue
        z = np.load(f"{A.OUT}/pairs/{c}.npz")
        lab = z["label"]
        sel = lab != 1 if kind == "tw" else np.ones(len(lab), bool)
        X.append(z["X"][sel])
        y.append(lab[sel] == 2)
    return h, kind, A.fit(np.vstack(X), np.concatenate(y))


def choose_q(curves, chroms, rec_min, crit="prec_known"):
    best = None
    for qi, q in enumerate(QS):
        m = A.metrics([curves[c][qi] for c in chroms])
        if m["recall"] is None or m["recall"] < rec_min - 1e-12 or m[crit] is None:
            continue
        if best is None or m[crit] > best[1] + 1e-12 or (abs(m[crit] - best[1]) <= 1e-12 and q > best[0]):
            best = (q, m[crit], qi)
    return best


def checks(P, ev_c, pc41):
    """Consistency with the frozen W24 code (module doc, item 2)."""
    from e3_anchor import candidates, select_evidence
    import simple_det as D
    out = {}
    w = json.load(open(f"{A.L.W}/evidence/chr18.json"))
    mine = json.load(open(f"{A.OUT}/evidence/chr18.json"))
    out["chr18_cage_equal_w24"] = [{k: v for k, v in x.items() if k != "mode_n"} for x in mine["cage"]] == w["cage"]
    out["chr18_quant_equal_w24"] = mine["quant"] == w["quant"]
    ow = sorted(tuple(x) for x in json.load(open(f"{A.L.W}/evidence/chr18.observed_rep23.json")))
    out["chr18_observed_equal_w24"] = ow == sorted(A.observed(["chr18"])["chr18"])
    import s3_features as S3
    out["chr18_seq_equal_data_human"] = S3.read_fasta(f"{A.ROOT}/data/human/chr18.fa").encode() == A.load_seq("chr18")
    ev18 = A.L.evidence("chr18", 1.0, 0)
    out["chr18_evidence_equal_w24_lib"] = all(np.array_equal(ev18[k][s], ev_c["chr18"][k][s]) for k in ("tss", "pas") for s in "+-")
    pairs_ok, r0_ok = True, True
    for c in TRAIN:
        ev = ev_c[c]
        Pc = P[c]
        cw = candidates(ev, A.DMAX, A.K_WIDE)
        g0, g1 = A.groups(Pc)
        mine_sp = [[(int(Pc["s"][i]), int(Pc["e"][i])) for i in range(a, b)] for a, b in zip(g0, g1)]
        pairs_ok &= mine_sp == [sp for _, _, sp in cw]
        sig = A.L.quant_signal(c)
        ref = select_evidence(candidates(ev, A.DMAX, A.K_R0), sig, "farthest", pc41[c])
        acc = A.method_units(Pc, "R0")
        r0_ok &= [(p["strand"], p["exons"][0][0], p["exons"][-1][1]) for p in ref] == \
                 [(A.STR[int(Pc["strand"][i])], int(Pc["s"][i]), int(Pc["e"][i])) for i in acc]
    out["pairs_equal_e3_candidates_k20"] = bool(pairs_ok)
    out["R0_equal_select_evidence_v41"] = bool(r0_ok)
    tib = D.read_tiberius("chr18")
    acc = A.method_units(P["chr18"], "R0", masked=A.masks(P["chr18"], tib))
    preds = [{"strand": A.STR[int(P["chr18"]["strand"][i])], "exons": [(int(P["chr18"]["s"][i]), int(P["chr18"]["e"][i]))]} for i in acc]
    r = A.L.score_unknown(A.L.load_truth("chr18"), preds, "chr18")
    e3 = json.load(open(f"{A.L.W}/dev/E3_chr18_f1.0_s0.json"))["grid"]["d50000_K3_ev-farthest+tib"]["evidence_only"]
    out["chr18_R0_tiberius_vs_w24_dev"] = {"auth": {k: r[k] for k in ("S50_sn", "S50_pr_known", "n_pred")},
                                           "w24": {k: e3[k] for k in ("S50_sn", "S50_pr_known", "n_pred")}}
    out["chr18_R0_tiberius_equal_w24_dev"] = all(abs(r[k] - e3[k]) < 1e-9 for k in ("S50_sn", "S50_pr_known", "n_pred"))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=6)
    a = ap.parse_args()
    t0 = time.time()
    os.makedirs(DEV, exist_ok=True)
    log = lambda *x: print(f"[{time.time() - t0:7.0f}s]", *x, flush=True)   # noqa: E731

    # ---------------------------------------------------------------- 1. inputs
    A.extract_seq(TRAIN)
    log("sequence ready")
    obs = A.observed(TRAIN)
    log("observed spans", {c: len(v) for c, v in obs.items()})
    PC.update(A.v41_pc(TRAIN))
    with Pool(a.procs) as pool:
        for c, n, sec in pool.imap_unordered(prep, TRAIN):
            log("pairs", c, n, f"{sec}s")
    P = {c: A.build_pairs(c, PC[c]) for c in TRAIN}
    ev_c = {c: A.evidence(c)[0] for c in TRAIN}
    truth = {c: A.load_truth(c)["ts"] for c in TRAIN}
    ocv = {c: A.obs_cover(obs[c]) for c in TRAIN}

    # ---------------------------------------------------------------- 2. checks
    chk = checks(P, ev_c, PC)
    log("checks", json.dumps(chk))
    assert all(v for k, v in chk.items() if isinstance(v, bool)), chk

    # ---------------------------------------------------------------- 3. LOCO
    with Pool(min(a.procs, 4)) as pool:
        fits = pool.map(fit_one, [(h, k) for h in TRAIN for k in ("tw", "tr")])
    oof = {k: {h: A.logit(m, P[h]["X"]) for h, kk, m in fits if kk == k} for k in ("tw", "tr")}
    log("LOCO fits done")

    acc = {c: {} for c in TRAIN}
    curves = {m: {} for m in ("R2", "R6", "model", "model_rest")}
    res_fixed = {m: {} for m in FIXED + ("P1",)}
    for c in TRAIN:
        Pc = P[c]
        for m in FIXED:
            acc[c][m] = A.method_units(Pc, m)
            res_fixed[m][c] = A.score_units(Pc, acc[c][m], truth[c], ocv[c])
        acc[c]["P1"] = A.method_units(Pc, "P1", score=oof["tw"][c])
        res_fixed["P1"][c] = A.score_units(Pc, acc[c]["P1"], truth[c], ocv[c])
        for m, kind, sc in (("R2", "R2", None), ("R6", "R6", None), ("model", "model", oof["tw"][c]), ("model_rest", "model", oof["tr"][c])):
            acc[c][m] = A.method_units(Pc, kind, score=sc, q=1.0)
            curves[m][c] = A.q_curve(Pc, acc[c][m], truth[c], ocv[c], QS)
        log("LOCO scored", c)
    pooled_fixed = {m: A.metrics(list(r.values())) for m, r in res_fixed.items()}
    r0_recall = pooled_fixed["R0"]["recall"]
    slack = {"R2": RULE_SLACK, "R6": RULE_SLACK, "model": MODEL_SLACK, "model_rest": MODEL_SLACK}
    rec_min, rec_min_model = r0_recall - RULE_SLACK, r0_recall - MODEL_SLACK
    qsel = {m: choose_q(curves[m], TRAIN, r0_recall - slack[m]) for m in ("R2", "R6", "model")}
    qsel["model_rest"] = choose_q(curves["model_rest"], TRAIN, rec_min_model, "prec_strict")
    q_found = {m: v is not None for m, v in qsel.items()}
    qsel = {m: v if v is not None else (1.0, None, len(QS) - 1) for m, v in qsel.items()}   # no eligible q: report q = 1
    nested = {m: {h: (choose_q(curves[m], [c for c in TRAIN if c != h], A.metrics([res_fixed["R0"][c] for c in TRAIN if c != h])["recall"] - slack[m],
                               "prec_strict" if m == "model_rest" else "prec_known") or (None,))[0] for h in TRAIN}
              for m in curves}
    at_q = {m: {c: curves[m][c][qsel[m][2]] for c in TRAIN} for m in curves}
    table = {m: A.metrics(list(r.values())) for m, r in {**res_fixed, **at_q}.items()}
    per_chrom = {m: {c: A.metrics([r[c]]) for c in TRAIN} for m, r in {**res_fixed, **at_q}.items()}
    elig = {m: table[m]["recall"] >= rec_min - 1e-12 for m in RULES}
    best_rule = max((m for m in RULES if elig[m]), key=lambda m: table[m]["prec_known"])
    curve_tab = {m: {str(q): A.metrics([curves[m][c][qi] for c in TRAIN]) for qi, q in enumerate(QS) if round(q * 100) % 5 == 0}
                 for m in curves}

    # preview on the training LOCO units (information only; q was chosen on the same units): paired 1 Mb block bootstrap
    lengths = {c: A.chrom_len(c) for c in TRAIN}
    bid = A.block_ids(TRAIN, lengths)

    def k_of(m, c, q):
        n = len(acc[c][m])
        return int(math.ceil(q * n - 1e-9)) if q < 1 else n
    Tb = {m: A.block_table(res_fixed[m], truth, {c: P[c]["s"][acc[c][m]] for c in TRAIN}, bid) for m in ("R0", "P1") + tuple(x for x in FIXED if x != "R0")}
    for m in ("R2", "R6", "model"):
        Tb[m] = A.block_table(at_q[m], truth, {c: P[c]["s"][acc[c][m][:k_of(m, c, qsel[m][0])]] for c in TRAIN}, bid)
    Wb = A.boot_weights(len(bid))
    rd = A.boot_metrics(Wb, Tb["model"])["recall"] - A.boot_metrics(Wb, Tb["R0"])["recall"]
    preview = {"model-R0": A.paired(Wb, Tb["model"], Tb["R0"]), f"model-{best_rule}": A.paired(Wb, Tb["model"], Tb[best_rule]),
               "P1-R0": A.paired(Wb, Tb["P1"], Tb["R0"]), "frac_resamples_model_minus_R0_recall_below_-0.03": float(np.mean(rd < -0.03)),
               "n_blocks": len(bid)}
    alt = {}
    for bound in (0.02, 0.01, 0.0):
        rm = r0_recall - bound
        alt[f"recall>=R0-{bound}"] = {m: (lambda b: {"q": b[0], **A.metrics([curves[m][c][b[2]] for c in TRAIN])} if b else None)(choose_q(curves[m], TRAIN, rm))
                                      for m in ("model", "R2")}
    preview["recall_bounds_compared"] = alt

    # pair-level diagnostics on the out-of-fold scores
    lab = np.concatenate([P[c]["label"] for c in TRAIN])
    s_tw = np.concatenate([oof["tw"][c] for c in TRAIN])
    s_tr = np.concatenate([oof["tr"][c] for c in TRAIN])
    tags = np.concatenate([P[c]["X"][:, A.FI["log_cage"]] for c in TRAIN])
    tw = lab != 1
    diag = {"n_pairs": int(len(lab)), "labels": {"true": int((lab == 2).sum()), "uncertain": int((lab == 1).sum()), "wrong": int((lab == 0).sum())},
            "auc_true_vs_wrong": {"model_oof": A.auc(s_tw[tw], lab[tw] == 2), "ssCAGE_tags_only": A.auc(tags[tw], lab[tw] == 2)},
            "auc_true_vs_rest": {"model_rest_oof": A.auc(s_tr, lab == 2), "ssCAGE_tags_only": A.auc(tags, lab == 2)},
            "calibration_oof_true_vs_wrong": A.calibration(A.sigmoid(s_tw[tw]), (lab[tw] == 2).astype(float))}
    labels_units = {m: {"true": table[m]["tp"], "uncertain": table[m]["unknown"], "wrong": table[m]["wrong"]} for m in table}

    # ---------------------------------------------------------------- 4. final models
    fin = {k: m for h, k, m in [fit_one((None, "tw")), fit_one((None, "tr"))]}
    coef = {k: {"intercept": fin[k]["w"][0], **{f: fin[k]["w"][i + 1] for i, f in enumerate(A.FEATURES)}} for k in fin}
    model = {"plan": "docs/auth/PLAN_AUTH.md section 5.3", "features": list(A.FEATURES), "logit_C": A.LOGIT_C,
             "true_vs_wrong": fin["tw"], "true_vs_rest": fin["tr"],
             "q": {"model": qsel["model"][0], "model_rest": qsel["model_rest"][0], "R2": qsel["R2"][0], "R6": qsel["R6"][0]},
             "q_found_under_recall_bound": q_found, "best_rule": best_rule, "R0_train_recall": r0_recall,
             "recall_min_model": rec_min_model, "recall_min_rules": rec_min, "q_constraint": Q_CONSTRAINT,
             "constants": {"DMAX": A.DMAX, "K_WIDE": A.K_WIDE, "K_R0": A.K_R0, "A_MAX": A.A_MAX, "PAS_HEX": A.PAS_HEX,
                           "WIN_ANTI": A.WIN_ANTI, "BIDIR": A.BIDIR, "WIN_GC": A.WIN_GC, "q_application": "top fraction per chromosome"},
             "train": list(TRAIN), "created": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    json.dump(model, open(f"{DEV}/MODEL.json", "w"), indent=1)
    dev = {"plan": "docs/auth/PLAN_AUTH.md section 7 step 3", "train": list(TRAIN), "checks": chk,
           "R0_train_recall": r0_recall, "recall_min_model": rec_min_model, "recall_min_rules": rec_min, "q_constraint": Q_CONSTRAINT, "q_selected": {m: v[0] for m, v in qsel.items()}, "q_found_under_recall_bound": q_found,
           "q_nested_leave_one_out": nested, "best_rule": best_rule, "rule_eligible": elig,
           "loco_pooled": table, "loco_block_bootstrap_preview": preview, "loco_per_chrom": per_chrom, "q_curves_every_0.05": curve_tab,
           "labels_of_units": labels_units, "pairs": diag, "coef_standardised": coef,
           "n_truth_train": sum(len(v) for v in truth.values()), "seconds": round(time.time() - t0, 1)}
    json.dump(dev, open(f"{DEV}/DEV.json", "w"), indent=1, default=float)
    log("q", dev["q_selected"], "best rule", best_rule, "recall_min model", round(rec_min_model, 4), "rules", round(rec_min, 4))
    for m in ("R0", "R1", "R2", "R3", "R4", "R5", "R6", "P1", "model", "model_rest"):
        t = table[m]
        log(f"{m:10s} recall {t['recall']:.3f}  prec_known {t['prec_known']:.3f}  strict {t['prec_strict']:.3f}  n {t['n_pred']}")


if __name__ == "__main__":
    main()
