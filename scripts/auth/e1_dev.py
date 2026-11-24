#!/usr/bin/env python3
"""AUTH-E1 development (docs/auth/e1/PLAN_E1.md sections 4-6), 14 training chromosomes only.

Groups: e1_labels.groups (W = 100).  Within a group every kept candidate end is ranked; y = supported (n_link >= 2).
Top-1 = the top-ranked end is supported.  Methods:
  E0  farthest of the K nearest kept ends, K in {1, 2, 3, 5, 10, 20} (chosen on the 13 other chromosomes in LOCO)
  E1  conditional logit on (log_len, end_rank)
  E2  conditional logit on (log_quant_signal, quant_nrep, log_len, log_quant_inside, log_cage_inside)
  E3  frozen AUTH model (outputs/AUTH/dev/MODEL.json, 'true vs wrong' score), argmax, no refit
  E4  conditional logit on the 18 AUTH features
Conditional logit: softmax within the group, loss = -sum_g log(sum_supported p) + ||w||^2 / (2C), C = 1, features
standardised on the fitting rows (L-BFGS).  Strongest baseline = best LOCO Top-1 among E0-E2.
  e1_dev.py -> outputs/AUTH_E1/dev/{DEV_E1.json, MODEL_E1.json}
"""
import json
import math
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import auth_lib as A                       # noqa: E402
import e1_labels as EL                     # noqa: E402

TRAIN = A.TRAIN
DEVD = f"{EL.OUT}/dev"
KS = (1, 2, 3, 5, 10, 20)
FEATS = {"E1": ("log_len", "end_rank"),
         "E2": ("log_quant_signal", "quant_nrep", "log_len", "log_quant_inside", "log_cage_inside"),
         "E4": tuple(A.FEATURES)}
C = 1.0


def load(c):
    z = np.load(f"{A.OUT}/pairs/{c}.npz")
    L = dict(np.load(f"{EL.OUT}/labels/{c}.npz"))
    gs = EL.groups(L, 100)
    return [{"X": z["X"][g], "rank": z["rank"][g], "y": L["link100"][g] >= EL.MIN_LINK, "s": int(z["s"][g[0]]),
             "chrom": c} for g in gs]


def fit_clogit(G, cols):
    ci = [A.FI[f] for f in cols]
    Xall = np.vstack([g["X"][:, ci] for g in G])
    mu, sd = Xall.mean(0), Xall.std(0)
    sd = np.where(sd > 0, sd, 1.0)
    Xs = [(g["X"][:, ci] - mu) / sd for g in G]
    Ys = [g["y"] for g in G]

    def f(w):
        loss, grad = 0.5 * w @ w / C, w / C
        for X, y in zip(Xs, Ys):
            s = X @ w
            m = s.max()
            e = np.exp(s - m)
            p = e / e.sum()
            q = e * y
            q = q / q.sum()
            loss -= math.log(e[y].sum()) - math.log(e.sum())
            grad -= X.T @ (q - p)
        return loss, grad
    r = minimize(f, np.zeros(len(ci)), jac=True, method="L-BFGS-B")
    return {"cols": list(cols), "mu": mu.tolist(), "sd": sd.tolist(), "w": r.x.tolist(), "converged": bool(r.success)}


def score_clogit(M, g):
    ci = [A.FI[f] for f in M["cols"]]
    return ((g["X"][:, ci] - np.array(M["mu"])) / np.array(M["sd"])) @ np.array(M["w"])


def pick_e0(g, K):
    ok = np.nonzero(g["rank"] < K)[0]
    if not len(ok):
        ok = np.arange(len(g["rank"]))
    return ok[np.argmax(g["rank"][ok])]


def top1(g, choice):
    return bool(g["y"][choice])


def rr(g, score):
    order = np.argsort(-score, kind="stable")
    return 1.0 / (1 + int(np.nonzero(g["y"][order])[0][0]))


def main():
    t0 = time.time()
    os.makedirs(DEVD, exist_ok=True)
    G = {c: load(c) for c in TRAIN}
    allg = [g for c in TRAIN for g in G[c]]
    auth = json.load(open(f"{A.OUT}/dev/MODEL.json"))["true_vs_wrong"]
    res = {m: {c: [] for c in TRAIN} for m in ("random", "E0", "E1", "E2", "E3", "E4")}
    mrr = {m: {c: [] for c in TRAIN} for m in ("E1", "E2", "E3", "E4")}
    k_loco = {}
    for h in TRAIN:
        tr = [g for c in TRAIN if c != h for g in G[c]]
        kbest = max(KS, key=lambda K: (np.mean([top1(g, pick_e0(g, K)) for g in tr]), -K))
        k_loco[h] = kbest
        fits = {m: fit_clogit(tr, FEATS[m]) for m in ("E1", "E2", "E4")}
        for g in G[h]:
            res["random"][h].append(float(g["y"].mean()))
            res["E0"][h].append(top1(g, pick_e0(g, kbest)))
            sc3 = A.logit(auth, g["X"])
            res["E3"][h].append(top1(g, int(np.argmax(sc3))))
            mrr["E3"][h].append(rr(g, sc3))
            for m in ("E1", "E2", "E4"):
                sc = score_clogit(fits[m], g)
                res[m][h].append(top1(g, int(np.argmax(sc))))
                mrr[m][h].append(rr(g, sc))
        print(h, len(G[h]), {m: round(float(np.mean(v[h])), 3) if v[h] else None for m, v in res.items()}, flush=True)
    pooled = {m: float(np.mean([x for c in TRAIN for x in v[c]])) for m, v in res.items()}
    pooled_mrr = {m: float(np.mean([x for c in TRAIN for x in v[c]])) for m, v in mrr.items()}
    best = max(("E0", "E1", "E2"), key=lambda m: pooled[m])
    # paired 1 Mb block bootstrap preview (LOCO predictions)
    lengths = {c: A.chrom_len(c) for c in TRAIN}
    bid = A.block_ids(TRAIN, lengths)
    B = np.zeros((len(bid), 6))
    cols = ("E0", "E1", "E2", "E3", "E4")
    for c in TRAIN:
        for i, g in enumerate(G[c]):
            b = bid[(c, (g["s"] - 1) // A.BLOCK)]
            B[b, 0] += 1
            for j, m in enumerate(cols):
                B[b, j + 1] += res[m][c][i]
    W = A.boot_weights(len(bid))
    S = W @ B
    acc = {m: S[:, j + 1] / S[:, 0] for j, m in enumerate(cols)}
    preview = {f"{m}-{best}": {"diff": pooled[m] - pooled[best], "ci95": np.nanpercentile(acc[m] - acc[best], [2.5, 97.5]).tolist()}
               for m in ("E3", "E4")}
    # final models on all 14 chromosomes
    kfin = max(KS, key=lambda K: (np.mean([top1(g, pick_e0(g, K)) for g in allg]), -K))
    fin = {m: fit_clogit(allg, FEATS[m]) for m in ("E1", "E2", "E4")}
    model = {"plan": "docs/auth/e1/PLAN_E1.md", "window": 100, "min_link": EL.MIN_LINK, "min_cov": EL.MIN_COV,
             "E0_K": kfin, "E1": fin["E1"], "E2": fin["E2"], "E4": fin["E4"], "E3": "outputs/AUTH/dev/MODEL.json true_vs_wrong",
             "strongest_baseline": best, "C": C, "train": list(TRAIN), "created": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    json.dump(model, open(f"{DEVD}/MODEL_E1.json", "w"), indent=1)
    dev = {"n_groups_train": len(allg), "n_candidates_mean": float(np.mean([len(g["y"]) for g in allg])),
           "n_supported_mean": float(np.mean([g["y"].sum() for g in allg])), "loco_top1": pooled, "loco_mrr": pooled_mrr,
           "loco_top1_per_chrom": {m: {c: float(np.mean(v[c])) for c in TRAIN} for m, v in res.items()},
           "E0_K_loco": k_loco, "E0_K_final": kfin, "strongest_baseline": best, "preview_block_bootstrap": preview,
           "coef_standardised": {m: dict(zip(fin[m]["cols"], fin[m]["w"])) for m in fin},
           "converged": {m: fin[m]["converged"] for m in fin}, "seconds": round(time.time() - t0, 1)}
    json.dump(dev, open(f"{DEVD}/DEV_E1.json", "w"), indent=1, default=float)
    print(json.dumps({k: dev[k] for k in ("n_groups_train", "loco_top1", "loco_mrr", "strongest_baseline", "E0_K_final", "preview_block_bootstrap")}, indent=1))


if __name__ == "__main__":
    main()
