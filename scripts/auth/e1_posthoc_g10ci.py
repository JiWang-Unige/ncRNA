#!/usr/bin/env python3
"""G10 confidence intervals (post hoc): signal-matched hard pairs on the test chromosomes, original and high-quality
labels; pairwise accuracy differences vs G8 with a bootstrap over TSS groups (2,000, seed 0).
  e1_posthoc_g10ci.py -> outputs/AUTH_E1/POSTHOC_G10_CI.json"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e1_posthoc as P   # noqa: E402

A = P.A


def pairs(G, min_cov=30):
    rows = []
    for gi, g in enumerate(G):
        if g.get("cov", 1e9) < min_cov:
            continue
        sig = np.expm1(g["X"][:, P.SIG])
        sc = P.scores(g)
        neg = np.nonzero(g["cnt"] == 0)[0]
        for i in np.nonzero(g["y"])[0]:
            for j in neg:
                a, b = max(sig[i], 1e-9), max(sig[j], 1e-9)
                if max(a, b) / min(a, b) <= 2 and abs(int(g["ends"][i]) - int(g["ends"][j])) > 200:
                    rows.append([gi] + [float(sc[m][i] > sc[m][j]) for m in ("G8", "E2", "E3", "E4")])
    return np.array(rows)


def ci(R):
    gids = np.unique(R[:, 0])
    rng = np.random.default_rng(0)
    by = {g: R[R[:, 0] == g, 1:] for g in gids}
    out = {"n_pairs": len(R), "n_groups": len(gids), "acc": dict(zip(("G8", "E2", "E3", "E4"), R[:, 1:].mean(0).round(4).tolist()))}
    diffs = {m: [] for m in ("E2", "E3", "E4")}
    for _ in range(2000):
        s = np.vstack([by[g] for g in rng.choice(gids, len(gids))])
        mu = s.mean(0)
        for k, m in ((1, "E2"), (2, "E3"), (3, "E4")):
            diffs[m].append(mu[k] - mu[0])
    out["minus_G8"] = {m: {"diff": float(R[:, 1 + k].mean() - R[:, 1].mean()), "ci95": np.percentile(v, [2.5, 97.5]).round(4).tolist()}
                       for k, (m, v) in enumerate(diffs.items(), 1)}
    return out


res = {}
Go = [g for c in A.TEST for g in P.orig_groups(c, True)]
res["test_orig"] = ci(pairs(Go))
Gh = []
for c in A.TEST:
    Gh += P.hq_groups(c, P.hq_reads(c))[0]
res["test_hq"] = ci(pairs(Gh))
json.dump(res, open(f"{P.EL.OUT}/POSTHOC_G10_CI.json", "w"), indent=1)
print(json.dumps(res, indent=1))
