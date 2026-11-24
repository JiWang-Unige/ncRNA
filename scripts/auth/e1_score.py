#!/usr/bin/env python3
"""AUTH-E1 one-time test (docs/auth/e1/PLAN_E1.md sections 6-7) on chr2, 6, 10, 16, 22 with everything frozen in
outputs/AUTH_E1/dev/MODEL_E1.json (E0 K, E1 / E2 / E4 conditional logits, strongest baseline) and the AUTH model (E3).
Verdict on Top-1 differences E3 - baseline and E4 - baseline (paired 1 Mb block bootstrap, 2,000, seed 0):
  GO    E3 or E4: diff >= +0.10 and CI low > 0
  STOP  E3 and E4: CI high < +0.05
  INCONCLUSIVE  test groups < 60, or neither.
Reported only: MRR; windows 50 / 200 bp; local coverage >= 30; the groups whose supported pair overlaps an AUTH truth lnc.
  e1_score.py --mode freeze | check | test | dev     (dev = the same code on the training chromosomes, in-sample)
"""
import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import auth_lib as A                       # noqa: E402
import e1_labels as EL                     # noqa: E402
import e1_dev as ED                        # noqa: E402

FROZEN, RESULT = f"{EL.OUT}/FROZEN", f"{EL.OUT}/RESULT_E1.json"
REQUIRED = ("scripts/auth/e1_labels.py", "scripts/auth/e1_dev.py", "scripts/auth/e1_score.py", "scripts/auth/auth_lib.py",
            "scripts/auth/auth_truth.py", "scripts/m1/w24/w24_lib.py", "outputs/AUTH/dev/MODEL.json",
            "outputs/AUTH_E1/dev/MODEL_E1.json", "outputs/AUTH_E1/dev/DEV_E1.json", "outputs/AUTH_E1/FEASIBILITY.json")
METHODS = ("E0", "E1", "E2", "E3", "E4")


def groups_of(chroms, w, label_dir, min_cov=EL.MIN_COV):
    G = []
    for c in chroms:
        z = np.load(f"{A.OUT}/pairs/{c}.npz")
        L = dict(np.load(f"{label_dir}/{c}.npz"))
        for g in EL.groups(L, w):
            if L[f"cov{w}"][g].max() < min_cov:
                continue
            G.append({"X": z["X"][g], "rank": z["rank"][g], "y": L[f"link{w}"][g] >= EL.MIN_LINK, "s": int(z["s"][g[0]]),
                      "e": z["e"][g], "ss": z["s"][g], "strand": A.STR[int(z["strand"][g[0]])], "chrom": c})
    return G


def evaluate(G, M, auth, chroms):
    hit = {m: [] for m in METHODS}
    mrr = {m: [] for m in METHODS if m != "E0"}
    for g in G:
        hit["E0"].append(ED.top1(g, ED.pick_e0(g, M["E0_K"])))
        for m in ("E1", "E2", "E3", "E4"):
            sc = A.logit(auth, g["X"]) if m == "E3" else ED.score_clogit(M[m], g)
            hit[m].append(ED.top1(g, int(np.argmax(sc))))
            mrr[m].append(ED.rr(g, sc))
    n = len(G)
    out = {"n_groups": n, "random": float(np.mean([g["y"].mean() for g in G])) if n else None,
           "top1": {m: float(np.mean(v)) if n else None for m, v in hit.items()},
           "mrr": {m: float(np.mean(v)) if n else None for m, v in mrr.items()}}
    if n:
        bid = A.block_ids(chroms, {c: A.chrom_len(c) for c in chroms})
        B = np.zeros((len(bid), 1 + len(METHODS)))
        for i, g in enumerate(G):
            b = bid[(g["chrom"], (g["s"] - 1) // A.BLOCK)]
            B[b, 0] += 1
            for j, m in enumerate(METHODS):
                B[b, j + 1] += hit[m][i]
        S = A.boot_weights(len(bid)) @ B
        acc = {m: S[:, j + 1] / S[:, 0] for j, m in enumerate(METHODS)}
        base = M["strongest_baseline"]
        out["paired_vs_baseline"] = {m: {"diff": out["top1"][m] - out["top1"][base],
                                         "ci95": np.nanpercentile(acc[m] - acc[base], [2.5, 97.5]).tolist()}
                                     for m in METHODS if m != base}
        out["ci95"] = {m: np.nanpercentile(acc[m], [2.5, 97.5]).tolist() for m in METHODS}
    return out


def verdict(ev, M):
    base = M["strongest_baseline"]
    if ev["n_groups"] < 60:
        return {"decision": "INCONCLUSIVE", "reason": "test groups < 60", "baseline": base}
    d = {m: ev["paired_vs_baseline"][m] for m in ("E3", "E4")}
    go = any(x["diff"] >= 0.10 and x["ci95"][0] > 0 for x in d.values())
    stop = all(x["ci95"][1] < 0.05 for x in d.values())
    dec = "GO" if go else "STOP" if stop else "INCONCLUSIVE"
    return {"decision": dec, "baseline": base, "E3_vs_baseline": d["E3"], "E4_vs_baseline": d["E4"]}


def lnc_subset(G):
    out = []
    for g in G:
        ts = A.load_truth(g["chrom"])["ts"]
        sup = np.nonzero(g["y"])[0]
        if any(t[1] == g["strand"] and t[2] <= int(g["e"][j]) and int(g["ss"][j]) <= t[3] for j in sup for t in ts):
            out.append(g)
    return out


def check_frozen():
    assert os.path.exists(FROZEN), "outputs/AUTH_E1/FROZEN missing"
    fz = json.load(open(FROZEN))
    assert all(p in fz for p in REQUIRED), "FROZEN incomplete"
    bad = [p for p, h in fz.items() if A.sha256(f"{A.ROOT}/{p}") != h]
    assert not bad, f"changed since the freeze: {bad}"
    assert not os.path.exists(RESULT), "the E1 test has already been read once"
    return fz


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True, choices=["freeze", "check", "test", "dev"])
    a = ap.parse_args()
    t0 = time.time()
    M = json.load(open(f"{EL.OUT}/dev/MODEL_E1.json"))
    auth = json.load(open(f"{A.OUT}/dev/MODEL.json"))["true_vs_wrong"]
    if a.mode == "freeze":
        assert not os.path.exists(FROZEN)
        fz = {p: A.sha256(f"{A.ROOT}/{p}") for p in REQUIRED}
        json.dump(fz, open(FROZEN, "w"), indent=1)
        print(json.dumps(fz, indent=1))
        return
    if a.mode == "dev":
        chroms, ldir = A.TRAIN, f"{EL.OUT}/labels"
    else:
        fz = check_frozen()
        print("FROZEN ok:", len(fz), flush=True)
        if a.mode == "check":
            return
        with open(f"{EL.OUT}/test_attempts.log", "a") as fh:
            fh.write(time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()) + " --mode test started\n")
        A.open_test()
        chroms, ldir = A.TEST, f"{EL.OUT}/labels_test"
        os.makedirs(ldir, exist_ok=False)
        for c in chroms:
            assert os.path.exists(f"{A.OUT}/pairs/{c}.npz"), c
            np.savez(f"{ldir}/{c}.npz", **EL.label_chrom(c))
    G = groups_of(chroms, 100, ldir)
    main_ev = evaluate(G, M, auth, chroms)
    out = {"plan": "docs/auth/e1/PLAN_E1.md", "mode": a.mode, "chroms": list(chroms), "main": main_ev, "verdict": verdict(main_ev, M),
           "window50": evaluate(groups_of(chroms, 50, ldir), M, auth, chroms),
           "window200": evaluate(groups_of(chroms, 200, ldir), M, auth, chroms),
           "cov30": evaluate(groups_of(chroms, 100, ldir, 30), M, auth, chroms),
           "lnc_subset": evaluate(lnc_subset(G), M, auth, chroms), "seconds": None}
    out["seconds"] = round(time.time() - t0, 1)
    path = RESULT if a.mode == "test" else f"{EL.OUT}/dev/SCORE_DEV_E1.json"
    json.dump(out, open(path, "w"), indent=1, default=float)
    print(json.dumps({"verdict": out["verdict"], "top1": main_ev["top1"], "n": main_ev["n_groups"], "random": main_ev["random"]}, indent=1, default=float))


if __name__ == "__main__":
    main()
