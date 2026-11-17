#!/usr/bin/env python3
"""M1 weeks 2-4, development on chr18 only (docs/m1/PLAN_M1_W24.md sections 3-4): evidence coverage of the WTC11 truth,
DNA-only v2.0 calls, the snapping baseline (DNA-only calls with ends snapped to the nearest evidence) and the
evidence-only pairing baseline (each ssCAGE TSS paired with the nearest downstream same-strand QuantSeq site; Tiberius
coding genes on the same strand removed).  Grids are explored here and chosen before the freeze.
  dev_baselines.py [--chrom chr18] -> outputs/M1/W24/dev/BASELINES_<chrom>.json
"""
import argparse
import json
import os

import numpy as np

import w24_lib as L

FRACS = (0.1, 0.3, 1.0)
SEEDS = (0, 1, 2)
DMAX = (2_000, 5_000, 10_000, 20_000, 50_000, 100_000)
SNAP_W = (500, 1_000, 2_000, 5_000, 10_000)


def overlaps(iv, s, e):
    return any(a <= e and s <= b for a, b in iv)


def tib_hit(tib, st, s, e):
    a = tib[st]
    return bool(len(a)) and bool(((a[:, 0] <= e) & (s <= a[:, 1])).any())


def pair(ev, tib, dmax):
    out = []
    for st in "+-":
        T, P = ev["tss"][st], ev["pas"][st]
        if not len(T) or not len(P):
            continue
        acc = []
        for i in np.argsort(-T[:, 1], kind="stable"):
            t = int(T[i, 0])
            if st == "+":
                j = np.searchsorted(P, t, side="right")
                if j >= len(P) or P[j] - t > dmax:
                    continue
                s, e = t, int(P[j])
            else:
                j = np.searchsorted(P, t, side="left") - 1
                if j < 0 or t - P[j] > dmax:
                    continue
                s, e = int(P[j]), t
            if e - s + 1 < 200 or overlaps(acc, s, e) or tib_hit(tib, st, s, e):
                continue
            acc.append((s, e))
            out.append({"strand": st, "exons": [(s, e)], "tss_count": int(T[i, 1])})
    return out


def snap(preds, ev, w):
    out = []
    for p in preds:
        st, s, e = p["strand"], p["exons"][0][0], p["exons"][-1][1]
        T = ev["tss"][st][:, 0] if len(ev["tss"][st]) else np.zeros(0, np.int64)
        P = ev["pas"][st]
        a, b = (s, e) if st == "+" else (e, s)          # a = 5' end, b = 3' end
        if len(T):
            k = int(np.abs(T - a).argmin())
            if abs(int(T[k]) - a) <= w:
                a = int(T[k])
        if len(P):
            k = int(np.abs(P - b).argmin())
            if abs(int(P[k]) - b) <= w:
                b = int(P[k])
        s2, e2 = (a, b) if st == "+" else (b, a)
        if e2 - s2 + 1 >= 200:
            out.append(dict(p, exons=[(s2, e2)]))
    return out


def mean_runs(rs):
    keys = ("S50_sn", "S50_pr", "L1_sn", "L1_pr", "n_pred")
    return {k: float(np.mean([r[k] for r in rs if r[k] is not None])) for k in keys}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", default="chr18")
    a = ap.parse_args()
    assert a.chrom == "chr18", "development is chr18 only"
    import simple_det as D
    truth = L.load_truth(a.chrom)
    tib = D.read_tiberius(a.chrom)
    ours = D.load_ours(a.chrom)
    model = [m for m in ours if m != "v1.9"][0]
    tau = ours[model]["tau"]
    dna = {r: [p for p in v if p["conf"] >= tau] for r, v in ours[model]["runs"].items()}
    res = {"chrom": a.chrom, "model": model, "tau": tau, "n_truth_lnc": len(truth["ts"]), "coverage": {}, "dna_only": {},
           "snap": {}, "pair": {}}
    res["dna_only"]["runs"] = {r: {k: v for k, v in L.score(truth, p).items() if k != "S50_ids"} for r, p in dna.items()}
    res["dna_only"]["mean"] = mean_runs(list(res["dna_only"]["runs"].values()))
    for f in FRACS:
        for sd in (SEEDS if f < 1 else (0,)):
            ev = L.evidence(a.chrom, f, sd)
            tag = f"f{f}_s{sd}"
            res["coverage"][tag] = {f"tol{t}": L.coverage(truth, ev, t, t) for t in (100, 500)}
            res["coverage"][tag]["n_tss"] = int(sum(len(v) for v in ev["tss"].values()))
            res["coverage"][tag]["n_pas"] = int(sum(len(v) for v in ev["pas"].values()))
            for w in SNAP_W:
                rs = [L.score(truth, snap(p, ev, w)) for p in dna.values()]
                res["snap"].setdefault(f"w{w}", {})[tag] = mean_runs(rs)
            for dm in DMAX:
                r = L.score(truth, pair(ev, tib, dm))
                res["pair"].setdefault(f"d{dm}", {})[tag] = {k: v for k, v in r.items() if k != "S50_ids"}
    os.makedirs(f"{L.W}/dev", exist_ok=True)
    json.dump(res, open(f"{L.W}/dev/BASELINES_{a.chrom}.json", "w"), indent=1)
    print(json.dumps({"dna_only": res["dna_only"]["mean"], "coverage_f1": res["coverage"]["f1.0_s0"]}, indent=1))
    for k in ("snap", "pair"):
        for g, v in res[k].items():
            x = v["f1.0_s0"]
            print(k, g, "f1.0", {kk: round(vv, 3) if isinstance(vv, float) else vv for kk, vv in x.items()})


if __name__ == "__main__":
    main()
