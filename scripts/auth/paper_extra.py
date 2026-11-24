#!/usr/bin/env python3
"""Pre-writing items from the ChatGPT round-3 checklist (docs/paper_auth/OUTLINE.md section 5), post hoc, no refit
of any frozen model.

1. Paired CIs (1 Mb block bootstrap, 2,000, seed 0, same resamples for both methods):
   a. G11 on the 261 E1 test groups: per group, does the top-ranked end of E3 (frozen AUTH) / E2 match an AUTH truth lnc
      locus at S50, and is it molecule-supported (>= 2 rep1 reads)?  Differences E3 - E2 for both rates.
   b. Region level (AUTH test chromosomes): full 18-feature model vs the 5-feature signal/density model, both trained on
      all candidates with the end chosen among legal candidates (g1_region variants *_all_legalsel); precision (known),
      recall, strict precision differences.
2. Locus cases, selection rule written before looking at any case:
   among E1 test groups where E3's top end matches S50 but is unsupported (< 2 reads) while E2's top end is supported,
   rank by local coverage (reads with 5' end within 100 bp of the TSS), descending; take the first 3 on different
   chromosomes.  Report every kept candidate end of the group (position, QuantSeq signal, linking reads, S50 match,
   rank by E2 / E3).
  paper_extra.py -> outputs/AUTH_E1/PAPER_EXTRA.json
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e1_posthoc as P                     # noqa: E402  (opens the test chromosomes; frozen models loaded there)
import g1_region as GR                     # noqa: E402

A, EL = P.A, P.EL


def s50_flags(g, ts):
    f = np.zeros(len(g["y"]), bool)
    for j in range(len(g["y"])):
        s, e = int(g["ss"][j]), int(g["ee"][j])
        for _, st, gs, ge in ts:
            if st != g["strand"]:
                continue
            ov = min(e, ge) - max(s, gs) + 1
            if ov > 0 and ov / (ge - gs + 1) >= 0.5 and ov / (e - s + 1) >= 0.5:
                f[j] = True
                break
    return f


def boot_diff(rows, chroms, cols):
    """rows: (chrom, start, values...) per group; paired block bootstrap of mean differences cols[0] - cols[1]."""
    bid = A.block_ids(chroms, {c: A.chrom_len(c) for c in chroms})
    B = np.zeros((len(bid), 1 + len(rows[0][2])))
    for c, s, v in rows:
        b = bid[(c, (s - 1) // A.BLOCK)]
        B[b, 0] += 1
        B[b, 1:] += v
    S = A.boot_weights(len(bid)) @ B
    return S[:, 1:] / S[:, [0]]


def main():
    res = {"doc": __doc__.split("  paper_extra.py")[0]}
    # ---- 1a G11 paired CIs
    G = [g for c in A.TEST for g in P.orig_groups(c, True)]
    ts = {c: A.load_truth(c)["ts"] for c in A.TEST}
    rows, cases = [], []
    for g in G:
        sc = P.scores(g)
        s50 = s50_flags(g, ts[g["chrom"]])
        k2, k3 = int(np.argmax(sc["E2"])), int(np.argmax(sc["E3"]))
        rows.append((g["chrom"], g["s"], np.array([s50[k3], s50[k2], g["y"][k3], g["y"][k2]], float)))
        if s50[k3] and not g["y"][k3] and g["y"][k2]:
            cases.append((g, s50, sc))
    R = boot_diff(rows, A.TEST, None)
    pt = np.mean([r[2] for r in rows], 0)
    res["G11_paired"] = {"n_groups": len(rows),
                         "s50_rate": {"E3": pt[0], "E2": pt[1], "diff": pt[0] - pt[1], "ci95": np.nanpercentile(R[:, 0] - R[:, 1], [2.5, 97.5]).tolist()},
                         "supported_rate": {"E3": pt[2], "E2": pt[3], "diff": pt[2] - pt[3], "ci95": np.nanpercentile(R[:, 2] - R[:, 3], [2.5, 97.5]).tolist()}}
    # ---- 2 locus cases (rule in the module doc)
    cases.sort(key=lambda t: -t[0]["cov"])
    picked, seen = [], set()
    for g, s50, sc in cases:
        if g["chrom"] in seen:
            continue
        seen.add(g["chrom"])
        r2 = np.argsort(np.argsort(-sc["E2"])) + 1
        r3 = np.argsort(np.argsort(-sc["E3"])) + 1
        picked.append({"chrom": g["chrom"], "strand": g["strand"], "tss": int(g["ss"][0]) if g["strand"] == "+" else int(g["ee"][0]),
                       "cov": g["cov"], "ends": [{"end": int(g["ends"][j]), "span": [int(g["ss"][j]), int(g["ee"][j])],
                                                  "quant_signal": float(np.expm1(g["X"][j, P.SIG])), "reads": int(g["cnt"][j]),
                                                  "s50": bool(s50[j]), "rank_E2": int(r2[j]), "rank_E3": int(r3[j])}
                                                 for j in range(len(g["y"]))]})
        if len(picked) == 3:
            break
    res["cases_eligible"] = len(cases)
    res["cases"] = picked
    # ---- 1b region: full vs 5-feature, legal selection
    pc = A.v41_pc(A.TRAIN + A.TEST)
    Pp = {c: A.build_pairs(c, pc[c]) for c in A.TRAIN + A.TEST}
    obs = A.observed(A.TEST)
    ocv = {c: A.obs_cover(obs[c]) for c in A.TEST}
    g1 = json.load(open(f"{A.OUT}/posthoc/G1_REGION.json"))["variants"]
    bid = A.block_ids(A.TEST, {c: A.chrom_len(c) for c in A.TEST})
    W = A.boot_weights(len(bid))
    T = {}
    for name, key in (("full", "full_all_legalsel"), ("signal_density", "signal_density_all_legalsel")):
        cols, q = GR.VARIANTS[name], g1[key]["q"]
        M = GR.fit(Pp, A.TRAIN, cols, False)
        acc = {c: A.method_units(Pp[c], "model", score=GR.score(M, Pp[c], cols, True), q=q) for c in A.TEST}
        rr = {c: A.score_units(Pp[c], acc[c], ts[c], ocv[c]) for c in A.TEST}
        T[name] = A.block_table(rr, ts, {c: Pp[c]["s"][acc[c]] for c in A.TEST}, bid)
        res.setdefault("region_check", {})[key] = A.metrics(list(rr.values()))
    res["region_full_minus_signal_density"] = A.paired(W, T["full"], T["signal_density"])
    json.dump(res, open(f"{EL.OUT}/PAPER_EXTRA.json", "w"), indent=1, default=float)
    print(json.dumps({k: res[k] for k in ("G11_paired", "cases_eligible", "region_full_minus_signal_density", "region_check")}, indent=1, default=float))
    for cs in picked:
        print(cs["chrom"], cs["strand"], cs["tss"], "cov", cs["cov"], [(e["end"], round(e["quant_signal"], 2), e["reads"], e["s50"], e["rank_E2"], e["rank_E3"]) for e in cs["ends"]])


if __name__ == "__main__":
    main()
