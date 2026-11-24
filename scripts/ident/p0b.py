#!/usr/bin/env python3
"""IDENT P0b-fix (docs/ident/P0a/DEFS_P0a.md v1.1, section 8): correctness of the method (ident_core) against an
independent exhaustive oracle (oracle_indep, which does not import ident_core).  Random instances are generated with
the oracle's semantics, not the method's.

F1/B1 path enumeration identical; k* and three states == oracle for delta in {0, 1, 2, unlimited}
F2/B2 every returned witness (z, x) verified directly: support, eta <= x <= M, zero off support, coverage,
      quantitative residuals, |Z| <= K, target in / out
F3/B3 CEGIS set: produced by the target, hits every oracle counterexample, diagnosis pool == oracle pool, cost ==
      oracle optimum (branch and bound over observations, own costs); INSUFFICIENT_POOL only if oracle says impossible
F4/F5/B4 ambiguity type under ends / short2 / pe / link == oracle (structural test; pairwise LP for
      information-limited in the ORIGINAL feasible region, current quantitative evidence kept); witnesses verified
B5 monotone in delta; certificate source == oracle (no-quant family / full family / parsimony); delta_break == oracle
F7/B6 + adversarial constructed cases with asserted outcomes (edge_cases)
UNRESOLVED / SKIPPED outcomes are counted, never treated as passes.
  p0b.py [--n 600] [--procs 8] [--seed 0] -> outputs/IDENT/P0b_fix/P0B_FIX.json
"""
import argparse
import json
import math
import os
import random
import sys
import time
from multiprocessing import get_context

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ident_core as C        # noqa: E402
import oracle_indep as O      # noqa: E402

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
OUT = f"{ROOT}/outputs/IDENT/{os.environ.get('P0B_TAG', 'P0b_fix')}"
DELTAS = (0, 1, 2, None)
TYPE_TAUS = ("ends", "short2", "pe", "link")
TYPE_DELTAS = (0, 1)


def to_ev(spec):
    loc = C.Locus(spec["n"], spec["edges"], spec["sources"], spec["sinks"])
    return C.Evidence(loc, hard=spec["hard"], quant=spec["quant"])


# ------------------------------------------------------------------ random instances (oracle semantics)
def random_instance(seed):
    rng = random.Random(seed)
    for _ in range(200):
        m = rng.randint(4, 7)
        edges = {(i, i + 1) for i in range(m - 1)}
        for i in range(m):
            for d in (2, 3):
                if i + d < m and rng.random() < 0.35:
                    edges.add((i, i + d))
        src_nodes = sorted(rng.sample(range(0, max(1, m // 2)), rng.randint(1, min(3, max(1, m // 2)))))
        snk_nodes = sorted(rng.sample(range(m // 2, m), rng.randint(1, min(3, m - m // 2))))
        spec = {"n": m, "edges": sorted(edges), "sources": {f"s{k}": (v, True) for k, v in enumerate(src_nodes)},
                "sinks": {f"t{k}": (v, True) for k, v in enumerate(snk_nodes)}, "hard": [], "quant": []}
        P = O.enum_paths(spec)
        if 2 <= len(P) <= 12:
            break
    else:
        return None
    truth = rng.sample(P, rng.randint(1, min(4, len(P))))
    x = {p: float(rng.randint(1, 20)) for p in truth}
    used_s, used_t = {p[0] for p in truth}, {p[2] for p in truth}
    spec["sources"] = {s: (v, s in used_s and rng.random() < 0.8) for s, (v, _) in spec["sources"].items()}
    spec["sinks"] = {t: (v, t in used_t and rng.random() < 0.8) for t, (v, _) in spec["sinks"].items()}
    hard = set()
    for p in truth:
        s, nodes, t = p
        for a in range(len(nodes) - 1):                    # junction-level short reads
            if rng.random() < 0.6:
                hard.add(("sub", nodes[a:a + 2]))
        if rng.random() < 0.35 and len(nodes) >= 3:        # a truncated long read
            a = rng.randint(0, len(nodes) - 2)
            b = rng.randint(a + 2, len(nodes))
            hard.add(("sub", nodes[a:b]))
        if rng.random() < 0.3:                              # a read pair from this transcript
            pe = sorted(O.tau_obs(p, "pe"))
            if pe:
                hard.add(rng.choice(pe))
        if rng.random() < 0.15:                             # a full-length molecule
            hard.add(("full", s, nodes, t))
    spec["hard"] = sorted(hard)
    if rng.random() < 0.5:
        y = {o: sum(O.yields(p, o) * x[p] for p in truth) for p in truth for o in O.tau_obs(p, "ends")}
        spec["quant"].append(("ends", {o: v for o, v in y.items() if v > 0}, 0.0, 1e-6))
    if rng.random() < 0.3:
        y = {o: sum(O.yields(p, o) * x[p] for p in truth) for p in truth for o in O.tau_obs(p, "short2") if len(o[1]) == 2}
        spec["quant"].append(("short2", {o: v for o, v in y.items() if v > 0}, 0.0, 1e-6))
    return spec


# ------------------------------------------------------------------ constructed cases (plain specs)
def spec(n, edges, sources, sinks, hard=(), quant=()):
    return {"n": n, "edges": sorted(edges), "sources": sources, "sinks": sinks, "hard": sorted(set(hard)), "quant": list(quant)}


def e1_case():
    """Two TSS x two TES through a shared middle exon; end marginals (2,1) / (2,1)."""
    return spec(5, {(0, 2), (1, 2), (2, 3), (2, 4)}, {"s1": (0, True), "s2": (1, True)}, {"t1": (3, True), "t2": (4, True)},
                hard={("sub", (0, 2)), ("sub", (1, 2)), ("sub", (2, 3)), ("sub", (2, 4))},
                quant=[("ends", {("src", "s1"): 2.0, ("src", "s2"): 1.0, ("snk", "t1"): 2.0, ("snk", "t2"): 1.0}, 0.0, 1e-6)])


def known_pair_case():
    """ChatGPT round-6 counterexample: x(s1,t1) = 1 known, total 4.  New end marginals determine x(s2,t2), so 'ends'
    must NOT be information-limited for (s2,t2)."""
    return spec(5, {(0, 2), (1, 2), (2, 3), (2, 4)}, {"s1": (0, True), "s2": (1, True)}, {"t1": (3, True), "t2": (4, True)},
                hard={("sub", (0, 2)), ("sub", (1, 2)), ("sub", (2, 3)), ("sub", (2, 4))},
                quant=[("link", {("link", "s1", "t1"): 1.0}, 0.0, 1e-6), ("long", {("sub", (2,)): 4.0}, 0.0, 1e-6)])


def full_case(split=False):
    """Two TSS on the same exon, an exon skip.  A full-length molecule (s1, 0-1-2, t1) forces exactly that transcript;
    the split version (subpath + link) can be met by two other transcripts."""
    hard = {("sub", (0, 1, 2)), ("link", "s1", "t1")} if split else {("full", "s1", (0, 1, 2), "t1")}
    return spec(3, {(0, 1), (1, 2), (0, 2)}, {"s1": (0, True), "s2": (0, False)}, {"t1": (2, True)}, hard=hard)


def paired_gap_case():
    """Mates 0-1 and 4-5 with exactly one node between them: only the short middle route fits."""
    return spec(6, {(0, 1), (1, 2), (2, 3), (3, 4), (2, 4), (4, 5)}, {"s": (0, True)}, {"t": (5, True)},
                hard={("paired", (0, 1), (4, 5), 1, 1)})


def paired_bubbles_case(split=False):
    """Two bubbles; mate 1 in the upper branch of bubble 1, mate 2 in the upper branch of bubble 2.  As one read pair
    they force upper-upper; as two independent reads lower/upper mixtures also explain them."""
    hard = {("sub", (1,)), ("sub", (4,))} if split else {("paired", (1,), (4,), 0, 2)}
    return spec(7, {(0, 1), (0, 2), (1, 3), (2, 3), (3, 4), (3, 5), (4, 6), (5, 6)}, {"s": (0, True)}, {"t": (6, True)}, hard=hard)


def joint_local_case():
    edges = {(0, 1), (0, 2), (1, 3), (2, 3), (3, 4), (4, 5), (4, 6), (5, 7), (6, 7)}
    return spec(8, edges, {"s": (0, True)}, {"t": (7, True)})


def read_pair_case():
    edges = {(0, 1), (1, 2), (1, 3), (2, 4), (3, 4), (4, 5), (1, 4)}
    return spec(6, edges, {"s": (0, True)}, {"t": (5, True)}, hard={("sub", (0, 1)), ("sub", (4, 5))})


def no_evidence_case():
    return spec(3, {(0, 1), (1, 2), (0, 2)}, {"s": (0, False)}, {"t": (2, False)})


def low_abundance_case(eps_rel=0.0):
    """One TSS, two TES (t2 open).  TSS count 1.5 with eta = 1: both transcripts cannot coexist unless the
    tolerance allows a count of 2."""
    return spec(3, {(0, 1), (0, 2)}, {"s": (0, True)}, {"t1": (1, True), "t2": (2, False)},
                quant=[("ends", {("src", "s"): 1.5}, eps_rel, 1e-6)])


def shared_chain_case():
    return spec(4, {(0, 1), (1, 2), (2, 3)}, {"s1": (0, True), "s2": (1, True)}, {"t": (3, True)}, hard={("sub", (1, 2)), ("sub", (2, 3))})


def unreachable_tss_case():
    return spec(4, {(0, 1), (1, 2), (2, 3)}, {"s1": (0, True), "s2": (3, True)}, {"t": (2, True)})


def empty_support_quant_case():
    """All boundaries open, no hard reads, but a positive end count: the empty support must be rejected."""
    return spec(3, {(0, 1), (1, 2)}, {"s": (0, False)}, {"t": (2, False)}, quant=[("ends", {("src", "s"): 3.0}, 0.0, 1e-6)])


def link_open_chain_case():
    return spec(5, {(0, 1), (0, 2), (1, 3), (2, 3), (3, 4)}, {"s": (0, True)}, {"t": (4, True)}, hard={("link", "s", "t")})


def symmetric_case():
    """Two symmetric alternatives: several equal-cost diagnoses."""
    return spec(4, {(0, 1), (0, 2), (1, 3), (2, 3)}, {"s": (0, True)}, {"t": (3, True)},
                quant=[("ends", {("src", "s"): 2.0, ("snk", "t"): 2.0}, 0.0, 1e-6)])


CONSTRUCTED = {"e1": e1_case, "known_pair": known_pair_case, "full": full_case, "full_split": lambda: full_case(True),
               "paired_gap": paired_gap_case, "paired_bubbles": paired_bubbles_case,
               "paired_bubbles_split": lambda: paired_bubbles_case(True), "joint_local": joint_local_case,
               "read_pair": read_pair_case, "no_evidence": no_evidence_case, "low_abundance": low_abundance_case,
               "low_abundance_tol": lambda: low_abundance_case(0.5), "shared_chain": shared_chain_case,
               "unreachable_tss": unreachable_tss_case, "empty_support_quant": empty_support_quant_case,
               "link_open_chain": link_open_chain_case, "symmetric": symmetric_case}


# ------------------------------------------------------------------ checks on one instance
def bump(d, k, v=1):
    d[k] = d.get(k, 0) + v


def check(args):
    seed, sp = args
    t0 = time.time()
    if sp is None:
        sp = random_instance(seed)
        if sp is None:
            return {"seed": seed, "skipped": True}
    ev = to_ev(sp)
    P = ev.locus.paths
    n = len(P)
    res = {"seed": seed, "n_paths": n, "fails": [], "scen": {}, "counts": {}}
    for k in C.SOLVER_LOG:
        C.SOLVER_LOG[k] = 0
    O.ORACLE_LOG["lp_unresolved"] = 0
    O.ORACLE_LOG["lp_retried"] = 0
    fails, counts, scen = res["fails"], res["counts"], res["scen"]
    if P != O.enum_paths(sp):                                                           # F1
        fails.append("F1 path enumeration differs")
        return res
    fam = O.family(sp, P)
    fam_nq = O.family(dict(sp, quant=[]), P)
    ks_or = min((len(Z) for Z in fam), default=None)
    prev, ts_by_d = None, {}
    for d in DELTAS:
        ts = C.three_state(ev, d)
        ts_by_d[d] = ts
        if ts["model"] != "OK":
            if ts["model"] == "INFEASIBLE_MODEL" and fam:
                fails.append(f"B1 INFEASIBLE_MODEL but oracle feasible (delta={d})")
            if ts["model"] == "UNRESOLVED":
                bump(counts, "UNRESOLVED_model")
            if ts["model"] == "INFEASIBLE_MODEL" and not fam:
                bump(counts, "INFEASIBLE_MODEL_agree")
            break
        if ts["kstar"] != ks_or:
            fails.append(f"B1 kstar {ts['kstar']} != oracle {ks_or}")
        K = ts["K"]
        ost = O.states(fam, n, K)
        states = {i: v["state"] for i, v in ts["paths"].items()}
        for i in range(n):
            if states[i] == "UNRESOLVED":
                bump(counts, "UNRESOLVED_state")
            elif states[i] != ost[i]:
                fails.append(f"B1 path {i} delta={d}: {states[i]} vs oracle {ost[i]}")
            for key in ("with", "without"):                                             # F2
                st, z, x = ts["paths"][i][key]
                if st == "feasible":
                    bad = O.check_witness(sp, P, z, x, K, must_in=i if key == "with" else None, must_out=i if key == "without" else None)
                    if bad:
                        fails.append(f"B2 witness {key} path {i} delta={d}: {bad[:3]}")
                    bump(counts, "witnesses_checked")
            bump(counts, states[i])
        if prev is not None:                                                            # B5
            for i in range(n):
                if states[i] == "NECESSARY" and prev[i] != "NECESSARY":
                    fails.append(f"B5 NECESSARY grew at delta={d} path {i}")
                if states[i] == "INCOMPATIBLE" and prev[i] != "INCOMPATIBLE":
                    fails.append(f"B5 INCOMPATIBLE grew at delta={d} path {i}")
        prev = states
        if d == 0:                                                                      # certificate source
            for i in range(n):
                if states[i] != "NECESSARY":
                    continue
                src, brk = C.certificate_source(ev, i, 0)
                if all(i in Z for Z in fam_nq) and fam_nq:
                    exp = "evidence"
                elif sp["quant"] and all(i in Z for Z in fam):
                    exp = "quantification"
                else:
                    exp = "parsimony"
                if src == "UNRESOLVED":
                    bump(counts, "UNRESOLVED_source")
                elif src != exp:
                    fails.append(f"B5 source {src} vs oracle {exp} path {i}")
                k0 = min((len(Z) for Z in fam if i not in Z), default=None)
                exp_brk = math.inf if k0 is None else k0 - ks_or
                if brk != "UNRESOLVED" and brk != exp_brk:
                    fails.append(f"B5 delta_break {brk} vs oracle {exp_brk} path {i}")
                bump(counts, f"source_{src}")
                if src == "parsimony":
                    bump(scen, "parsimony_necessity")
        for i in range(n):                                                              # B3, B4
            if states[i] != "POSSIBLE":
                continue
            Fm = [Z for Z in fam if (K is None or len(Z) <= K) and i not in Z]
            cex = [[P[j] for j in Z] for Z in Fm]
            pool, cost = C.diagnosis_pool(P[i], C.PROTOCOLS, ev.locus)
            pool_or = set().union(*(O.tau_obs(P[i], tau, sp) for tau in C.PROTOCOLS))
            if set(pool) != pool_or:
                fails.append(f"B3 pool differs path {i}")
            if any(abs(cost[o] - O.cost_of(o, C.PROTOCOLS)) > 1e-12 for o in pool):
                fails.append(f"B3 costs differ path {i}")
            dg = C.diagnose(ev, i, K)
            opt = O.optimal_cover(sorted(pool_or), C.PROTOCOLS, cex)
            if dg["status"] == "SUFFICIENT":
                if not all(O.yields(P[i], o) > 0 for o in dg["obs"]):
                    fails.append(f"B3 obs not produced by target path {i}")
                if not all(any(not any(O.yields(p, o) for p in Z) for o in dg["obs"]) for Z in cex):
                    fails.append(f"B3 not sufficient path {i}")
                if opt[0] in ("impossible", "unresolved"):
                    if opt[0] == "unresolved":
                        bump(counts, "UNRESOLVED_cover_oracle")
                    else:
                        fails.append(f"B3 oracle impossible but method sufficient path {i}")
                elif abs(opt[0] - dg["cost"]) > 1e-9:
                    fails.append(f"B3 cost {dg['cost']} vs oracle optimum {opt[0]} path {i}")
                else:
                    bump(counts, "B3_optimum_verified")
                kinds = {o[0] for o in dg["obs"]}
                if len(cex) >= 2 and len(dg["obs"]) >= 2:
                    bump(scen, "multi_counterexample_cover")
                if dg["obs"] and all(o[0] == "sub" and len(o[1]) <= 2 for o in dg["obs"]) and len(dg["obs"]) >= 2:
                    bump(scen, "joint_local_sufficient")
                if "paired" in kinds:
                    bump(scen, "diagnosis_uses_read_pair")
                if "full" in kinds:
                    bump(scen, "diagnosis_needs_full_length")
                link = ("link", P[i][0], P[i][2])
                if not all(not any(O.yields(p, link) for p in Z) for Z in cex):
                    bump(scen, "link_insufficient_chain_needed")
            elif dg["status"] == "INSUFFICIENT_POOL":
                if opt[0] != "impossible":
                    fails.append(f"B3 claimed insufficient but oracle cover exists path {i}")
            else:
                bump(counts, "UNRESOLVED_diagnosis")
            if d not in TYPE_DELTAS:
                continue
            for tau in TYPE_TAUS:
                at, w = C.ambiguity_type(ev, i, K, tau)
                bump(counts, f"{tau}_{at}")
                if at == "UNRESOLVED":
                    continue
                if O.structural(P, Fm, i, tau):
                    exp = "STRUCTURAL"
                else:
                    il = O.info_limited(sp, P, fam, i, K, tau)
                    exp = None if il == "SKIPPED" else ("INFORMATION_LIMITED" if il else "QUANT_DISTINGUISHABLE")
                if exp is None:
                    bump(counts, "type_oracle_skipped")
                    if at == "STRUCTURAL":
                        fails.append(f"B4 STRUCTURAL but oracle not structural path {i} tau={tau}")
                elif at != exp:
                    fails.append(f"B4 type {at} vs oracle {exp} path {i} tau={tau} delta={d}")
                else:
                    bump(counts, "B4_type_verified")
                if at == "INFORMATION_LIMITED":
                    bump(scen, f"info_limited_{tau}")
                    bad = O.check_witness(sp, P, w["z1"], w["x1"], K, must_in=i) + O.check_witness(sp, P, w["z2"], w["x2"], K, must_out=i)
                    space = sorted({o for p in P for o in O.tau_obs(p, tau)})
                    for o in space:
                        a = sum(O.yields(P[j], o) * w["x1"][j] for j in range(n))
                        b = sum(O.yields(P[j], o) * w["x2"][j] for j in range(n))
                        if abs(a - b) > 1e-6:
                            bad.append(f"tau obs differ {o}")
                            break
                    if bad:
                        fails.append(f"B4 info-limited witness path {i} tau={tau}: {bad[:3]}")
    for i in range(n):                                                                  # evidence-forced at all delta
        if all(i in Z for Z in fam_nq) and fam_nq:
            for d, ts in ts_by_d.items():
                if ts["model"] == "OK" and ts["paths"][i]["state"] not in ("NECESSARY", "UNRESOLVED"):
                    fails.append(f"B5 evidence-forced path {i} not NECESSARY at delta={d}")
    kinds = {o[0] for o in sp["hard"]}
    for k in ("full", "paired"):
        if k in kinds:
            bump(scen, f"instances_with_{k}_evidence")
    res["solver"] = dict(C.SOLVER_LOG, oracle_lp_retried=O.ORACLE_LOG["lp_retried"])
    if O.ORACLE_LOG["lp_unresolved"]:
        res["oracle_lp_unresolved"] = O.ORACLE_LOG["lp_unresolved"]
        fails.append(f"ORACLE lp unresolved x{O.ORACLE_LOG['lp_unresolved']} (comparison not trustworthy)")
    res["seconds"] = round(time.time() - t0, 2)
    return res


# ------------------------------------------------------------------ asserted outcomes (B6 + adversarial)
def _idx(P, path):
    return P.index(path)


def edge_cases():
    out = {}
    # B6: conflicts, timeouts, open boundary, truncated read, path cap
    loc = C.Locus(3, {(0, 1), (1, 2)}, {"s": (0, True)}, {"t": (2, True)})
    out["conflict_unproducible_obs"] = C.three_state(C.Evidence(loc, hard={("sub", (0, 2))}))["model"] == "INFEASIBLE_MODEL"
    ev2 = C.Evidence(loc, quant=[("ends", {("src", "s"): 5.0, ("snk", "t"): 9.0}, 0.0, 1e-6)])
    out["conflict_quant"] = C.three_state(ev2)["model"] == "INFEASIBLE_MODEL"
    out["unreachable_observed_tss"] = C.three_state(to_ev(unreachable_tss_case()))["model"] == "INFEASIBLE_MODEL"
    loc3 = C.Locus(3, {(0, 1), (1, 2)}, {"s1": (0, True), "s2": (1, False)}, {"t": (2, True)})
    ts3 = C.three_state(C.Evidence(loc3), None)
    out["open_boundary_not_forced"] = all(v["state"] != "NECESSARY" for i, v in ts3["paths"].items() if loc3.paths[i][0] == "s2")
    loc4 = C.Locus(4, {(0, 1), (1, 2), (2, 3), (0, 2)}, {"s": (0, True)}, {"t": (3, True)})
    ts4 = C.three_state(C.Evidence(loc4, hard={("sub", (1, 2))}), 0)
    out["truncated_read_as_subpath"] = ts4["paths"][_idx(loc4.paths, ("s", (0, 1, 2, 3), "t"))]["state"] == "NECESSARY"
    m = 16
    edges = {(i, i + 1) for i in range(m - 1)} | {(i, i + 2) for i in range(m - 2)}
    loc5 = C.Locus(m, edges, {"s": (0, True)}, {"t": (m - 1, True)})
    ev5 = C.Evidence(loc5, hard={("sub", (i, i + 1)) for i in range(m - 1)})
    r5 = C.three_state(ev5, 0, time_limit=1e-4)
    out["timeout_unresolved"] = r5["model"] == "UNRESOLVED" or any(v["state"] == "UNRESOLVED" for v in r5.get("paths", {}).values())
    out["timeout_certificate_source_unresolved"] = C.certificate_source(ev5, 0, 0, time_limit=1e-4)[0] == "UNRESOLVED"
    big = {(2 * k, 2 * k + 1) for k in range(14)} | {(2 * k, 2 * k + 2) for k in range(14)} | {(2 * k + 1, 2 * k + 2) for k in range(14)}
    locb = C.Locus(29, big, {"s": (0, True)}, {"t": (28, True)})
    out["path_cap_unresolved"] = locb.too_many and C.three_state(C.Evidence(locb))["model"] == "UNRESOLVED"
    # E1: parsimony-made necessity, information-limited pairing under ends
    ev = to_ev(e1_case()); P = ev.locus.paths
    j = [k for k, p in enumerate(P) if (p[0], p[2]) == ("s2", "t2")][0]
    t0, t1 = C.three_state(ev, 0), C.three_state(ev, 1)
    out["e1_delta0_necessary"] = t0["paths"][j]["state"] == "NECESSARY"
    out["e1_delta1_possible"] = t1["paths"][j]["state"] == "POSSIBLE"
    out["e1_source_parsimony"] = C.certificate_source(ev, j, 0)[0] == "parsimony"
    out["e1_ends_information_limited"] = C.ambiguity_type(ev, j, t1["K"], "ends")[0] == "INFORMATION_LIMITED"
    out["e1_link_structural"] = C.ambiguity_type(ev, j, t1["K"], "link")[0] == "STRUCTURAL"
    # ChatGPT round-6 counterexample: known pair -> ends not information-limited
    ev = to_ev(known_pair_case()); P = ev.locus.paths
    j = [k for k, p in enumerate(P) if (p[0], p[2]) == ("s2", "t2")][0]
    t1 = C.three_state(ev, 1)
    out["known_pair_possible"] = t1["model"] == "OK" and t1["paths"][j]["state"] == "POSSIBLE"
    out["known_pair_ends_not_info_limited"] = C.ambiguity_type(ev, j, t1["K"], "ends")[0] == "QUANT_DISTINGUISHABLE"
    # full-length molecule vs split constraints
    ev = to_ev(full_case()); P = ev.locus.paths
    j = _idx(P, ("s1", (0, 1, 2), "t1"))
    out["full_forces_exact_transcript"] = C.three_state(ev, None)["paths"][j]["state"] == "NECESSARY" and C.certificate_source(ev, j, 0)[0] == "evidence"
    ev = to_ev(full_case(True))
    out["split_full_not_forced"] = C.three_state(ev, None)["paths"][j]["state"] == "POSSIBLE"
    # read pairs
    ev = to_ev(paired_gap_case()); P = ev.locus.paths
    ts = C.three_state(ev, None)
    out["paired_gap_selects_route"] = (ts["paths"][_idx(P, ("s", (0, 1, 2, 4, 5), "t"))]["state"] == "NECESSARY"
                                      and ts["paths"][_idx(P, ("s", (0, 1, 2, 3, 4, 5), "t"))]["state"] == "POSSIBLE")
    ev = to_ev(paired_bubbles_case()); P = ev.locus.paths
    j = _idx(P, ("s", (0, 1, 3, 4, 6), "t"))
    out["paired_same_molecule_forces"] = C.three_state(ev, None)["paths"][j]["state"] == "NECESSARY"
    ev = to_ev(paired_bubbles_case(True))
    out["paired_split_not_forced"] = C.three_state(ev, None)["paths"][j]["state"] == "POSSIBLE"
    ev = to_ev(read_pair_case()); P = ev.locus.paths
    t = C.three_state(ev, 0)
    out["read_pair_middle_ambiguous"] = len(P) == 3 and all(v["state"] == "POSSIBLE" for v in t["paths"].values())
    j = _idx(P, ("s", (0, 1, 2, 4, 5), "t"))
    dg = C.diagnose(ev, j, t["K"])
    out["read_pair_middle_resolved_by_one_junction"] = dg["status"] == "SUFFICIENT" and dg["cost"] == 1.0
    # joint local sufficiency
    ev = to_ev(joint_local_case()); P = ev.locus.paths
    j = _idx(P, ("s", (0, 1, 3, 4, 5, 7), "t"))
    t = C.three_state(ev, 0)
    dg = C.diagnose(ev, j, t["K"], protocols=("short2",))
    out["joint_local_two_junctions_short2"] = dg["status"] == "SUFFICIENT" and len(dg["obs"]) == 2 and abs(dg["cost"] - 2.0) < 1e-9 \
        and all(o[0] == "sub" and len(o[1]) == 2 for o in dg["obs"])
    out["joint_local_short2_alone_structural"] = C.ambiguity_type(ev, j, t["K"], "short2")[0] == "STRUCTURAL"
    dg = C.diagnose(ev, j, t["K"])          # all protocols: one read pair spanning both bubbles is cheaper (1.5)
    out["joint_local_read_pair_cheapest"] = dg["status"] == "SUFFICIENT" and len(dg["obs"]) == 1 and dg["obs"][0][0] == "paired" \
        and abs(dg["cost"] - 1.5) < 1e-9
    # empty support / abundance floor / tolerance
    t = C.three_state(to_ev(no_evidence_case()), 0)
    out["no_evidence_kstar0_all_incompatible"] = t["kstar"] == 0 and all(v["state"] == "INCOMPATIBLE" for v in t["paths"].values())
    out["no_evidence_unlimited_all_possible"] = all(v["state"] == "POSSIBLE" for v in C.three_state(to_ev(no_evidence_case()), None)["paths"].values())
    out["empty_support_rejected_by_quant"] = C.three_state(to_ev(empty_support_quant_case()), 0)["kstar"] == 1
    ev = to_ev(low_abundance_case()); P = ev.locus.paths
    j = _idx(P, ("s", (0, 2), "t2"))
    out["abundance_floor_excludes"] = C.three_state(ev, None)["paths"][j]["state"] == "INCOMPATIBLE"
    out["tolerance_allows"] = C.three_state(to_ev(low_abundance_case(0.5)), None)["paths"][j]["state"] == "POSSIBLE"
    t = C.three_state(to_ev(shared_chain_case()), 0)
    out["shared_chain_both_necessary"] = all(v["state"] == "NECESSARY" for v in t["paths"].values())
    # permanent regression: instances where HiGHS presolve wrongly reported infeasible (P0b-fix, 2026-11-24)
    for s in (329, 391, 414):
        out[f"regression_presolve_seed{s}"] = not check((s, None))["fails"]
    # class state: two paths sharing the chain 1-2 through different TSS; the class is necessary, neither path is
    ev = to_ev(spec(3, {(0, 1), (1, 2)}, {"s1": (0, False), "s2": (1, False)}, {"t": (2, True)}, hard={("sub", (1, 2))}))
    P = ev.locus.paths
    cls = [k for k, p in enumerate(P) if p[1][-2:] == (1, 2)]
    t = C.three_state(ev, None)
    out["class_necessary_members_possible"] = (C.class_state(ev, cls, None) == "NECESSARY"
                                               and all(t["paths"][k]["state"] == "POSSIBLE" for k in cls))
    # full-length observation not offered for an open-ended path
    loc = C.Locus(3, {(0, 1), (1, 2)}, {"s": (0, False)}, {"t": (2, True)})
    out["no_full_for_open_end"] = not any(o[0] == "full" for o in C.protocol_obs(loc.paths[0], "long", loc))
    out["prune_unreachable"] = C.Locus(4, {(0, 1), (1, 2), (0, 3)}, {"s": (0, True)}, {"t": (2, True)}).paths == [("s", (0, 1, 2), "t")]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=600)
    ap.add_argument("--procs", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    ec = edge_cases()
    print("edge cases", json.dumps(ec), flush=True)
    jobs = [(name, f()) for name, f in CONSTRUCTED.items()] + [(a.seed * 100000 + k, None) for k in range(a.n)]
    rs = []
    with get_context("spawn").Pool(a.procs) as pool, open(f"{OUT}/instances.jsonl", "w") as fh:
        for r in pool.imap_unordered(check, jobs, chunksize=1):
            rs.append(r)
            fh.write(json.dumps(r, default=str) + "\n"); fh.flush()
    done = [r for r in rs if not r.get("skipped")]
    fails = [(r["seed"], f) for r in done for f in r["fails"]]
    counts, scen = {}, {}
    for r in done:
        for k, v in r["counts"].items():
            bump(counts, k, v)
        for k, v in r["scen"].items():
            bump(scen, k, v)
    solver = {}
    for r in done:
        for k, v in r.get("solver", {}).items():
            bump(solver, k, v)
    summary = {"solver_rechecks": solver, "instances": len(done), "constructed": len(CONSTRUCTED), "skipped": len(rs) - len(done), "edge_cases": ec,
               "edge_cases_all_pass": all(ec.values()), "n_fail": len(fails), "fails": fails[:200],
               "constructed_fails": [(r["seed"], r["fails"]) for r in done if isinstance(r["seed"], str) and r["fails"]],
               "state_counts": counts, "scenarios": scen, "paths_total": sum(r["n_paths"] for r in done),
               "max_seconds": max(r.get("seconds", 0) for r in done), "seconds": round(time.time() - t0, 1)}
    json.dump(summary, open(f"{OUT}/P0B_FIX.json", "w"), indent=1, default=str)
    import scipy
    summary["versions"] = {"scipy": scipy.__version__, "numpy": np.__version__}
    json.dump(summary, open(f"{OUT}/P0B_FIX.json", "w"), indent=1, default=str)
    print(json.dumps({k: summary[k] for k in ("versions", "solver_rechecks", "instances", "skipped", "edge_cases_all_pass", "n_fail", "constructed_fails",
                                              "state_counts", "scenarios", "paths_total", "max_seconds", "seconds")}, indent=1, default=str))
    if fails:
        print("FIRST FAILS", fails[:30])


if __name__ == "__main__":
    main()
