#!/usr/bin/env python3
"""IDENT P0b (docs/ident/P0a/DEFS_P0a.md v1, section 8): correctness of the method against exhaustive enumeration.

B1 three states == oracle for delta in {0, 1, 2, unlimited}
B2 every witness (with / without) is a feasible support under an independent checker
B3 the CEGIS sufficient set hits every oracle counterexample and its cost equals the optimum from all oracle
   counterexamples (independent cover solve); INSUFFICIENT_POOL claims verified by an unhittable oracle counterexample
B4 every INFORMATION_LIMITED witness: both supports feasible (coverage, boundaries, size <= K), path in / out,
   normalised tau vectors equal
B5 NECESSARY and INCOMPATIBLE sets never grow with delta; 'evidence'-sourced necessities are NECESSARY at every delta
B6 constructed edge cases: conflicting evidence -> INFEASIBLE_MODEL; time limit -> UNRESOLVED; open boundaries not
   forced; truncated reads as subpaths
Required scenario cases (counted on random graphs + constructed): multiple counterexamples needing different
observations; parsimony-made necessity; ends linked but chain open; information-limited under some protocol;
two local observations jointly sufficient.
  p0b.py [--n 600] [--procs 8] [--seed 0] -> outputs/IDENT/P0b/P0B.json
"""
import argparse
import json
import math
import os
import random
import sys
import time
from multiprocessing import Pool, get_context

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ident_core as C        # noqa: E402

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
OUT = f"{ROOT}/outputs/IDENT/P0b"
DELTAS = (0, 1, 2, None)


# ------------------------------------------------------------------ instances
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
        sources = {f"s{k}": (v, True) for k, v in enumerate(src_nodes)}
        sinks = {f"t{k}": (v, True) for k, v in enumerate(snk_nodes)}
        loc = C.Locus(m, edges, sources, sinks)
        n = len(loc.paths)
        if 2 <= n <= 12:
            break
    else:
        return None
    k = rng.randint(1, min(4, n))
    truth = rng.sample(range(n), k)
    x = {j: float(rng.randint(1, 20)) for j in truth}
    P = loc.paths
    used_s = {P[j][0] for j in truth}
    used_t = {P[j][2] for j in truth}
    # boundaries: unused ones are open; some used ones become open (truncated evidence)
    for s in list(sources):
        sources[s] = (sources[s][0], s in used_s and rng.random() < 0.8)
    for t in list(sinks):
        sinks[t] = (sinks[t][0], t in used_t and rng.random() < 0.8)
    loc = C.Locus(m, edges, sources, sinks)
    P = loc.paths
    truth = [P.index(p) for p in [C.Locus(m, edges, sources, sinks).paths[j] for j in truth]]
    hard = set()
    for j in truth:
        nodes = P[j][1]
        for a in range(len(nodes) - 1):                    # junction-level short reads
            if rng.random() < 0.7:
                hard.add(("sub", nodes[a:a + 2]))
        if rng.random() < 0.4 and len(nodes) >= 3:         # a truncated long read
            a = rng.randint(0, len(nodes) - 2)
            b = rng.randint(a + 2, len(nodes))
            hard.add(("sub", nodes[a:b]))
        if rng.random() < 0.15:                            # a full-length molecule
            hard.add(("sub", nodes)); hard.add(("link", P[j][0], P[j][2]))
    quant = []
    xv = np.zeros(len(P))
    for j in truth:
        xv[j] = x[j]
    if rng.random() < 0.5:
        y = {o: float(sum(C.produces(P[j], o) * xv[j] for j in range(len(P)))) for o in C.protocol_space(loc, "ends")}
        quant.append(("ends", {o: v for o, v in y.items() if v > 0}, 0.0, 1e-6))
    if rng.random() < 0.3:
        y = {o: float(sum(C.produces(P[j], o) * xv[j] for j in range(len(P)))) for o in C.protocol_space(loc, "short2") if len(o[1]) == 2}
        quant.append(("short2", {o: v for o, v in y.items() if v > 0}, 0.0, 1e-6))
    ev = C.Evidence(loc, hard=hard, quant=quant)
    return ev, truth


def e1_case():
    """Two TSS x two TES through a shared middle exon; end marginals (2,1) / (2,1): parsimony-made necessity and
    information-limited pairing under 'ends'."""
    # nodes: 0=A1 1=A2 2=C 3=D1 4=D2
    loc = C.Locus(5, {(0, 2), (1, 2), (2, 3), (2, 4)}, {"s1": (0, True), "s2": (1, True)}, {"t1": (3, True), "t2": (4, True)})
    y = {("src", "s1"): 2.0, ("src", "s2"): 1.0, ("snk", "t1"): 2.0, ("snk", "t2"): 1.0}
    hard = {("sub", (0, 2)), ("sub", (1, 2)), ("sub", (2, 3)), ("sub", (2, 4))}
    return C.Evidence(loc, hard=hard, quant=[("ends", y, 0.0, 1e-6)])


def joint_local_case():
    """Two independent bubbles joined by a linker exon: p = 0-1-3-4-5-7.  Each competitor differs in one bubble only,
    so one junction per bubble is needed; the only single observation spanning both is a 4-node subpath (long reads)."""
    edges = {(0, 1), (0, 2), (1, 3), (2, 3), (3, 4), (4, 5), (4, 6), (5, 7), (6, 7)}
    loc = C.Locus(8, edges, {"s": (0, True)}, {"t": (7, True)})
    return C.Evidence(loc)


def read_pair_case():
    """Read pair: both mates observed (0-1 and 4-5), several middle routes compatible with both."""
    edges = {(0, 1), (1, 2), (1, 3), (2, 4), (3, 4), (4, 5), (1, 4)}
    loc = C.Locus(6, edges, {"s": (0, True)}, {"t": (5, True)})
    return C.Evidence(loc, hard={("sub", (0, 1)), ("sub", (4, 5))})


# ------------------------------------------------------------------ checks on one instance
def check(args):
    seed, inst = args
    t0 = time.time()
    if inst is None:
        built = random_instance(seed)
        if built is None:
            return {"seed": seed, "skipped": True}
        ev, truth = built
    else:
        ev, truth = inst, None
    P, n = ev.locus.paths, len(ev.locus.paths)
    fam = C.oracle_family(ev, True)
    res = {"seed": seed, "n_paths": n, "fails": [], "scen": {}, "counts": {}}
    prev = None
    nec_evidence = set()
    ts_by_d = {}
    for d in DELTAS:
        ts = C.three_state(ev, d)
        ts_by_d[d] = ts
        kstar = min((len(Z) for Z in fam), default=None)
        if ts["model"] != "OK":
            if fam:
                res["fails"].append(f"B1 model {ts['model']} but oracle feasible (delta={d})")
            break
        if ts["kstar"] != kstar:
            res["fails"].append(f"B1 kstar {ts['kstar']} != oracle {kstar}")
        K = ts["K"]
        ost = C.oracle_states(fam, n, K)
        states = {i: v["state"] for i, v in ts["paths"].items()}
        for i in range(n):
            if states[i] != ost[i]:
                res["fails"].append(f"B1 path {i} delta={d}: {states[i]} vs oracle {ost[i]}")
            # B2 witnesses
            for key in ("with", "without"):
                st, z, x = ts["paths"][i][key]
                if st == "feasible":
                    Z = tuple(np.nonzero(z)[0])
                    ok, _ = C.support_feasible(ev, Z, True)
                    good = ok and (K is None or len(Z) <= K) and ((i in Z) == (key == "with"))
                    xs = np.asarray(x)
                    good = good and all(xs[j] >= C.ETA - 1e-6 for j in Z) and all(abs(xs[j]) < 1e-6 for j in range(n) if j not in Z)
                    if not good:
                        res["fails"].append(f"B2 witness {key} path {i} delta={d}")
            res["counts"][states[i]] = res["counts"].get(states[i], 0) + 1
        # B5 monotonicity
        if prev is not None:
            for i in range(n):
                if states[i] == "NECESSARY" and prev[i] != "NECESSARY":
                    res["fails"].append(f"B5 NECESSARY grew at delta={d} path {i}")
                if states[i] == "INCOMPATIBLE" and prev[i] != "INCOMPATIBLE":
                    res["fails"].append(f"B5 INCOMPATIBLE grew at delta={d} path {i}")
        prev = states
        # certificate sources (delta = 0 only)
        if d == 0:
            for i in range(n):
                if states[i] == "NECESSARY":
                    src, brk = C.certificate_source(ev, i, 0)
                    res["counts"][f"source_{src}"] = res["counts"].get(f"source_{src}", 0) + 1
                    if src == "evidence":
                        nec_evidence.add(i)
                    if src == "parsimony":
                        res["scen"]["parsimony_necessity"] = res["scen"].get("parsimony_necessity", 0) + 1
            # B3 diagnosis and B4 information-limited
            for i in range(n):
                if states[i] != "POSSIBLE":
                    continue
                Fm = [Z for Z in fam if len(Z) <= K and i not in Z]
                cex_all = [[P[j] for j in Z] for Z in Fm]
                dg = C.diagnose(ev, i, K)
                pool = sorted({o for tau in C.PROTOCOLS for o in C.protocol_obs(P[i], tau)})
                cost = {o: C.obs_cost(o, C.PROTOCOLS) for o in pool}
                if dg["status"] == "SUFFICIENT":
                    if not all(any(C.hits(Z, o) for o in dg["obs"]) for Z in cex_all):
                        res["fails"].append(f"B3 not sufficient path {i}")
                    opt = C.min_cost_cover(pool, cost, cex_all)
                    if opt is None or abs(sum(cost[o] for o in opt) - dg["cost"]) > 1e-9:
                        res["fails"].append(f"B3 cost {dg['cost']} vs oracle optimum path {i}")
                    if len(dg["obs"]) >= 2 and len({o for o in dg["obs"]}) >= 2:
                        res["scen"]["multi_counterexample_cover"] = res["scen"].get("multi_counterexample_cover", 0) + 1
                    kinds = {o[0] for o in dg["obs"]}
                    if all(o[0] == "sub" and len(o[1]) <= 2 for o in dg["obs"]) and len(dg["obs"]) >= 2:
                        single = any(all(C.hits(Z, o) for Z in cex_all) for o in pool if o[0] == "sub" and len(o[1]) <= 2)
                        if not single:
                            res["scen"]["joint_local_sufficient"] = res["scen"].get("joint_local_sufficient", 0) + 1
                    if kinds == {"sub"} and any(o[0] == "link" for o in pool):
                        # ends-link alone insufficient but chain observation needed
                        link = [o for o in pool if o[0] == "link"][0]
                        if not all(C.hits(Z, link) for Z in cex_all):
                            res["scen"]["link_insufficient_chain_needed"] = res["scen"].get("link_insufficient_chain_needed", 0) + 1
                elif dg["status"] == "INSUFFICIENT_POOL":
                    if all(any(C.hits(Z, o) for o in pool) for Z in cex_all):
                        res["fails"].append(f"B3 claimed insufficient but every counterexample hittable path {i}")
                for tau in ("ends", "short2"):
                    at, w = C.ambiguity_type(ev, i, K, tau)
                    res["counts"][f"{tau}_{at}"] = res["counts"].get(f"{tau}_{at}", 0) + 1
                    if at == "INFORMATION_LIMITED":
                        res["scen"][f"info_limited_{tau}"] = res["scen"].get(f"info_limited_{tau}", 0) + 1
                        z1, x1, z2, x2, B = w["z1"], w["x1"], w["z2"], w["x2"], w["B"]
                        Z1, Z2 = tuple(np.nonzero(z1)[0]), tuple(np.nonzero(z2)[0])
                        ok = C.support_feasible(ev, Z1, False)[0] and C.support_feasible(ev, Z2, False)[0]
                        ok = ok and len(Z1) <= K and len(Z2) <= K and i in Z1 and i not in Z2
                        ok = ok and np.allclose(B @ x1, B @ x2, atol=1e-6) and abs((B @ x1).sum() - 1) < 1e-6 and abs((B @ x2).sum() - 1) < 1e-6
                        if not ok:
                            res["fails"].append(f"B4 info-limited witness path {i} tau={tau}")
    for i in nec_evidence:
        for d, ts in ts_by_d.items():
            if ts["model"] == "OK" and ts["paths"][i]["state"] != "NECESSARY":
                res["fails"].append(f"B5 evidence-forced path {i} not NECESSARY at delta={d}")
    res["seconds"] = round(time.time() - t0, 2)
    return res


def edge_cases():
    out = {}
    # conflicting evidence: an observation no candidate can produce
    loc = C.Locus(3, {(0, 1), (1, 2)}, {"s": (0, True)}, {"t": (2, True)})
    out["conflict_unproducible_obs"] = C.three_state(C.Evidence(loc, hard={("sub", (0, 2))}))["model"] == "INFEASIBLE_MODEL"
    # conflicting quantities: end marginals that disagree
    loc2 = C.Locus(3, {(0, 1), (1, 2)}, {"s": (0, True)}, {"t": (2, True)})
    ev2 = C.Evidence(loc2, quant=[("ends", {("src", "s"): 5.0, ("snk", "t"): 9.0}, 0.0, 1e-6)])
    out["conflict_quant"] = C.three_state(ev2)["model"] == "INFEASIBLE_MODEL"
    # open boundary not forced
    loc3 = C.Locus(3, {(0, 1), (1, 2)}, {"s1": (0, True), "s2": (1, False)}, {"t": (2, True)})
    ts3 = C.three_state(C.Evidence(loc3), None)
    st = {loc3.paths[i]: v["state"] for i, v in ts3["paths"].items()}
    out["open_boundary_not_forced"] = all(v != "NECESSARY" for p, v in st.items() if p[0] == "s2")
    # truncated read: subpath not touching any boundary
    loc4 = C.Locus(4, {(0, 1), (1, 2), (2, 3), (0, 2)}, {"s": (0, True)}, {"t": (3, True)})
    ts4 = C.three_state(C.Evidence(loc4, hard={("sub", (1, 2))}), 0)
    st4 = {loc4.paths[i]: v["state"] for i, v in ts4["paths"].items()}
    out["truncated_read_as_subpath"] = st4.get(("s", (0, 1, 2, 3), "t")) == "NECESSARY"
    # time limit -> UNRESOLVED: a wide graph with a tiny time limit
    m = 16
    edges = {(i, i + 1) for i in range(m - 1)} | {(i, i + 2) for i in range(m - 2)}
    loc5 = C.Locus(m, edges, {"s": (0, True)}, {"t": (m - 1, True)})
    hard5 = {("sub", (i, i + 1)) for i in range(m - 1)}
    r5 = C.three_state(C.Evidence(loc5, hard=hard5), 0, time_limit=1e-4)
    out["timeout_unresolved"] = r5["model"] == "UNRESOLVED" or any(v["state"] == "UNRESOLVED" for v in r5.get("paths", {}).values())
    # E1-type case
    ev = e1_case()
    P = ev.locus.paths
    idx = {(p[0], p[2]): i for i, p in enumerate(P)}
    t0, t1 = C.three_state(ev, 0), C.three_state(ev, 1)
    out["e1_delta0_s2t2_necessary"] = t0["paths"][idx[("s2", "t2")]]["state"] == "NECESSARY"
    out["e1_delta1_s2t2_possible"] = t1["paths"][idx[("s2", "t2")]]["state"] == "POSSIBLE"
    out["e1_source_parsimony"] = C.certificate_source(ev, idx[("s2", "t2")], 0)[0] == "parsimony"
    K1 = t1["K"]
    at_ends, _ = C.ambiguity_type(ev, idx[("s2", "t2")], K1, "ends")
    at_link, dl = C.ambiguity_type(ev, idx[("s2", "t2")], K1, "link")
    at_s3, _ = C.ambiguity_type(ev, idx[("s2", "t2")], K1, "short3")
    out["e1_ends_information_limited"] = at_ends == "INFORMATION_LIMITED"
    out["e1_link_structural"] = at_link == "STRUCTURAL"
    out["e1_short3_structural"] = at_s3 == "STRUCTURAL"
    # joint local sufficiency: two short2 junctions together, no single short observation suffices
    ev = joint_local_case(); P = ev.locus.paths
    i = P.index(("s", (0, 1, 3, 4, 5, 7), "t"))
    t = C.three_state(ev, 0)
    dg = C.diagnose(ev, i, t["K"])
    out["joint_local_possible"] = t["paths"][i]["state"] == "POSSIBLE"
    out["joint_local_two_junctions"] = (dg["status"] == "SUFFICIENT" and len(dg["obs"]) == 2
                                        and all(o[0] == "sub" and len(o[1]) == 2 for o in dg["obs"]) and abs(dg["cost"] - 2.0) < 1e-9)
    out["joint_local_short2_alone_structural"] = C.ambiguity_type(ev, i, t["K"], "short2")[0] == "STRUCTURAL"
    # read pair compatible with several middle routes
    ev = read_pair_case(); P = ev.locus.paths
    t = C.three_state(ev, 0)
    out["read_pair_middle_ambiguous"] = len(P) == 3 and all(v["state"] == "POSSIBLE" for v in t["paths"].values())
    i = P.index(("s", (0, 1, 2, 4, 5), "t"))
    dg = C.diagnose(ev, i, t["K"])
    out["read_pair_middle_resolved_by_junction"] = (dg["status"] == "SUFFICIENT" and dg["cost"] == 1.0
                                                    and all(o[0] == "sub" and 2 in o[1] for o in dg["obs"]))
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
    jobs = [(a.seed * 100000 + k, None) for k in range(a.n)] + [(-1, e1_case()), (-2, joint_local_case()), (-3, read_pair_case())]
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
            counts[k] = counts.get(k, 0) + v
        for k, v in r["scen"].items():
            scen[k] = scen.get(k, 0) + v
    summary = {"instances": len(done), "skipped": len(rs) - len(done), "edge_cases": ec, "edge_cases_all_pass": all(ec.values()),
               "n_fail": len(fails), "fails": fails[:200], "state_counts": counts, "scenarios": scen,
               "paths_total": sum(r["n_paths"] for r in done), "max_seconds": max(r["seconds"] for r in done),
               "seconds": round(time.time() - t0, 1)}
    json.dump(summary, open(f"{OUT}/P0B.json", "w"), indent=1, default=str)
    print(json.dumps({k: summary[k] for k in ("instances", "skipped", "edge_cases_all_pass", "n_fail", "state_counts", "scenarios", "paths_total", "max_seconds", "seconds")}, indent=1, default=str))
    if fails:
        print("FIRST FAILS", fails[:20])


if __name__ == "__main__":
    main()
