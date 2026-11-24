#!/usr/bin/env python3
"""Independent exhaustive oracle for IDENT P0b-fix (DEFS_P0a.md v1.1).  Deliberately does NOT import ident_core:
own path enumeration, own observation semantics, own protocol observation pools, own costs, own cover optimum
(subset enumeration), own feasibility LPs (HiGHS interior point; the method uses HiGHS branch-and-bound).

A locus spec is plain data: {"n", "edges", "sources": {name: (node, observed)}, "sinks": {...}, "hard": [...],
"quant": [(tau, {obs: y}, eps_rel, eps_abs), ...]}.
"""
import itertools

import numpy as np
from scipy.optimize import linprog

ETA, BIGM = 1.0, 1e4
GAP = (0, 2)          # read-pair gap (nodes) and mate length <= 2, as in the definitions
MATE = 2
COST = {"ends": 0.5, "short2": 1.0, "short3": 1.5, "pe": 1.5, "long": 3.0, "link": 5.0}
TOL = 1e-6
ORACLE_LOG = {"lp_unresolved": 0, "lp_retried": 0}


# ------------------------------------------------------------------ paths and observation semantics
def enum_paths(spec):
    nxt = {}
    for u, v in spec["edges"]:
        nxt.setdefault(u, set()).add(v)
    ends_at = {}
    for t, (v, _) in spec["sinks"].items():
        ends_at.setdefault(v, set()).add(t)
    found = []

    def walk(s, trail):
        for t in ends_at.get(trail[-1], ()):
            found.append((s, tuple(trail), t))
        for w in sorted(nxt.get(trail[-1], ())):
            walk(s, trail + [w])
    for s, (v0, _) in spec["sources"].items():
        walk(s, [v0])
    return sorted(set(found))


def _starts(nodes, piece):
    piece = list(piece)
    return [k for k in range(0, len(nodes) - len(piece) + 1) if list(nodes[k:k + len(piece)]) == piece]


def yields(path, o):
    """How many times transcript `path` yields observation o (independent re-implementation)."""
    s, nodes, t = path
    kind = o[0]
    if kind == "src":
        return 1 if o[1] == s else 0
    if kind == "snk":
        return 1 if o[1] == t else 0
    if kind == "link":
        return 1 if (o[1], o[2]) == (s, t) else 0
    if kind == "full":
        return 1 if (o[1], tuple(o[2]), o[3]) == (s, tuple(nodes), t) else 0
    if kind == "sub":
        return len(_starts(nodes, o[1]))
    if kind == "paired":
        a, b, lo, hi = o[1], o[2], o[3], o[4]
        c = 0
        for i in _starts(nodes, a):
            for j in _starts(nodes, b):
                gap = j - i - len(a)
                if lo <= gap <= hi:
                    c += 1
        return c
    raise ValueError(o)


def tau_obs(path, tau, spec=None):
    s, nodes, t = path
    L = len(nodes)
    windows = {tuple(nodes[i:i + w]) for w in range(1, L + 1) for i in range(0, L - w + 1)}
    if tau == "ends":
        return {("src", s), ("snk", t)}
    if tau == "link":
        return {("link", s, t)}
    if tau == "short2":
        return {("sub", w) for w in windows if len(w) <= 2}
    if tau == "short3":
        return {("sub", w) for w in windows if len(w) <= 3}
    if tau == "long":
        full_ok = spec is None or (spec["sources"][s][1] and spec["sinks"][t][1])
        return {("sub", w) for w in windows} | ({("full", s, tuple(nodes), t)} if full_ok else set())
    if tau == "pe":
        out = set()
        for i in range(L):
            for la in range(1, MATE + 1):
                for gap in range(GAP[0], GAP[1] + 1):
                    for lb in range(1, MATE + 1):
                        j = i + la + gap
                        if i + la <= L and j + lb <= L:
                            out.add(("paired", tuple(nodes[i:i + la]), tuple(nodes[j:j + lb]), GAP[0], GAP[1]))
        return out
    raise ValueError(tau)


def cost_of(o, protocols):
    c = []
    if o[0] in ("src", "snk") and "ends" in protocols:
        c.append(COST["ends"])
    if o[0] == "link" and "link" in protocols:
        c.append(COST["link"])
    if o[0] == "paired" and "pe" in protocols:
        c.append(COST["pe"])
    if o[0] == "full" and "long" in protocols:
        c.append(COST["long"])
    if o[0] == "sub":
        L = len(o[1])
        if L <= 2 and "short2" in protocols:
            c.append(COST["short2"])
        if L <= 3 and "short3" in protocols:
            c.append(COST["short3"])
        if "long" in protocols:
            c.append(COST["long"])
    return min(c) if c else float("inf")


def must_cover(spec):
    b = {("src", s) for s, (_, ob) in spec["sources"].items() if ob} | {("snk", t) for t, (_, ob) in spec["sinks"].items() if ob}
    return sorted(set(map(tuple, spec["hard"])) | b)


# ------------------------------------------------------------------ feasibility of a support
def _quant_rows(spec, P, Z, offset, nvar):
    A, lo, hi = [], [], []
    for tau, y, er, ea in spec["quant"]:
        for o, yo in y.items():
            r = np.zeros(nvar)
            for k, j in enumerate(Z):
                r[offset + k] = yields(P[j], o)
            A.append(r); lo.append(yo - (er * yo + ea)); hi.append(yo + (er * yo + ea))
    return A, lo, hi


def _lp(A, lo, hi, nvar, A_eq=None):
    """Feasibility of lo <= A v <= hi (and A_eq v = 0), ETA <= v <= BIGM, by interior point."""
    if nvar == 0:
        ok = all(l <= TOL and h >= -TOL for l, h in zip(lo, hi))
        return ok, np.zeros(0)
    A_ub, b_ub = [], []
    for r, l, h in zip(A, lo, hi):
        A_ub.append(r); b_ub.append(h)
        A_ub.append(-r); b_ub.append(-l)
    res = linprog(np.zeros(nvar), A_ub=np.array(A_ub) if A_ub else None, b_ub=np.array(b_ub) if b_ub else None,
                  A_eq=np.array(A_eq) if A_eq else None, b_eq=np.zeros(len(A_eq)) if A_eq else None,
                  bounds=[(ETA, BIGM)] * nvar, method="highs-ipm")
    if res.status not in (0, 2):                       # interior point did not conclude: retry with dual simplex
        res = linprog(np.zeros(nvar), A_ub=np.array(A_ub) if A_ub else None, b_ub=np.array(b_ub) if b_ub else None,
                      A_eq=np.array(A_eq) if A_eq else None, b_eq=np.zeros(len(A_eq)) if A_eq else None,
                      bounds=[(ETA, BIGM)] * nvar, method="highs-ds")
        ORACLE_LOG["lp_retried"] += 1
        if res.status not in (0, 2):
            ORACLE_LOG["lp_unresolved"] += 1
    return res.status == 0, (res.x if res.status == 0 else None)


def support_ok(spec, P, Z):
    Zp = [P[j] for j in Z]
    for o in must_cover(spec):
        if not any(yields(p, o) for p in Zp):
            return False
    A, lo, hi = _quant_rows(spec, P, list(Z), 0, len(Z))
    return _lp(A, lo, hi, len(Z))[0]


def family(spec, P):
    return [frozenset(Z) for r in range(0, len(P) + 1) for Z in itertools.combinations(range(len(P)), r) if support_ok(spec, P, Z)]


def states(fam, n, K):
    F = [Z for Z in fam if K is None or len(Z) <= K]
    if not F:
        return None
    return {i: ("NECESSARY" if all(i in Z for Z in F) else "INCOMPATIBLE" if not any(i in Z for Z in F) else "POSSIBLE")
            for i in range(n)}


def check_witness(spec, P, z, x, K, must_in=None, must_out=None):
    """Directly verify a returned (z, x): support from z, x zero off support, ETA..BIGM on it, coverage, quantitative
    residuals, |Z| <= K, and the target path in / out.  -> list of violated items (empty = valid)."""
    bad = []
    z, x = np.asarray(z), np.asarray(x, dtype=float)
    Z = [j for j in range(len(P)) if z[j] == 1]
    if any(v not in (0, 1) for v in z):
        bad.append("z not binary")
    for j in range(len(P)):
        if j in Z and not (ETA - 1e-7 <= x[j] <= BIGM + 1e-7):
            bad.append(f"x[{j}] out of [eta, M]")
        if j not in Z and abs(x[j]) > 1e-9:
            bad.append(f"x[{j}] nonzero off support")
    for o in must_cover(spec):
        if not any(yields(P[j], o) for j in Z):
            bad.append(f"uncovered {o}")
    for tau, y, er, ea in spec["quant"]:
        for o, yo in y.items():
            pred = sum(yields(P[j], o) * x[j] for j in range(len(P)))
            if abs(pred - yo) > er * yo + ea + 1e-7:
                bad.append(f"quant residual {o}")
    if K is not None and len(Z) > K:
        bad.append("support > K")
    if must_in is not None and must_in not in Z:
        bad.append("target missing")
    if must_out is not None and must_out in Z:
        bad.append("target present")
    return bad


# ------------------------------------------------------------------ diagnosis oracles
def optimal_cover(pool, protocols, cex_paths, max_size=7):
    """Minimum total cost of a subset of pool hitting every counterexample, by subset enumeration over distinct hit
    masks.  -> (cost, chosen) | ('impossible', None) | ('unresolved', None)."""
    if not cex_paths:
        return 0.0, []
    best_for_mask = {}
    for o in pool:
        mask = 0
        for k, Z in enumerate(cex_paths):
            if not any(yields(p, o) for p in Z):
                mask |= 1 << k
        c = cost_of(o, protocols)
        if mask and (mask not in best_for_mask or c < best_for_mask[mask][0]):
            best_for_mask[mask] = (c, o)
    full = (1 << len(cex_paths)) - 1
    union = 0
    for m in best_for_mask:
        union |= m
    if union != full:
        return "impossible", None
    items = sorted(best_for_mask.items(), key=lambda kv: kv[1][0])
    # drop dominated masks (subset of another mask with cost <=)
    keep = [(m, c, o) for m, (c, o) in items if not any(m2 != m and (m | m2) == m2 and c2 <= c for m2, (c2, _) in items)]
    # exact branch and bound: take the lowest uncovered counterexample, branch on every observation that hits it
    best = [float("inf"), None]
    calls = [0]

    def rec(covered, cost, chosen):
        calls[0] += 1
        if calls[0] > 2_000_000:
            raise RuntimeError
        if cost >= best[0] - 1e-12:
            return
        if covered == full:
            best[0], best[1] = cost, list(chosen)
            return
        k = 0
        while covered >> k & 1:
            k += 1
        for m, c, o in keep:
            if m >> k & 1:
                chosen.append(o)
                rec(covered | m, cost + c, chosen)
                chosen.pop()
    try:
        rec(0, 0.0, [])
    except RuntimeError:
        return "unresolved", None
    return best[0], best[1]


def structural(P, F_minus, i, tau):
    """Some tau observation of path i excludes each counterexample (so a tau-only sufficient set exists)."""
    pool = tau_obs(P[i], tau)
    return all(any(not any(yields(P[j], o) for j in Z) for o in pool) for Z in F_minus)


def info_limited(spec, P, fam, i, K, tau, cap=3000):
    """Exists Z1 (with i), Z2 (without i), both feasible and |Z| <= K, with x1, x2 meeting the current quantitative
    evidence and B_tau x1 = B_tau x2?  -> True / False / 'SKIPPED' (too many pairs)."""
    F = [Z for Z in fam if K is None or len(Z) <= K]
    with_i = [sorted(Z) for Z in F if i in Z]
    without = [sorted(Z) for Z in F if i not in Z]
    if len(with_i) * len(without) > cap:
        return "SKIPPED"
    space = sorted({o for p in P for o in tau_obs(p, tau)})
    for Z1 in with_i:
        for Z2 in without:
            nv = len(Z1) + len(Z2)
            A, lo, hi = _quant_rows(spec, P, Z1, 0, nv)
            A2, lo2, hi2 = _quant_rows(spec, P, Z2, len(Z1), nv)
            A_eq = []
            for o in space:
                r = np.zeros(nv)
                for k, j in enumerate(Z1):
                    r[k] = yields(P[j], o)
                for k, j in enumerate(Z2):
                    r[len(Z1) + k] = -yields(P[j], o)
                if r.any():
                    A_eq.append(r)
            if _lp(A + A2, lo + lo2, hi + hi2, nv, A_eq)[0]:
                return True
    return False
