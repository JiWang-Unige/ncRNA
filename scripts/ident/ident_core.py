#!/usr/bin/env python3
"""IDENT core (docs/ident/P0a/DEFS_P0a.md v1.1): evidence-constrained transcript identification on one locus.

Objects: a DAG of exon segments; virtual sources (TSS clusters) and sinks (TES clusters), each 'observed' or 'open';
candidate transcripts = all source->sink paths.  Observations:
  ('sub', nodes)                    contiguous subpath (short read, truncated long read)
  ('paired', m1, m2, gmin, gmax)    read pair: both mates in the same transcript, m2 after m1, gmin..gmax nodes between
  ('full', s, nodes, t)             full-length molecule: exactly one transcript
  ('link', s, t)                    end-to-end link (TIF-seq-like)
  ('src', s) / ('snk', t)           end presence (CAGE / 3'-end seq)
Hard observations must be produced by a selected path.  Optional quantitative evidence per protocol:
|B_tau x - y| <= eps_rel*y + eps_abs.  Supports: z_p in {0,1}, eta*z_p <= x_p <= M*z_p; parsimony sum z <= K = k* + delta.

Method = MILP (scipy HiGHS branch-and-bound).  The P0b oracle lives in oracle_indep.py and does not import this file.
"""
import math
import time

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, linprog, milp

PROTOCOLS = ("ends", "short2", "short3", "pe", "long", "link")
DEFAULT_COST = {"ends": 0.5, "short2": 1.0, "short3": 1.5, "pe": 1.5, "long": 3.0, "link": 5.0}
ETA, BIGM = 1.0, 1e4
PE_MATE_MAX, PE_GAP = 2, (0, 2)
MAX_PATHS = 5000


# ------------------------------------------------------------------ locus, paths, observations
class Locus:
    def __init__(self, n_nodes, edges, sources, sinks, max_paths=MAX_PATHS):
        """sources / sinks: {name: (node, observed_bool)}; edges: iterable of (u, v) with u < v."""
        self.n = n_nodes
        self.edges = sorted(set(edges))
        self.succ = {v: sorted(w for u, w in self.edges if u == v) for v in range(n_nodes)}
        self.sources, self.sinks = dict(sources), dict(sinks)
        self.max_paths = max_paths
        self.too_many = False
        self.paths = self._enumerate()

    def _enumerate(self):
        sink_at = {}
        for t, (v, _) in self.sinks.items():
            sink_at.setdefault(v, []).append(t)
        reach = set(sink_at)                                   # nodes that can reach some sink
        changed = True
        while changed:
            changed = False
            for u, v in self.edges:
                if v in reach and u not in reach:
                    reach.add(u); changed = True
        out = set()
        for s, (v0, _) in sorted(self.sources.items()):
            if v0 not in reach:
                continue
            stack = [(v0, (v0,))]
            while stack:
                v, nodes = stack.pop()
                for t in sink_at.get(v, []):
                    out.add((s, nodes, t))
                    if len(out) > self.max_paths:
                        self.too_many = True
                        return []
                for w in self.succ[v]:
                    if w in reach:
                        stack.append((w, nodes + (w,)))
        return sorted(out)


def _occurrences(nodes, sub):
    L = len(sub)
    return [i for i in range(len(nodes) - L + 1) if nodes[i:i + L] == sub]


def produces(path, o):
    """Number of ways path produces observation o."""
    s, nodes, t = path
    k = o[0]
    if k == "src":
        return int(s == o[1])
    if k == "snk":
        return int(t == o[1])
    if k == "link":
        return int(s == o[1] and t == o[2])
    if k == "full":
        return int(path == (o[1], o[2], o[3]))
    if k == "sub":
        return len(_occurrences(nodes, o[1]))
    if k == "paired":
        m1, m2, gmin, gmax = o[1], o[2], o[3], o[4]
        return sum(1 for i in _occurrences(nodes, m1) for j in _occurrences(nodes, m2)
                   if gmin <= j - (i + len(m1)) <= gmax)
    raise ValueError(o)


def _pe_obs(nodes):
    out = set()
    n = len(nodes)
    for i in range(n):
        for a in range(1, PE_MATE_MAX + 1):
            if i + a > n:
                break
            for g in range(PE_GAP[0], PE_GAP[1] + 1):
                j = i + a + g
                for b in range(1, PE_MATE_MAX + 1):
                    if j + b <= n:
                        out.add(("paired", nodes[i:i + a], nodes[j:j + b], PE_GAP[0], PE_GAP[1]))
    return out


def protocol_obs(path, tau, locus=None):
    """Observations of protocol tau that a path produces (the candidate pool for diagnosis).  A full-length observation
    is offered only for paths whose source and sink are both observed boundaries (when the locus is given)."""
    s, nodes, t = path
    if tau == "ends":
        return [("src", s), ("snk", t)]
    if tau == "link":
        return [("link", s, t)]
    if tau == "pe":
        return sorted(_pe_obs(nodes))
    subs = [("sub", nodes[i:j]) for i in range(len(nodes)) for j in range(i + 1, len(nodes) + 1)]
    if tau == "short2":
        return [o for o in subs if len(o[1]) <= 2]
    if tau == "short3":
        return [o for o in subs if len(o[1]) <= 3]
    if tau == "long":
        full_ok = locus is None or (locus.sources[s][1] and locus.sinks[t][1])
        return subs + ([("full", s, nodes, t)] if full_ok else [])
    raise ValueError(tau)


def protocol_space(locus, tau):
    return sorted({o for p in locus.paths for o in protocol_obs(p, tau, locus)})


class Evidence:
    """hard: list of observations; quant: list of (tau, {obs: y}, eps_rel, eps_abs)."""

    def __init__(self, locus, hard=(), quant=()):
        self.locus = locus
        self.hard = sorted(set(hard))
        bnd = [("src", s) for s, (_, ob) in locus.sources.items() if ob] + [("snk", t) for t, (_, ob) in locus.sinks.items() if ob]
        self.boundary = sorted(set(bnd))
        self.quant = list(quant)

    def coverage_obs(self):
        return sorted(set(self.hard) | set(self.boundary))


# ------------------------------------------------------------------ MILP method
SOLVER_LOG = {"infeasible_rechecked": 0, "presolve_false_infeasible": 0, "presolve_recheck_unresolved": 0,
              "presolve_confirmed_infeasible": 0, "witness_check_failed": 0, "both_infeasible": 0}


def _milp(c, cons, bnd, integ, time_limit):
    """HiGHS MILP; every 'infeasible' is re-checked with presolve off (P0b-fix found presolve declaring feasible
    models infeasible).  If the re-check finds a solution, that solution is used (it is then polished and verified)."""
    r = milp(c, constraints=cons, bounds=bnd, integrality=integ, options={"time_limit": time_limit})
    if r.status == 2:
        SOLVER_LOG["infeasible_rechecked"] += 1
        r2 = milp(c, constraints=cons, bounds=bnd, integrality=integ, options={"time_limit": time_limit, "presolve": False})
        if r2.status == 0 and r2.x is not None:
            SOLVER_LOG["presolve_false_infeasible"] += 1     # counted as false only once a witness exists (verified later)
            return r2
        if r2.status != 2:
            SOLVER_LOG["presolve_recheck_unresolved"] += 1
            return r2
        SOLVER_LOG["presolve_confirmed_infeasible"] += 1
    return r


class Model:
    def __init__(self, ev, use_quant=True, time_limit=30.0):
        self.ev, self.P = ev, ev.locus.paths
        self.n = len(self.P)
        self.use_quant = use_quant and bool(ev.quant)
        self.time_limit = time_limit
        self.cov = ev.coverage_obs()

    def quant_rows(self, offset, width):
        """Rows |B x - y| <= tol on x placed at columns offset .. offset+n of a width-wide vector."""
        A, lo, hi = [], [], []
        if not self.use_quant:
            return A, lo, hi
        for tau, y, er, ea in self.ev.quant:
            for o, yo in y.items():
                r = np.zeros(width)
                for i, p in enumerate(self.P):
                    r[offset + i] = produces(p, o)
                tol = er * yo + ea
                A.append(r); lo.append(yo - tol); hi.append(yo + tol)
        return A, lo, hi

    def struct_rows(self, zoff, xoff, width, K=None, extra_cov=()):
        """eta z <= x <= M z, coverage, parsimony; None if some observation no candidate can produce."""
        n = self.n
        A, lo, hi = [], [], []
        for i in range(n):
            r = np.zeros(width); r[zoff + i] = -ETA; r[xoff + i] = 1.0; A.append(r); lo.append(0.0); hi.append(np.inf)
            r = np.zeros(width); r[zoff + i] = -BIGM; r[xoff + i] = 1.0; A.append(r); lo.append(-np.inf); hi.append(0.0)
        for o in list(self.cov) + list(extra_cov):
            r = np.zeros(width)
            for i, p in enumerate(self.P):
                if produces(p, o):
                    r[zoff + i] = 1.0
            if not r.any():
                return None
            A.append(r); lo.append(1.0); hi.append(np.inf)
        if K is not None:
            r = np.zeros(width); r[zoff:zoff + n] = 1.0; A.append(r); lo.append(-np.inf); hi.append(K)
        return A, lo, hi

    def _base(self, K=None, fix=None, extra_cov=(), objective="none"):
        n = self.n
        sr = self.struct_rows(0, n, 2 * n, K, extra_cov)
        if sr is None:
            return None
        A, lo, hi = sr
        qa, ql, qh = self.quant_rows(n, 2 * n)
        A, lo, hi = A + qa, lo + ql, hi + qh
        lb, ub = np.zeros(2 * n), np.r_[np.ones(n), np.full(n, BIGM)]
        for i, v in (fix or {}).items():
            lb[i] = ub[i] = v
        c = np.r_[np.ones(n), np.zeros(n)] if objective == "minz" else np.zeros(2 * n)
        if not A:
            A, lo, hi = [np.zeros(2 * n)], [0.0], [0.0]
        return c, LinearConstraint(np.array(A), lo, hi), Bounds(lb, ub)

    def solve(self, **kw):
        """-> ('feasible', z, x) | ('infeasible', None, None) | ('unresolved', None, None)."""
        if self.ev.locus.too_many:
            return "unresolved", None, None
        b = self._base(**kw)
        if b is None:
            return "infeasible", None, None
        c, cons, bnd = b
        integ = np.r_[np.ones(self.n), np.zeros(self.n)]
        r = _milp(c, cons, bnd, integ, self.time_limit)
        if r.status == 0 and r.x is not None:
            z = np.round(r.x[:self.n]).astype(int)
            x = self._polish(z)
            if x is None or not self.verify(z, x, **{k: kw.get(k) for k in ("K", "fix", "extra_cov")}):
                SOLVER_LOG["witness_check_failed"] += 1
                return "unresolved", None, None
            return "feasible", z, x
        if r.status == 2:
            return "infeasible", None, None
        return "unresolved", None, None

    def verify(self, z, x, K=None, fix=None, extra_cov=None):
        """Direct check of a returned witness against every constraint of the model it came from."""
        Z = [j for j in range(self.n) if z[j] == 1]
        if any(v not in (0, 1) for v in z):
            return False
        for i, v in (fix or {}).items():
            if z[i] != v:
                return False
        if K is not None and len(Z) > K:
            return False
        for j in range(self.n):
            if (j in Z and not (ETA - 1e-7 <= x[j] <= BIGM + 1e-7)) or (j not in Z and abs(x[j]) > 1e-9):
                return False
        for o in list(self.cov) + list(extra_cov or ()):
            if not any(produces(self.P[j], o) for j in Z):
                return False
        if self.use_quant:
            for tau, y, er, ea in self.ev.quant:
                for o, yo in y.items():
                    if abs(sum(produces(self.P[j], o) * x[j] for j in Z) - yo) > er * yo + ea + 1e-7:
                        return False
        return True

    def _polish(self, z):
        """Exact abundances on the chosen support: x_j = 0 off support, eta <= x_j <= M on it, quantitative rows by LP.
        Removes MILP integrality / big-M residues (x ~ 1e-6 on unselected paths)."""
        Z = list(np.nonzero(z)[0])
        x = np.zeros(self.n)
        rows = []
        if self.use_quant:
            for tau, y, er, ea in self.ev.quant:
                for o, yo in y.items():
                    rows.append(([float(produces(self.P[j], o)) for j in Z], yo, er * yo + ea))
        if not Z:
            return x if all(abs(yo) <= tol for _, yo, tol in rows) else None
        if not rows:
            x[Z] = ETA
            return x
        A_ub = [a for a, _, _ in rows] + [[-v for v in a] for a, _, _ in rows]
        b_ub = [yo + tol for _, yo, tol in rows] + [-(yo - tol) for _, yo, tol in rows]
        r = linprog(np.zeros(len(Z)), A_ub=np.array(A_ub), b_ub=np.array(b_ub), bounds=[(ETA, BIGM)] * len(Z), method="highs")
        if r.status != 0:
            return None
        x[Z] = r.x
        return x

    def kstar(self, fix=None):
        st, z, x = self.solve(fix=fix, objective="minz")
        return (int(z.sum()), z, x) if st == "feasible" else (st, None, None)


def three_state(ev, delta=0, use_quant=True, time_limit=30.0, deadline=None):
    """Per candidate path: NECESSARY / POSSIBLE / INCOMPATIBLE (+ witnesses), or model-level INFEASIBLE_MODEL / UNRESOLVED.
    deadline: absolute time.time() budget for the whole locus."""
    if ev.locus.too_many:
        return {"model": "UNRESOLVED"}
    if not ev.locus.paths:
        return {"model": "INFEASIBLE_MODEL" if ev.coverage_obs() or any(y for _, y, _, _ in ev.quant) else "OK",
                "kstar": 0, "K": None if delta is None else delta, "paths": {}}
    m = Model(ev, use_quant, time_limit)
    ks = m.kstar()
    if ks[0] == "infeasible":
        return {"model": "INFEASIBLE_MODEL"}
    if ks[0] == "unresolved":
        return {"model": "UNRESOLVED"}
    k = ks[0]
    K = None if delta is None else k + delta
    out = {"model": "OK", "kstar": k, "K": K, "paths": {}}
    for i in range(m.n):
        if deadline is not None and time.time() > deadline:
            return {"model": "UNRESOLVED"}
        w1 = m.solve(K=K, fix={i: 1})
        w0 = m.solve(K=K, fix={i: 0})
        if "unresolved" in (w1[0], w0[0]):
            st = "UNRESOLVED"
        elif w0[0] == "infeasible" and w1[0] == "infeasible":   # contradicts k*: solver inconsistency
            SOLVER_LOG["both_infeasible"] += 1
            st = "UNRESOLVED"
        elif w0[0] == "infeasible":
            st = "NECESSARY"
        elif w1[0] == "infeasible":
            st = "INCOMPATIBLE"
        else:
            st = "POSSIBLE"
        out["paths"][i] = {"state": st, "with": w1, "without": w0}
    return out


def certificate_source(ev, i, delta=0, time_limit=30.0):
    """For a path NECESSARY at delta: 'evidence' / 'quantification' / 'parsimony' / 'UNRESOLVED'; and delta_break
    (int, math.inf if no solution avoids the path, or 'UNRESOLVED')."""
    def nec(use_quant, K):
        st = Model(ev, use_quant, time_limit).solve(K=K, fix={i: 0})[0]
        return None if st == "unresolved" else st == "infeasible"
    ks = Model(ev, True, time_limit).kstar()
    if not isinstance(ks[0], int):
        return "UNRESOLVED", "UNRESOLVED"
    k = ks[0]
    src = None
    for use_quant, K, name in ((False, None, "evidence"), (True, None, "quantification"), (True, k + delta, "parsimony")):
        if name == "quantification" and not ev.quant:
            continue
        r = nec(use_quant, K)
        if r is None:
            src = "UNRESOLVED"
            break
        if r:
            src = name
            break
    k0 = Model(ev, True, time_limit).kstar(fix={i: 0})[0]
    brk = (k0 - k) if isinstance(k0, int) else (math.inf if k0 == "infeasible" else "UNRESOLVED")
    return src, brk


# ------------------------------------------------------------------ sufficient-observation diagnosis (CEGIS)
def hits(paths_in_Z, o):
    return not any(produces(p, o) for p in paths_in_Z)


def min_cost_cover(pool, cost, counterexamples, time_limit=30.0):
    """Minimum-cost subset of pool hitting every counterexample.  -> ('ok', obs) | ('impossible', None) | ('unresolved', None)."""
    if not counterexamples:
        return "ok", []
    H = np.array([[1.0 if hits(Z, o) else 0.0 for o in pool] for Z in counterexamples])
    if (H.sum(1) == 0).any():
        return "impossible", None
    c = np.array([cost[o] for o in pool])
    r = _milp(c, LinearConstraint(H, 1.0, np.inf), Bounds(0, 1), np.ones(len(pool)), time_limit)
    if r.status == 0 and r.x is not None:
        return "ok", [o for o, u in zip(pool, r.x) if u > 0.5]
    return ("impossible", None) if r.status == 2 else ("unresolved", None)


def obs_cost(o, protocols, cost_table=None):
    """Cheapest protocol (among those allowed) that can produce o, under cost_table (default DEFAULT_COST)."""
    cost_table = cost_table or DEFAULT_COST
    best = math.inf
    for tau in protocols:
        if (o[0] == "link" and tau == "link") or (o[0] in ("src", "snk") and tau == "ends") \
                or (o[0] == "paired" and tau == "pe") or (o[0] == "full" and tau == "long"):
            best = min(best, cost_table[tau])
        elif o[0] == "sub" and tau in ("short2", "short3", "long"):
            L = len(o[1])
            if tau == "long" or (tau == "short2" and L <= 2) or (tau == "short3" and L <= 3):
                best = min(best, cost_table[tau])
    return best


def diagnosis_pool(path, protocols, locus=None, cost_table=None):
    pool = sorted({o for tau in protocols for o in protocol_obs(path, tau, locus)})
    return pool, {o: obs_cost(o, protocols, cost_table) for o in pool}


def diagnose(ev, i, K, protocols=PROTOCOLS, time_limit=30.0, max_iter=500, cost_table=None, use_quant=True, members=None):
    """Minimum-cost conditionally sufficient observation set for path i (or for the class `members`, all absent in a
    counterexample; pool = observations produced by every member), K fixed.  -> dict."""
    m = Model(ev, use_quant, time_limit)
    members = list(members) if members is not None else [i]
    pools = [set(diagnosis_pool(m.P[j], protocols, ev.locus, cost_table)[0]) for j in members]
    pool = sorted(set.intersection(*pools))
    cost = {o: obs_cost(o, protocols, cost_table) for o in pool}
    cex, O = [], []
    for it in range(max_iter):
        st, z, x = m.solve(K=K, fix={j: 0 for j in members}, extra_cov=O)
        if st == "unresolved":
            return {"status": "UNRESOLVED"}
        if st == "infeasible":
            return {"status": "SUFFICIENT", "obs": O, "cost": float(sum(cost[o] for o in O)), "n_counterexamples": len(cex), "iterations": it}
        Z = [m.P[j] for j in np.nonzero(z)[0]]
        cex.append(Z)
        cs, O = min_cost_cover(pool, cost, cex, time_limit)
        if cs == "unresolved":
            return {"status": "UNRESOLVED"}
        if cs == "impossible":
            return {"status": "INSUFFICIENT_POOL", "counterexample": Z, "n_counterexamples": len(cex)}
    return {"status": "UNRESOLVED"}


def info_limited(ev, i, K, tau, time_limit=30.0):
    """In the ORIGINAL feasible region (hard, boundaries, current quantitative evidence, eta, K): two solutions, one with
    path i and one without, with identical expected tau observations B_tau x1 = B_tau x2 (same abundance scale).
    -> witness dict | None (no such pair) | 'UNRESOLVED'."""
    m = Model(ev, True, time_limit)
    n = m.n
    space = protocol_space(ev.locus, tau)
    B = np.array([[float(produces(p, o)) for p in m.P] for o in space])
    N = 4 * n                                                   # z1, x1, z2, x2
    A, lo, hi = [], [], []
    for zoff, xoff in ((0, n), (2 * n, 3 * n)):
        sr = m.struct_rows(zoff, xoff, N, K)
        if sr is None:
            return None
        qa, ql, qh = m.quant_rows(xoff, N)
        A += sr[0] + qa; lo += sr[1] + ql; hi += sr[2] + qh
    for k in range(len(space)):
        r = np.zeros(N); r[n:2 * n] = B[k]; r[3 * n:] = -B[k]; A.append(r); lo.append(0.0); hi.append(0.0)
    lb = np.zeros(N)
    ub = np.r_[np.ones(n), np.full(n, BIGM), np.ones(n), np.full(n, BIGM)]
    lb[i] = ub[i] = 1
    lb[2 * n + i] = ub[2 * n + i] = 0
    integ = np.r_[np.ones(n), np.zeros(n), np.ones(n), np.zeros(n)]
    r = _milp(np.zeros(N), LinearConstraint(np.array(A), lo, hi), Bounds(lb, ub), integ, time_limit)
    if r.status != 0 or r.x is None:
        return None if r.status == 2 else "UNRESOLVED"
    z1, z2 = np.round(r.x[:n]).astype(int), np.round(r.x[2 * n:3 * n]).astype(int)
    pol = _polish_pair(m, z1, z2, B)
    if pol is None:
        return "UNRESOLVED"
    x1, x2 = pol
    return {"z1": z1, "x1": x1, "z2": z2, "x2": x2, "space": space, "B": B}


def _polish_pair(m, z1, z2, B):
    """Exact x1, x2 on fixed supports with the quantitative rows and B x1 = B x2 (removes MILP residues)."""
    n = m.n
    Z1, Z2 = list(np.nonzero(z1)[0]), list(np.nonzero(z2)[0])
    nv = len(Z1) + len(Z2)
    if nv == 0:
        return np.zeros(n), np.zeros(n)
    A_eq, b_eq, A_ub, b_ub = [], [], [], []
    for k in range(B.shape[0]):
        A_eq.append([B[k, j] for j in Z1] + [-B[k, j] for j in Z2]); b_eq.append(0.0)
    if m.use_quant:
        for blk, Zb in ((0, Z1), (1, Z2)):
            for tau, y, er, ea in m.ev.quant:
                for o, yo in y.items():
                    a = [0.0] * nv
                    for k, j in enumerate(Zb):
                        a[blk * len(Z1) + k] = float(produces(m.P[j], o))
                    tol = er * yo + ea
                    A_ub.append(a); b_ub.append(yo + tol)
                    A_ub.append([-v for v in a]); b_ub.append(-(yo - tol))
    r = linprog(np.zeros(nv), A_ub=np.array(A_ub) if A_ub else None, b_ub=np.array(b_ub) if b_ub else None,
                A_eq=np.array(A_eq), b_eq=np.array(b_eq), bounds=[(ETA, BIGM)] * nv, method="highs")
    if r.status != 0:
        return None
    x1, x2 = np.zeros(n), np.zeros(n)
    x1[Z1] = r.x[:len(Z1)]
    x2[Z2] = r.x[len(Z1):]
    return x1, x2


def ambiguity_type(ev, i, K, tau, time_limit=30.0):
    """STRUCTURAL / INFORMATION_LIMITED / QUANT_DISTINGUISHABLE / UNRESOLVED for a POSSIBLE path under protocol tau."""
    d = diagnose(ev, i, K, protocols=(tau,), time_limit=time_limit)
    if d["status"] == "SUFFICIENT":
        return "STRUCTURAL", d
    if d["status"] == "UNRESOLVED":
        return "UNRESOLVED", d
    w = info_limited(ev, i, K, tau, time_limit)
    if w == "UNRESOLVED":
        return "UNRESOLVED", None
    return ("INFORMATION_LIMITED", w) if w is not None else ("QUANT_DISTINGUISHABLE", d)


def class_state(ev, members, K, use_quant=True, time_limit=30.0):
    """Three state of an equivalence class (e.g. all paths sharing an intron chain): NECESSARY iff no feasible solution
    has every member absent; INCOMPATIBLE iff no feasible solution has some member present."""
    m = Model(ev, use_quant, time_limit)
    absent = m.solve(K=K, fix={j: 0 for j in members})[0]
    if absent == "unresolved":
        return "UNRESOLVED"
    if absent == "infeasible":
        return "NECESSARY"
    b = m._base(K=K, extra_cov=())
    if b is None:
        return "UNRESOLVED"
    c, cons, bnd = b
    n = m.n
    r = np.zeros(2 * n)
    r[list(members)] = 1.0
    A = np.vstack([cons.A, r])
    lo, hi = np.r_[cons.lb, 1.0], np.r_[cons.ub, np.inf]
    res = _milp(c, LinearConstraint(A, lo, hi), bnd, np.r_[np.ones(n), np.zeros(n)], time_limit)
    if res.status == 0:
        return "POSSIBLE"
    return "INCOMPATIBLE" if res.status == 2 else "UNRESOLVED"
