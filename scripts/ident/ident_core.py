#!/usr/bin/env python3
"""IDENT core (docs/ident/P0a/DEFS_P0a.md v1): evidence-constrained transcript identification on one locus.

Objects: a DAG of exon segments; virtual sources (TSS clusters) and sinks (TES clusters), each 'observed' or 'open';
candidate transcripts = all source->sink paths.  Observations: ('sub', nodes) contiguous subpath, ('link', s, t)
end-to-end link, ('src', s) / ('snk', t) end presence.  Hard observations must be produced by a selected path.
Optional quantitative evidence per protocol: |B_tau x - y| <= eps_rel*y + eps_abs.
Supports: z_p in {0,1}, eta*z_p <= x_p <= M*z_p; parsimony sum z <= K = k* + delta.

Method = MILP (scipy HiGHS).  Oracle (P0b only) = enumerate every support subset, LP feasibility per subset.
"""
import itertools
import math

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, linprog, milp

PROTOCOLS = ("ends", "short2", "short3", "long", "link")
DEFAULT_COST = {"ends": 0.5, "short2": 1.0, "short3": 1.5, "long": 3.0, "link": 5.0}
ETA, BIGM = 1.0, 1e4


# ------------------------------------------------------------------ locus, paths, observations
class Locus:
    def __init__(self, n_nodes, edges, sources, sinks):
        """sources / sinks: {name: (node, observed_bool)}; edges: iterable of (u, v) with u < v."""
        self.n = n_nodes
        self.succ = {v: sorted(w for u, w in edges if u == v) for v in range(n_nodes)}
        self.sources, self.sinks = dict(sources), dict(sinks)
        self.paths = self._enumerate()

    def _enumerate(self):
        sink_at = {}
        for t, (v, _) in self.sinks.items():
            sink_at.setdefault(v, []).append(t)
        out = []
        for s, (v0, _) in sorted(self.sources.items()):
            stack = [(v0, (v0,))]
            while stack:
                v, nodes = stack.pop()
                for t in sink_at.get(v, []):
                    out.append((s, nodes, t))
                for w in self.succ[v]:
                    stack.append((w, nodes + (w,)))
        return sorted(set(out))


def produces(path, o):
    """Number of ways path produces observation o (0/1 in a DAG)."""
    s, nodes, t = path
    k = o[0]
    if k == "src":
        return int(s == o[1])
    if k == "snk":
        return int(t == o[1])
    if k == "link":
        return int(s == o[1] and t == o[2])
    sub = o[1]
    L = len(sub)
    return sum(1 for i in range(len(nodes) - L + 1) if nodes[i:i + L] == sub)


def protocol_obs(path, tau):
    """Observations of protocol tau that a path produces (the candidate pool for diagnosis)."""
    s, nodes, t = path
    if tau == "ends":
        return [("src", s), ("snk", t)]
    if tau == "link":
        return [("link", s, t)]
    subs = [("sub", nodes[i:j]) for i in range(len(nodes)) for j in range(i + 1, len(nodes) + 1)]
    if tau == "short2":
        return [o for o in subs if len(o[1]) <= 2]
    if tau == "short3":
        return [o for o in subs if len(o[1]) <= 3]
    if tau == "long":
        return subs
    raise ValueError(tau)


def protocol_space(locus, tau):
    return sorted({o for p in locus.paths for o in protocol_obs(p, tau)})


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
class Model:
    def __init__(self, ev, use_quant=True, time_limit=30.0):
        self.ev, self.P = ev, ev.locus.paths
        self.n = len(self.P)
        self.use_quant = use_quant and bool(ev.quant)
        self.time_limit = time_limit
        self.cov = ev.coverage_obs()

    def _base(self, K=None, fix=None, extra_cov=(), objective="none"):
        n = self.n
        A, lo, hi = [], [], []

        def row(zc=None, xc=None):
            r = np.zeros(2 * n)
            if zc:
                for i, v in zc.items():
                    r[i] += v
            if xc:
                for i, v in xc.items():
                    r[n + i] += v
            return r
        for i in range(n):                                    # eta z <= x <= M z
            A.append(row({i: -ETA}, {i: 1.0})); lo.append(0.0); hi.append(np.inf)
            A.append(row({i: -BIGM}, {i: 1.0})); lo.append(-np.inf); hi.append(0.0)
        for o in list(self.cov) + list(extra_cov):          # coverage
            prod = {i: 1.0 for i, p in enumerate(self.P) if produces(p, o)}
            if not prod:
                return None                                   # observation no candidate can produce
            A.append(row(prod)); lo.append(1.0); hi.append(np.inf)
        if self.use_quant:
            for tau, y, er, ea in self.ev.quant:
                for o, yo in y.items():
                    xc = {i: float(produces(p, o)) for i, p in enumerate(self.P) if produces(p, o)}
                    tol = er * yo + ea
                    A.append(row(None, xc)); lo.append(yo - tol); hi.append(yo + tol)
        if K is not None:
            A.append(row({i: 1.0 for i in range(n)})); lo.append(-np.inf); hi.append(K)
        lb, ub = np.zeros(2 * n), np.r_[np.ones(n), np.full(n, BIGM)]
        for i, v in (fix or {}).items():
            lb[i] = ub[i] = v
        c = np.r_[np.ones(n), np.zeros(n)] if objective == "minz" else np.zeros(2 * n)
        return c, LinearConstraint(np.array(A), lo, hi), Bounds(lb, ub)

    def solve(self, **kw):
        """-> ('feasible', z, x) | ('infeasible', None, None) | ('unresolved', None, None)."""
        b = self._base(**kw)
        if b is None:
            return "infeasible", None, None
        c, cons, bnd = b
        integ = np.r_[np.ones(self.n), np.zeros(self.n)]
        r = milp(c, constraints=cons, bounds=bnd, integrality=integ, options={"time_limit": self.time_limit})
        if r.status == 0 and r.x is not None:
            z = np.round(r.x[:self.n]).astype(int)
            x = self._polish(z)
            if x is None:                                     # MILP tolerance accepted a support the exact LP rejects
                return "unresolved", None, None
            return "feasible", z, x

        if r.status == 2:
            return "infeasible", None, None
        return "unresolved", None, None

    def _polish(self, z):
        """Exact abundances on the chosen support: x_j = 0 off support, eta <= x_j <= M on it, quantitative rows by LP.
        Removes MILP integrality / big-M residues (x ~ 1e-6 on unselected paths)."""
        Z = list(np.nonzero(z)[0])
        x = np.zeros(self.n)
        if not Z:
            return x
        if not self.use_quant:
            x[Z] = ETA
            return x
        A_ub, b_ub = [], []
        for tau, y, er, ea in self.ev.quant:
            for o, yo in y.items():
                a = [float(produces(self.P[j], o)) for j in Z]
                tol = er * yo + ea
                A_ub.append(a); b_ub.append(yo + tol)
                A_ub.append([-v for v in a]); b_ub.append(-(yo - tol))
        r = linprog(np.zeros(len(Z)), A_ub=np.array(A_ub), b_ub=np.array(b_ub), bounds=[(ETA, BIGM)] * len(Z), method="highs")
        if r.status != 0:
            return None
        x[Z] = r.x
        return x

    def kstar(self, fix=None):
        st, z, x = self.solve(fix=fix, objective="minz")
        return (int(z.sum()), z, x) if st == "feasible" else (st, None, None)


def three_state(ev, delta=0, use_quant=True, time_limit=30.0):
    """Per candidate path: NECESSARY / POSSIBLE / INCOMPATIBLE (+ witnesses), or model-level INFEASIBLE_MODEL / UNRESOLVED."""
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
        w1 = m.solve(K=K, fix={i: 1})
        w0 = m.solve(K=K, fix={i: 0})
        if "unresolved" in (w1[0], w0[0]):
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
    """For a path NECESSARY at delta: evidence-forced / quantification-forced / parsimony-dependent; and delta_break."""
    def nec(use_quant, K):
        m = Model(ev, use_quant, time_limit)
        return m.solve(K=K, fix={i: 0})[0] == "infeasible"
    if nec(False, None):
        src = "evidence"
    elif ev.quant and nec(True, None):
        src = "quantification"
    elif nec(True, Model(ev, True, time_limit).kstar()[0] + delta):
        src = "parsimony"
    else:
        src = None
    m = Model(ev, True, time_limit)
    k = m.kstar()[0]
    k0 = m.kstar(fix={i: 0})[0]
    brk = (k0 - k) if isinstance(k0, int) else math.inf
    return src, brk


# ------------------------------------------------------------------ sufficient-observation diagnosis (CEGIS)
def hits(paths_in_Z, o):
    return not any(produces(p, o) for p in paths_in_Z)


def min_cost_cover(pool, cost, counterexamples):
    """Minimum-cost subset of pool hitting every counterexample (list of path lists); None if impossible."""
    if not counterexamples:
        return []
    H = np.array([[1.0 if hits(Z, o) else 0.0 for o in pool] for Z in counterexamples])
    if (H.sum(1) == 0).any():
        return None
    c = np.array([cost[o] for o in pool])
    r = milp(c, constraints=LinearConstraint(H, 1.0, np.inf), bounds=Bounds(0, 1), integrality=np.ones(len(pool)))
    if r.status != 0:
        return None
    return [o for o, u in zip(pool, r.x) if u > 0.5]


def obs_cost(o, protocols):
    """Cheapest protocol (among those allowed) that can produce o."""
    best = math.inf
    for tau in protocols:
        if o[0] == "link" and tau == "link":
            best = min(best, DEFAULT_COST[tau])
        elif o[0] in ("src", "snk") and tau == "ends":
            best = min(best, DEFAULT_COST[tau])
        elif o[0] == "sub" and tau in ("short2", "short3", "long"):
            L = len(o[1])
            if tau == "long" or (tau == "short2" and L <= 2) or (tau == "short3" and L <= 3):
                best = min(best, DEFAULT_COST[tau])
    return best


def diagnose(ev, i, K, protocols=PROTOCOLS, time_limit=30.0, max_iter=500):
    """Minimum-cost conditionally sufficient observation set for path i (K fixed).  -> dict."""
    m = Model(ev, True, time_limit)
    p = m.P[i]
    pool = sorted({o for tau in protocols for o in protocol_obs(p, tau)})
    cost = {o: obs_cost(o, protocols) for o in pool}
    cex, O = [], []
    for it in range(max_iter):
        st, z, x = m.solve(K=K, fix={i: 0}, extra_cov=O)
        if st == "unresolved":
            return {"status": "UNRESOLVED"}
        if st == "infeasible":
            return {"status": "SUFFICIENT", "obs": O, "cost": float(sum(cost[o] for o in O)), "n_counterexamples": len(cex), "iterations": it}
        Z = [m.P[j] for j in np.nonzero(z)[0]]
        cex.append(Z)
        O = min_cost_cover(pool, cost, cex)
        if O is None:
            return {"status": "INSUFFICIENT_POOL", "counterexample": Z, "n_counterexamples": len(cex)}
    return {"status": "UNRESOLVED"}


def info_limited(ev, i, K, tau, time_limit=30.0):
    """Two feasible supports (with / without path i), each normalised to total tau-mass 1, with identical tau
    observation vectors?  Current quantitative evidence is not used (scale-free comparison).  -> witness or None."""
    loc = ev.locus
    P, n = loc.paths, len(loc.paths)
    space = protocol_space(loc, tau)
    B = np.array([[float(produces(p, o)) for p in P] for o in space])
    cov = ev.coverage_obs()
    eta = 1e-3
    A, lo, hi = [], [], []
    N = 4 * n                                                   # z1, x1, z2, x2

    def blank():
        return np.zeros(N)
    for blk in (0, 2):
        zo, xo = blk * n, (blk + 1) * n
        for j in range(n):
            r = blank(); r[xo + j] = 1; r[zo + j] = -eta; A.append(r); lo.append(0); hi.append(np.inf)
            r = blank(); r[xo + j] = 1; r[zo + j] = -1.0; A.append(r); lo.append(-np.inf); hi.append(0)
        for o in cov:
            r = blank()
            for j, p in enumerate(P):
                if produces(p, o):
                    r[zo + j] = 1
            if not r.any():
                return None
            A.append(r); lo.append(1); hi.append(np.inf)
        if K is not None:
            r = blank(); r[zo:zo + n] = 1; A.append(r); lo.append(-np.inf); hi.append(K)
        r = blank(); r[xo:xo + n] = B.sum(0); A.append(r); lo.append(1); hi.append(1)     # total tau mass = 1
    for k in range(len(space)):                                                         # identical tau vectors
        r = blank(); r[n:2 * n] = B[k]; r[3 * n:4 * n] = -B[k]; A.append(r); lo.append(0); hi.append(0)
    lb, ub = np.zeros(N), np.r_[np.ones(n), np.ones(n), np.ones(n), np.ones(n)]
    lb[i] = ub[i] = 1                 # path i in solution 1
    lb[2 * n + i] = ub[2 * n + i] = 0  # path i not in solution 2
    integ = np.r_[np.ones(n), np.zeros(n), np.ones(n), np.zeros(n)]
    r = milp(np.zeros(N), constraints=LinearConstraint(np.array(A), lo, hi), bounds=Bounds(lb, ub), integrality=integ,
             options={"time_limit": time_limit})
    if r.status != 0 or r.x is None:
        return None if r.status == 2 else "UNRESOLVED"
    z1, x1 = np.round(r.x[:n]).astype(int), r.x[n:2 * n]
    z2, x2 = np.round(r.x[2 * n:3 * n]).astype(int), r.x[3 * n:]
    return {"z1": z1, "x1": x1, "z2": z2, "x2": x2, "space": space, "B": B}


def ambiguity_type(ev, i, K, tau, time_limit=30.0):
    d = diagnose(ev, i, K, protocols=(tau,), time_limit=time_limit)
    if d["status"] == "SUFFICIENT":
        return "STRUCTURAL", d
    if d["status"] == "UNRESOLVED":
        return "UNRESOLVED", d
    w = info_limited(ev, i, K, tau, time_limit)
    if w == "UNRESOLVED":
        return "UNRESOLVED", None
    return ("INFORMATION_LIMITED", w) if w is not None else ("SAMPLING_LIMITED", d)


# ------------------------------------------------------------------ oracle (P0b)
def support_feasible(ev, Z_idx, use_quant=True):
    """Independent check: coverage + boundaries + (LP) quantitative feasibility for a support set."""
    P = ev.locus.paths
    Zp = [P[j] for j in Z_idx]
    for o in ev.coverage_obs():
        if not any(produces(p, o) for p in Zp):
            return False, None
    if not (use_quant and ev.quant):
        return True, None
    if not Z_idx:
        return False, None
    A_ub, b_ub = [], []
    for tau, y, er, ea in ev.quant:
        for o, yo in y.items():
            a = [float(produces(P[j], o)) for j in Z_idx]
            tol = er * yo + ea
            A_ub.append(a); b_ub.append(yo + tol)
            A_ub.append([-v for v in a]); b_ub.append(-(yo - tol))
    r = linprog(np.zeros(len(Z_idx)), A_ub=np.array(A_ub), b_ub=np.array(b_ub), bounds=[(ETA, BIGM)] * len(Z_idx), method="highs")
    return (r.status == 0), (r.x if r.status == 0 else None)


def oracle_family(ev, use_quant=True):
    n = len(ev.locus.paths)
    fam = []
    for r in range(1, n + 1):
        for Z in itertools.combinations(range(n), r):
            ok, _ = support_feasible(ev, Z, use_quant)
            if ok:
                fam.append(frozenset(Z))
    return fam


def oracle_states(fam, n, K):
    F = [Z for Z in fam if K is None or len(Z) <= K]
    if not F:
        return None
    return {i: ("NECESSARY" if all(i in Z for Z in F) else "INCOMPATIBLE" if all(i not in Z for Z in F) else "POSSIBLE")
            for i in range(n)}
