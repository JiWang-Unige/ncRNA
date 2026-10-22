"""G2-PH grammar (protocol v1.7 section 14): per-base phase-type semi-Markov model, coding x lncRNA, both strands.

States: IG, and per strand s and type k
  exon roles by genomic position S (single), A (leftmost), M (middle), Z (rightmost), each a chain of NE stages
  (self-loop + advance per stage -> negative-binomial-like lengths); intron: L (first base), NBR parallel branches of
  NI stages (mixture of phase-type lengths), R (last base).
Paths: IG -> S -> IG | IG -> A -> L -> I -> R -> (M -> L -> I -> R)* -> Z -> IG ; a gene never changes k or s.
Full model: NE = 3, NBR = 2, NI = 3 (81 states); ablation G2-PH-1: NE = NBR = NI = 1 (29 states).
Emissions: the G1 base state's per-base potential (kEs / kLs / kIs / kRs / IG) + G1 motif gates.
Boundary bonuses (as G1B): gene entry IG -> S0 / A0 and exit S_last / Z_last -> IG get coef x TSS / TES logit.
Parameters are free log-potentials defined on the + strand; the - strand is the mirror image (A <-> Z, L <-> R,
stage order reversed), so log Z is invariant to reverse complementation by construction.

Inference: numba float64 sparse forward-backward over windows (parallel over a batch) returning log Z, posteriors
collapsed to the 17 G1 states, expected transition counts and expected boundary-edge use per position.  Training
uses the surrogate sum(x * dpost) + sum(T * dC) + sum(bonus * dE), whose gradient equals grad(log Z).
"""
import numpy as np
import torch
from torch import nn

from . import g1

NEG = g1.NEG
ROLES = ("S", "A", "M", "Z")
TRANS = ("s_in", "e5_in", "s_out", "e3_out", "donor", "acc_mid", "acc_last")

try:
    from numba import njit, prange
except ImportError:  # pragma: no cover
    njit = None


class Space:
    """Expanded state space for a stage configuration."""

    def __init__(self, NE=3, NBR=2, NI=3):
        self.NE, self.NBR, self.NI = NE, NBR, NI
        names, base = ["IG"], [g1.IDX["IG"]]
        for s in g1.STRANDS:
            for k in g1.TYPES:
                for r in ROLES:
                    for i in range(NE):
                        names.append((k, s, "E", r, i))
                        base.append(g1.IDX[f"{k}E{s}"])
                names.append((k, s, "L", None, 0))
                base.append(g1.IDX[f"{k}L{s}"])
                for b in range(NBR):
                    for i in range(NI):
                        names.append((k, s, "I", b, i))
                        base.append(g1.IDX[f"{k}I{s}"])
                names.append((k, s, "R", None, 0))
                base.append(g1.IDX[f"{k}R{s}"])
        self.names = names
        self.idx = {n: i for i, n in enumerate(names)}
        self.base = np.array(base, np.int64)
        self.n = len(names)
        self.mirror = np.array([self.idx[self._mirror(n)] for n in names], np.int64)

    def _mirror(self, n):
        if n == "IG":
            return n
        k, s, part, r, i = n
        s2 = "-" if s == "+" else "+"
        if part == "E":
            return (k, s2, "E", {"S": "S", "A": "Z", "M": "M", "Z": "A"}[r], self.NE - 1 - i)
        if part == "I":
            return (k, s2, "I", r, self.NI - 1 - i)
        return (k, s2, "R" if part == "L" else "L", None, 0)

    def plus_edges(self):
        """+ strand edges: list of (from, to, param key); param key = (group, k, ...) resolved by G2PHParams."""
        E, I, NE, NI = self.idx, [], self.NE, self.NI
        out = []
        for k in g1.TYPES:
            ex = lambda r, i: E[(k, "+", "E", r, i)]
            Lb, Rb = E[(k, "+", "L", None, 0)], E[(k, "+", "R", None, 0)]
            out += [(0, ex("S", 0), ("t", k, "s_in")), (0, ex("A", 0), ("t", k, "e5_in")),
                    (ex("S", NE - 1), 0, ("t", k, "s_out")), (ex("Z", NE - 1), 0, ("t", k, "e3_out")),
                    (ex("A", NE - 1), Lb, ("t", k, "donor")), (ex("M", NE - 1), Lb, ("t", k, "donor")),
                    (Rb, ex("M", 0), ("t", k, "acc_mid")), (Rb, ex("Z", 0), ("t", k, "acc_last"))]
            for ri, r in enumerate(ROLES):
                for i in range(NE):
                    out.append((ex(r, i), ex(r, i), ("es", k, ri, i)))
                    if i < NE - 1:
                        out.append((ex(r, i), ex(r, i + 1), ("ea", k, ri, i)))
            for b in range(self.NBR):
                ib = lambda i: E[(k, "+", "I", b, i)]
                out.append((Lb, ib(0), ("br", k, b)))
                for i in range(NI):
                    out.append((ib(i), ib(i), ("is", k, b, i)))
                    out.append((ib(i), ib(i + 1) if i < NI - 1 else Rb, ("ia", k, b, i)))
        return out

    def all_edges(self):
        """(from, to, key) for both strands; - strand edges are mirrors (reversed) of + strand edges."""
        pe = self.plus_edges()
        inv = self.mirror  # mirror is an involution
        return [(0, 0, ("ig",))] + pe + [(int(inv[v]), int(inv[u]), key) for u, v, key in pe]

    def boundary_edges(self):
        """the 16 position-dependent edges: (from, to, boundary channel, coef index (0 start / 1 end), type index,
        side 'entry' (bonus at t) / 'exit' (bonus at t-1))."""
        E, NE, out = self.idx, self.NE, []
        for ki, k in enumerate(g1.TYPES):
            for s in g1.STRANDS:
                (chL, cfL), (chR, cfR) = ((0, 0), (1, 1)) if s == "+" else ((3, 1), (2, 0))
                for r in ("S", "A"):
                    out.append((0, E[(k, s, "E", r, 0)], chL, cfL, ki, "entry"))
                for r in ("S", "Z"):
                    out.append((E[(k, s, "E", r, NE - 1)], 0, chR, cfR, ki, "exit"))
        return out


class G2PHParams(nn.Module):
    def __init__(self, space):
        super().__init__()
        self.space = space
        nk = len(g1.TYPES)
        self.trans = nn.Parameter(torch.zeros(nk, len(TRANS)))
        self.e_stay = nn.Parameter(torch.zeros(nk, len(ROLES), space.NE))
        self.e_adv = nn.Parameter(torch.zeros(nk, len(ROLES), max(space.NE - 1, 1)))
        self.branch = nn.Parameter(torch.zeros(nk, space.NBR))
        i_stay = torch.zeros(nk, space.NBR, space.NI)
        if space.NBR > 1:  # break the symmetry between the short and the long branch
            i_stay[:, 1:] = 1.0
        self.i_stay = nn.Parameter(i_stay)
        self.i_adv = nn.Parameter(torch.zeros(nk, space.NBR, space.NI))
        self.coef = nn.Parameter(torch.ones(nk, 2))
        self._edges = space.all_edges()

    def _value(self, key):
        if key[0] == "ig":                       # IG self-loop: the 0 reference, as in G1
            return self.trans.new_zeros(())
        g, k = key[0], g1.TYPES.index(key[1])
        if g == "t":
            return self.trans[k, TRANS.index(key[2])]
        if g == "es":
            return self.e_stay[k, key[2], key[3]]
        if g == "ea":
            return self.e_adv[k, key[2], key[3]]
        if g == "br":
            return self.branch[k, key[2]]
        if g == "is":
            return self.i_stay[k, key[2], key[3]]
        return self.i_adv[k, key[2], key[3]]

    def forward(self):
        """[n, n] transition log-potentials (NEG = forbidden), differentiable."""
        n = self.space.n
        frm = torch.tensor([u for u, _, _ in self._edges])
        to = torch.tensor([v for _, v, _ in self._edges])
        vals = torch.stack([self._value(key) for _, _, key in self._edges])
        T = torch.full((n, n), NEG, dtype=vals.dtype, device=vals.device)
        return T.index_put((frm.to(vals.device), to.to(vals.device)), vals)

    def bonus(self, bnd, clamp=10.0):
        """bnd [..., L, 4] logits -> per-step bonus [..., L, 16] (row t = step t-1 -> t; row 0 unused = 0)."""
        b = bnd.clamp(-clamp, clamp)
        cols = []
        for _, _, ch, cf, ki, side in self.space.boundary_edges():
            v = self.coef[ki, cf].to(b.dtype) * b[..., ch]
            if side == "entry":
                v = torch.cat([torch.zeros_like(v[..., :1]), v[..., 1:]], -1)
            else:
                v = torch.cat([torch.zeros_like(v[..., :1]), v[..., :-1]], -1)
            cols.append(v)
        return torch.stack(cols, -1)


# ---- numba forward-backward ---------------------------------------------------------------------------------

def _tables(space, T):
    T = np.asarray(T, np.float64)
    n = space.n
    preds = [np.flatnonzero(T[:, j] > NEG / 2) for j in range(n)]
    succs = [np.flatnonzero(T[i, :] > NEG / 2) for i in range(n)]
    mk = lambda lst: (np.array([len(p) for p in lst], np.int64),
                      np.array([np.pad(p, (0, max(map(len, lst)) - len(p))) for p in lst], np.int64))
    pn, pi = mk(preds)
    sn, si = mk(succs)
    eid = np.full((n, n), -1, np.int64)
    for e, (u, v, *_r) in enumerate(space.boundary_edges()):
        eid[u, v] = e
    return T, pi, pn, si, sn, eid


if njit is not None:
    @njit(cache=True)
    def _lse(a, b):
        if a == -np.inf:
            return b
        if b == -np.inf:
            return a
        if a > b:
            return a + np.log1p(np.exp(b - a))
        return b + np.log1p(np.exp(a - b))

    @njit(cache=True)
    def _fb_one(x, T, pi, pn, si, sn, eid, ev, base, nbase, C, Eout, post):
        L, n = x.shape
        a = np.full((L, n), -np.inf)
        b = np.full((L, n), -np.inf)
        for j in range(n):
            a[0, j] = x[0, j]
        for t in range(1, L):
            for j in range(n):
                s = -np.inf
                for q in range(pn[j]):
                    i = pi[j, q]
                    c = a[t - 1, i] + T[i, j]
                    if eid[i, j] >= 0:
                        c += ev[t, eid[i, j]]
                    s = _lse(s, c)
                a[t, j] = s + x[t, j]
        for j in range(n):
            b[L - 1, j] = 0.0
        for t in range(L - 2, -1, -1):
            for i in range(n):
                s = -np.inf
                for q in range(sn[i]):
                    j = si[i, q]
                    c = T[i, j] + x[t + 1, j] + b[t + 1, j]
                    if eid[i, j] >= 0:
                        c += ev[t + 1, eid[i, j]]
                    s = _lse(s, c)
                b[t, i] = s
        z = -np.inf
        for j in range(n):
            z = _lse(z, a[L - 1, j])
        for t in range(L):
            for j in range(n):
                p = np.exp(a[t, j] + b[t, j] - z)
                post[t, base[j]] += p
        for t in range(1, L):
            for j in range(n):
                bj = x[t, j] + b[t, j] - z
                for q in range(pn[j]):
                    i = pi[j, q]
                    c = a[t - 1, i] + T[i, j] + bj
                    e = eid[i, j]
                    if e >= 0:
                        c += ev[t, e]
                    w = np.exp(c)
                    C[i, j] += w
                    if e >= 0:
                        Eout[t, e] += w
        return z

    @njit(cache=True, parallel=True)
    def _fb_batch(X, T, pi, pn, si, sn, eid, EV, base, nbase):
        B, L, n = X.shape
        ne = EV.shape[2]
        Z = np.empty(B)
        P = np.zeros((B, L, nbase))
        Cs = np.zeros((B, n, n))
        Es = np.zeros((B, L, ne))
        for bi in prange(B):
            Z[bi] = _fb_one(X[bi], T, pi, pn, si, sn, eid, EV[bi], base, nbase, Cs[bi], Es[bi], P[bi])
        return Z, P, Cs, Es


def expand(space, x17, gates, allowed=None):
    """x17 [B, L, 17] potentials, gates [B, L, 17], allowed bool [B, L, 17] -> [B, L, n] expanded log-potentials."""
    x = np.asarray(x17, np.float64)[..., space.base] + np.asarray(gates, np.float64)[..., space.base]
    if allowed is not None:
        x = np.where(np.asarray(allowed)[..., space.base], x, NEG)
    return np.ascontiguousarray(x)


def forward_backward(space, X, T, EV):
    """X [B, L, n] expanded potentials, T [n, n], EV [B, L, 16] bonuses -> (logZ [B], post17 [B, L, 17],
    expected transition counts [B, n, n], expected boundary-edge use [B, L, 16])."""
    T, pi, pn, si, sn, eid = _tables(space, T)
    return _fb_batch(np.ascontiguousarray(X, np.float64), T, pi, pn, si, sn, eid,
                     np.ascontiguousarray(EV, np.float64), space.base, g1.S)


def viterbi(space, X, T, EV):
    """MAP expanded path for one sequence (X [L, n], EV [L, 16]) -> (base-state path int8 [L], score)."""
    T, pi, pn, si, sn, eid = _tables(space, T)
    path, score = g1._viterbi_nb(np.ascontiguousarray(X, np.float64), T, pi, pn, np.ones(space.n, np.bool_), eid,
                                 np.ascontiguousarray(EV, np.float64))
    return space.base[path.astype(np.int64)].astype(np.int8), score


def feasible(space, gates, allowed):
    """bool [B]: does the observed (allowed) partition contain at least one legal path?  Depends only on the masks
    and the grammar structure, so it is computed with zero potentials and zero transition scores."""
    T0 = np.where(G2PHParams(space)().detach().numpy() > NEG / 2, 0.0, NEG)   # only legal edges are visited
    ok = (np.asarray(gates)[..., space.base] > NEG / 2) & np.asarray(allowed)[..., space.base]
    X = np.where(ok, 0.0, -1e9)                     # log(#paths) can reach 1e5, so mask far below NEG
    B, L = X.shape[:2]
    z, *_ = forward_backward(space, X, T0, np.zeros((B, L, 16)))
    return z > -1e8


def nll_surrogate(space, prm, x17, bnd, gates, allowed, feas):
    """Partial-observation CRF loss for a batch (protocol v1.7): x17 [B, L, 17], bnd [B, L, 4] torch (with grad);
    gates / allowed numpy; feas bool [B] from feasible().  Returns (loss [B] on x17's device, valued
    log Z - log Z_obs with the exact gradient; 0 for infeasible windows) and (z, z_obs)."""
    dev = x17.device
    T = prm()
    bon = prm.bonus(bnd.double())
    Tn, bn = T.detach().cpu().numpy(), bon.detach().cpu().numpy()
    xn = x17.detach().double().cpu().numpy()
    z, p, c, e = forward_backward(space, expand(space, xn, gates), Tn, bn)
    zo, po, co, eo = forward_backward(space, expand(space, xn, gates, allowed), Tn, bn)
    sur = ((x17.double() * torch.as_tensor(p - po, device=dev)).sum((1, 2))
           + (bon * torch.as_tensor(e - eo, device=bon.device)).sum((1, 2)).to(dev)
           + (T[None] * torch.as_tensor(c - co, device=T.device)).sum((1, 2)).to(dev))
    val = torch.as_tensor(z - zo, device=dev)
    loss = sur + (val - sur.detach())
    keep = torch.as_tensor(np.asarray(feas), device=dev)
    return torch.where(keep, loss, torch.zeros_like(loss)), (z, zo)
