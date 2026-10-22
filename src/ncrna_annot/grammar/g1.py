"""G1 grammar (Phase 1 protocol v1, D1): coding x lncRNA union HMM-CRF over bases, both strands.

States (17): IG, then per strand s in (+, -) and type k in (m = mRNA, l = lncRNA): E, L, I, R
  E exon, L / R = left / right intron boundary base in genomic left-to-right order, I intron interior.
Paths (left to right):  IG -> E -> (L -> I* -> R -> E)* -> IG ; a gene never switches type or strand.
Splice motifs gate the boundary bases (genomic orientation):
  + strand: L at i needs seq[i:i+2] in {GT, GC} (donor); R at j needs seq[j-1:j+1] == AG (acceptor)
  - strand: L at i needs seq[i:i+2] == CT (rc AG);      R at j needs seq[j-1:j+1] in {AC, GC} (rc donor)
Strand parameters are tied in mirror image so that log Z is invariant to reverse complementation
(minus-strand start = plus end, minus L-entry = plus acceptor, minus R-entry = plus donor).

Training: log_partition() (differentiable two-level log-space scan; posteriors = d logZ / d emissions).
Decoding: viterbi() (numba, CPU, whole chromosome).

G1B (protocol v1.6 section 13): optional position-dependent bonuses on the gene-entry (IG -> kE s) and gene-exit
(kE s -> IG) edges from a boundary head (TSS / TES logits per strand), see boundary_edges(); all inference
functions take `edges=(frm, to, vals)` with vals[..., t-1, e] = bonus of edge e on the step (t-1 -> t).
"""
import numpy as np
import torch
from torch import nn

NEG = -1e5
STRANDS, TYPES, PARTS = ("+", "-"), ("m", "l"), ("E", "L", "I", "R")
STATES = ["IG"] + [f"{k}{p}{s}" for s in STRANDS for k in TYPES for p in PARTS]
IDX = {n: i for i, n in enumerate(STATES)}
S = len(STATES)
BASE = {"A": 0, "C": 1, "G": 2, "T": 3}
MOTIF_L = {"+": {"GT", "GC"}, "-": {"CT"}}
MOTIF_R = {"+": {"AG"}, "-": {"AC", "GC"}}
# per type, the free transition parameters (IG self-loop is the 0 reference)
PARAMS = ("start", "end", "e_stay", "i_stay", "donor", "acceptor")


def encode(seq):
    """DNA string -> int8 codes (A0 C1 G2 T3, other 4)."""
    a = np.frombuffer(seq.upper().encode(), dtype=np.uint8)
    out = np.full(len(a), 4, np.int8)
    for ch, v in BASE.items():
        out[a == ord(ch)] = v
    return out


def _pairs(codes):
    """dinucleotide code at (i, i+1) for i in 0..L-2 as 2-char strings is slow; use ints 0..24."""
    c = codes.astype(np.int16)
    return c[:-1] * 5 + c[1:]


def _dinuc_set(strings):
    return {BASE[a] * 5 + BASE[b] for a, b in strings}


def gate_mask(codes):
    """[L, S] additive mask: NEG where a state's motif condition fails at that base, else 0."""
    L = len(codes)
    g = np.zeros((L, S), np.float32)
    d = np.full(L + 1, -1, np.int16)
    d[:L - 1] = _pairs(codes)  # d[i] = dinucleotide (i, i+1)
    for s in STRANDS:
        okL = np.isin(d[:L], list(_dinuc_set(MOTIF_L[s])))
        prev = np.concatenate([[-1], d[:L - 1]])  # dinucleotide (i-1, i)
        okR = np.isin(prev, list(_dinuc_set(MOTIF_R[s])))
        for k in TYPES:
            g[~okL, IDX[f"{k}L{s}"]] = NEG
            g[~okR, IDX[f"{k}R{s}"]] = NEG
    return g


class G1Transitions(nn.Module):
    """Tied transition log-potentials -> [S, S] matrix (NEG = forbidden)."""

    def __init__(self):
        super().__init__()
        self.p = nn.ParameterDict({k: nn.Parameter(torch.zeros(len(PARAMS))) for k in TYPES})

    def forward(self):
        T = torch.full((S, S), NEG, dtype=self.p["m"].dtype, device=self.p["m"].device)
        T = T.clone()
        T[0, 0] = 0.0
        for k in TYPES:
            st, en, es, is_, do, ac = self.p[k]
            for s in STRANDS:
                E, Lb, I, R = (IDX[f"{k}{p}{s}"] for p in PARTS)
                plus = s == "+"
                T[0, E] = st if plus else en      # gene starts (left end)
                T[E, 0] = en if plus else st      # gene ends (right end)
                T[E, E] = es
                T[E, Lb] = do if plus else ac     # left intron boundary
                T[Lb, I] = 0.0
                T[Lb, R] = 0.0                    # 2-bp intron (L immediately followed by R)
                T[I, I] = is_
                T[I, R] = ac if plus else do      # right intron boundary
                T[R, E] = 0.0
        return T


BOUNDARY = ("TSS+", "TES+", "TSS-", "TES-")


class BoundaryCoef(nn.Module):
    """G1B: per-type weights of the boundary logits on gene entry / exit (tied in mirror image)."""

    def __init__(self, init=1.0):
        super().__init__()
        self.c = nn.ParameterDict({k: nn.Parameter(torch.full((2,), float(init))) for k in TYPES})  # (start, end)


def boundary_edges(bnd, coef, clamp=10.0):
    """bnd [..., L, 4] boundary logits (BOUNDARY order) -> (frm, to, vals [..., L-1, 8]).

    + gene starts at t: IG -> kE+ on step t, bonus c_start * TSS+[t];  ends at t-1: kE+ -> IG, c_end * TES+[t-1]
    - gene left end at t: IG -> kE-, c_end * TES-[t];                 right end at t-1: kE- -> IG, c_start * TSS-[t-1]
    """
    b = bnd.clamp(-clamp, clamp)
    tss_p, tes_p, tss_m, tes_m = (b[..., i] for i in range(4))
    frm, to, vals = [], [], []
    for k in TYPES:
        cs, ce = coef.c[k][0], coef.c[k][1]
        Ep, Em = IDX[f"{k}E+"], IDX[f"{k}E-"]
        for f_, t_, v in ((0, Ep, cs * tss_p[..., 1:]), (Ep, 0, ce * tes_p[..., :-1]),
                          (0, Em, ce * tes_m[..., 1:]), (Em, 0, cs * tss_m[..., :-1])):
            frm.append(f_)
            to.append(t_)
            vals.append(v)
    return frm, to, torch.stack(vals, -1)


def _edge_table(frm, to, vals, n, L):
    """numpy edge-id table [n, n] (-1 = none) and per-step values [L, nE] (row t = step t-1 -> t)."""
    eid = np.full((n, n), -1, np.int64)
    if frm is None:
        return eid, np.zeros((L, 1))
    for e, (i, j) in enumerate(zip(frm, to)):
        eid[i, j] = e
    ev = np.zeros((L, len(frm)))
    ev[1:] = np.asarray(vals, np.float64)
    return eid, ev


def _lmm(A, B):
    """log-semiring matmul over the last two dims."""
    return torch.logsumexp(A.unsqueeze(-1) + B.unsqueeze(-3), dim=-2)


def _add_edges(W, frm, to, vals):
    add = torch.zeros_like(W)
    add[:, :, frm, to] = vals.to(W.dtype)
    return W + add


def log_partition(emit, T, gates, allowed=None, block=64, edges=None):
    """log Z for a batch.  emit [B, L, S], T [S, S], gates [B, L, S] (0 / NEG), allowed bool [B, L, S] or None.
    edges: optional (frm, to, vals [B, L-1, nE]) position-dependent edge bonuses (G1B).

    Any state may begin or end a window except that motif gates still apply.  Exact (no pruning).
    """
    x = emit + gates
    if allowed is not None:
        x = x.masked_fill(~allowed, NEG)
    Bt, L, _ = x.shape
    alpha = x[:, 0]                                         # [B, S]
    if L == 1:
        return torch.logsumexp(alpha, -1)
    W = T[None, None] + x[:, 1:, None, :]                   # [B, L-1, S, S]: from i (row) to j (col)
    if edges is not None:
        W = _add_edges(W, *edges)
    n = W.shape[1]
    nb = -(-n // block)
    pad = nb * block - n
    if pad:
        eye = torch.full((S, S), NEG, dtype=W.dtype, device=W.device)
        eye.fill_diagonal_(0.0)
        W = torch.cat([W, eye.expand(Bt, pad, S, S)], 1)
    W = W.view(Bt, nb, block, S, S)
    P = W[:, :, 0]
    for k in range(1, block):
        P = _lmm(P, W[:, :, k])
    for b in range(nb):
        alpha = torch.logsumexp(alpha.unsqueeze(-1) + P[:, b], dim=-2)
    return torch.logsumexp(alpha, -1)


def marginals(emit, T, gates, allowed=None, edges=None):
    """Posterior state marginals [B, L, S] via autograd (d logZ / d emit)."""
    e = emit.detach().requires_grad_(True)
    if edges is not None:
        edges = (edges[0], edges[1], edges[2].detach())
    z = log_partition(e, T.detach(), gates, allowed, edges=edges)
    (g,) = torch.autograd.grad(z.sum(), e)
    return g


try:
    from numba import njit
except ImportError:  # pragma: no cover
    njit = None


def _viterbi_py(x, T, start_ok, eid=None, ev=None):
    L, n = x.shape
    bp = np.zeros((L, n), np.int8)
    v = x[0] + np.where(start_ok, 0.0, NEG)
    preds = [np.flatnonzero(T[:, j] > NEG / 2) for j in range(n)]
    for t in range(1, L):
        nv = np.empty(n)
        for j in range(n):
            c = v[preds[j]] + T[preds[j], j]
            if eid is not None:
                e = eid[preds[j], j]
                c = c + np.where(e >= 0, ev[t, np.maximum(e, 0)], 0.0)
            a = int(np.argmax(c))
            bp[t, j] = preds[j][a]
            nv[j] = c[a] + x[t, j]
        v = nv
    path = np.empty(L, np.int8)
    path[-1] = int(np.argmax(v))
    for t in range(L - 1, 0, -1):
        path[t - 1] = bp[t, path[t]]
    return path, float(v.max())


if njit is not None:
    @njit(cache=True)
    def _viterbi_nb(x, T, pred_idx, pred_n, start_ok, eid, ev):
        L, n = x.shape
        bp = np.zeros((L, n), np.int8)
        v = np.empty(n)
        for j in range(n):
            v[j] = x[0, j] + (0.0 if start_ok[j] else NEG)
        nv = np.empty(n)
        for t in range(1, L):
            for j in range(n):
                best, arg = -1e300, 0
                for q in range(pred_n[j]):
                    i = pred_idx[j, q]
                    c = v[i] + T[i, j]
                    if eid[i, j] >= 0:
                        c += ev[t, eid[i, j]]
                    if c > best:
                        best, arg = c, i
                bp[t, j] = arg
                nv[j] = best + x[t, j]
            for j in range(n):
                v[j] = nv[j]
        path = np.empty(L, np.int8)
        m, a = -1e300, 0
        for j in range(n):
            if v[j] > m:
                m, a = v[j], j
        path[L - 1] = a
        for t in range(L - 1, 0, -1):
            path[t - 1] = bp[t, path[t]]
        return path, m


def viterbi(emit, T, gates, allowed=None, edges=None):
    """MAP state path for one sequence.  emit [L, S] (np/torch), T [S, S], gates [L, S].
    edges: optional (frm, to, vals [L-1, nE]) as in log_partition (unbatched, numpy or torch)."""
    x = np.asarray(emit, np.float64) + np.asarray(gates, np.float64)
    if allowed is not None:
        x = np.where(np.asarray(allowed), x, NEG)
    T = np.asarray(T, np.float64)
    start_ok = np.ones(S, np.bool_)
    frm, to, vals = edges if edges is not None else (None, None, None)
    eid, ev = _edge_table(frm, to, None if vals is None else _np(vals), S, len(x))
    if njit is None:
        return _viterbi_py(x, T, start_ok, eid, ev)
    preds = [np.flatnonzero(T[:, j] > NEG / 2) for j in range(S)]
    pred_n = np.array([len(p) for p in preds], np.int64)
    pred_idx = np.zeros((S, pred_n.max()), np.int64)
    for j, p in enumerate(preds):
        pred_idx[j, :len(p)] = p
    return _viterbi_nb(np.ascontiguousarray(x), T, pred_idx, pred_n, start_ok, eid, ev)


def _np(a):
    return a.detach().cpu().numpy() if torch.is_tensor(a) else np.asarray(a)


if njit is not None:
    @njit(cache=True)
    def _lse_acc(acc, v):
        if v == -np.inf:
            return acc
        if acc == -np.inf:
            return v
        if acc > v:
            return acc + np.log1p(np.exp(v - acc))
        return v + np.log1p(np.exp(acc - v))

    @njit(cache=True)
    def _fb_nb(x, T, pred_idx, pred_n, succ_idx, succ_n, eid, ev):
        L, n = x.shape
        a = np.full((L, n), -np.inf)
        b = np.full((L, n), -np.inf)
        for j in range(n):
            a[0, j] = x[0, j]
        for t in range(1, L):
            for j in range(n):
                s = -np.inf
                for q in range(pred_n[j]):
                    i = pred_idx[j, q]
                    c = a[t - 1, i] + T[i, j]
                    if eid[i, j] >= 0:
                        c += ev[t, eid[i, j]]
                    s = _lse_acc(s, c)
                a[t, j] = s + x[t, j]
        for j in range(n):
            b[L - 1, j] = 0.0
        for t in range(L - 2, -1, -1):
            for i in range(n):
                s = -np.inf
                for q in range(succ_n[i]):
                    j = succ_idx[i, q]
                    c = T[i, j] + x[t + 1, j] + b[t + 1, j]
                    if eid[i, j] >= 0:
                        c += ev[t + 1, eid[i, j]]
                    s = _lse_acc(s, c)
                b[t, i] = s
        z = -np.inf
        for j in range(n):
            z = _lse_acc(z, a[L - 1, j])
        post = np.exp(a + b - z)
        return post, z


def marginals_np(emit, T, gates, allowed=None, edges=None):
    """Exact posteriors [L, S] (float64) and log Z of one sequence by a sparse CPU forward-backward (numba).
    Same model as marginals() / log_partition() (incl. G1B edges), without the float32 block-scan drift."""
    x = np.asarray(emit, np.float64) + np.asarray(gates, np.float64)
    if allowed is not None:
        x = np.where(np.asarray(allowed), x, NEG)
    T = np.asarray(T, np.float64)
    frm, to, vals = edges if edges is not None else (None, None, None)
    eid, ev = _edge_table(frm, to, None if vals is None else _np(vals), S, len(x))
    preds = [np.flatnonzero(T[:, j] > NEG / 2) for j in range(S)]
    succs = [np.flatnonzero(T[i, :] > NEG / 2) for i in range(S)]
    mk = lambda lst: (np.array([len(p) for p in lst], np.int64),
                      np.array([np.pad(p, (0, max(map(len, lst)) - len(p))) for p in lst], np.int64))
    pred_n, pred_idx = mk(preds)
    succ_n, succ_idx = mk(succs)
    return _fb_nb(np.ascontiguousarray(x), T, pred_idx, pred_n, succ_idx, succ_n, eid, ev)


def path_to_genes(path, offset=0):
    """State path -> list of genes (type, strand, [(start, end) 1-based closed exons])."""
    genes, cur = [], None
    for t, st in enumerate(np.asarray(path)):
        name = STATES[st]
        if name == "IG":
            if cur:
                genes.append(cur)
                cur = None
            continue
        k, part, s = name[0], name[1], name[2]
        if cur is None or (cur["type"], cur["strand"]) != (k, s):
            if cur:
                genes.append(cur)
            cur = {"type": k, "strand": s, "exons": []}
        if part == "E":
            if cur["exons"] and cur["exons"][-1][1] == offset + t:  # extend (1-based end == t)
                cur["exons"][-1][1] = offset + t + 1
            else:
                cur["exons"].append([offset + t + 1, offset + t + 1])
    if cur:
        genes.append(cur)
    return [g for g in genes if g["exons"]]


# ---- decode-time minimum durations (protocol v1.3 section 10, "S3") ----------------------------------

def expand_duration(T, min_exon=20, min_intron=50):
    """Expanded state space with minimum exon / intron lengths (intron length counts L and R bases).

    Returns (T_exp [Se, Se], orig [Se] original-state index of each expanded state).
    Chains reuse the original self-loop score for internal steps; entry/exit scores are G1's.
    Runs that touch the window/chromosome ends may be shorter (any state may start or end a path).
    """
    T = np.asarray(T, np.float64)
    orig, slot = [0], {}
    for s in STRANDS:
        for k in TYPES:
            E, Lb, I, R = (IDX[f"{k}{p}{s}"] for p in PARTS)
            slot[(k, s, "E")] = list(range(len(orig), len(orig) + min_exon))
            orig += [E] * min_exon
            slot[(k, s, "L")] = [len(orig)]
            orig.append(Lb)
            n_i = max(min_intron - 2, 1)
            slot[(k, s, "I")] = list(range(len(orig), len(orig) + n_i))
            orig += [I] * n_i
            slot[(k, s, "R")] = [len(orig)]
            orig.append(R)
    Se = len(orig)
    Te = np.full((Se, Se), NEG)
    Te[0, 0] = T[0, 0]
    for s in STRANDS:
        for k in TYPES:
            E, Lb, I, R = (IDX[f"{k}{p}{s}"] for p in PARTS)
            es, ls, is_, rs = slot[(k, s, "E")], slot[(k, s, "L")][0], slot[(k, s, "I")], slot[(k, s, "R")][0]
            Te[0, es[0]] = T[0, E]                       # gene start
            for a, b in zip(es, es[1:]):
                Te[a, b] = T[E, E]
            Te[es[-1], es[-1]] = T[E, E]
            Te[es[-1], 0] = T[E, 0]                      # gene end (after >= min_exon bases)
            Te[es[-1], ls] = T[E, Lb]
            Te[ls, is_[0]] = T[Lb, I]
            for a, b in zip(is_, is_[1:]):
                Te[a, b] = T[I, I]
            Te[is_[-1], is_[-1]] = T[I, I]
            Te[is_[-1], rs] = T[I, R]
            Te[rs, es[0]] = T[R, E]
    return Te, np.array(orig, np.int64)


if njit is not None:
    @njit(cache=True)
    def _viterbi_dur_nb(x, Te, orig, pred_idx, pred_n, eid, ev):
        L = x.shape[0]
        n = Te.shape[0]
        bp = np.zeros((L, n), np.int16)
        v = np.empty(n)
        for j in range(n):
            v[j] = x[0, orig[j]]
        nv = np.empty(n)
        for t in range(1, L):
            for j in range(n):
                best, arg = -1e300, 0
                for q in range(pred_n[j]):
                    i = pred_idx[j, q]
                    c = v[i] + Te[i, j]
                    if eid[i, j] >= 0:
                        c += ev[t, eid[i, j]]
                    if c > best:
                        best, arg = c, i
                bp[t, j] = arg
                nv[j] = best + x[t, orig[j]]
            for j in range(n):
                v[j] = nv[j]
        m, a = -1e300, 0
        for j in range(n):
            if v[j] > m:
                m, a = v[j], j
        path = np.empty(L, np.int64)
        path[L - 1] = a
        for t in range(L - 1, 0, -1):
            path[t - 1] = bp[t, path[t]]
        out = np.empty(L, np.int8)
        for t in range(L):
            out[t] = orig[path[t]]
        return out, m


def viterbi_duration(emit, T, gates, allowed=None, min_exon=20, min_intron=50, edges=None):
    """MAP path under G1 with minimum exon/intron lengths; returns original-state path and score.
    edges (G1B): bonuses on IG -> kE s / kE s -> IG are applied to the expanded gene entry / exit edges."""
    x = np.asarray(emit, np.float64) + np.asarray(gates, np.float64)
    if allowed is not None:
        x = np.where(np.asarray(allowed), x, NEG)
    Te, orig = expand_duration(T, min_exon, min_intron)
    n = len(orig)
    if edges is None:
        eid, ev = _edge_table(None, None, None, n, len(x))
    else:
        frm, to, vals = edges
        efrm, eto = [], []
        for i, j in zip(frm, to):  # map original entry / exit edges onto the expanded chain ends
            if i == 0:
                efrm.append(0)
                eto.append(int(np.flatnonzero((orig == j) & (Te[0] > NEG / 2))[0]))
            else:
                cand = np.flatnonzero((orig == i) & (Te[:, 0] > NEG / 2))
                efrm.append(int(cand[0]))
                eto.append(0)
        eid, ev = _edge_table(efrm, eto, _np(vals), n, len(x))
    preds = [np.flatnonzero(Te[:, j] > NEG / 2) for j in range(n)]
    pred_n = np.array([len(p) for p in preds], np.int64)
    pred_idx = np.zeros((n, pred_n.max()), np.int64)
    for j, p in enumerate(preds):
        pred_idx[j, :len(p)] = p
    return _viterbi_dur_nb(np.ascontiguousarray(x), Te, orig, pred_idx, pred_n, eid, ev)
