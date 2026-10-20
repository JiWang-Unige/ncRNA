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


def _lmm(A, B):
    """log-semiring matmul over the last two dims."""
    return torch.logsumexp(A.unsqueeze(-1) + B.unsqueeze(-3), dim=-2)


def log_partition(emit, T, gates, allowed=None, block=128):
    """log Z for a batch.  emit [B, L, S], T [S, S], gates [B, L, S] (0 / NEG), allowed bool [B, L, S] or None.

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


def marginals(emit, T, gates, allowed=None):
    """Posterior state marginals [B, L, S] via autograd (d logZ / d emit)."""
    e = emit.detach().requires_grad_(True)
    z = log_partition(e, T.detach(), gates, allowed)
    (g,) = torch.autograd.grad(z.sum(), e)
    return g


try:
    from numba import njit
except ImportError:  # pragma: no cover
    njit = None


def _viterbi_py(x, T, start_ok):
    L, n = x.shape
    bp = np.zeros((L, n), np.int8)
    v = x[0] + np.where(start_ok, 0.0, NEG)
    preds = [np.flatnonzero(T[:, j] > NEG / 2) for j in range(n)]
    for t in range(1, L):
        nv = np.empty(n)
        for j in range(n):
            c = v[preds[j]] + T[preds[j], j]
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
    def _viterbi_nb(x, T, pred_idx, pred_n, start_ok):
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


def viterbi(emit, T, gates, allowed=None):
    """MAP state path for one sequence.  emit [L, S] (np/torch), T [S, S], gates [L, S]."""
    x = np.asarray(emit, np.float64) + np.asarray(gates, np.float64)
    if allowed is not None:
        x = np.where(np.asarray(allowed), x, NEG)
    T = np.asarray(T, np.float64)
    start_ok = np.ones(S, np.bool_)
    if njit is None:
        return _viterbi_py(x, T, start_ok)
    preds = [np.flatnonzero(T[:, j] > NEG / 2) for j in range(S)]
    pred_n = np.array([len(p) for p in preds], np.int64)
    pred_idx = np.zeros((S, pred_n.max()), np.int64)
    for j, p in enumerate(preds):
        pred_idx[j, :len(p)] = p
    return _viterbi_nb(np.ascontiguousarray(x), T, pred_idx, pred_n, start_ok)


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
