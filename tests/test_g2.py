"""G2 (protocol v1.6 section 13): exactness of the semi-Markov lattice against brute-force enumeration."""
import math

import numpy as np
import pytest
import torch

from ncrna_annot.grammar import g1, g2


@pytest.fixture(autouse=True)
def _float64():
    old = torch.get_default_dtype()
    torch.set_default_dtype(torch.float64)
    yield
    torch.set_default_dtype(old)


# ---- independent brute force ----------------------------------------------------------------------------------

def node_kinds(L, cands):
    kinds = {}
    for ch, arr in enumerate(cands):
        off, kd, s = g2.CAND_NODE[ch]
        for p in arr:
            x = int(p) + off
            if 0 < x < L:
                kinds.setdefault(x, set()).add((kd, s))
    return [0] + sorted(kinds) + [L], kinds


def start_ok(m, y, kinds):
    if m == 0:
        return True
    ks = kinds.get(m, set())
    if y == "IG":
        return ("GR", "+") in ks or ("GR", "-") in ks
    c, s = y[1], y[2]
    return ({"S": "GL", "A": "GL", "M": "IR", "Z": "IR", "N": "IL"}[c], s) in ks


def end_ok(n, y, kinds, L):
    if n == L:
        return True
    ks = kinds.get(n, set())
    if y == "IG":
        return ("GL", "+") in ks or ("GL", "-") in ks
    c, s = y[1], y[2]
    return ({"S": "GR", "Z": "GR", "A": "IL", "M": "IL", "N": "IR"}[c], s) in ks


def follows(prev, y):
    if prev is None:
        return True
    return g2.TS[g2.LIDX[prev], g2.LIDX[y]] == 0


def seg_score(a, b, y, L, x, gates, bnd, prm, h, allowed, min_len):
    """torch scalar or None (invalid)."""
    cl, cr = a != 0, b != L
    ln = b - a
    if y == "IG":
        if allowed is not None and not all(allowed[t, g1.IDX["IG"]] for t in range(a, b)):
            return None
        return sum(x[t, g1.IDX["IG"]] for t in range(a, b)) + prm.slope_ig * ln
    k, c, s = y[0], y[1], y[2]
    ki = g1.TYPES.index(k)
    sc = torch.zeros(())
    if c == "N":
        if cl and cr and ln < 2:
            return None
        for t in range(a, b):
            st = "L" if (t == a and cl) else ("R" if (t == b - 1 and cr) else "I")
            idx = g1.IDX[f"{k}{st}{s}"]
            if st in "LR" and gates[t, idx] <= g1.NEG / 2:
                return None
            if allowed is not None:
                if st in "LR" or (cl and cr):
                    if not allowed[t, idx]:
                        return None
                elif not any(allowed[t, g1.IDX[f"{k}{p}{s}"]] for p in "LIR"):
                    return None
            sc = sc + x[t, idx]
    else:
        E = g1.IDX[f"{k}E{s}"]
        if allowed is not None and not all(allowed[t, E] for t in range(a, b)):
            return None
        sc = sc + sum(x[t, E] for t in range(a, b))
    r = g2.ROLES.index(g2.role(c, s))
    sc = sc + prm.slope[ki, r] * ln
    (chL, cfL), (chR, cfR) = g2.end_channels(s)
    if c in "SA" and cl:
        sc = sc + prm.coef[ki, cfL] * bnd[a, chL]
    if c in "SZ" and cr:
        sc = sc + prm.coef[ki, cfR] * bnd[b - 1, chR]
    tl, tr = g2.trans_names(c, s)
    if tl and cl:
        sc = sc + prm.trans[ki, g2.TRANS.index(tl)]
    if tr and cr:
        sc = sc + prm.trans[ki, g2.TRANS.index(tr)]
    if prm.level >= 2 and cl and cr:
        sc = sc + prm.psi[ki, r, min(int(math.floor(math.log2(ln))), g2.NB - 1)]
    if prm.level >= 3 and c != "N":
        M = prm.mlp
        u = h @ prm.proj.weight.T
        body = u[a:b, :M].mean(0)
        left = u[max(a - 50, 0):a, M:2 * M].sum(0) / max(a - max(a - 50, 0), 1)
        right = u[b:min(b + 50, L), 2 * M:].sum(0) / max(min(b + 50, L) - b, 1)
        z = body + left + right + prm.u_len * math.log(ln) + prm.b1
        col = g2.EXON_LABELS.index(y)
        sc = sc + torch.nn.functional.gelu(z) @ prm.out.weight[col] + prm.out.bias[col]
    if min_len is not None and cl and cr and ln < (min_len[1] if c == "N" else min_len[0]):
        return None
    return sc


def brute(L, cands, x, gates, bnd, prm, h=None, allowed=None, min_len=None):
    nodes, kinds = node_kinds(L, cands)
    paths, cache = [], {}

    def score(m, n, y):
        if (m, n, y) not in cache:
            cache[(m, n, y)] = seg_score(m, n, y, L, x, gates, bnd, prm, h, allowed, min_len)
        return cache[(m, n, y)]

    def rec(m, prev, acc, segs):
        for y in g2.LABELS:
            if not follows(prev, y) or not start_ok(m, y, kinds):
                continue
            for n in nodes:
                if n <= m or not end_ok(n, y, kinds, L):
                    continue
                sc = score(m, n, y)
                if sc is None:
                    continue
                if n == L:
                    paths.append((acc + sc, segs + [(m, n, y)]))
                else:
                    rec(n, y, acc + sc, segs + [(m, n, y)])

    rec(0, None, torch.zeros(()), [])
    scores = torch.stack([p[0] for p in paths])
    return torch.logsumexp(scores, 0), paths


# ---- fixtures -------------------------------------------------------------------------------------------------

def setup(seed, level, L=12, obs=False):
    rng = np.random.default_rng(seed)
    codes = rng.integers(0, 4, L).astype(np.int8)
    gates = torch.tensor(g1.gate_mask(codes), dtype=torch.float64)
    x = torch.tensor(rng.normal(size=(L, g1.S)), requires_grad=True)
    bnd = torch.tensor(rng.normal(size=(L, 4)), requires_grad=True)
    H = 6
    h = torch.tensor(rng.normal(size=(L, H)), requires_grad=True)
    prm = g2.G2Params(level=level, hidden=H, mlp=3).double()
    with torch.no_grad():
        for p in prm.parameters():
            p.copy_(torch.tensor(rng.normal(size=p.shape)) * 0.5)
    # candidates: splice at motif-valid bases where possible, plus a few ends
    motif = lambda st: [t for t in range(L) if gates[t, g1.IDX[st]] > g1.NEG / 2]
    pick = lambda xs, n: sorted(rng.choice(xs, min(n, len(xs)), replace=False).tolist()) if xs else []
    cands = [pick(motif("lL+"), 2), pick(motif("lR+"), 2), pick(motif("lL-"), 1), pick(motif("lR-"), 1),
             pick(list(range(1, L - 1)), 2), pick(list(range(1, L - 1)), 2), pick(list(range(1, L - 1)), 1),
             pick(list(range(1, L - 1)), 1)]
    allowed = None
    if obs:
        allowed = rng.random((L, g1.S)) < 0.7
        allowed[:, g1.IDX["IG"]] = True
    return codes, gates, x, bnd, h, prm, cands, allowed


def lattice_logz(L, cands, x, gates, bnd, prm, h, allowed=None, min_len=None, K=64):
    lat = g2.collate([g2.build_lattice(L, cands, K=K)])
    al = None if allowed is None else torch.tensor(allowed)[None]
    SEG, IG = g2.segment_scores(x[None], gates[None], bnd[None], lat, prm, h=h[None] if h is not None else None,
                                allowed=al, min_len=min_len)
    return SEG, IG, lat


# ---- tests ----------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("level", [1, 2, 3])
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_logz_and_gradients_match_bruteforce(level, seed):
    L = 12
    codes, gates, x, bnd, h, prm, cands, _ = setup(seed, level, L)
    zb, paths = brute(L, cands, x, gates, bnd, prm, h)
    SEG, IG, lat = lattice_logz(L, cands, x, gates, bnd, prm, h)
    z, sur = g2.log_partition(SEG, IG, lat)
    assert abs(float(z[0]) - float(zb)) < 1e-8, (float(z[0]), float(zb), len(paths))
    leaves = [x, bnd, h] + list(prm.parameters())
    gb = torch.autograd.grad(zb, leaves, allow_unused=True)
    gs = torch.autograd.grad(sur.sum(), leaves, allow_unused=True)
    for a, b in zip(gb, gs):
        a = torch.zeros(()) if a is None else a
        b = torch.zeros(()) if b is None else b
        assert torch.allclose(a, b, atol=1e-8), (a, b)


@pytest.mark.parametrize("seed", [3, 4])
def test_observed_partition_and_nll(seed):
    L = 12
    codes, gates, x, bnd, h, prm, cands, allowed = setup(seed, 2, L, obs=True)
    zb, _ = brute(L, cands, x, gates, bnd, prm, h)
    zob, _ = brute(L, cands, x, gates, bnd, prm, h, allowed=allowed)
    SEG, IG, lat = lattice_logz(L, cands, x, gates, bnd, prm, h)
    SEGo, IGo, _ = lattice_logz(L, cands, x, gates, bnd, prm, h, allowed=allowed)
    zo, _ = g2.log_partition(SEGo, IGo, lat)
    assert abs(float(zo[0]) - float(zob)) < 1e-8
    loss, feas = g2.nll(SEG, IG, SEGo, IGo, lat)
    assert bool(feas[0]) and abs(float(loss[0]) - float(zb - zob)) < 1e-8
    leaves = [x, bnd] + list(prm.parameters())
    gb = torch.autograd.grad(zb - zob, leaves, allow_unused=True)
    gl = torch.autograd.grad(loss.sum(), leaves, allow_unused=True)
    for a, b in zip(gb, gl):
        a = torch.zeros(()) if a is None else a
        b = torch.zeros(()) if b is None else b
        assert torch.allclose(a, b, atol=1e-8)


@pytest.mark.parametrize("seed", [5, 6, 7])
def test_viterbi_matches_bruteforce(seed):
    L = 12
    codes, gates, x, bnd, h, prm, cands, _ = setup(seed, 3, L)
    for min_len in (None, (2, 3)):
        _, paths = brute(L, cands, x, gates, bnd, prm, h, min_len=min_len)
        best = max(paths, key=lambda p: float(p[0]))
        SEG, IG, lat = lattice_logz(L, cands, x, gates, bnd, prm, h, min_len=min_len)
        segs, score = g2.viterbi(SEG, IG, lat)
        assert abs(float(score[0]) - float(best[0])) < 1e-8
        assert segs[0] == best[1]


def test_batch_padding_and_k_cap():
    L = 12
    out = []
    for seed in (0, 8):
        codes, gates, x, bnd, h, prm, cands, _ = setup(seed, 1, L)
        out.append((codes, gates, x, bnd, cands, prm))
    prm = out[0][5]
    lats = [g2.build_lattice(L, o[4], K=64) for o in out]
    lat = g2.collate(lats)
    x = torch.stack([o[2] for o in out]).detach()
    gates = torch.stack([o[1] for o in out])
    bnd = torch.stack([o[3] for o in out]).detach()
    SEG, IG = g2.segment_scores(x, gates, bnd, lat, prm)
    z, _ = g2.log_partition(SEG, IG, lat)
    for b in range(2):
        zb, _ = brute(L, out[b][4], x[b], gates[b], bnd[b], prm)
        assert abs(float(z[b]) - float(zb)) < 1e-8
    # a look-back cap can only remove paths
    lat1 = g2.collate([g2.build_lattice(L, out[0][4], K=1)])
    SEG1, IG1 = g2.segment_scores(x[:1], gates[:1], bnd[:1], lat1, prm)
    z1, _ = g2.log_partition(SEG1, IG1, lat1)
    assert float(z1[0]) <= float(z[0]) + 1e-9


def test_segments_to_genes():
    segs = [(0, 3, "IG"), (3, 5, "lA+"), (5, 8, "lN+"), (8, 10, "lZ+"), (10, 11, "IG"), (11, 12, "mS-")]
    assert g2.segments_to_genes(segs, offset=100) == [
        {"type": "l", "strand": "+", "exons": [[104, 105], [109, 110]]},
        {"type": "m", "strand": "-", "exons": [[112, 112]]}]
