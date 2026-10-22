"""P1.0 tests for the G1 grammar: exactness against brute-force path enumeration, legality, symmetry."""
import itertools

import numpy as np
import pytest
import torch

from ncrna_annot.grammar import g1



@pytest.fixture(autouse=True)
def _float64():
    """exact comparisons in double precision, restored afterwards (other test modules use float32)."""
    old = torch.get_default_dtype()
    torch.set_default_dtype(torch.float64)
    yield
    torch.set_default_dtype(old)

RC = str.maketrans("ACGT", "TGCA")


def rand_trans(seed):
    m = g1.G1Transitions()
    with torch.no_grad():
        g = torch.Generator().manual_seed(seed)
        for k in g1.TYPES:
            m.p[k].copy_(torch.randn(len(g1.PARAMS), generator=g))
    return m().double().detach()


def brute(seq, emit, T, allowed=None):
    """All S^L paths: log Z, posterior marginals, best path and score."""
    L = len(seq)
    x = emit + g1.gate_mask(g1.encode(seq)).astype(np.float64)
    if allowed is not None:
        x = np.where(allowed, x, g1.NEG)
    P = np.array(list(itertools.product(range(g1.S), repeat=L)))
    sc = x[np.arange(L), P].sum(1) + T[P[:, :-1], P[:, 1:]].sum(1)
    z = np.logaddexp.reduce(sc)
    w = np.exp(sc - z)
    marg = np.zeros((L, g1.S))
    for t in range(L):
        np.add.at(marg[t], P[:, t], w)
    b = int(np.argmax(sc))
    return z, marg, P[b], sc[b]


SEQS = ["AGTAG", "GTCAG", "CTAAC", "ACTGC", "GGTAG"]


@pytest.mark.parametrize("seq", SEQS)
@pytest.mark.parametrize("block", [1, 2, 3, 128])
def test_logz_matches_bruteforce(seq, block):
    rng = np.random.default_rng(hash(seq) % 1000 + block)
    emit = rng.normal(size=(len(seq), g1.S)) * 2
    T = rand_trans(block)
    z_b, _, _, _ = brute(seq, emit, T.numpy())
    gates = torch.tensor(g1.gate_mask(g1.encode(seq)))[None].double()
    z = g1.log_partition(torch.tensor(emit)[None], T, gates, block=block)
    assert abs(float(z) - z_b) < 1e-8


@pytest.mark.parametrize("seq", SEQS)
def test_marginals_and_viterbi_match_bruteforce(seq):
    rng = np.random.default_rng(7)
    emit = rng.normal(size=(len(seq), g1.S)) * 2
    T = rand_trans(3)
    z_b, marg_b, path_b, sc_b = brute(seq, emit, T.numpy())
    gates = g1.gate_mask(g1.encode(seq))
    marg = g1.marginals(torch.tensor(emit)[None], T, torch.tensor(gates)[None].double())[0].numpy()
    assert np.abs(marg - marg_b).max() < 1e-8
    path, sc = g1.viterbi(emit, T.numpy(), gates)
    assert abs(sc - sc_b) < 1e-8
    assert list(path) == list(path_b)


def test_observed_partition_matches_bruteforce():
    seq = "AGTAG"
    rng = np.random.default_rng(11)
    emit = rng.normal(size=(5, g1.S))
    allowed = rng.random((5, g1.S)) < 0.5
    allowed[:, 0] = True  # keep at least one feasible path
    T = rand_trans(5)
    z_b, _, _, _ = brute(seq, emit, T.numpy(), allowed)
    gates = torch.tensor(g1.gate_mask(g1.encode(seq)))[None].double()
    z = g1.log_partition(torch.tensor(emit)[None], T, gates, torch.tensor(allowed)[None])
    assert abs(float(z) - z_b) < 1e-8


def mirror_index():
    """state index under reverse complement: strand swap, L<->R."""
    m = {0: 0}
    for k in g1.TYPES:
        for p, q in (("E", "E"), ("I", "I"), ("L", "R"), ("R", "L")):
            m[g1.IDX[f"{k}{p}+"]] = g1.IDX[f"{k}{q}-"]
            m[g1.IDX[f"{k}{p}-"]] = g1.IDX[f"{k}{q}+"]
    return np.array([m[i] for i in range(g1.S)])


def test_reverse_complement_invariance():
    rng = np.random.default_rng(3)
    seq = "".join(rng.choice(list("ACGT"), 300))
    rc = seq.translate(RC)[::-1]
    emit = rng.normal(size=(300, g1.S))
    mi = mirror_index()
    emit_rc = np.empty_like(emit)
    emit_rc[:, mi] = emit[::-1]
    T = rand_trans(9)
    z1 = g1.log_partition(torch.tensor(emit)[None], T, torch.tensor(g1.gate_mask(g1.encode(seq)))[None].double())
    z2 = g1.log_partition(torch.tensor(emit_rc)[None], T, torch.tensor(g1.gate_mask(g1.encode(rc)))[None].double())
    assert abs(float(z1) - float(z2)) < 1e-7


def test_viterbi_paths_are_legal_on_long_random_input():
    rng = np.random.default_rng(5)
    seq = "".join(rng.choice(list("ACGT"), 20000))
    emit = rng.normal(size=(20000, g1.S)) * 3
    T = rand_trans(1).numpy()
    path, _ = g1.viterbi(emit, T, g1.gate_mask(g1.encode(seq)))
    names = [g1.STATES[i] for i in path]
    for t in range(1, len(path)):
        assert T[path[t - 1], path[t]] > g1.NEG / 2, (t, names[t - 1], names[t])
        a, b = names[t - 1], names[t]
        if a != "IG" and b != "IG":
            assert a[0] == b[0] and a[2] == b[2], "type/strand switch inside a gene"
        if b[1:2] == "L":
            assert seq[t:t + 2] in g1.MOTIF_L[b[2]]
        if b[1:2] == "R":
            assert seq[t - 1:t + 1] in g1.MOTIF_R[b[2]]


def test_path_to_genes():
    st = g1.IDX
    path = [0, st["lE+"], st["lE+"], st["lL+"], st["lI+"], st["lR+"], st["lE+"], 0, st["mE-"], 0]
    genes = g1.path_to_genes(np.array(path), offset=100)
    assert genes == [{"type": "l", "strand": "+", "exons": [[102, 103], [107, 107]]},
                     {"type": "m", "strand": "-", "exons": [[109, 109]]}]


def _runs_ok(names, min_exon, min_intron):
    """duration constraints on runs fully inside the sequence."""
    L = len(names)
    t = 0
    while t < L:
        n = names[t]
        if n != "IG" and n[1] == "E":
            u = t
            while u < L and names[u] == n:
                u += 1
            if t > 0 and u < L and u - t < min_exon:
                return False
            t = u
        elif n != "IG" and n[1] == "L":
            u = t
            while u < L and names[u] != "IG" and names[u][1] != "E":
                u += 1
            if t > 0 and u < L and u - t < min_intron:
                return False
            t = u
        else:
            t += 1
    return True


@pytest.mark.parametrize("seq", ["AGTAGTAGTA", "GTAAGGTCAG"])
def test_duration_viterbi_matches_constrained_bruteforce(seq):
    rng = np.random.default_rng(13)
    L = len(seq)
    sub = [0] + [g1.IDX[f"l{p}+"] for p in g1.PARTS]
    allowed = np.zeros((L, g1.S), bool)
    allowed[:, sub] = True
    emit = rng.normal(size=(L, g1.S)) * 2
    T = rand_trans(4).numpy()
    gates = g1.gate_mask(g1.encode(seq))
    x = np.where(allowed, emit + gates, g1.NEG)
    best, best_path = -np.inf, None
    for P in itertools.product(sub, repeat=L):
        if any(T[a, b] < g1.NEG / 2 for a, b in zip(P, P[1:])) or any(x[t, p] < g1.NEG / 2 for t, p in enumerate(P)):
            continue
        if not _runs_ok([g1.STATES[p] for p in P], 2, 4):
            continue
        sc = sum(x[t, p] for t, p in enumerate(P)) + sum(T[a, b] for a, b in zip(P, P[1:]))
        if sc > best:
            best, best_path = sc, P
    path, sc = g1.viterbi_duration(emit, T, gates, allowed, min_exon=2, min_intron=4)
    assert abs(sc - best) < 1e-8
    assert _runs_ok([g1.STATES[p] for p in path], 2, 4)


# ---- G1B (protocol v1.6 section 13): position-dependent boundary bonuses on gene entry / exit ----------------

def rand_boundary(L, seed):
    rng = np.random.default_rng(seed)
    coef = g1.BoundaryCoef()
    with torch.no_grad():
        for k in g1.TYPES:
            coef.c[k].copy_(torch.tensor(rng.normal(size=2)))
    bnd = torch.tensor(rng.normal(size=(L, 4)) * 2)
    return coef, bnd


def edge_score(P, frm, to, vals):
    """bonus of each path [N, L] given edges; vals [L-1, nE]."""
    out = np.zeros(len(P))
    for e, (i, j) in enumerate(zip(frm, to)):
        hit = (P[:, :-1] == i) & (P[:, 1:] == j)
        out += (hit * vals[None, :, e]).sum(1)
    return out


def brute_edges(seq, emit, T, edges, allowed=None):
    L = len(seq)
    x = emit + g1.gate_mask(g1.encode(seq)).astype(np.float64)
    if allowed is not None:
        x = np.where(allowed, x, g1.NEG)
    P = np.array(list(itertools.product(range(g1.S), repeat=L)))
    frm, to, vals = edges
    sc = x[np.arange(L), P].sum(1) + T[P[:, :-1], P[:, 1:]].sum(1) + edge_score(P, frm, to, vals.detach().numpy())
    z = np.logaddexp.reduce(sc)
    w = np.exp(sc - z)
    marg = np.zeros((L, g1.S))
    for t in range(L):
        np.add.at(marg[t], P[:, t], w)
    b = int(np.argmax(sc))
    return z, marg, P[b], sc[b]


@pytest.mark.parametrize("seq", SEQS[:3])
def test_g1b_logz_marginals_viterbi_match_bruteforce(seq):
    rng = np.random.default_rng(21)
    L = len(seq)
    emit = rng.normal(size=(L, g1.S)) * 2
    T = rand_trans(6)
    coef, bnd = rand_boundary(L, 22)
    edges = g1.boundary_edges(bnd, coef)
    z_b, marg_b, path_b, sc_b = brute_edges(seq, emit, T.numpy(), edges)
    gates = torch.tensor(g1.gate_mask(g1.encode(seq)))[None].double()
    bedges = g1.boundary_edges(bnd[None], coef)
    for block in (1, 2, 64):
        z = g1.log_partition(torch.tensor(emit)[None], T, gates, block=block, edges=bedges)
        assert abs(float(z) - z_b) < 1e-8
    marg = g1.marginals(torch.tensor(emit)[None], T, gates, edges=bedges)[0].numpy()
    assert np.abs(marg - marg_b).max() < 1e-8
    path, sc = g1.viterbi(emit, T.numpy(), gates[0].numpy(), edges=edges)
    assert abs(sc - sc_b) < 1e-8
    assert list(path) == list(path_b)


def test_g1b_observed_partition_and_gradient():
    seq = "AGTAG"
    rng = np.random.default_rng(23)
    emit = rng.normal(size=(5, g1.S))
    allowed = rng.random((5, g1.S)) < 0.6
    allowed[:, 0] = True
    T = rand_trans(7)
    coef, bnd = rand_boundary(5, 24)
    z_b, _, _, _ = brute_edges(seq, emit, T.numpy(), g1.boundary_edges(bnd, coef), allowed)
    gates = torch.tensor(g1.gate_mask(g1.encode(seq)))[None].double()
    b = bnd[None].clone().requires_grad_(True)
    z = g1.log_partition(torch.tensor(emit)[None], T, gates, torch.tensor(allowed)[None], edges=g1.boundary_edges(b, coef))
    assert abs(float(z) - z_b) < 1e-8
    z.sum().backward()
    assert b.grad is not None and torch.isfinite(b.grad).all()
    assert all(coef.c[k].grad is not None for k in g1.TYPES)


def test_g1b_reverse_complement_invariance():
    rng = np.random.default_rng(25)
    L = 300
    seq = "".join(rng.choice(list("ACGT"), L))
    rc = seq.translate(RC)[::-1]
    emit = rng.normal(size=(L, g1.S))
    mi = mirror_index()
    emit_rc = np.empty_like(emit)
    emit_rc[:, mi] = emit[::-1]
    coef, bnd = rand_boundary(L, 26)
    bnd_rc = bnd.flip(0)[:, [2, 3, 0, 1]]  # TSS+ <-> TSS-, TES+ <-> TES-
    T = rand_trans(10)
    z1 = g1.log_partition(torch.tensor(emit)[None], T, torch.tensor(g1.gate_mask(g1.encode(seq)))[None].double(),
                          edges=g1.boundary_edges(bnd[None], coef))
    z2 = g1.log_partition(torch.tensor(emit_rc)[None], T, torch.tensor(g1.gate_mask(g1.encode(rc)))[None].double(),
                          edges=g1.boundary_edges(bnd_rc[None], coef))
    assert abs(float(z1) - float(z2)) < 1e-7


@pytest.mark.parametrize("seq", ["AGTAGTAGTA", "GTAAGGTCAG"])
def test_g1b_duration_viterbi_matches_constrained_bruteforce(seq):
    rng = np.random.default_rng(27)
    L = len(seq)
    sub = [0] + [g1.IDX[f"l{p}+"] for p in g1.PARTS]
    allowed = np.zeros((L, g1.S), bool)
    allowed[:, sub] = True
    emit = rng.normal(size=(L, g1.S)) * 2
    T = rand_trans(8).numpy()
    coef, bnd = rand_boundary(L, 28)
    frm, to, vals = g1.boundary_edges(bnd * 3, coef)
    V = vals.detach().numpy()
    gates = g1.gate_mask(g1.encode(seq))
    x = np.where(allowed, emit + gates, g1.NEG)
    eb = {(i, j): e for e, (i, j) in enumerate(zip(frm, to))}
    best = -np.inf
    for P in itertools.product(sub, repeat=L):
        if any(T[a, b] < g1.NEG / 2 for a, b in zip(P, P[1:])) or any(x[t, p] < g1.NEG / 2 for t, p in enumerate(P)):
            continue
        if not _runs_ok([g1.STATES[p] for p in P], 2, 4):
            continue
        sc = sum(x[t, p] for t, p in enumerate(P)) + sum(T[a, b] for a, b in zip(P, P[1:]))
        sc += sum(V[t, eb[(a, b)]] for t, (a, b) in enumerate(zip(P, P[1:])) if (a, b) in eb)
        best = max(best, sc)
    path, sc = g1.viterbi_duration(emit, T, gates, allowed, min_exon=2, min_intron=4, edges=(frm, to, vals))
    assert abs(sc - best) < 1e-8
    assert _runs_ok([g1.STATES[p] for p in path], 2, 4)


@pytest.mark.parametrize("seq", SEQS[:3])
def test_marginals_np_match_bruteforce(seq):
    rng = np.random.default_rng(31)
    L = len(seq)
    emit = rng.normal(size=(L, g1.S)) * 2
    T = rand_trans(11)
    gates = g1.gate_mask(g1.encode(seq))
    z_b, marg_b, _, _ = brute(seq, emit, T.numpy())
    post, z = g1.marginals_np(emit, T.numpy(), gates)
    assert abs(z - z_b) < 1e-8 and np.abs(post - marg_b).max() < 1e-8
    coef, bnd = rand_boundary(L, 32)
    edges = g1.boundary_edges(bnd, coef)
    z_b, marg_b, _, _ = brute_edges(seq, emit, T.numpy(), edges)
    post, z = g1.marginals_np(emit, T.numpy(), gates, edges=edges)
    assert abs(z - z_b) < 1e-8 and np.abs(post - marg_b).max() < 1e-8
    allowed = rng.random((L, g1.S)) < 0.6
    allowed[:, 0] = True
    z_b, marg_b, _, _ = brute_edges(seq, emit, T.numpy(), edges, allowed)
    post, z = g1.marginals_np(emit, T.numpy(), gates, allowed, edges=edges)
    assert abs(z - z_b) < 1e-8 and np.abs(post - marg_b).max() < 1e-8


def test_marginals_np_match_torch_float64_on_long_input():
    rng = np.random.default_rng(33)
    L = 2000
    seq = "".join(rng.choice(list("ACGT"), L))
    emit = rng.normal(size=(L, g1.S)) * 2
    T = rand_trans(12)
    gates = g1.gate_mask(g1.encode(seq))
    coef, bnd = rand_boundary(L, 34)
    post, z = g1.marginals_np(emit, T.numpy(), gates, edges=g1.boundary_edges(bnd, coef))
    m = g1.marginals(torch.tensor(emit)[None], T, torch.tensor(gates)[None].double(),
                     edges=g1.boundary_edges(bnd[None], coef))[0].numpy()
    assert np.abs(post - m).max() < 1e-8 and np.abs(post.sum(1) - 1).max() < 1e-10
