"""G2-PH (protocol v1.7 section 14): exactness and symmetry tests."""
import itertools

import numpy as np
import pytest
import torch

from ncrna_annot.grammar import g1, g2ph

RC = str.maketrans("ACGT", "TGCA")


@pytest.fixture(autouse=True)
def _float64():
    old = torch.get_default_dtype()
    torch.set_default_dtype(torch.float64)
    yield
    torch.set_default_dtype(old)


def rand_params(space, seed):
    prm = g2ph.G2PHParams(space).double()
    g = np.random.default_rng(seed)
    with torch.no_grad():
        for p in prm.parameters():
            p.copy_(torch.tensor(g.normal(size=p.shape)))
    return prm


def mirror17():
    m = {0: 0}
    for k in g1.TYPES:
        for p, q in (("E", "E"), ("I", "I"), ("L", "R"), ("R", "L")):
            m[g1.IDX[f"{k}{p}+"]] = g1.IDX[f"{k}{q}-"]
            m[g1.IDX[f"{k}{p}-"]] = g1.IDX[f"{k}{q}+"]
    return np.array([m[i] for i in range(g1.S)])


def brute(space, X, T, EV, states):
    """enumerate paths over `states` (all others are impossible) -> z, post17, C, E, best path, best score."""
    L = X.shape[0]
    P = np.array(list(itertools.product(states, repeat=L)))
    sc = X[np.arange(L), P].sum(1) + T[P[:, :-1], P[:, 1:]].sum(1)
    bes = space.boundary_edges()
    for e, (u, v, *_r) in enumerate(bes):
        hit = (P[:, :-1] == u) & (P[:, 1:] == v)
        sc += (hit * EV[None, 1:, e]).sum(1)
    z = np.logaddexp.reduce(sc)
    w = np.exp(sc - z)
    post = np.zeros((L, g1.S))
    for t in range(L):
        np.add.at(post[t], space.base[P[:, t]], w)
    C = np.zeros((space.n, space.n))
    for t in range(1, L):
        np.add.at(C, (P[:, t - 1], P[:, t]), w)
    E = np.zeros((L, len(bes)))
    for e, (u, v, *_r) in enumerate(bes):
        hit = (P[:, :-1] == u) & (P[:, 1:] == v)
        E[1:, e] = (hit * w[:, None]).sum(0)
    b = int(np.argmax(sc))
    return z, post, C, E, P[b], sc[b]


@pytest.mark.parametrize("cfg,seq", [((2, 1, 2), "AGTAGT"), ((1, 2, 1), "GTAAGTC"), ((2, 1, 2), "CTGACG")])
def test_restricted_bruteforce(cfg, seq):
    space = g2ph.Space(*cfg)
    L = len(seq)
    rng = np.random.default_rng(sum(map(ord, seq)))
    k, s = "l", "+" if seq[0] != "C" else "-"
    allowed = np.zeros((L, g1.S), bool)
    allowed[:, [0] + [g1.IDX[f"{k}{p}{s}"] for p in "ELIR"]] = True
    x17 = rng.normal(size=(L, g1.S)) * 2
    gates = g1.gate_mask(g1.encode(seq))
    prm = rand_params(space, 3)
    bnd = torch.tensor(rng.normal(size=(L, 4)) * 2)
    T = prm().detach().numpy()
    EV = prm.bonus(bnd).detach().numpy()
    X = g2ph.expand(space, x17[None], gates[None], allowed[None])[0]
    states = [0] + [i for i, n in enumerate(space.names) if n != "IG" and n[0] == k and n[1] == s]
    z_b, post_b, C_b, E_b, path_b, sc_b = brute(space, X, T, EV, states)
    z, post, C, E = g2ph.forward_backward(space, X[None], T, EV[None])
    assert abs(z[0] - z_b) < 1e-8
    assert np.abs(post[0] - post_b).max() < 1e-8
    assert np.abs(C[0] - C_b).max() < 1e-8
    assert np.abs(E[0] - E_b).max() < 1e-8
    path, score = g2ph.viterbi(space, X, T, EV)
    assert abs(score - sc_b) < 1e-8
    assert list(path) == list(space.base[path_b])


def dense_logz(space, prm, x17, bnd, gates, allowed=None):
    """independent torch float64 forward over the full expanded space (differentiable)."""
    T = prm()
    bon = prm.bonus(bnd)
    xe = x17[:, torch.as_tensor(space.base)] + torch.as_tensor(gates[:, space.base])
    if allowed is not None:
        xe = torch.where(torch.as_tensor(allowed[:, space.base]), xe, torch.full_like(xe, g1.NEG))
    bes = space.boundary_edges()
    frm = torch.tensor([u for u, *_ in bes])
    to = torch.tensor([v for _, v, *_ in bes])
    a = xe[0]
    for t in range(1, x17.shape[0]):
        Tt = T.index_put((frm, to), bon[t], accumulate=True)
        a = torch.logsumexp(a[:, None] + Tt, 0) + xe[t]
    return torch.logsumexp(a, 0)


@pytest.mark.parametrize("cfg", [(3, 2, 3), (1, 1, 1)])
def test_full_space_matches_dense_torch_and_gradients(cfg):
    space = g2ph.Space(*cfg)
    rng = np.random.default_rng(7)
    L = 40
    seq = "".join(rng.choice(list("ACGT"), L))
    gates = g1.gate_mask(g1.encode(seq))
    allowed = rng.random((L, g1.S)) < 0.8
    allowed[:, 0] = True
    prm = rand_params(space, 8)
    x17 = torch.tensor(rng.normal(size=(L, g1.S)), requires_grad=True)
    bnd = torch.tensor(rng.normal(size=(L, 4)) * 2, requires_grad=True)
    zd = dense_logz(space, prm, x17, bnd, gates)
    zod = dense_logz(space, prm, x17, bnd, gates, allowed)
    feas = g2ph.feasible(space, gates[None], allowed[None])
    assert bool(feas[0])
    loss, (z, zo) = g2ph.nll_surrogate(space, prm, x17[None], bnd[None], gates[None], allowed[None], feas)
    assert abs(z[0] - float(zd)) < 1e-8 and abs(zo[0] - float(zod)) < 1e-8
    assert abs(float(loss[0]) - float(zd - zod)) < 1e-8
    leaves = [x17, bnd] + list(prm.parameters())
    gd = torch.autograd.grad(zd - zod, leaves, allow_unused=True)
    gs = torch.autograd.grad(loss.sum(), leaves, allow_unused=True)
    for a, b in zip(gd, gs):
        a = torch.zeros(()) if a is None else a
        b = torch.zeros(()) if b is None else b
        assert torch.allclose(a, b, atol=1e-7), (a, b)


def test_reverse_complement_invariance():
    space = g2ph.Space(3, 2, 3)
    rng = np.random.default_rng(9)
    L = 300
    seq = "".join(rng.choice(list("ACGT"), L))
    rc = seq.translate(RC)[::-1]
    x17 = rng.normal(size=(L, g1.S))
    x_rc = np.empty_like(x17)
    x_rc[:, mirror17()] = x17[::-1]
    bnd = torch.tensor(rng.normal(size=(L, 4)) * 2)
    bnd_rc = bnd.flip(0)[:, [2, 3, 0, 1]]
    prm = rand_params(space, 10)
    T = prm().detach().numpy()
    z1, *_ = g2ph.forward_backward(space, g2ph.expand(space, x17[None], g1.gate_mask(g1.encode(seq))[None]), T,
                                   prm.bonus(bnd).detach().numpy()[None])
    z2, *_ = g2ph.forward_backward(space, g2ph.expand(space, x_rc[None], g1.gate_mask(g1.encode(rc))[None]), T,
                                   prm.bonus(bnd_rc).detach().numpy()[None])
    assert abs(z1[0] - z2[0]) < 1e-7


def test_posteriors_sum_to_one_and_feasibility():
    space = g2ph.Space()
    rng = np.random.default_rng(11)
    L = 2000
    seq = "".join(rng.choice(list("ACGT"), L))
    gates = g1.gate_mask(g1.encode(seq))[None]
    prm = rand_params(space, 12)
    z, post, C, E = g2ph.forward_backward(space, g2ph.expand(space, rng.normal(size=(1, L, g1.S)) * 3, gates),
                                          prm().detach().numpy(), prm.bonus(torch.tensor(rng.normal(size=(1, L, 4)))).detach().numpy())
    assert np.abs(post[0].sum(1) - 1).max() < 1e-9
    assert abs(C[0].sum() - (L - 1)) < 1e-6
    ok = np.zeros((2, 50, g1.S), bool)
    ok[:, :, 0] = True
    ok[1, 20, :] = False                    # no state allowed at one base
    assert g2ph.feasible(space, np.zeros((2, 50, g1.S)), ok).tolist() == [True, False]


def test_viterbi_genes_are_legal():
    space = g2ph.Space()
    rng = np.random.default_rng(13)
    L = 5000
    seq = "".join(rng.choice(list("ACGT"), L))
    prm = rand_params(space, 14)
    X = g2ph.expand(space, rng.normal(size=(1, L, g1.S)) * 3, g1.gate_mask(g1.encode(seq))[None])[0]
    path, _ = g2ph.viterbi(space, X, prm().detach().numpy(), np.zeros((L, 16)))
    names = [g1.STATES[i] for i in path]
    for t in range(1, L):
        a, b = names[t - 1], names[t]
        if a != "IG" and b != "IG":
            assert a[0] == b[0] and a[2] == b[2]
        if b[1:2] == "L":
            assert seq[t:t + 2] in g1.MOTIF_L[b[2]]
        if b[1:2] == "R":
            assert seq[t - 1:t + 1] in g1.MOTIF_R[b[2]]
    assert len(g1.path_to_genes(path)) > 0
