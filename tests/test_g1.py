"""P1.0 tests for the G1 grammar: exactness against brute-force path enumeration, legality, symmetry."""
import itertools

import numpy as np
import pytest
import torch

from ncrna_annot.grammar import g1

torch.set_default_dtype(torch.float64)
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
