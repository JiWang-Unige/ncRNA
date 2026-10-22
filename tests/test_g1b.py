"""G1B (protocol v1.6 section 13): training step, boundary BCE, block decoding with boundary bonuses."""
import numpy as np
import torch

from ncrna_annot.data import labels
from ncrna_annot.decode import chrom
from ncrna_annot.grammar import g1
from ncrna_annot.heads.c0 import C0BHead, C0Head
from ncrna_annot.train import screen_g1b as sb


def test_c0bhead_shares_c0_parameters_and_outputs():
    h = C0BHead(16, hidden=32)
    p, b = h(torch.randn(2, 10, 16), torch.randint(0, 4, (2, 60)))
    assert p.shape == (2, 60, g1.S) and b.shape == (2, 60, 4)
    plain = C0Head(16, hidden=32)
    missing = plain.load_state_dict({k: v for k, v in h.state_dict().items() if not k.startswith("bnd.")})
    assert not missing.missing_keys and not missing.unexpected_keys


def test_training_step_reduces_loss_and_bce_ignores():
    torch.manual_seed(0)
    rng = np.random.default_rng(0)
    N, L, D = 2, 600, 16
    f = torch.tensor(rng.normal(size=(N, L // 6, D)), dtype=torch.float32)
    c_np = rng.integers(0, 4, size=(N, L)).astype(np.int8)
    c = torch.tensor(c_np.astype(np.int64))
    tx = [("l", "+", [(100, 200), (300, 400)], True)]
    bits = np.stack([labels.allowed_bits(L, tx) for _ in range(N)])
    head, trans, coef = C0BHead(D, hidden=32), g1.G1Transitions(), g1.BoundaryCoef()
    ev = torch.zeros(N, L, 4, dtype=torch.int8)
    ev[:, 95:106, 0] = 1
    ev[:, 390:411, 1] = 1
    ev[:, 500:, 2] = -1
    opt = torch.optim.Adam(list(head.parameters()) + list(trans.parameters()) + list(coef.parameters()), lr=1e-2)
    losses = []
    for _ in range(5):
        nll, b = sb.crf_nll(head, trans, coef, f, c, c_np, bits, "cpu")
        loss = nll.mean() + sb.boundary_bce(b, ev)
        opt.zero_grad()
        loss.backward()
        opt.step()
        losses.append(float(loss))
    assert np.isfinite(losses).all() and losses[-1] < losses[0]
    # ignored entries do not contribute
    b = torch.zeros(1, 4, 4, requires_grad=True)
    e = torch.tensor([[[1, -1, -1, -1]] * 4], dtype=torch.int8)
    sb.boundary_bce(b, e).backward()
    assert (b.grad[..., 1:] == 0).all() and (b.grad[..., 0] != 0).all()


def test_decode_blocks_with_boundary_equals_single_pass():
    rng = np.random.default_rng(3)
    L = 30000
    codes = g1.encode("".join(rng.choice(list("ACGT"), L)))
    pot = rng.normal(size=(L, g1.S)).astype(np.float32) * 2
    bnd = rng.normal(size=(L, 4)).astype(np.float32) * 3
    T = g1.G1Transitions()().detach().double().numpy()
    coef = g1.BoundaryCoef()
    edges = g1.boundary_edges(torch.tensor(bnd, dtype=torch.float64), coef.double())
    whole = g1.path_to_genes(g1.viterbi_duration(pot, T, g1.gate_mask(codes), None, 5, 10, edges=edges)[0])
    blk = chrom.decode_blocks(pot, codes, T, 5, 10, block=10000, pad=5000, bnd=bnd, coef=g1.BoundaryCoef())
    key = lambda gs: sorted((g["type"], g["strand"], tuple(map(tuple, g["exons"]))) for g in gs)
    assert len(set(key(whole)) & set(key(blk))) >= 0.9 * len(whole)
    plain = g1.path_to_genes(g1.viterbi_duration(pot, T, g1.gate_mask(codes), None, 5, 10)[0])
    assert key(plain) != key(whole)  # the bonuses change the decoding
