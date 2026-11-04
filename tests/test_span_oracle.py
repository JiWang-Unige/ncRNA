"""S1b (PLAN_S1B.md): the span-oracle lnc mask of splice_diag.py."""
import os
import sys

import numpy as np
import torch

from ncrna_annot.grammar import g1

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts", "phase1"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts", "phase1", "splice_diag"))
import splice_diag as sd  # noqa: E402

TRUTH = {
    "g1": {"strand": "+", "cls": "lnc", "tx": {"t1": [(11, 20), (41, 60)], "t2": [(5, 15)]}},
    "g2": {"strand": "-", "cls": "lnc", "tx": {"t3": [(70, 90)]}},
    "g3": {"strand": "+", "cls": "mrna", "tx": {"t4": [(100, 140)]}},
}


def test_span_mask_bounds():
    m = sd.span_mask(TRUTH, 150)
    assert m.shape == (150, 2)
    plus, minus = m[:, 0] == 0, m[:, 1] == 0
    assert np.flatnonzero(plus).tolist() == list(range(4, 60))     # 1-based 5..60 over both isoforms
    assert np.flatnonzero(minus).tolist() == list(range(69, 90))   # 1-based 70..90
    assert (m[~plus, 0] == g1.NEG).all() and (m[~minus, 1] == g1.NEG).all()  # mRNA gene g3 opens nothing


def test_masked_marginals_put_no_lnc_mass_outside_spans():
    rng = np.random.default_rng(0)
    L = 150
    codes = rng.integers(0, 4, L).astype(np.int8)
    emit = rng.normal(0, 1, (L, g1.S))
    T = g1.G1Transitions()().double().detach().numpy()
    m = sd.span_mask(TRUTH, L)
    x = emit.copy()
    for k, s in enumerate(g1.STRANDS):
        for part in g1.PARTS:
            x[:, g1.IDX[f"l{part}{s}"]] += m[:, k]
    post, _ = g1.marginals_np(x, T, g1.gate_mask(codes))
    for k, s in enumerate(g1.STRANDS):
        lnc = post[:, [g1.IDX[f"l{part}{s}"] for part in g1.PARTS]].sum(1)
        assert lnc[m[:, k] < 0].max() < 1e-12
    np.testing.assert_allclose(post.sum(1), 1.0, atol=1e-9)
    # inside the + span lnc states keep probability mass
    assert post[4:60, [g1.IDX[f"l{p}+"] for p in g1.PARTS]].sum() > 0.1
    torch.set_default_dtype(torch.float32)
