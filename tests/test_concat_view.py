"""v3.0 Concat view: rows, slices and index lists match np.concatenate, also over an Overlay part."""
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, f"{ROOT}/scripts/phase1/v30")
sys.path.insert(0, f"{ROOT}/scripts/phase1/v20")
from concat_view import Concat  # noqa: E402
from feat_overlay import Overlay  # noqa: E402


def test_concat_matches_numpy():
    rng = np.random.default_rng(0)
    base = rng.normal(size=(7, 3, 2)).astype(np.float16)
    new = rng.normal(size=(2, 3, 2)).astype(np.float16)
    ov = Overlay(base, new, [1, 5])
    ref_ov = base.copy()
    ref_ov[1], ref_ov[5] = new[0], new[1]
    extra = rng.normal(size=(4, 3, 2)).astype(np.float16)
    c = Concat([ov, extra])
    ref = np.concatenate([ref_ov, extra])
    assert len(c) == 11 and c.shape == ref.shape and c.dtype == ref.dtype
    for i in range(11):
        assert np.array_equal(c[i], ref[i])
    assert np.array_equal(c[-1], ref[-1])
    assert np.array_equal(c[3:9], ref[3:9])
    idx = rng.permutation(11)[:6]
    assert np.array_equal(c[idx], ref[idx])
    assert np.array_equal(c[list(idx)], ref[list(idx)])
