"""v2.0 feature overlay (scripts/phase1/v20/feat_overlay.py)."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts", "phase1", "v20"))
from feat_overlay import Overlay  # noqa: E402


def test_overlay_rows():
    base = np.arange(10 * 3 * 2, dtype=np.float16).reshape(10, 3, 2)
    new = -np.ones((2, 3, 2), np.float16) * np.array([1, 2], np.float16)[:, None, None]
    ov = Overlay(base, new, [7, 2])
    assert ov.shape == base.shape and len(ov) == 10
    assert np.array_equal(ov[7], new[0]) and np.array_equal(ov[2], new[1]) and np.array_equal(ov[3], base[3])
    got = ov[np.array([1, 2, 7, 9])]
    assert np.array_equal(got, np.stack([base[1], new[1], new[0], base[9]]))
    assert np.array_equal(ov[0:3], np.stack([base[0], base[1], new[1]]))
