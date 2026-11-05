"""ag_colab_job.py window tiling: the kept centres tile the chromosome once, and each kept centre is the genome."""
import importlib.util
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location("ag_job", f"{ROOT}/scripts/phase1/s2ag/ag_colab_job.py")
J = importlib.util.module_from_spec(spec)
spec.loader.exec_module(J)


def runner():
    r = object.__new__(J.Runner)
    r.np, r.letters = np, np.frombuffer(b"ACGTN", np.uint8)
    return r


def test_tiling_and_sequence():
    rng = np.random.default_rng(0)
    R = runner()
    for L in (J.STEP * 2 + 17, J.STEP * 3, 1000):
        codes = rng.integers(0, 5, L).astype(np.uint8)
        cover = np.zeros(L, int)
        for k in range(-(-L // J.STEP)):
            seq = R.seq(codes, k)
            assert len(seq) == J.W
            lo, hi = k * J.STEP, min((k + 1) * J.STEP, L)
            cover[lo:hi] += 1
            assert seq[J.H:J.H + hi - lo] == R.letters[codes[lo:hi]].tobytes().decode()
            s = k * J.STEP - J.H
            if s < 0:
                assert set(seq[:-s]) == {"N"}
            if s + J.W > L:
                assert set(seq[L - s:]) == {"N"}
        assert (cover == 1).all()
