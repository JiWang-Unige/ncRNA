"""G2-0 candidate generation (protocol v1.6 section 13)."""
import numpy as np

from ncrna_annot.decode import g2cand


def test_peaks_one_per_plateau_and_radius():
    x = np.zeros(400)
    x[100:110] = 1.0          # plateau -> one peak
    x[200] = 0.5
    x[215] = 0.6              # within 25 bp of 200: only the larger survives
    pv = g2cand.peak_values(x)
    pk = np.flatnonzero(pv >= 0)
    assert sum(90 <= p <= 120 for p in pk) == 1
    assert sum(180 <= p <= 240 for p in pk) == 1


def test_choose_theta_reaches_target_recall():
    rng = np.random.default_rng(0)
    vals = {c: rng.random(1000) for c in range(8)}
    th = g2cand.choose_theta(vals)
    for c in range(8):
        r = np.mean(vals[c] >= th[c])
        t = g2cand.TARGET["splice" if c < 4 else "ends"]
        assert r >= t and r < t + 0.01


def test_candidates_and_true_values_consistent():
    rng = np.random.default_rng(1)
    L = 3000
    sc = rng.random((L, 8)).astype(np.float32) * 0.1
    sc[500, 0] = 0.9
    sc[1500:1521, 4] = 0.95
    th = np.full(8, 0.5)
    c = g2cand.candidates(sc, th)
    assert 500 in c[0] and len(c[0]) == 1
    assert len(c[4]) == 1 and abs(c[4][0] - 1510) <= 10
    v = g2cand.true_event_values(sc, np.array([500, 1540, 1700]), np.array([0, 4, 4]))
    assert v[0] == np.float32(0.9) and v[1] >= 0.5 and v[2] < 0.5   # 1540 within 50 bp of the peak, 1700 not
