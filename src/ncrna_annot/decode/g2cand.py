"""G2-0 candidate events from a frozen G1B model (protocol v1.6 section 13).

Candidate channels (genomic orientation, the base where the event sits):
  0 L+  1 R+  2 L-  3 R-      splice: first / last intron base (G1 L / R states); + donor = L+, - donor = R-
  4 left+ 5 right+ 6 left- 7 right-   gene ends: + TSS = left+, + TES = right+, - TES = left-, - TSS = right-
Scores: splice = max over types k of the G1B posterior of the k{L,R}{s} state (0 off-motif, the states are gated);
        ends = sigmoid of the G1B boundary logits.
Candidates: splice positions with score >= theta; ends = peaks of the 21-bp box-smoothed probability
  (local maxima within +-25 bp) with smoothed value >= theta.
Thresholds (set on training windows): splice recall of exact canonical events >= 99%; ends recall within
  +-50 bp >= 95%.
"""
import numpy as np

from ..grammar import g1

CHANNELS = ("L+", "R+", "L-", "R-", "left+", "right+", "left-", "right-")
# exact event channel (prep_events.EV order: donor+, acceptor+, TSS+, TES+, donor-, acceptor-, TSS-, TES-) -> candidate channel
EV_TO_CAND = {0: 0, 1: 1, 2: 4, 3: 5, 4: 3, 5: 2, 6: 7, 7: 6}
SMOOTH, RADIUS, TOL = 21, 25, 50
TARGET = {"splice": 0.99, "ends": 0.95}
_BND = {"left+": 0, "right+": 1, "left-": 3, "right-": 2}   # index into g1.BOUNDARY (TSS+, TES+, TSS-, TES-)


def scores(post, bnd):
    """post [L, 17] G1B posteriors, bnd [L, 4] boundary logits -> float32 [L, 8] candidate scores."""
    post, bnd = np.asarray(post, np.float32), np.asarray(bnd, np.float32)
    out = np.zeros((len(post), 8), np.float32)
    for j, (part, s) in enumerate((("L", "+"), ("R", "+"), ("L", "-"), ("R", "-"))):
        out[:, j] = np.maximum(post[:, g1.IDX[f"m{part}{s}"]], post[:, g1.IDX[f"l{part}{s}"]])
    for j, name in enumerate(CHANNELS[4:], start=4):
        out[:, j] = 1.0 / (1.0 + np.exp(-np.clip(bnd[:, _BND[name]], -30, 30)))
    return out


def smooth(x, w=SMOOTH):
    c = np.concatenate([[0.0], np.cumsum(np.asarray(x, np.float64))])
    h = w // 2
    i = np.arange(len(x))
    lo, hi = np.maximum(i - h, 0), np.minimum(i + h + 1, len(x))
    return (c[hi] - c[lo]) / (hi - lo)


def peak_values(x, radius=RADIUS):
    """smoothed values at local maxima (the leftmost of a plateau), -1 elsewhere."""
    s = smooth(x)
    L = len(s)
    pad = np.concatenate([np.full(radius, -np.inf), s, np.full(radius, -np.inf)])
    win = np.lib.stride_tricks.sliding_window_view(pad, 2 * radius + 1)
    is_max = s >= win.max(1)
    left = win[:, :radius].max(1)
    is_max &= s > left                     # strict w.r.t. the left side: one peak per plateau
    out = np.full(L, -1.0)
    out[is_max] = s[is_max]
    return out


def candidates(sc, theta):
    """sc [L, 8] scores, theta [8] -> list of 8 sorted int arrays of candidate positions."""
    out = []
    for j in range(8):
        if j < 4:
            out.append(np.flatnonzero((sc[:, j] >= theta[j]) & (sc[:, j] > 0)))
        else:
            out.append(np.flatnonzero(peak_values(sc[:, j]) >= theta[j]))
    return out


def true_event_values(sc, pos, ch, tol=TOL):
    """best achievable candidate value for each exact event (pos, candidate channel ch) in one window:
    splice -> score at the base; ends -> max peak value within +-tol (0 if none)."""
    vals = np.zeros(len(pos))
    pk = {}
    for i, (p, c) in enumerate(zip(pos, ch)):
        if c < 4:
            vals[i] = sc[p, c]
        else:
            if c not in pk:
                pk[c] = peak_values(sc[:, c])
            w = pk[c][max(p - tol, 0):p + tol + 1]
            vals[i] = max(float(w.max()), 0.0)
    return vals


def choose_theta(values_by_channel):
    """largest thresholds reaching the target recall per channel: values_by_channel[c] = best values of true events."""
    theta = np.zeros(8)
    for c in range(8):
        v = np.sort(np.asarray(values_by_channel[c], np.float64))
        if not len(v):
            continue
        target = TARGET["splice" if c < 4 else "ends"]
        k = int(np.floor((1 - target) * len(v)))        # number of events allowed to fall below theta
        theta[c] = max(v[k], 1e-6)
    return theta
