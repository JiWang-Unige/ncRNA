"""S1 (PLAN_S1.md): site conventions and the weighted AUROC of s1_analyze.py."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts", "phase1", "s1"))
import s1_analyze as s1  # noqa: E402
from ncrna_annot.grammar import g1  # noqa: E402


def test_canonical_mask_matches_intron_motifs():
    # exon(0-4) GT...AG intron(1-based 6..15) exon ; and the reverse complement for the minus strand
    seq = "AAAAA" + "GTCCCCCCAG" + "TTTTT"
    codes = g1.encode(seq)
    m = s1.canonical_mask(codes)
    (lc, lp), (rc, rp) = s1.junction_sites(np.array([6]), np.array([15]), np.array(["+"]))
    assert m[lp[0], lc[0]] and m[rp[0], rc[0]]          # don+ at 0-based 4, acc+ at 0-based 15
    assert (lp[0], rp[0]) == (4, 15) and (lc[0], rc[0]) == (1, 0)
    rcseq = seq[::-1].translate(str.maketrans("ACGT", "TGCA"))
    m2 = s1.canonical_mask(g1.encode(rcseq))
    (lc, lp), (rc, rp) = s1.junction_sites(np.array([6]), np.array([15]), np.array(["-"]))
    assert m2[lp[0], lc[0]] and m2[rp[0], rc[0]]        # acc- left, don- right
    assert m.sum() >= 2 and m2.sum() >= 2


def test_weighted_auroc_equals_bruteforce():
    rng = np.random.default_rng(1)
    pos, neg = rng.normal(1, 1, 50).round(1), rng.normal(0, 1, 80).round(1)
    wp, wn = rng.integers(0, 3, 50).astype(float), rng.integers(0, 3, 80).astype(float)
    order = np.argsort(neg, kind="stable")
    lo, hi = s1.wauroc_parts(pos, neg[order])
    got = s1.wauroc(lo, hi, wp, wn[order])
    num = sum(wp[i] * wn[j] * (1.0 if pos[i] > neg[j] else 0.5 if pos[i] == neg[j] else 0.0) for i in range(50) for j in range(80))
    assert abs(got - num / (wp.sum() * wn.sum())) < 1e-12
    grp = rng.integers(0, 3, 50)
    per = s1.wauroc(lo, hi, wp, wn[order], grp, 3)
    for k in range(3):
        sel = grp == k
        ref = s1.wauroc(lo[sel], hi[sel], wp[sel], wn[order])
        assert abs(per[k] - ref) < 1e-12
