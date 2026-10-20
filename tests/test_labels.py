import numpy as np

from ncrna_annot.data import labels as L
from ncrna_annot.grammar import g1


def names(bits):
    return {g1.STATES[i] for i in range(g1.S) if (int(bits) >> i) & 1}


def test_precedence_and_strands():
    tx = [("m", "+", [(10, 20), (40, 50)]), ("l", "+", [(15, 25), (45, 60)]), ("l", "-", [(70, 75), (90, 95)])]
    c = L.strand_classes(100, tx, "+")
    assert list(c[10:15]) == [2] * 5 and list(c[15:25]) == [4] * 10
    assert list(c[25:40]) == [3] * 15 and list(c[40:45]) == [2] * 5 and list(c[45:60]) == [4] * 15
    b = L.allowed_bits(100, tx)
    assert names(b[0]) == {"IG"}
    assert names(b[12]) == {"mE+"}
    assert names(b[30]) == {"lL+", "lI+", "lR+"}
    assert names(b[72]) == {"lE-"}
    assert names(b[80]) == {"lL-", "lI-", "lR-"}


def test_both_strands_union_and_unknown():
    tx = [("m", "+", [(10, 20)]), ("l", "-", [(15, 25)])]
    unk = np.zeros(40, bool)
    unk[30:35] = True
    b = L.allowed_bits(40, tx, unk)
    assert names(b[17]) == {"mE+", "lE-"}
    assert names(b[32]) == {"IG"} | {f"l{p}{s}" for s in "+-" for p in "ELIR"}
    assert names(b[37]) == {"IG"}
    m = L.bits_to_bool(b)
    assert m.shape == (40, g1.S) and m[17, g1.IDX["mE+"]] and not m[17, 0]


def test_canonical_only_and_alt_isoform_unknown():
    canon = ("l", "+", [(10, 20), (40, 50)], True)
    alt = ("l", "+", [(10, 20), (28, 32), (40, 50)], False)   # alternative middle exon
    b = L.allowed_bits(60, [canon, alt])
    assert names(b[15]) == {"lE+"}
    assert names(b[25]) == {"lL+", "lI+", "lR+"}                # intron in both isoforms
    assert names(b[30]) == {"IG"} | {f"l{p}+" for p in "ELIR"}  # alt exon -> unknown for this gene


def test_repair_makes_observation_feasible():
    seq = "A" * 10 + "CCCCCCCCCC" + "A" * 10 + "CCCCCCCCCC" + "A" * 20   # no GT/AG anywhere: splicing impossible
    codes = g1.encode(seq)
    b = L.allowed_bits(len(seq), [("l", "+", [(10, 20), (30, 40)], True)])
    T = g1.G1Transitions()().detach().double().numpy()
    _, sc = g1.viterbi(np.zeros((len(seq), g1.S)), T, g1.gate_mask(codes), L.bits_to_bool(b))
    assert sc < g1.NEG / 2                                       # infeasible before repair
    rb, relaxed, nviol = L.repair(b, codes, flank=5)
    _, sc2 = g1.viterbi(np.zeros((len(seq), g1.S)), T, g1.gate_mask(codes), L.bits_to_bool(rb))
    assert sc2 > g1.NEG / 2 and relaxed > 0 and nviol >= 1
    assert (rb & b == b).all()                                   # repair only adds allowed states
