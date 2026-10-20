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
