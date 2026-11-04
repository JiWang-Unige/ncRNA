"""S4 (PLAN_BATCH2.md C): chain reading and interval mapping of s4_synteny.py against hand-computed coordinates."""
import io
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts", "phase1", "s4"))
import s4_synteny as S  # noqa: E402

CHAIN = """chain 100 tchr 100 + 10 30 qchr 200 - 50 70 1
5 3 2
12

chain 50 tchr 100 + 40 50 qplus 300 + 100 110 2
10

chain 10 other 100 + 0 10 qchr 200 + 0 10 3
10
"""


def test_reverse_chain_mapping():
    ch = S.read_chains(io.StringIO(CHAIN), "tchr")
    assert ch["t0"].tolist() == [10, 18, 40] and ch["size"].tolist() == [5, 12, 10]
    # t [12, 14) -> reverse-strand q [52, 54) -> plus-strand [200 - 54, 200 - 52)
    assert S.map_interval(ch, 12, 14) == [(0, -1, 146, 148)]
    # t [14, 20): 1 base in block 1 (q 54), gap t 15-17, 2 bases in block 2 (q 57, 58)
    assert S.map_interval(ch, 14, 20) == [(0, -1, 145, 146), (0, -1, 141, 143)]


def test_plus_chain_and_lift():
    ch = S.read_chains(io.StringIO(CHAIN), "tchr")
    assert S.map_interval(ch, 42, 45) == [(1, 1, 102, 105)]
    # exons 1-based (43, 45) and (47, 48) -> all 5 bases on qplus, strand kept
    assert S.lift(ch, [(43, 45), (47, 48)], "+") == (1, "+", 102, 108)
    # mostly unaligned exon -> not lifted
    assert S.lift(ch, [(60, 99)], "+") is None
    # reverse chain flips the strand
    assert S.lift(ch, [(13, 14)], "+")[1] == "-"
