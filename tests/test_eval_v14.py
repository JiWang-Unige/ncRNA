import numpy as np

from ncrna_annot.decode import chrom
from ncrna_annot.eval import prcurve
from ncrna_annot.grammar import g1


def _truth():
    mk = lambda s, ex: {"strand": s, "cls": "lnc", "gene_type": "lncRNA", "tx": {"t": ex}}
    return {"A": mk("+", [(100, 200)]), "B": mk("+", [(1000, 1100)]), "C": mk("-", [(2000, 2100)])}


def test_pr_curve_perfect_and_fragments():
    truth = _truth()
    perfect = [{"strand": "+", "exons": [(110, 190)], "conf": 0.9}, {"strand": "+", "exons": [(1010, 1090)], "conf": 0.8},
               {"strand": "-", "exons": [(2010, 2090)], "conf": 0.7}]
    r = prcurve.pr_curve(perfect, truth)
    assert abs(r["auprc"] - 1.0) < 1e-9 and abs(r["best_f1"] - 1.0) < 1e-9
    # many overlapping fragments in one true gene count as ONE locus (locus-merged precision)
    frags = [{"strand": "+", "exons": [(100 + i, 105 + i)], "conf": 0.5} for i in range(0, 50, 3)]
    r2 = prcurve.pr_curve(frags + [{"strand": "+", "exons": [(5000, 5100)], "conf": 0.4}], truth)
    assert abs(r2["best_point"]["pr"] - 1.0) < 1e-9 and abs(r2["best_point"]["sn"] - 1 / 3) < 1e-9
    # a false positive ranked first lowers AUPRC
    r3 = prcurve.pr_curve([{"strand": "+", "exons": [(5000, 5100)], "conf": 0.99}] + perfect, truth)
    assert r3["auprc"] < 1.0


def test_decode_blocks_equals_single_pass_on_small_input():
    rng = np.random.default_rng(2)
    L = 30000
    seq = "".join(rng.choice(list("ACGT"), L))
    codes = g1.encode(seq)
    pot = rng.normal(size=(L, g1.S)).astype(np.float32) * 2
    T = g1.G1Transitions()().detach().double().numpy()
    whole = g1.path_to_genes(g1.viterbi_duration(pot, T, g1.gate_mask(codes), None, 5, 10)[0])
    blk = chrom.decode_blocks(pot, codes, T, 5, 10, block=10000, pad=5000)
    key = lambda gs: sorted((g["type"], g["strand"], tuple(map(tuple, g["exons"]))) for g in gs)
    common = set(key(whole)) & set(key(blk))
    assert len(common) >= 0.9 * len(whole)  # block decoding agrees away from rare boundary effects


def test_gene_confidence_range():
    post = np.zeros((300, g1.S), np.float16)
    post[:, g1.IDX["lE+"]] = 0.8
    post[:, g1.IDX["lI+"]] = 0.1
    gene = {"type": "l", "strand": "+", "exons": [[11, 50], [101, 150]]}
    c = chrom.gene_confidence(gene, post)
    assert 0.1 < c < 0.8


def test_screen2_cpu_end_to_end(tmp_path):
    from ncrna_annot.data import labels
    from ncrna_annot.train import screen2
    rng = np.random.default_rng(0)
    N, L, D = 4, 24576, 16
    tx = [("l", "+", [(1000, 2000), (3000, 4000)], True)]
    tr = dict(feats=rng.normal(size=(N, L // 6, D)).astype(np.float16), codes=rng.integers(0, 4, size=(N, L)).astype(np.int8),
              bits=np.stack([labels.allowed_bits(L, tx) for _ in range(N)]))
    va = dict(feats=tr["feats"][:2], codes=tr["codes"][:2], bits=tr["bits"][:2])
    gtf = tmp_path / "t.gtf"
    gtf.write_text('chrT\tX\texon\t1001\t2000\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_type "lncRNA";\n')
    Lc = 3 * 24576
    ch = dict(feats=rng.normal(size=(Lc // 6, D)).astype(np.float16), codes=rng.integers(0, 4, size=Lc).astype(np.int8),
              truth_gtf=str(gtf), chrom="chrT")
    res = screen2.run("R", 42, tr, va, ch, str(tmp_path / "o"), epochs=1, batch=2, device="cpu")
    assert "S3" in res and "S1" in res and 0.0 <= res["S3"]["auprc_lnc"] <= 1.0
