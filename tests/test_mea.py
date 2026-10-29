"""v1.9 development (chr18 only): MEA decoding under the G1 grammar, thresholded lnc L1 grid, weight averaging of
the last epochs, and the eval-only entry points (G1B and arm R) with posteriors saved for CPU decoders."""
import numpy as np
import torch

from ncrna_annot.decode import chrom
from ncrna_annot.grammar import g1
from ncrna_annot.train import screen2, screen_g1b

from test_formal_arms import D, _chrom, _data


def _gene_seq(rng, L=1200):
    """random sequence with a + strand lncRNA gene: exons [100, 200), [300, 400) 0-based, intron GT ... AG."""
    s = list(rng.choice(list("ACGT"), L))
    s[200:202], s[298:300] = "GT", "AG"
    path = np.zeros(L, np.int8)
    path[100:200] = g1.IDX["lE+"]
    path[200] = g1.IDX["lL+"]
    path[201:299] = g1.IDX["lI+"]
    path[299] = g1.IDX["lR+"]
    path[300:400] = g1.IDX["lE+"]
    return g1.encode("".join(s)), path


def test_mea_recovers_a_clean_path_and_obeys_the_grammar():
    rng = np.random.default_rng(0)
    codes, path = _gene_seq(rng)
    post = np.full((len(path), g1.S), 0.1 / g1.S)
    post[np.arange(len(path)), path] += 0.9
    T = g1.G1Transitions()().detach().double().numpy()
    got, _ = g1.mea_duration(post, T, g1.gate_mask(codes))
    assert (got == path).all()
    for kw in ({"log_gain": True}, {"level": "exon"}, {"level": "exon", "log_gain": True}):
        assert (g1.mea_duration(post, T, g1.gate_mask(codes), **kw)[0] == path).all(), kw
    assert g1.path_to_genes(got) == [{"type": "l", "strand": "+", "exons": [[101, 200], [301, 400]]}]
    # a 10-bp "exon" in the posteriors cannot be called (minimum exon 20)
    short = np.full((300, g1.S), 0.01)
    short[:, 0] = 0.6
    short[100:110, g1.IDX["lE+"]] = 0.99
    got2, _ = g1.mea_duration(short, T, np.zeros((300, g1.S)))
    assert all(sum(b - a + 1 for a, b in g["exons"]) >= 20 for g in g1.path_to_genes(got2))


def test_mea_gamma_calls_split_genes():
    """A gene whose posterior is split between exon and intron states (IG 0.4, exon 0.3, intron 0.3 over 400 bp):
    per-state scores with gamma 1 (like S1's argmax) call nothing; gamma 2 calls it."""
    L = 1000
    post = np.zeros((L, g1.S))
    post[:, 0] = 1.0
    post[300:700, 0] = 0.4
    post[300:700, g1.IDX["lE+"]] = 0.3
    post[300:700, g1.IDX["lI+"]] = 0.3
    T = g1.G1Transitions()().detach().double().numpy()
    n = lambda **kw: int((g1.mea_duration(post, T, np.zeros((L, g1.S)), **kw)[0] != 0).sum())
    assert n(gamma=1.0) == 0 and n(gamma=2.0) >= 380 and n(gamma=4.0) >= n(gamma=2.0)


def test_exon_level_needs_exon_evidence_but_joins_supported_exons():
    """An intron-only posterior stretch (a tile inside a long intron): state-level MEA invents two minimum-length
    exons around it, exon-level MEA calls nothing; two supported exons around a supported intron stay one gene."""
    T = g1.G1Transitions()().detach().double().numpy()
    L = 6000
    post = np.zeros((L, g1.S))
    post[:, 0] = 1.0
    post[1000:5000, 0] = 0.0
    post[1000:5000, g1.IDX["mI-"]] = 1.0
    gates = np.zeros((L, g1.S))
    assert len(g1.path_to_genes(g1.mea_duration(post, T, gates)[0])) == 1  # spurious 2-exon gene
    assert g1.path_to_genes(g1.mea_duration(post, T, gates, level="exon")[0]) == []
    post[900:1000, 0], post[900:1000, g1.IDX["mE-"]] = 0.0, 1.0
    post[5000:5100, 0], post[5000:5100, g1.IDX["mE-"]] = 0.0, 1.0
    genes = g1.path_to_genes(g1.mea_duration(post, T, gates, level="exon")[0])
    assert [len(g["exons"]) for g in genes] == [2] and genes[0]["exons"][0][0] == 901


def test_decode_blocks_mea_equals_single_pass():
    rng = np.random.default_rng(2)
    L = 30000
    codes = g1.encode("".join(rng.choice(list("ACGT"), L)))
    post = rng.dirichlet(np.full(g1.S, 0.5), size=L).astype(np.float32)
    T = g1.G1Transitions()().detach().double().numpy()
    whole = g1.path_to_genes(g1.mea_duration(post, T, g1.gate_mask(codes), 2.0, 5, 10)[0])
    blk = chrom.decode_blocks(post, codes, T, 5, 10, block=10000, pad=5000, mode="mea", gamma=2.0)
    key = lambda gs: sorted((g["type"], g["strand"], tuple(map(tuple, g["exons"]))) for g in gs)
    assert len(set(key(whole)) & set(key(blk))) >= 0.9 * len(whole)


def test_load_ckpt_averages_weights(tmp_path):
    a = {"head": {"w": torch.ones(3)}, "trans": {"p": torch.zeros(2)}, "epoch": 10}
    b = {"head": {"w": torch.full((3,), 3.0)}, "trans": {"p": torch.full((2,), 2.0)}, "epoch": 11}
    torch.save(a, tmp_path / "ep10.pt")
    torch.save(b, tmp_path / "ep11.pt")
    avg = screen2.load_ckpt([str(tmp_path / "ep10.pt"), str(tmp_path / "ep11.pt")], "cpu")
    assert torch.equal(avg["head"]["w"], torch.full((3,), 2.0)) and torch.equal(avg["trans"]["p"], torch.ones(2))
    assert avg["epoch"] == [10, 11]
    assert screen2.load_ckpt(str(tmp_path / "ep10.pt"), "cpu")["epoch"] == 10
    assert screen2.mea_specs("1,x2,pv1,pvx0.5") == [("MEA_g1", 1.0, "state", False), ("MEAX_g2", 2.0, "exon", False),
                                                    ("PV_g1", 1.0, "state", True), ("PVX_g0.5", 0.5, "exon", True)]
    assert screen2.mea_specs((2.0,)) == [("MEA_g2", 2.0, "state", False)]


def test_eval_only_mea_grid_average_and_saved_posteriors(tmp_path):
    rng = np.random.default_rng(5)
    tr, va, ch = _data(rng, 4), _data(rng, 2), _chrom(rng, tmp_path)
    q = lambda *_: None
    screen_g1b.run(1, tr, va, None, str(tmp_path / "g"), epochs=2, batch=2, device="cpu", log=q, keep_epochs=True)
    cks = [str(tmp_path / "g" / f"ep0{e}.pt") for e in (1, 2)]
    po = str(tmp_path / "post.npy")
    ev = screen_g1b.evaluate_checkpoint(cks, D, ch, str(tmp_path / "e"), device="cpu", log=q, with_s1=True,
                                        mea_gammas="1,x4", conf_grid=screen2.CONF_GRID, post_out=po)
    assert ev["best_epoch"] == [1, 2] and {"S3", "S1", "MEA_g1", "MEAX_g4"} <= set(ev)
    assert 0 <= ev["MEAX_g4"]["auprc_lnc_exonconf"] <= 1
    grid = ev["MEA_g1"]["conf_grid"]
    assert [g["tau"] for g in grid] == list(screen2.CONF_GRID) and all(a["n"] >= b["n"] for a, b in zip(grid, grid[1:]))
    post = np.load(po)
    assert post.shape == (len(ch["codes"]), g1.S) and post.dtype == np.float16
    # arm R (AB0) eval-only
    screen2.run("R", 0, tr, va, ch, str(tmp_path / "r"), epochs=1, batch=2, device="cpu", log=q, keep_epochs=True)
    er = screen2.evaluate_checkpoint([str(tmp_path / "r" / "ep01.pt")], D, ch, str(tmp_path / "re"), device="cpu",
                                     log=q, with_s1=True, mea_gammas="2,pv1", conf_grid=(0.5,), exact_post=True)
    assert {"S3", "S1", "MEA_g2", "PV_g1"} <= set(er) and len(er["MEA_g2"]["conf_grid"]) == 1


def test_exact_chromosome_marginals_match_autograd_and_sum_to_one():
    rng = np.random.default_rng(8)
    L = 30000
    codes = rng.integers(0, 4, L).astype(np.int8)
    pot = torch.tensor(rng.normal(size=(L, g1.S)), dtype=torch.float32)
    bnd = torch.tensor(rng.normal(size=(L, 4)), dtype=torch.float32)
    T, coef = g1.G1Transitions()().detach(), g1.BoundaryCoef()
    fn = lambda s, e: (pot[s:e], bnd[s:e])
    a = chrom.chromosome_marginals(fn, L, T, codes, "cpu", coef=coef).astype(np.float64)
    b = chrom.chromosome_marginals(fn, L, T, codes, "cpu", coef=coef, exact=True).astype(np.float64)
    assert np.abs(b.sum(1) - 1).max() < 5e-3  # float16 storage only
    assert np.abs(a - b).max() < 0.2 and np.abs(a - b).mean() < 0.01
