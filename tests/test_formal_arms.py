"""Protocol v1.8 arms on CPU with random features: G1B lnc-only training + S1 decoding (AB3 / AB1), and the R arm
(AB0, screen2.run) pause / resume reproducing an uninterrupted run, with every epoch's weights kept (keep_epochs)."""
import json

import numpy as np
import torch

from ncrna_annot.data import labels
from ncrna_annot.train import screen2, screen_g1b

L, D = 600, 16
TX = [("l", "+", [(100, 200), (300, 400)]), ("m", "-", [(450, 500)])]
GTF = ('chrT\tX\texon\t101\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_type "lncRNA";\n'
       'chrT\tX\texon\t301\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_type "lncRNA";\n'
       'chrT\tX\texon\t1451\t1500\t.\t-\t.\tgene_id "G2"; transcript_id "T2"; gene_type "protein_coding";\n')


def _data(rng, n):
    return dict(feats=rng.normal(size=(n, L // 6, D)).astype(np.float16),
                codes=rng.integers(0, 4, size=(n, L)).astype(np.int8),
                bits=np.stack([labels.allowed_bits(L, TX) for _ in range(n)]),
                events=np.zeros((n, L, 8), np.int8))


def _chrom(rng, tmp_path):
    gtf = tmp_path / "truth.gtf"
    gtf.write_text(GTF)
    return dict(feats=rng.normal(size=(30000 // 6, D)).astype(np.float16),
                codes=rng.integers(0, 4, size=30000).astype(np.int8), truth_gtf=str(gtf), chrom="chrT")


def test_g1b_lnc_only_and_s1(tmp_path):
    rng = np.random.default_rng(5)
    tr, va, ch = _data(rng, 4), _data(rng, 2), _chrom(rng, tmp_path)
    res = screen_g1b.run(1, tr, va, ch, str(tmp_path / "ab3"), epochs=1, batch=2, device="cpu", log=lambda *_: None,
                         lnc_only=True)
    assert res["lnc_only"] is True and res["S3"]["n_pred_mrna"] == 0  # mRNA states are impossible
    ev = screen_g1b.evaluate_checkpoint(str(tmp_path / "ab3" / "best.pt"), D, ch, str(tmp_path / "ab1"), device="cpu",
                                        log=lambda *_: None, with_s1=True)
    assert "S1" in ev and "S3" in ev and np.isfinite(ev["S1"]["auprc_lnc"])


def test_r_arm_pause_resume(tmp_path):
    rng = np.random.default_rng(6)
    tr, va, ch = _data(rng, 6), _data(rng, 2), _chrom(rng, tmp_path)
    q = lambda *_: None
    a = screen2.run("R", 0, tr, va, ch, str(tmp_path / "a"), epochs=3, batch=2, device="cpu", log=q, resume=True)
    n = {"i": 0}
    def pause():
        n["i"] += 1
        return n["i"] == 5
    try:
        screen2.run("R", 0, tr, va, ch, str(tmp_path / "b"), epochs=3, batch=2, device="cpu", log=q, pause=pause, resume=True,
                    keep_epochs=True)
        raise AssertionError("expected Paused")
    except screen2.Paused:
        pass
    b = screen2.run("R", 0, tr, va, ch, str(tmp_path / "b"), epochs=3, batch=2, device="cpu", log=q, resume=True,
                    keep_epochs=True)
    assert a["best_epoch"] == b["best_epoch"] and abs(a["best_val_nll"] - b["best_val_nll"]) < 1e-9
    ra = [json.loads(l)["val_nll"] for l in open(tmp_path / "a" / "metrics.jsonl")]
    rb = [json.loads(l)["val_nll"] for l in open(tmp_path / "b" / "metrics.jsonl")]
    assert ra == rb and "S1" in a
    ck = {k: torch.load(tmp_path / "b" / f"{k}.pt") for k in ["ep01", "ep02", "ep03", "best", "last"]}
    same = lambda x, y: all(torch.equal(x["head"][k], y["head"][k]) for k in x["head"])
    assert same(ck["ep03"], ck["last"]) and same(ck[f"ep{b['best_epoch']:02d}"], ck["best"])
    assert not (tmp_path / "a" / "ep01.pt").exists()


def test_gradient_accumulation_matches_full_batch(tmp_path):
    """batch 8 as 2 x 4 (shared-GPU memory): same validation NLL trajectory as the full batch, for G1B and R."""
    rng = np.random.default_rng(7)
    tr, va, ch = _data(rng, 8), _data(rng, 2), _chrom(rng, tmp_path)
    q = lambda *_: None
    for arm in ("G1B", "R"):
        runs = []
        for acc in (1, 2):
            out = str(tmp_path / f"{arm}{acc}")
            if arm == "G1B":
                screen_g1b.run(0, tr, va, None, out, epochs=2, batch=4, device="cpu", log=q, accum=acc)
            else:
                screen2.run("R", 0, tr, va, ch, out, epochs=2, batch=4, device="cpu", log=q, accum=acc)
            runs.append([json.loads(l)["val_nll"] for l in open(f"{out}/metrics.jsonl")])
        assert np.allclose(runs[0], runs[1], rtol=1e-5, atol=1e-6), (arm, runs)
