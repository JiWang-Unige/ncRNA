"""End-to-end CPU smoke of the Phase 1 training loop with random features (no backbone)."""
import json

import numpy as np

from ncrna_annot.data import labels
from ncrna_annot.decode import genes as dec
from ncrna_annot.grammar import g1
from ncrna_annot.train import loop


def test_train_loop_cpu(tmp_path):
    rng = np.random.default_rng(0)
    N, L, D = 6, 600, 32
    feats = rng.normal(size=(N, L // 6, D)).astype(np.float16)
    codes = rng.integers(0, 4, size=(N, L)).astype(np.int8)
    tx = [("l", "+", [(100, 200), (300, 400)]), ("m", "-", [(450, 500)])]
    bits = np.stack([labels.allowed_bits(L, tx) for _ in range(N)])
    gtf = tmp_path / "truth.gtf"
    gtf.write_text(
        'chrT\tX\texon\t101\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_type "lncRNA";\n'
        'chrT\tX\texon\t301\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_type "lncRNA";\n'
        'chrT\tX\texon\t1451\t1500\t.\t-\t.\tgene_id "G2"; transcript_id "T2"; gene_type "protein_coding";\n')
    val = {"feats": rng.normal(size=(3000 // 6, D)).astype(np.float16),
           "codes": rng.integers(0, 4, size=3000).astype(np.int8), "truth_gtf": str(gtf), "chrom": "chrT"}
    best = loop.train(feats, codes, bits, val, str(tmp_path / "run"), epochs=2, batch=4, device="cpu")
    recs = [json.loads(l) for l in open(tmp_path / "run" / "metrics.jsonl")]
    assert len(recs) == 2 and all(np.isfinite(r["loss_per_base"]) for r in recs)
    assert (tmp_path / "run" / "best.pt").exists() and best[1] in (1, 2)
    assert (tmp_path / "run" / "val_ep1.gff3").read_text().startswith("##gff-version 3")


def test_chromosome_potentials_tiling_matches_direct():
    import torch
    from ncrna_annot.heads.c0 import C0Head
    rng = np.random.default_rng(1)
    head = C0Head(16, hidden=32)
    ntok = 9000
    feats = rng.normal(size=(ntok, 16)).astype(np.float32)
    codes = rng.integers(0, 4, size=ntok * 6).astype(np.int8)
    pot = dec.chromosome_potentials(head, feats, codes, device="cpu")
    assert pot.shape == (ntok * 6, g1.S) and np.isfinite(pot).all()
    # the DNA conv is local and the GRU sees >=256 tokens of context: interior chunk equals a direct call
    with torch.no_grad():
        direct = head(torch.tensor(feats[3000:7096])[None], torch.tensor(codes[18000:42576].astype(np.int64))[None])[0].numpy()
    assert np.abs(pot[18000 + 6 * 1000:18000 + 6 * 1100] - direct[6 * 1000:6 * 1100]).max() < 0.5
