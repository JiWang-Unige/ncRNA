#!/usr/bin/env python3
"""v1.5 diagnostic step 1: decode a whole chromosome on CPU with a P1-S2 checkpoint (S3, lnc >= 200 nt) -> GFF3.

usage: diag_decode.py <best.pt> <feats.npy> <codes.npy> <chrom> <out.gff3> [lnc_only 0|1]
"""
import os
import sys
import time

import numpy as np
import torch

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ncrna_annot.decode import chrom as C  # noqa: E402
from ncrna_annot.decode import genes as dec  # noqa: E402
from ncrna_annot.grammar import g1  # noqa: E402
from ncrna_annot.heads.c0 import C0Head  # noqa: E402
from ncrna_annot.train.screen import M_STATES  # noqa: E402


def main():
    ck, fpath, cpath, chrom, out = sys.argv[1:6]
    lnc_only = len(sys.argv) > 6 and sys.argv[6] == "1"
    torch.set_num_threads(4)
    t0 = time.time()
    feats = np.load(fpath, mmap_mode="r")
    codes = np.load(cpath, mmap_mode="r")
    head, trans = C0Head(feats.shape[-1]), g1.G1Transitions()
    st = torch.load(ck, map_location="cpu")
    head.load_state_dict(st["head"])
    trans.load_state_dict(st["trans"])
    L = len(feats) * 6
    pot = dec.chromosome_potentials(head, feats, codes[:L], device="cpu")
    if lnc_only:
        pot[:, M_STATES] = g1.NEG
    print("potentials s", round(time.time() - t0), flush=True)
    genes = C.decode_blocks(pot, codes[:L], trans().detach().double().numpy(), 20, 50)
    genes = [g for g in genes if g["type"] != "l" or sum(b - a + 1 for a, b in g["exons"]) >= 200]
    dec.write_gff3(genes, chrom, out)
    print("genes", len(genes), "lnc", sum(g["type"] == "l" for g in genes), "s", round(time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
