#!/usr/bin/env python3
"""G2-0 candidate scores on CPU (protocol v1.6 section 13; used while the Kaggle GPU quota is exhausted).

usage: g2c_cpu.py <g1b best.pt> <out_dir> <split> [shard i/n] [threads]
Writes <out_dir>/<split>_cand_scores.npy (float32 [n, 24576, 8]: posteriors of true sites go below the float16 range, decode/g2cand.py channels); with a shard,
writes <split>_cand_scores.part<i>.npy holding only that shard's rows (merge with --merge).
Posteriors by the exact float64 numba forward-backward (g1.marginals_np), same model as the GPU kernel.
"""
import os
import sys
import time

import numpy as np
import torch

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ncrna_annot.decode import g2cand  # noqa: E402
from ncrna_annot.grammar import g1  # noqa: E402
from ncrna_annot.heads.c0 import C0BHead  # noqa: E402

ROOT = os.path.expanduser("~/ncrna-annot")


def main():
    ck, out, sp = sys.argv[1:4]
    i, n = map(int, sys.argv[4].split("/")) if len(sys.argv) > 4 else (0, 1)
    torch.set_num_threads(int(sys.argv[5]) if len(sys.argv) > 5 else 4)
    os.makedirs(out, exist_ok=True)
    st = torch.load(ck, map_location="cpu")
    head, trans, coef = C0BHead(512), g1.G1Transitions(), g1.BoundaryCoef()
    head.load_state_dict(st["head"])
    trans.load_state_dict(st["trans"])
    coef.load_state_dict(st["coef"])
    head.eval()
    T = trans().detach().double().numpy()
    feats = np.load(f"{ROOT}/data/p05_feat/feat_{sp}.npy", mmap_mode="r")
    codes = np.load(f"{ROOT}/data/p1s_v12/{sp}_codes.npy", mmap_mode="r")
    rows = np.arange(len(codes))[i::n]
    res = np.zeros((len(rows), codes.shape[1], 8), np.float32)
    t0 = time.time()
    for r, w in enumerate(rows):
        with torch.no_grad():
            p, b = head(torch.as_tensor(np.asarray(feats[w:w + 1], np.float32)),
                        torch.as_tensor(np.asarray(codes[w:w + 1]).astype(np.int64)))
            edges = g1.boundary_edges(b[0].double(), coef.double())
        c = np.asarray(codes[w])
        post, _ = g1.marginals_np(p[0].double().numpy(), T, g1.gate_mask(c), edges=edges)
        res[r] = g2cand.scores(post, b[0].numpy())
        if r % 50 == 0:
            print(sp, f"shard {i}/{n}", r, len(rows), round(time.time() - t0), flush=True)
    name = f"{sp}_cand_scores" + (f".part{i}" if n > 1 else "")
    np.save(f"{out}/{name}.npy", res)
    print("DONE", name, round(time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
