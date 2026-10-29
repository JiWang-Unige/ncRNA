#!/usr/bin/env python3
"""v1.9 development (chr18 only, CPU): decode saved chr18 posteriors with S1 / MEA and score them like the formal runs.

  dev_decode.py --post outputs/PHASE1/FORMAL/<run>/chr18_post.npy --mea 1,x1,x2,pv1 [--s1] --out <json>
([pv][x]gamma as in g1.mea_duration / screen2.mea_specs)
Posteriors come from train_formal_ont.py --eval-only ... --dump-post (same model, same tiling as the run's own
evaluation); --post a.npy,b.npy,... decodes their per-base mean (posterior ensemble, saved as <out>.post.npy).  The grammar's allowed transitions do not depend on the trained parameters, so the default
G1Transitions() gives the structure for mea_duration.  chr20 / T2T chr20 are never read here.
"""
import argparse
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--post", required=True)
    ap.add_argument("--mea", default="1")
    ap.add_argument("--s1", action="store_true")
    ap.add_argument("--codes", default=f"{ROOT}/data/p1/chr18_codes.npy")
    ap.add_argument("--truth", default=f"{ROOT}/outputs/PHASE1/DIAG_L2/chr18_v50.gtf")
    ap.add_argument("--chrom", default="chr18")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    assert a.chrom == "chr18", "development decoders run on chr18 only"
    from ncrna_annot.decode import chrom
    from ncrna_annot.eval import ruler
    from ncrna_annot.grammar import g1
    from ncrna_annot.train import screen2 as s2
    paths = a.post.split(",")
    post = np.load(paths[0], mmap_mode="r")
    if len(paths) > 1:  # posterior ensemble: per-base mean of several checkpoints' posteriors (float16 on disk)
        parts = [post] + [np.load(p, mmap_mode="r") for p in paths[1:]]
        assert all(p.shape == post.shape for p in parts)
        post = np.lib.format.open_memmap(a.out + ".post.npy", "w+", np.float16, post.shape)
        for i in range(0, len(post), 1 << 22):
            post[i:i + (1 << 22)] = sum(np.asarray(p[i:i + (1 << 22)], np.float32) for p in parts) / len(parts)
    codes = np.load(a.codes, mmap_mode="r")[:len(post)]
    T = g1.G1Transitions()().detach().double().numpy()
    truth = ruler.read_gtf(a.truth, a.chrom)
    res = {"post": a.post}
    modes = ([("S1", "argmax", 1.0, "state", False)] if a.s1 else []) + [(n, "mea", g, lv, lg) for n, g, lv, lg in s2.mea_specs(a.mea)]
    for name, mode, gm, lv, lg in modes:
        t0 = time.time()
        arr = np.asarray(post, np.float32) if mode == "argmax" else post
        genes = chrom.decode_blocks(arr, codes, T, 20, 50, mode=mode, gamma=gm, level=lv, log_gain=lg)
        res[name] = s2.summarize(genes, post, truth, s2.CONF_GRID)
        res[name]["decode_s"] = round(time.time() - t0)
        print(json.dumps({name: {k: v for k, v in res[name].items() if k != "curve"}}), flush=True)
        json.dump(res, open(a.out, "w"), indent=1)


if __name__ == "__main__":
    main()
