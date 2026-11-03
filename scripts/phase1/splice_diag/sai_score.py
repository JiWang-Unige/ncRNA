#!/usr/bin/env python3
"""chr18 splice diagnostic, step 1: SpliceAI scores for every base of chr18 on both strands, and how well they
recover the chr18 truth splice sites (SpliceAI was trained on chr18: an in-training check of the port, not a
measure of generalisation).

  sai_score.py --codes data/p1/chr18_codes.npy --truth outputs/PHASE1/DIAG_L2/chr18_v50.gtf --chrom chr18 \
      --out outputs/PHASE1/SPLICE_DIAG/sai_chr18.npy

Writes <out> float16 [L, 4] = (acceptor +, donor +, acceptor -, donor -) at genomic 0-based positions
(acceptor = first exon base, donor = last exon base, in transcript orientation) and <out>.json (timings, port check,
truth-site statistics).  For a truth intron (s, e) (1-based, inclusive) the sites are at 0-based s - 2 (left: donor
on +, acceptor on -) and e (right: acceptor on +, donor on -).
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import torch

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sai_torch  # noqa: E402
from ncrna_annot.eval import ruler  # noqa: E402

COL = {("+", "L"): 1, ("+", "R"): 0, ("-", "L"): 2, ("-", "R"): 3}  # (strand, intron side) -> column (donor / acceptor)


def truth_sites(truth, cls):
    """{column: sorted 0-based site positions} over all transcripts of the class's genes."""
    out = {c: set() for c in range(4)}
    for g in truth.values():
        if g["cls"] != cls:
            continue
        for ex in g["tx"].values():
            for s, e in ruler.introns(ruler._merge(ex)):
                out[COL[(g["strand"], "L")]].add(s - 2)
                out[COL[(g["strand"], "R")]].add(e)
    return {c: np.array(sorted(v), np.int64) for c, v in out.items()}


def site_stats(sc, sites):
    res = {}
    for c, pos in sites.items():
        if not len(pos):
            continue
        v = sc[pos, c].astype(np.float32)
        k = len(pos)
        top = np.argpartition(-sc[:, c].astype(np.float32), k)[:k]
        res[["acc+", "don+", "acc-", "don-"][c]] = {
            "n": k, "median": float(np.median(v)), "frac_ge_0.5": float((v >= 0.5).mean()), "frac_ge_0.1": float((v >= 0.1).mean()),
            "topk_acc": float(np.isin(top, pos).mean()),
            "offset_mean": {str(o): float(sc[np.clip(pos + o, 0, len(sc) - 1), c].astype(np.float32).mean()) for o in (-2, -1, 0, 1, 2)}}
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--codes", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--export", default="ext/spliceai/export")
    a = ap.parse_args()
    torch.backends.cudnn.allow_tf32 = False
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.cuda.set_per_process_memory_fraction(float(os.environ.get("NCRNA_GPU_FRAC", "0.30")), 0)
    exp = f"{ROOT}/{a.export}"
    meta = {"port_max_abs_diff": sai_torch.check(exp), "started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    assert meta["port_max_abs_diff"] < 1e-4, meta
    codes = np.load(f"{ROOT}/{a.codes}", mmap_mode="r")
    models = sai_torch.load_models(exp)
    t0 = time.time()
    sc = np.zeros((len(codes), 4), np.float16)
    sc[:, 0:2] = sai_torch.score(models, codes, "+")
    meta["plus_s"] = round(time.time() - t0)
    sc[:, 2:4] = sai_torch.score(models, codes, "-")
    meta["score_s"] = round(time.time() - t0)
    out = f"{ROOT}/{a.out}"
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.save(out, sc)
    truth = ruler.read_gtf(f"{ROOT}/{a.truth}", a.chrom)
    meta["sites"] = {cls: site_stats(sc, truth_sites(truth, cls)) for cls in ("mrna", "lnc")}
    json.dump(meta, open(out + ".json", "w"), indent=1)
    print(json.dumps(meta), flush=True)


if __name__ == "__main__":
    main()
