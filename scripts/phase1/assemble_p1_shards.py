#!/usr/bin/env python3
"""Merge feature shards (extract_p1_shard.py) into the P1.1 cache data/p1_feat/ on ont, with checks.

  assemble_p1_shards.py --shards DIR [--feat-root ~/ncrna-annot/data/p1_feat] [--codes-dir ~/ncrna-annot/data/p1]
                        [--ref-dir DIR2] [--check-only]

- verifies each shard's sha256 against its .json
- writes rows into the same memmaps extract_p1.py produces (train.npy [N, 4096, 4096]; <chrom>.npy [len//6, 4096])
- keeps <name>.cover.json (which windows / tiles are filled) so partial caches are explicit
- --ref-dir: for every shard with a same-stem shard in DIR2 (e.g. ont bf16 vs Colab fp16), prints cosine stats;
  protocol v1.4 consistency rule: median cosine >= 0.99, else stop and investigate before any read-out.
"""
import argparse
import glob
import hashlib
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ncrna_annot.features import generanno as G  # noqa: E402

ROOT = os.path.expanduser("~/ncrna-annot")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 24), b""):
            h.update(blk)
    return h.hexdigest()


def cosine(a, b):
    a = np.asarray(a, np.float32).reshape(-1, a.shape[-1])
    b = np.asarray(b, np.float32).reshape(-1, b.shape[-1])
    c = (a * b).sum(1) / (np.sqrt((a * a).sum(1) * (b * b).sum(1)) + 1e-12)
    return {"cos_median": float(np.median(c)), "cos_p01": float(np.percentile(c, 1)), "cos_min": float(c.min())}


def open_arr(path, shape):
    return np.lib.format.open_memmap(path, "r+" if os.path.exists(path) else "w+", np.float16, shape)


def cover(path):
    return json.load(open(path)) if os.path.exists(path) else []


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shards", required=True)
    ap.add_argument("--feat-root", default=f"{ROOT}/data/p1_feat")
    ap.add_argument("--codes-dir", default=f"{ROOT}/data/p1")
    ap.add_argument("--ref-dir", default=None)
    ap.add_argument("--check-only", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.feat_root, exist_ok=True)
    ok = True
    for mp in sorted(glob.glob(os.path.join(a.shards, "*.json"))):
        m = json.load(open(mp))
        npy = mp[:-5] + ".npy"
        if not os.path.exists(npy):
            print("MISSING", npy); ok = False; continue
        digest = sha256(npy)
        if digest != m["sha256"]:
            print("SHA256 MISMATCH", npy); ok = False; continue
        rep = {"stem": m["stem"], "kind": m["kind"], "compute": m.get("compute"), "gpu": m.get("gpu"), "sha_ok": True}
        if a.ref_dir:
            rp = os.path.join(a.ref_dir, m["stem"] + ".npy")
            if os.path.exists(rp):
                rep["ref"] = cosine(np.load(npy, mmap_mode="r"), np.load(rp, mmap_mode="r"))
                rep["ref_pass"] = rep["ref"]["cos_median"] >= 0.99
                ok = ok and rep["ref_pass"]
        if not a.check_only:
            sh = np.load(npy, mmap_mode="r")
            if m["kind"] == "train":
                wname = m.get("name", "train")
                n = m.get("n_total") or len(np.load(os.path.join(a.codes_dir, f"{wname}_codes.npy"), mmap_mode="r"))
                arr = open_arr(os.path.join(a.feat_root, f"{wname}.npy"), (n, G.NT, G.DIM))
                ga, gb = m["index_range"]
                arr[ga:gb] = sh
                arr.flush()
                cp = os.path.join(a.feat_root, f"{wname}.cover.json")
                cv = cover(cp); cv.append([ga, gb]); json.dump(sorted(cv), open(cp, "w"))
                rep["wrote"] = [ga, gb]
            else:
                name, L = m["name"], m["chrom_len"]
                arr = open_arr(os.path.join(a.feat_root, f"{name}.npy"), (L // G.TOK, G.DIM))
                for k, s, lo, hi, r0, r1 in m["tiles"]:
                    arr[s // G.TOK + lo:s // G.TOK + hi] = sh[r0:r1]
                arr.flush()
                cp = os.path.join(a.feat_root, f"{name}.cover.json")
                cv = cover(cp); cv.append(m["tile_range"]); json.dump(sorted(cv), open(cp, "w"))
                rep["wrote_tiles"] = m["tile_range"]
        print(json.dumps(rep), flush=True)
    print("ASSEMBLE_OK" if ok else "ASSEMBLE_FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
