#!/usr/bin/env python3
"""P1.1 feature cache on ont (protocol v1 section 2): GENERanno-1.2b-cds layers 17+23, fp16, bf16 compute.

data/p1_feat/train.npy [3000, 4096, 4096]; chr18.npy, chr20.npy [len//6, 4096]; meta.json.
Resumable: each array has a .done marker per window / tile; rerun skips finished parts.
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ncrna_annot.features import generanno as G  # noqa: E402

ROOT = os.path.expanduser("~/ncrna-annot")
D, O = f"{ROOT}/data/p1", f"{ROOT}/data/p1_feat"


def open_arr(path, shape):
    return np.lib.format.open_memmap(path, "r+" if os.path.exists(path) else "w+", np.float16, shape)


def main():
    os.makedirs(O, exist_ok=True)
    tok, model = G.load()
    t0 = time.time()
    codes = np.load(f"{D}/train_codes.npy", mmap_mode="r")
    arr = open_arr(f"{O}/train.npy", (len(codes), G.NT, G.DIM))
    done = f"{O}/train.done"
    start = int(open(done).read()) if os.path.exists(done) else 0
    for i in range(start, len(codes)):
        arr[i] = G.features(tok, model, [codes[i]])[0].cpu().numpy()
        if (i + 1) % 100 == 0:
            arr.flush()
            open(done, "w").write(str(i + 1))
            print(json.dumps({"train": i + 1, "t": round(time.time() - t0)}), flush=True)
    arr.flush()
    open(done, "w").write(str(len(codes)))
    for name in ("chr18", "chr20"):
        c = np.load(f"{D}/{name}_codes.npy", mmap_mode="r")
        starts = G.chromosome_tiles(len(c))
        arr = open_arr(f"{O}/{name}.npy", (len(c) // G.TOK, G.DIM))
        done = f"{O}/{name}.done"
        k0 = int(open(done).read()) if os.path.exists(done) else 0
        for k in range(k0, len(starts)):
            s = starts[k]
            lo, hi = G.kept_tokens(s, starts)
            arr[s // G.TOK + lo:s // G.TOK + hi] = G.features(tok, model, [c[s:s + G.WIN]])[0, lo:hi].cpu().numpy()
            if (k + 1) % 200 == 0:
                arr.flush()
                open(done, "w").write(str(k + 1))
                print(json.dumps({name: k + 1, "of": len(starts), "t": round(time.time() - t0)}), flush=True)
        arr.flush()
        open(done, "w").write(str(len(starts)))
    json.dump({"layers": G.LAYERS, "dim": G.DIM, "dtype": "float16", "compute": str(G.half_dtype()),
               "elapsed_s": round(time.time() - t0)}, open(f"{O}/meta.json", "w"), indent=1)


if __name__ == "__main__":
    main()
