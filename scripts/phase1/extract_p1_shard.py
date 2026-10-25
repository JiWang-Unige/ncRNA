#!/usr/bin/env python3
"""P1.1 features as self-describing shards (runs on Colab / Kaggle / ont GPU):
GENERanno-1.2b-cds hidden layers 17+23 (4,096-d), stored fp16, computed bf16 (Ampere+) or fp16 (T4).

  train windows   extract_p1_shard.py --kind train --codes train_codes.npy --range A:B --out DIR [--index-offset I]
  chromosome      extract_p1_shard.py --kind chrom --name chr18 --codes chr18_codes.npy --tiles K0:K1 --out DIR
                  [--chrom-len L --bp-offset O]   (codes may be a slice of the chromosome starting at bp O)

Each shard = <stem>.npy (fp16) + <stem>.json (sha256, shape, tile index, GPU, compute dtype, timings,
cosine vs --ref if given).  Stems carry global indices (train_00000_00020, chr18_t00000_00010), so a
shard made on any machine drops into data/p1_feat/ via assemble_p1_shards.py.  Protocol v1 section 2.
"""
import argparse
import hashlib
import json
import os
import sys
import time


def parse():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kind", choices=["train", "chrom"], required=True)
    ap.add_argument("--codes", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--src", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
    ap.add_argument("--range", default=None, help="train: A:B window indices within --codes")
    ap.add_argument("--index-offset", type=int, default=0, help="train: global index of codes[0]")
    ap.add_argument("--name", default=None, help="chrom: chromosome name; train: window-set name (default train)")
    ap.add_argument("--n-total", type=int, default=None, help="train: size of the full window set (default len(codes) when --index-offset is 0, else read by the assembler from <codes-dir>/<name>_codes.npy)")
    ap.add_argument("--tiles", default=None, help="chrom: K0:K1 tile indices (tiles of the full chromosome)")
    ap.add_argument("--chrom-len", type=int, default=None, help="chrom: full chromosome length (bp)")
    ap.add_argument("--bp-offset", type=int, default=0, help="chrom: bp position of codes[0] in the chromosome")
    ap.add_argument("--batch", type=int, default=1)
    ap.add_argument("--ref", default=None, help="shard .npy from another machine to compare against (cosine)")
    ap.add_argument("--cpu-stub", action="store_true", help="tests only: random features, no model")
    return ap.parse_args()


def cosine_stats(a, b):
    import numpy as np
    a = np.asarray(a, np.float32).reshape(-1, a.shape[-1])
    b = np.asarray(b, np.float32).reshape(-1, b.shape[-1])
    num = (a * b).sum(1)
    den = np.sqrt((a * a).sum(1) * (b * b).sum(1)) + 1e-12
    c = num / den
    return {"cos_median": float(np.median(c)), "cos_p01": float(np.percentile(c, 1)), "cos_min": float(c.min()),
            "n_tokens": int(c.size)}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 24), b""):
            h.update(blk)
    return h.hexdigest()


def main():
    a = parse()
    sys.path.insert(0, a.src)
    import numpy as np
    from ncrna_annot.features import generanno as G
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    if a.cpu_stub:  # tests: deterministic random features, no torch / model
        tok, model, compute, gpu = None, None, "stub", "cpu"
        rng = np.random.default_rng(0)

        class _Arr(np.ndarray):
            def cpu(self):
                return self

            def numpy(self):
                return np.asarray(self)

        def feats(batch):
            return rng.standard_normal((len(batch), G.NT, G.DIM), dtype=np.float32).astype(np.float16).view(_Arr)
    else:
        import torch
        frac = float(os.environ.get("NCRNA_GPU_FRAC", "0") or 0)
        if frac > 0:  # share the card with the user's other (non-TE) jobs: cap this process
            torch.cuda.set_per_process_memory_fraction(frac, 0)
        tok, model = G.load()
        compute, gpu = str(G.half_dtype()), torch.cuda.get_device_name(0)

        def feats(batch):
            return G.features(tok, model, batch)
    load_s = round(time.time() - t0, 1)
    codes = np.load(a.codes, mmap_mode="r")
    meta = {"kind": a.kind, "layers": list(G.LAYERS), "dim": G.DIM, "dtype": "float16", "compute": compute,
            "gpu": gpu, "model_load_s": load_s, "host": os.uname().nodename}
    if a.kind == "train":
        A, B = map(int, a.range.split(":"))
        ga, gb = A + a.index_offset, B + a.index_offset
        wname = a.name or "train"
        stem = f"{wname}_{ga:05d}_{gb:05d}"
        arr = np.empty((B - A, G.NT, G.DIM), np.float16)
        t1 = time.time()
        for i in range(A, B, a.batch):
            j = min(i + a.batch, B)
            arr[i - A:j - A] = feats([codes[k] for k in range(i, j)]).cpu().numpy()
            print(json.dumps({"win": j - A, "of": B - A, "t": round(time.time() - t0)}), flush=True)
        meta.update(name=wname, n_total=a.n_total or (len(codes) if a.index_offset == 0 else None), index_range=[ga, gb],
                    shape=list(arr.shape), window_s=round((time.time() - t1) / (B - A), 2))
    else:
        K0, K1 = map(int, a.tiles.split(":"))
        L = a.chrom_len or (len(codes) + a.bp_offset)
        starts = G.chromosome_tiles(L)
        tiles, rows = [], 0
        for k in range(K0, K1):
            s = starts[k]
            lo, hi = G.kept_tokens(s, starts)
            tiles.append([k, int(s), lo, hi, rows, rows + hi - lo])
            rows += hi - lo
        arr = np.empty((rows, G.DIM), np.float16)
        t1 = time.time()
        for n, (k, s, lo, hi, r0, r1) in enumerate(tiles):
            c = codes[s - a.bp_offset:s - a.bp_offset + G.WIN]
            assert len(c) == G.WIN, (k, s, a.bp_offset, len(codes))
            arr[r0:r1] = feats([c])[0, lo:hi].cpu().numpy()
            print(json.dumps({"tile": n + 1, "of": len(tiles), "t": round(time.time() - t0)}), flush=True)
        stem = f"{a.name}_t{K0:05d}_{K1:05d}"
        meta.update(name=a.name, chrom_len=L, n_tiles_total=len(starts), tile_range=[K0, K1],
                    tiles=tiles, shape=list(arr.shape), tile_s=round((time.time() - t1) / max(1, len(tiles)), 2))
    path = os.path.join(a.out, stem + ".npy")
    np.save(path, arr)
    meta.update(stem=stem, sha256=sha256(path), bytes=os.path.getsize(path), elapsed_s=round(time.time() - t0, 1))
    if a.ref:
        ref = np.load(a.ref, mmap_mode="r")
        assert ref.shape == arr.shape, (ref.shape, arr.shape)
        meta["ref"] = {"path": a.ref, **cosine_stats(arr, ref)}
    json.dump(meta, open(os.path.join(a.out, stem + ".json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in meta.items() if k != "tiles"}), flush=True)
    print("total_s", meta["elapsed_s"], flush=True)


if __name__ == "__main__":
    main()
