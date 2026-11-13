#!/usr/bin/env python3
"""PLAN_END.md section 1 (stage 0): raw boundary logits (g1.BOUNDARY order TSS+, TES+, TSS-, TES-) of the v2.0 main
models (uniform average of epochs 10-12) on chr18, forward in 8 x 24,576 bp pieces as screen_g1b.evaluate_chromosome.
No decoding.  -> outputs/PHASE1/ENDS/<name>/<run>_bnd.npy (float16 [L, 4]).
Run through scripts/phase1/run_formal_ont.sh --script ends/dump_bnd.py.

  dump_bnd.py --name stage0 [--runs v20-g1b-42-ke,v20-g1b-123-ke,v20-g1b-456-ke]
"""
import argparse
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
TOK, WIN = 6, 24576


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--runs", default="v20-g1b-42-ke,v20-g1b-123-ke,v20-g1b-456-ke")
    a = ap.parse_args()
    import torch
    from ncrna_annot.heads.c0 import C0BHead
    from ncrna_annot.train.screen2 import load_ckpt
    torch.cuda.set_per_process_memory_fraction(float(os.environ.get("NCRNA_GPU_FRAC", "0.30")), 0)
    dev = torch.device("cuda")
    out = f"{ROOT}/outputs/PHASE1/ENDS/{a.name}"
    os.makedirs(out, exist_ok=True)
    cf = np.load(f"{ROOT}/data/p1_feat/chr18.npy", mmap_mode="r")
    cc = np.load(f"{ROOT}/data/p1/chr18_codes.npy", mmap_mode="r")
    L = len(cf) * TOK
    for run in a.runs.split(","):
        dst = f"{out}/{run}_bnd.npy"
        if os.path.exists(dst + ".done"):
            continue
        t0 = time.time()
        head = C0BHead(cf.shape[-1]).to(dev)
        ck = load_ckpt([f"{ROOT}/outputs/PHASE1/FORMAL/{run}/ep{e}.pt" for e in (10, 11, 12)], dev)
        head.load_state_dict(ck["head"])
        head.eval()
        res = np.lib.format.open_memmap(dst, mode="w+", dtype=np.float16, shape=(L, 4))
        with torch.no_grad():
            for s in range(0, L, WIN * 8):
                e = min(s + WIN * 8, L)
                f = torch.as_tensor(np.asarray(cf[s // TOK:e // TOK], np.float32), device=dev)[None]
                c = torch.as_tensor(np.asarray(cc[s:e], np.int64), device=dev)[None]
                _, b = head(f, c)
                res[s:e] = b[0].float().cpu().numpy().astype(np.float16)
        res.flush()
        open(dst + ".done", "w").write(json.dumps({"epochs": ck["epoch"], "L": L, "s": round(time.time() - t0)}))
        print(json.dumps({"run": run, "L": L, "s": round(time.time() - t0)}), flush=True)
        del head
        torch.cuda.empty_cache()


if __name__ == "__main__":
    main()
