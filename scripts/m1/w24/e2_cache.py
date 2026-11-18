#!/usr/bin/env python3
"""M1 weeks 2-4, E2 (evidence in the model): cache, per v2.0 run and per 24,576-bp tile, the head's emission potentials
[WIN, S] and boundary logits [WIN, 4] (float32, exactly what p15_read.posteriors feeds to the exact forward-backward),
plus the run's transitions and boundary coefficients.  GPU once per run (through the GPU gate); every evidence condition
is then CPU only (e2_decode.py).
  e2_cache.py --target chr18 [--runs ...] -> outputs/M1/W24/cache/<target>/<run>/{emit.npy, bnd.npy, meta.json}
"""
import argparse
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1"))
import p15_read as R  # noqa: E402

CFG = f"{ROOT}/configs/v20.json"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", default="chr18")
    ap.add_argument("--runs", default="")
    a = ap.parse_args()
    import torch
    torch.cuda.set_per_process_memory_fraction(float(os.environ.get("NCRNA_GPU_FRAC", "0.30")), 0)
    from ncrna_annot.features import generanno as G
    from ncrna_annot.grammar import g1
    from ncrna_annot.heads.c0 import C0BHead
    from ncrna_annot.train import screen2 as s2
    from ncrna_annot.train import screen_g1b as sg
    cfg = json.load(open(CFG))
    R.check_code(cfg)
    feats, codes, tg = R.load_target(cfg, a.target)
    L = len(feats) * sg.TOK
    starts = G.chromosome_tiles(L)
    runs = [r for r in a.runs.split(",") if r] or cfg["targets"][a.target]["runs"]
    for name in runs:
        out = f"{ROOT}/outputs/M1/W24/cache/{a.target}/{name}"
        if os.path.exists(f"{out}/meta.json"):
            print(json.dumps({"skip": name}), flush=True)
            continue
        if R.gate_busy():
            sys.exit(75)
        os.makedirs(out, exist_ok=True)
        run = cfg["runs"][name]
        assert run["arm"] == "G1B"
        bad = [p for p, h in zip(run["ckpts"], run["ckpt_sha256"]) if R.sha256(f"{ROOT}/{p}") != h]
        assert not bad, bad
        ck = s2.load_ckpt([f"{ROOT}/{p}" for p in run["ckpts"]], "cuda")
        head, coef = C0BHead(feats.shape[-1]).to("cuda"), g1.BoundaryCoef().to("cuda")
        coef.load_state_dict(ck["coef"])
        head.load_state_dict(ck["head"])
        head.eval()
        trans = g1.G1Transitions().to("cuda")
        trans.load_state_dict(ck["trans"])
        T = trans().detach().double().cpu().numpy()
        em = np.lib.format.open_memmap(f"{out}/emit.npy", "w+", np.float32, (len(starts), G.WIN, g1.S))
        bd = np.lib.format.open_memmap(f"{out}/bnd.npy", "w+", np.float32, (len(starts), G.WIN, 4))
        t0 = time.time()
        for k, s in enumerate(starts):
            f = torch.as_tensor(np.asarray(feats[s // sg.TOK:(s + G.WIN) // sg.TOK], np.float32), device="cuda")[None]
            c = torch.as_tensor(np.asarray(codes[s:s + G.WIN], np.int64), device="cuda")[None]
            with torch.no_grad():
                p, b = head(f, c)
                p = sg._mask_lnc_only(p, run["lnc_only"])[0]
            em[k] = p.float().cpu().numpy()
            bd[k] = b[0].float().cpu().numpy()
        em.flush()
        bd.flush()
        np.save(f"{out}/T.npy", T)
        json.dump({"run": name, "target": a.target, "L": L, "starts": starts, "coef": {k: v.detach().cpu().tolist()
                   for k, v in coef.c.items()}, "lnc_only": run["lnc_only"], "s": round(time.time() - t0)},
                  open(f"{out}/meta.json", "w"))
        print(json.dumps({"done": name, "s": round(time.time() - t0), "tiles": len(starts)}), flush=True)


if __name__ == "__main__":
    main()
