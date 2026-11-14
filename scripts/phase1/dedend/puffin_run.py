#!/usr/bin/env python3
"""PLAN_DEDEND: Puffin (Dudnyk et al. 2024) FANTOM_CAGE prediction over whole chromosomes, both strands.

The model class below is the forward part of puffin.py (github.com/jzhoulab/puffin, commit 9bc81eb) without the
selene_sdk dependency: one-hot A,C,G,T in this order, any other base = 0.25 in every channel (as
selene_sdk.sequences.sequence_to_encoding), 325 bp context on each side.  Weights resources/puffin.pth, loaded with
strict=True so that a key mismatch cannot pass silently.  Output channels: 0 = FANTOM_CAGE + strand, 9 = FANTOM_CAGE
- strand (puffin.py targeti_map / targeti_rev_map).

  puffin_run.py --name puffin --chroms chr8,chr9,chr18,chr20   -> outputs/PHASE1/DEDEND/puffin/<chrom>.npy  float16 [L, 2]
"""
import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
PUFFIN = os.path.expanduser(os.environ.get("PUFFIN_DIR", "~/tools/dedend/puffin"))
PAD, CHUNK = 325, 2_000_000


def gate_busy():
    r = subprocess.run(["bash", "-c", f". {ROOT}/scripts/phase1/gpu_gate.sh; gpu_why_busy"], capture_output=True, text=True)
    why = r.stdout.strip() if r.returncode == 0 else ""
    return "" if why.startswith("only ") else why


def read_fasta(path):
    with open(path) as fh:
        return "".join(l.strip() for l in fh if not l.startswith(">"))


def onehot(s):
    lut = np.full(256, 4, np.uint8)
    for i, b in enumerate("ACGT"):
        lut[ord(b)] = i
        lut[ord(b.lower())] = i
    idx = lut[np.frombuffer(s.encode(), np.uint8)]
    x = np.full((len(s), 4), 0.25, np.float32)
    k = idx < 4
    x[k] = 0.0
    x[np.flatnonzero(k), idx[k]] = 1.0
    return x


def build(dev):
    import torch
    from torch import nn
    from torch_fftconv import FFTConv1d

    class Puffin(nn.Module):
        def __init__(self):
            super().__init__()
            self.conv = nn.Conv1d(4, 10, kernel_size=51, padding=25)
            self.conv_inr = nn.Conv1d(4, 10, kernel_size=15, padding=7)
            self.conv_sim = nn.Conv1d(4, 32, kernel_size=3, padding=1)
            self.activation = nn.Softplus()
            self.softplus = nn.Softplus()
            self.deconv = FFTConv1d(10 * 2, 10, kernel_size=601, padding=300)
            self.deconv_sim = FFTConv1d(64, 10, kernel_size=601, padding=300)
            self.deconv_inr = nn.ConvTranspose1d(20, 10, kernel_size=15, padding=7)
            self.scaler = nn.Parameter(torch.ones(1))
            self.scaler2 = nn.Parameter(torch.ones(1))

        def forward(self, x):
            y = torch.cat([self.conv(x), self.conv(x.flip([1, 2])).flip([2])], 1)
            y_sim = torch.cat([self.conv_sim(x), self.conv_sim(x.flip([1, 2])).flip([2])], 1)
            y_inr = torch.cat([self.conv_inr(x), self.conv_inr(x.flip([1, 2])).flip([2])], 1)
            return self.softplus(self.deconv(self.activation(y)) + self.deconv_inr(self.activation(y_inr))
                                 + self.deconv_sim(self.activation(y_sim)))

    m = Puffin()
    sd = torch.load(f"{PUFFIN}/resources/puffin.pth", map_location="cpu")
    missing, unexpected = m.load_state_dict(sd, strict=False)
    assert not missing, missing   # every parameter of the forward model must come from the file
    print(json.dumps({"puffin_weights": len(sd), "unused_keys": sorted(unexpected)}), flush=True)
    return m.to(dev).eval()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--chroms", default="chr8,chr9,chr18,chr20")
    a = ap.parse_args()
    import torch
    torch.cuda.set_per_process_memory_fraction(float(os.environ.get("NCRNA_GPU_FRAC", "0.30")), 0)
    dev = torch.device("cuda")
    out = f"{ROOT}/outputs/PHASE1/DEDEND/{a.name}"
    os.makedirs(out, exist_ok=True)
    model = build(dev)
    for c in a.chroms.split(","):
        dst = f"{out}/{c}.npy"
        if os.path.exists(dst + ".done"):
            continue
        seq = read_fasta(f"{ROOT}/data/human/{c}.fa")
        L = len(seq)
        res = np.lib.format.open_memmap(dst, mode="r+" if os.path.exists(dst) else "w+", dtype=np.float16, shape=(L, 2))
        prog = f"{dst}.progress"
        s0 = int(open(prog).read()) if os.path.exists(prog) else 0
        t0 = time.time()
        for s in range(s0, L, CHUNK):
            why = gate_busy()
            if why:
                res.flush()
                print(json.dumps({"yield_gpu": why, "chrom": c, "at": s}), flush=True)
                sys.exit(75)
            e = min(s + CHUNK, L)
            a0, b0 = max(0, s - PAD), min(L, e + PAD)
            x = onehot(seq[a0:b0])
            x = np.concatenate([np.full((PAD - (s - a0), 4), 0.25, np.float32), x,
                                np.full((PAD - (b0 - e), 4), 0.25, np.float32)])
            with torch.no_grad():
                y = model(torch.as_tensor(x.T[None], device=dev))[0].float().cpu().numpy()
            y = y[:, PAD:PAD + (e - s)]
            res[s:e, 0] = y[0].astype(np.float16)
            res[s:e, 1] = y[9].astype(np.float16)
            res.flush()
            open(prog, "w").write(str(e))
        open(dst + ".done", "w").write(json.dumps({"L": L, "s": round(time.time() - t0)}))
        print(json.dumps({"chrom": c, "L": L, "s": round(time.time() - t0)}), flush=True)


if __name__ == "__main__":
    main()
