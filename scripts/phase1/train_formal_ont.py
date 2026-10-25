#!/usr/bin/env python3
"""Phase 1 formal-feature runs on ont (P1.1 cache: GENERanno L17+L23, 4,096-d): G1B training + chr18 evaluation.

  train_formal_ont.py --name NAME --seed S [--labels v50uev] [--epochs 12] [--smoke] [--eval-only]
Data (ont): data/p1_feat/{train,val,chr18}.npy (memmaps), data/p1v12/{train,val}_{codes,bits_<labels>,events}.npy,
chr18 truth outputs/PHASE1/DIAG_L2/chr18_v50.gtf.  Output: outputs/PHASE1/FORMAL/<NAME>/ (metrics.jsonl, last.pt,
best.pt, result.json, chr18.gff3).  Settings other than features / window set are those of the screening arm
(protocol v1.6 section 13): 12 epochs, AdamW 1e-3 cosine to 1e-5, batch 8, checkpoint by validation NLL.
GPU sharing (rule 2026-10-26): every 25 batches the GPU gate (scripts/phase1/gpu_gate.sh) is checked; if a TE or
another user's job is on the card the run stops with exit code 75 after the last full epoch was saved, and
run_formal_ont.sh waits and restarts it.  The process memory is capped by NCRNA_GPU_FRAC (default 0.45).
chr20 / T2T chr20 are never read here.
"""
import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np

ROOT = os.path.expanduser("~/ncrna-annot")
sys.path.insert(0, os.path.join(ROOT, "src"))


def gate_busy():
    r = subprocess.run(["bash", "-c", f". {ROOT}/scripts/phase1/gpu_gate.sh; gpu_why_busy"], capture_output=True, text=True)
    return r.stdout.strip() if r.returncode == 0 else ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--labels", default="v50uev", choices=["v50uev", "v50f", "v41uev", "v50uevlnc"])
    ap.add_argument("--epochs", type=int, default=12)
    ap.add_argument("--smoke", action="store_true", help="64 train / 8 val windows, 1 epoch, first 3 Mb of chr18")
    ap.add_argument("--eval-only", action="store_true")
    a = ap.parse_args()
    import torch
    frac = float(os.environ.get("NCRNA_GPU_FRAC", "0.45"))
    torch.cuda.set_per_process_memory_fraction(frac, 0)
    from ncrna_annot.train import screen_g1b as sg
    out = f"{ROOT}/outputs/PHASE1/FORMAL/{a.name}"
    os.makedirs(out, exist_ok=True)
    F, D = f"{ROOT}/data/p1_feat", f"{ROOT}/data/p1v12"
    mm = lambda p: np.load(p, mmap_mode="r")
    tr = dict(feats=mm(f"{F}/train.npy"), codes=mm(f"{D}/train_codes.npy"), bits=mm(f"{D}/train_bits_{a.labels}.npy"),
              events=mm(f"{D}/train_events.npy"))
    va = dict(feats=mm(f"{F}/val.npy"), codes=mm(f"{D}/val_codes.npy"), bits=mm(f"{D}/val_bits_{a.labels}.npy"),
              events=mm(f"{D}/val_events.npy"))
    ch = dict(feats=mm(f"{F}/chr18.npy"), codes=mm(f"{ROOT}/data/p1/chr18_codes.npy"),
              truth_gtf=f"{ROOT}/outputs/PHASE1/DIAG_L2/chr18_v50.gtf", chrom="chr18", gff_out=f"{out}/chr18.gff3")
    for k in ("train", "val", "chr18"):  # refuse to train on a partially filled cache
        cov = json.load(open(f"{F}/{k}.cover.json"))
        need = len(tr["codes"]) if k == "train" else len(va["codes"]) if k == "val" else None
        if need is not None:
            got = sum(b - a_ for a_, b in cov)
            assert got == need, f"{k} cache covers {got} of {need} windows"
    if a.smoke:
        tr = {k: np.asarray(v[:64]) for k, v in tr.items()}
        va = {k: np.asarray(v[:8]) for k, v in va.items()}
        ch.update(feats=ch["feats"][:500000], codes=ch["codes"][:3000000])
        a.epochs = 1
    meta = dict(name=a.name, seed=a.seed, labels=a.labels, epochs=a.epochs, smoke=a.smoke, gpu_frac=frac,
                features="GENERanno-1.2b-cds L17+L23 fp16 (data/p1_feat)", started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    json.dump(meta, open(f"{out}/run_meta.json", "w"), indent=1)
    n = {"i": 0}

    def pause():
        n["i"] += 1
        if n["i"] % 25:
            return False
        why = gate_busy()
        if why:
            print(json.dumps({"pause": why, "t": time.strftime("%H:%M:%SZ", time.gmtime())}), flush=True)
        return bool(why)

    try:
        if a.eval_only:
            if gate_busy():
                sys.exit(75)
            sg.evaluate_checkpoint(f"{out}/best.pt", 4096, ch, out, device="cuda")
        else:
            sg.run(a.seed, tr, va, ch, out, epochs=a.epochs, device="cuda",
                   pause=pause, resume=True)
    except sg.Paused:
        sys.exit(75)
    print("RUN_DONE", flush=True)


if __name__ == "__main__":
    main()
