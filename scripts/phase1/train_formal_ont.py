#!/usr/bin/env python3
"""Phase 1 formal-feature runs on ont (P1.1 cache: GENERanno L17+L23, 4,096-d): G1B training + chr18 evaluation.

  train_formal_ont.py --name NAME --seed S [--arm G1B|R|G2PH] [--labels v50uev] [--lnc-only] [--epochs 12] [--smoke]
  train_formal_ont.py --name ab1-s1-42 --eval-only --from g1b-42 --with-s1      (AB1: S1 decoding of a main model)
Arms (protocol v1.8 section 15): G1B main / AB2 (--labels v50f) / AB4 (--labels v41uev) / AB3 (--labels v50uevlnc
--lnc-only); AB0 = --arm R (G1 without boundary head, screen2); AB6 = --arm G2PH (screen_g2ph, numba CPU FB).
Data (ont): data/p1_feat/{train,val,chr18}.npy (memmaps), data/p1v12/{train,val}_{codes,bits_<labels>,events}.npy,
chr18 truth outputs/PHASE1/DIAG_L2/chr18_v50.gtf.  Output: outputs/PHASE1/FORMAL/<NAME>/ (metrics.jsonl, last.pt,
best.pt, result.json, chr18.gff3).  Settings other than features / window set are those of the screening arm
(protocol v1.6 section 13): 12 epochs, AdamW 1e-3 cosine to 1e-5, batch 8, checkpoint by validation NLL.
GPU sharing (rule 2026-10-26): every 5 batches the GPU gate (scripts/phase1/gpu_gate.sh) is checked; if a TE or
another user's job is on the card the run stops with exit code 75 after the last full epoch was saved, and
run_formal_ont.sh waits and restarts it; a CUDA out-of-memory error (another job took the card) is treated the same
way.  Memory: NCRNA_GPU_FRAC (default 0.30) and batch 8 as 2 x 4 gradient accumulation (--accum).
chr20 / T2T chr20 are never read here.
"""
import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
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
    ap.add_argument("--arm", default="G1B", choices=["G1B", "R", "G2PH"])
    ap.add_argument("--lnc-only", action="store_true", help="AB3: mRNA states masked in training and decoding")
    ap.add_argument("--from", dest="src", default=None, help="eval-only: take best.pt from this run name")
    ap.add_argument("--with-s1", action="store_true", help="eval-only: also decode S1 (posterior argmax), AB1")
    ap.add_argument("--accum", type=int, default=int(os.environ.get("NCRNA_ACCUM", "2")),
                    help="batch 8 as accum x (8/accum) micro-batches (shared GPU)")
    ap.add_argument("--online-features", action="store_true", default=os.environ.get("NCRNA_ONLINE") == "1",
                    help="compute backbone features on the fly (win3090: no feature cache on disk)")
    a = ap.parse_args()
    import torch
    frac = float(os.environ.get("NCRNA_GPU_FRAC", "0.30"))
    torch.cuda.set_per_process_memory_fraction(frac, 0)
    from ncrna_annot.train import screen_g1b as sg
    from ncrna_annot.train import screen2 as s2
    from ncrna_annot.train import screen_g2ph as sp
    out = f"{ROOT}/outputs/PHASE1/FORMAL/{a.name}"
    os.makedirs(out, exist_ok=True)
    F, D = f"{ROOT}/data/p1_feat", f"{ROOT}/data/p1v12"
    mm = lambda p: np.load(p, mmap_mode="r")
    tr = dict(feats=mm(f"{F}/train.npy") if os.path.exists(f"{F}/train.npy") else None, codes=mm(f"{D}/train_codes.npy"), bits=mm(f"{D}/train_bits_{a.labels}.npy"),
              events=mm(f"{D}/train_events.npy"))
    va = dict(feats=mm(f"{F}/val.npy") if os.path.exists(f"{F}/val.npy") else None, codes=mm(f"{D}/val_codes.npy"), bits=mm(f"{D}/val_bits_{a.labels}.npy"),
              events=mm(f"{D}/val_events.npy"))
    ch = dict(feats=mm(f"{F}/chr18.npy") if os.path.exists(f"{F}/chr18.npy") else None, codes=mm(f"{ROOT}/data/p1/chr18_codes.npy"),
              truth_gtf=f"{ROOT}/outputs/PHASE1/DIAG_L2/chr18_v50.gtf", chrom="chr18", gff_out=f"{out}/chr18.gff3")
    if a.online_features:  # same features as data/p1_feat (GENERanno L17+L23, bf16 compute, fp16), computed per use
        from ncrna_annot.features import generanno as G
        from ncrna_annot.features import online
        tok, model = G.load()
        fn = online.backbone_fn(tok, model, batch=int(os.environ.get("NCRNA_BB_BATCH", "2")))
        tr["feats"] = online.OnlineWindows(tr["codes"], fn)
        va["feats"] = online.OnlineWindows(va["codes"], fn, cache_all=True)
        ch["feats"] = online.OnlineChromosome(ch["codes"], fn)
    for k in ([] if a.online_features else ("train", "val", "chr18")):  # refuse a partially filled cache
        cov = json.load(open(f"{F}/{k}.cover.json"))
        need = len(tr["codes"]) if k == "train" else len(va["codes"]) if k == "val" else None
        if need is not None:
            got = sum(b - a_ for a_, b in cov)
            assert got == need, f"{k} cache covers {got} of {need} windows"
    if a.smoke:
        tr = {k: np.asarray(v[:64]) for k, v in tr.items()}
        va = {k: np.asarray(v[:8]) for k, v in va.items()}
        ch.update(codes=ch["codes"][:3000000])
        ch["feats"] = (online.OnlineChromosome(ch["codes"], fn) if a.online_features else ch["feats"][:500000])
        a.epochs = 1
    meta = dict(name=a.name, host=os.uname().nodename, online_features=a.online_features, accum=a.accum, arm=a.arm, seed=a.seed, labels=a.labels, lnc_only=a.lnc_only, eval_only=a.eval_only,
                src=a.src, with_s1=a.with_s1, epochs=a.epochs, smoke=a.smoke, gpu_frac=frac,
                features="GENERanno-1.2b-cds L17+L23 fp16 (data/p1_feat)", started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    json.dump(meta, open(f"{out}/run_meta.json", "w"), indent=1)
    n = {"i": 0}

    def pause():
        n["i"] += 1
        if n["i"] % 5:
            return False
        why = gate_busy()
        if why:
            print(json.dumps({"pause": why, "t": time.strftime("%H:%M:%SZ", time.gmtime())}), flush=True)
        return bool(why)

    try:
        if a.eval_only:
            if gate_busy():
                sys.exit(75)
            ck = f"{ROOT}/outputs/PHASE1/FORMAL/{a.src or a.name}/best.pt"
            if a.arm == "G2PH":
                sp.evaluate_checkpoint(ck, "G2PH", 4096, ch, out, device="cuda")
            else:
                sg.evaluate_checkpoint(ck, 4096, ch, out, device="cuda", lnc_only=a.lnc_only, with_s1=a.with_s1)
        elif a.arm == "G1B":
            sg.run(a.seed, tr, va, ch, out, epochs=a.epochs, device="cuda", pause=pause, resume=True, lnc_only=a.lnc_only,
                   accum=a.accum)
        elif a.arm == "R":
            s2.run("R", a.seed, tr, va, ch, out, lnc_only=a.lnc_only, epochs=a.epochs, device="cuda", pause=pause, resume=True,
                   accum=a.accum)
        else:  # G2PH: training (GPU head + numba CPU forward-backward), then evaluation in the same run
            if not os.path.exists(f"{out}/train_summary.json"):
                sp.run("G2PH", a.seed, tr, va, out, epochs=a.epochs, device="cuda", pause=pause, resume=True, accum=a.accum)
            if gate_busy():
                sys.exit(75)
            sp.evaluate_checkpoint(f"{out}/best.pt", "G2PH", 4096, ch, out, device="cuda")
    except sg.Paused:
        sys.exit(75)
    except torch.OutOfMemoryError as e:  # another job (TE / coding) took the card: yield, resume from last.pt
        print(json.dumps({"yield_oom": str(e)[:200]}), flush=True)
        sys.exit(75)
    print("RUN_DONE", flush=True)


if __name__ == "__main__":
    main()
