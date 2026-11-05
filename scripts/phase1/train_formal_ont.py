#!/usr/bin/env python3
"""Phase 1 formal-feature runs on ont (P1.1 cache: GENERanno L17+L23, 4,096-d): G1B training + chr18 evaluation.

  train_formal_ont.py --name NAME --seed S [--arm G1B|R|G2PH] [--labels v50uev] [--lnc-only] [--epochs 12] [--smoke]
  train_formal_ont.py --name ab1-s1-42 --eval-only --from g1b-42 --with-s1      (AB1: S1 decoding of a main model)
  train_formal_ont.py --name dev-x --eval-only --from ab2-f-123 --ckpt avg10-12 --with-s1 --mea 1,x1 --dump-post
      (v1.9 development on chr18: weight average of epochs 10-12, MEA decoding, posteriors saved for CPU decoders)
Arms (protocol v1.8 section 15): G1B main / AB2 (--labels v50f) / AB4 (--labels v41uev) / AB3 (--labels v50uevlnc
--lnc-only); AB0 = --arm R (G1 without boundary head, screen2); AB6 = --arm G2PH (screen_g2ph, numba CPU FB).
Data (ont): data/p1_feat/{train,val,chr18}.npy (memmaps), data/p1v12/{train,val}_{codes,bits_<labels>,events}.npy,
chr18 truth outputs/PHASE1/DIAG_L2/chr18_v50.gtf.  Output: outputs/PHASE1/FORMAL/<NAME>/ (metrics.jsonl, last.pt,
best.pt, ep<NN>.pt (every epoch, 5.7 MB each), result.json, chr18.gff3).  Settings other than features / window set are those of the screening arm
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
import re
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
    ap.add_argument("--ckpt", default="best",
                    help="eval-only: checkpoint of the source run: best, last (final epoch), epNN (diagnostics only) or "
                         "avgNN-MM (v1.9 development: uniform weight average of epochs NN..MM)")
    ap.add_argument("--mea", default="",
                    help="eval-only (v1.9 development): comma list of [pv][x]gamma (screen2.mea_specs), each decoded with g1.mea_duration; "
                         "also adds the thresholded lnc L1 grid (conf_grid) to every decoder")
    ap.add_argument("--exact-post", action="store_true",
                    help="eval-only (v1.9 development): posteriors by the float64 CPU forward-backward (no float32 drift)")
    ap.add_argument("--dump-post", action="store_true",
                    help="eval-only (v1.9 development): save the chr18 posteriors as <out>/chr18_post.npy (float16, 2.7 GB)")
    ap.add_argument("--lnc-entry-bias", default="",
                    help="eval-only, G1B (diagnostic): comma list of per-gene log priors b on IG -> lncRNA entry, "
                         "each decoded as an extra S3_lncbias<b> block in result.json")
    ap.add_argument("--accum", type=int, default=int(os.environ.get("NCRNA_ACCUM", "2")),
                    help="batch 8 as accum x (8/accum) micro-batches (shared GPU)")
    ap.add_argument("--data", default="v12", choices=["v12", "v20", "v30h6", "v30hm"],
                    help="v20 (protocol section 17): training windows data/p1v20 (chr9 windows replaced), features = "
                         "data/p1_feat/train.npy with the replaced rows from data/p1_feat_v20/train_new.npy; val / chr18 unchanged. "
                         "v30h6 / v30hm (section 18): the v20 windows followed by the 3,000 new human (data/p1v30/human_new) "
                         "or mouse (data/p1v30/mouse) windows, features from data/p1_feat_v30/<set>.npy")
    ap.add_argument("--online-features", action="store_true", default=os.environ.get("NCRNA_ONLINE") == "1",
                    help="compute backbone features on the fly (win3090: no feature cache on disk)")
    a = ap.parse_args()
    if not re.fullmatch(r"best|last|ep\d\d|avg\d\d-\d\d", a.ckpt):
        ap.error(f"--ckpt {a.ckpt}: expected best, last, epNN or avgNN-MM")
    import torch
    if os.environ.get("NCRNA_THREADS"):
        torch.set_num_threads(int(os.environ["NCRNA_THREADS"]))
    frac = float(os.environ.get("NCRNA_GPU_FRAC", "0.30"))
    torch.cuda.set_per_process_memory_fraction(frac, 0)
    from ncrna_annot.train import screen_g1b as sg
    from ncrna_annot.train import screen2 as s2
    from ncrna_annot.train import screen_g2ph as sp
    out = f"{ROOT}/outputs/PHASE1/FORMAL/{a.name}"
    os.makedirs(out, exist_ok=True)
    F, D = f"{ROOT}/data/p1_feat", f"{ROOT}/data/p1v12"
    mm = lambda p: np.load(p, mmap_mode="r")
    DT = f"{ROOT}/data/p1v20" if a.data in ("v20", "v30h6", "v30hm") else D
    tr = dict(feats=mm(f"{F}/train.npy") if os.path.exists(f"{F}/train.npy") else None, codes=mm(f"{DT}/train_codes.npy"), bits=mm(f"{DT}/train_bits_{a.labels}.npy"),
              events=mm(f"{DT}/train_events.npy"))
    if a.data in ("v20", "v30h6", "v30hm"):
        sys.path.insert(0, f"{ROOT}/scripts/phase1/v20")
        from feat_overlay import Overlay
        rows = [r["row"] for r in json.load(open(f"{DT}/replaced.json"))["replaced"]]
        new = mm(f"{ROOT}/data/p1_feat_v20/train_new.npy")
        cov = json.load(open(f"{ROOT}/data/p1_feat_v20/train_new.cover.json"))
        assert sum(b - a_ for a_, b in cov) == len(rows) == len(new), "v20 new-window feature cache incomplete"
        tr["feats"] = Overlay(tr["feats"], new, rows)
    n_base = len(tr["codes"])   # windows covered by data/p1_feat/train.npy (checked below)
    if a.data in ("v30h6", "v30hm"):
        sys.path.insert(0, f"{ROOT}/scripts/phase1/v30")
        from concat_view import Concat
        part = "human_new" if a.data == "v30h6" else "mouse"
        P3 = f"{ROOT}/data/p1v30/{part}"
        extra = dict(feats=mm(f"{ROOT}/data/p1_feat_v30/{part}.npy"), codes=mm(f"{P3}/train_codes.npy"),
                     bits=mm(f"{P3}/train_bits_{a.labels}.npy"), events=mm(f"{P3}/train_events.npy"))
        cov = json.load(open(f"{ROOT}/data/p1_feat_v30/{part}.cover.json"))
        assert sum(b - a_ for a_, b in cov) == len(extra["codes"]) == len(extra["feats"]) == len(extra["bits"]) == len(extra["events"]), \
            f"v30 {part} cache or labels incomplete"
        tr = {k: Concat([tr[k], extra[k]]) for k in tr}
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
        va["feats"] = online.OnlineWindows(va["codes"], fn, cache_all=os.environ.get("NCRNA_VAL_CACHE", "1") == "1")
        ch["feats"] = online.OnlineChromosome(ch["codes"], fn)
    for k in ([] if a.online_features else ("train", "val", "chr18")):  # refuse a partially filled cache
        cov = json.load(open(f"{F}/{k}.cover.json"))
        need = n_base if k == "train" else len(va["codes"]) if k == "val" else None
        if need is not None:
            got = sum(b - a_ for a_, b in cov)
            assert got == need, f"{k} cache covers {got} of {need} windows"
    if a.smoke:
        tr = {k: np.asarray(v[:64]) for k, v in tr.items()}
        va = {k: np.asarray(v[:8]) for k, v in va.items()}
        ch.update(codes=ch["codes"][:3000000])
        ch["feats"] = (online.OnlineChromosome(ch["codes"], fn) if a.online_features else ch["feats"][:500000])
        a.epochs = 1
    meta = dict(name=a.name, data=a.data, host=os.uname().nodename, online_features=a.online_features, accum=a.accum, arm=a.arm, seed=a.seed, labels=a.labels, lnc_only=a.lnc_only, eval_only=a.eval_only,
                src=a.src, ckpt=a.ckpt, lnc_entry_bias=a.lnc_entry_bias, with_s1=a.with_s1, epochs=a.epochs, smoke=a.smoke, gpu_frac=frac,
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
            src = f"{ROOT}/outputs/PHASE1/FORMAL/{a.src or a.name}"
            ck = f"{src}/{a.ckpt}.pt"
            m = re.fullmatch(r"avg(\d\d)-(\d\d)", a.ckpt)
            if m:
                ck = [f"{src}/ep{e:02d}.pt" for e in range(int(m.group(1)), int(m.group(2)) + 1)]
            dev = dict(mea_gammas=a.mea, conf_grid=s2.CONF_GRID if a.mea else (),
                       post_out=f"{out}/chr18_post.npy" if a.dump_post else None, exact_post=a.exact_post)
            if a.arm == "G2PH":
                sp.evaluate_checkpoint(ck, "G2PH", 4096, ch, out, device="cuda")
            elif a.arm == "R":
                s2.evaluate_checkpoint(ck, 4096, ch, out, device="cuda", lnc_only=a.lnc_only, with_s1=a.with_s1, **dev)
            else:
                sg.evaluate_checkpoint(ck, 4096, ch, out, device="cuda", lnc_only=a.lnc_only, with_s1=a.with_s1,
                                       lnc_entry_bias=[float(x) for x in a.lnc_entry_bias.split(",") if x], **dev)
        elif a.arm == "G1B":
            sg.run(a.seed, tr, va, ch, out, epochs=a.epochs, device="cuda", pause=pause, resume=True, lnc_only=a.lnc_only,
                   accum=a.accum, keep_epochs=True)
        elif a.arm == "R":
            s2.run("R", a.seed, tr, va, ch, out, lnc_only=a.lnc_only, epochs=a.epochs, device="cuda", pause=pause, resume=True,
                   accum=a.accum, keep_epochs=True)
        else:  # G2PH: training (GPU head + numba CPU forward-backward), then evaluation in the same run
            if not os.path.exists(f"{out}/train_summary.json"):
                sp.run("G2PH", a.seed, tr, va, out, epochs=a.epochs, device="cuda", pause=pause, resume=True, accum=a.accum,
                       keep_epochs=True)
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
