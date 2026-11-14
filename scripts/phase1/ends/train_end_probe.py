#!/usr/bin/env python3
"""PLAN_END.md section 2 (stage 1): end probe = C0 hidden layers + 8 class-specific end logits, BCE only (no CRF).

v2.0 training windows (data/p1v20 codes, data/p1_feat/train.npy with the v2.0 replaced rows), labels
data/p1ends/train_<arm>.npy, fixed epochs, the last epoch is the model (no validation selection: the validation
windows are on chr18).  Afterwards the model is run over chr18 in 8 x 24,576 bp pieces (as screen_g1b.evaluate_chromosome)
and the logits are saved as <out>/chr18_end.npy (float16 [L, 8], channels prep_end_labels.CH8).
Run through scripts/phase1/run_formal_ont.sh --script ends/train_end_probe.py (GPU gate, lock; exit 75 = yielded,
resumed from the last finished epoch).

  train_end_probe.py --name end-ann-42 --arm ann [--seed 42 --epochs 4]
"""
import argparse
import json
import math
import os
import subprocess
import sys
import time

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1/v20"))
TOK, WIN = 6, 24576


def gate_busy():
    r = subprocess.run(["bash", "-c", f". {ROOT}/scripts/phase1/gpu_gate.sh; gpu_why_busy"], capture_output=True, text=True)
    return r.stdout.strip() if r.returncode == 0 else ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--arm", required=True, choices=["ann", "ev"])
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--epochs", type=int, default=4)
    ap.add_argument("--batch", type=int, default=8)
    ap.add_argument("--accum", type=int, default=2)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--chrom", default="chr18", help="chromosome for the forward pass")
    ap.add_argument("--forward-only", action="store_true",
                    help="no training: load ep<epochs>.pt of an existing run and run the forward on --chrom (PLAN_REP20)")
    ap.add_argument("--smoke", action="store_true", help="16 windows, 1 epoch, first 2 Mb of chr18")
    a = ap.parse_args()
    import torch
    from torch import nn
    from feat_overlay import Overlay
    from ncrna_annot.heads.c0 import C0Head
    torch.cuda.set_per_process_memory_fraction(float(os.environ.get("NCRNA_GPU_FRAC", "0.30")), 0)
    dev = torch.device("cuda")
    out = f"{ROOT}/outputs/PHASE1/ENDS/{a.name}"
    os.makedirs(out, exist_ok=True)
    log = open(f"{out}/train.log" if not a.forward_only else f"{out}/forward_{a.chrom}.log", "a")

    def say(**kw):
        kw["t"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        print(json.dumps(kw), flush=True)
        log.write(json.dumps(kw) + "\n")
        log.flush()

    mm = lambda p: np.load(p, mmap_mode="r")  # noqa: E731
    rows = [r["row"] for r in json.load(open(f"{ROOT}/data/p1v20/replaced.json"))["replaced"]]
    feats = Overlay(mm(f"{ROOT}/data/p1_feat/train.npy"), mm(f"{ROOT}/data/p1_feat_v20/train_new.npy"), rows)
    codes = mm(f"{ROOT}/data/p1v20/train_codes.npy")
    lab = mm(f"{ROOT}/data/p1ends/train_{a.arm}.npy")
    assert len(feats) == len(codes) == len(lab), (len(feats), len(codes), len(lab))
    N = 16 if a.smoke else len(codes)
    epochs = 1 if a.smoke else a.epochs

    class EndProbe(C0Head):
        def __init__(self, d):
            super().__init__(d)
            self.out = nn.Identity()          # the 17-state output layer of C0 is not used
            self.end = nn.Linear(256, 8)

        def forward(self, f, c):
            return self.end(self.hidden(f, c))

    torch.manual_seed(a.seed)
    rng = np.random.default_rng(a.seed)
    model = EndProbe(feats.shape[-1]).to(dev)
    opt = torch.optim.AdamW(model.parameters(), lr=a.lr, weight_decay=1e-2)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=epochs * math.ceil(N / a.batch), eta_min=1e-5)
    start = 0
    last = f"{out}/last.pt"
    if a.forward_only:
        st = torch.load(f"{out}/ep{epochs:02d}.pt", map_location=dev, weights_only=False)
        model.load_state_dict(st["model"])
        start = epochs
        say(forward_only=a.chrom, weights=f"ep{epochs:02d}.pt", epoch=st["epoch"])
    elif os.path.exists(last):
        st = torch.load(last, map_location=dev, weights_only=False)
        model.load_state_dict(st["model"]); opt.load_state_dict(st["opt"]); sched.load_state_dict(st["sched"])
        rng.bit_generator.state = st["rng"]; torch.set_rng_state(st["torch_rng"].cpu())
        start = st["epoch"]
        say(resume_from_epoch=start)
    meta = dict(vars(a), n_windows=N, epochs_run=epochs, host=os.uname().nodename)
    if not a.forward_only:
        json.dump(meta, open(f"{out}/run_meta.json", "w"), indent=1)
    t0 = time.time()
    for ep in range(start, epochs):
        model.train()
        order, tot, nb = rng.permutation(N), 0.0, 0
        for i in range(0, N, a.batch):
            if nb % 25 == 0:
                why = gate_busy()
                if why and not why.startswith("only "):   # our own allocation counts as used memory: never yield for it
                    say(yield_gpu=why, epoch=ep + 1, batch=nb)
                    sys.exit(75)
            idx = np.sort(order[i:i + a.batch])
            opt.zero_grad()
            bl = 0.0
            for part in np.array_split(idx, a.accum):
                if not len(part):
                    continue
                f = torch.as_tensor(np.asarray(feats[part], np.float32), device=dev)
                c = torch.as_tensor(np.asarray(codes[part]).astype(np.int64), device=dev)
                y = torch.as_tensor(np.asarray(lab[part]), device=dev)
                z = model(f, c)
                keep = y >= 0
                loss = nn.functional.binary_cross_entropy_with_logits(z[keep], y[keep].float()) * (len(part) / len(idx))
                loss.backward()
                bl += float(loss)
            opt.step()
            sched.step()
            tot += bl
            nb += 1
        say(epoch=ep + 1, train_bce=tot / nb, s=round(time.time() - t0))
        torch.save({"model": model.state_dict(), "epoch": ep + 1}, f"{out}/ep{ep + 1:02d}.pt")
        torch.save({"model": model.state_dict(), "opt": opt.state_dict(), "sched": sched.state_dict(),
                    "rng": rng.bit_generator.state, "torch_rng": torch.get_rng_state(), "epoch": ep + 1}, last)

    # forward on --chrom with the last epoch
    C, CU = a.chrom, a.chrom.upper()
    if os.path.exists(f"{out}/{C}_end.npy") and os.path.exists(f"{out}/{CU}_DONE"):
        say(forward=C, status="already done")
        return
    model.eval()
    cf = mm(f"{ROOT}/data/p1_feat/{C}.npy")
    cc = mm(f"{ROOT}/data/p1/{C}_codes.npy")
    L = len(cf) * TOK
    if a.smoke:
        L = min(WIN * 8 * 10, L)
    res = np.lib.format.open_memmap(f"{out}/{C}_end.npy", mode="w+", dtype=np.float16, shape=(L, 8))
    t1 = time.time()
    with torch.no_grad():
        for s in range(0, L, WIN * 8):
            e = min(s + WIN * 8, L)
            f = torch.as_tensor(np.asarray(cf[s // TOK:e // TOK], np.float32), device=dev)[None]
            c = torch.as_tensor(np.asarray(cc[s:e], np.int64), device=dev)[None]
            res[s:e] = model(f, c)[0].float().cpu().numpy().astype(np.float16)
    res.flush()
    open(f"{out}/{CU}_DONE", "w").write(f"{L}\n")
    say(forward_done=C, L=L, s=round(time.time() - t1))


if __name__ == "__main__":
    main()
