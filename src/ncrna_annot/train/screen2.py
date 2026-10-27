"""P1-S2 (protocol v1.4 section 11): 12-epoch cosine training, checkpoint by validation NLL on the 300 chr18
windows, then ONE whole-chr18 evaluation: S3 decoding (+ S1 for arm R), gene confidence, gene-level PR/AUPRC.
"""
import json
import math
import os
import time

import numpy as np
import torch

from ..data import labels
from ..decode import chrom, scratch
from ..eval import prcurve, ruler
from ..features import generanno as G
from ..grammar import g1
from ..heads.c0 import C0Head
from .screen import M_STATES

TOK = 6


class Paused(Exception):
    """Raised by the training loops when the pause hook asks to yield the GPU (ont: a TE or other user's job)."""


def _emit(head, f, c, lnc_only):
    e = head(f, c)
    if lnc_only:
        e = e.clone()
        e[..., M_STATES] = g1.NEG
    return e


def val_nll(head, trans, va, lnc_only, device):
    head.eval()
    tot, n = 0.0, 0
    with torch.no_grad():
        T = trans()
        for i in range(0, len(va["codes"]), 4):
            f = torch.as_tensor(np.asarray(va["feats"][i:i + 4], np.float32), device=device)
            c_np = np.asarray(va["codes"][i:i + 4])
            c = torch.as_tensor(c_np.astype(np.int64), device=device)
            gates = torch.as_tensor(np.stack([g1.gate_mask(x) for x in c_np]), device=device)
            al = torch.as_tensor(labels.bits_to_bool(np.asarray(va["bits"][i:i + 4])), device=device)
            e = _emit(head, f, c, lnc_only)
            tot += float(((g1.log_partition(e, T, gates) - g1.log_partition(e, T, gates, al)) / c.shape[1]).sum())
            n += len(c_np)
    return tot / n


def evaluate_chromosome(head, trans, ch, lnc_only, device, with_s1=False):
    """ch = dict(feats [ntok, K] memmap, codes int8 [L], truth_gtf, chrom)."""
    head.eval()
    feats, codes = ch["feats"], ch["codes"]
    L = len(feats) * TOK
    T = trans().detach()
    Tn = T.double().cpu().numpy()

    def emit_fn(s, e):
        f = torch.as_tensor(np.asarray(feats[s // TOK:e // TOK], np.float32), device=device)[None]
        c = torch.as_tensor(np.asarray(codes[s:e], np.int64), device=device)[None]
        with torch.no_grad():
            return _emit(head, f, c, lnc_only)[0]

    t0 = time.time()
    pot = scratch.zeros((L, g1.S), np.float32)
    for s in range(0, L, G.WIN * 8):
        e = min(s + G.WIN * 8, L)
        pot[s:e] = emit_fn(s, e).float().cpu().numpy()
    post = chrom.chromosome_marginals(emit_fn, L, T, codes[:L], device)
    truth = ruler.read_gtf(ch["truth_gtf"], ch["chrom"])
    out = {"eval_s": None}
    modes = [("S3", "s3", pot)] + ([("S1", "argmax", scratch.astype(post, np.float32))] if with_s1 else [])
    for name, mode, arr in modes:
        genes = chrom.decode_blocks(arr, codes[:L], Tn, 20, 50, mode=mode)
        out[name] = summarize(genes, post, truth)
    out["eval_s"] = round(time.time() - t0)
    return out


def summarize(genes, post, truth):
    """lnc >= 200 nt filter, gene confidence, gene-level PR / AUPRC and ruler L1 / L2 (protocol v1.4)."""
    genes = [g for g in genes if g["type"] != "l" or sum(b - a + 1 for a, b in g["exons"]) >= 200]
    for g in genes:
        g["conf"] = chrom.gene_confidence(g, post)
    lnc = [{"strand": g["strand"], "exons": [tuple(x) for x in g["exons"]], "conf": g["conf"]}
           for g in genes if g["type"] == "l"]
    pr = prcurve.pr_curve(lnc, truth)
    preds = [{"id": str(i), "strand": g["strand"], "cls": "lnc" if g["type"] == "l" else "mrna",
              "exons": [tuple(x) for x in g["exons"]]} for i, g in enumerate(genes)]
    idx = ruler.ExonIndex(truth)
    sl = ruler.score_class("lnc", truth, preds, idx, 0)
    sm = ruler.score_class("mrna", truth, preds, idx, 0)
    return {"auprc_lnc": pr["auprc"], "best_locus_f1": pr["best_f1"], "best_point": pr["best_point"],
            "n_truth_lnc": pr["n_truth"], "n_pred_lnc": pr["n_pred"],
            "lnc_L1": {"sn": sl["L1_locus"]["sn"], "pr": sl["L1_locus"]["pr_strict"]},
            "lnc_L2": {"sn": sl["L2_intron_chain"]["sn"], "pr": sl["L2_intron_chain"]["pr_strict"]},
            "mrna_L1": {"sn": sm["L1_locus"]["sn"], "pr": sm["L1_locus"]["pr_strict"]},
            "mrna_L2": {"sn": sm["L2_intron_chain"]["sn"], "pr": sm["L2_intron_chain"]["pr_strict"]},
            "n_pred_mrna": sum(1 for g in genes if g["type"] == "m"), "curve": pr["curve"]}


def run(arm, seed, tr, va, ch, out, lnc_only=False, epochs=12, batch=8, lr=1e-3, device="cuda", log=print,
        pause=None, resume=False, accum=1):
    """Arms R / S2 / S4 (G1 without boundary head).
    pause: optional callable checked before every batch; True aborts the epoch with Paused.  resume: keep last.pt
    after every epoch and continue from it (see screen_g1b.run)."""
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    os.makedirs(out, exist_ok=True)
    N = len(tr["codes"])
    head = C0Head(tr["feats"].shape[-1]).to(device)
    trans = g1.G1Transitions().to(device)
    params = list(head.parameters()) + list(trans.parameters())
    opt = torch.optim.AdamW(params, lr=lr, weight_decay=1e-2)
    steps = epochs * math.ceil(N / batch)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=steps, eta_min=1e-5)
    best, t0, start = (math.inf, 0), time.time(), 0
    last = f"{out}/last.pt"
    if resume and os.path.exists(last):
        st = torch.load(last, map_location=device, weights_only=False)
        head.load_state_dict(st["head"]); trans.load_state_dict(st["trans"])
        opt.load_state_dict(st["opt"]); sched.load_state_dict(st["sched"])
        rng.bit_generator.state = st["rng"]; torch.set_rng_state(st["torch_rng"].cpu())  # map_location moved it to the GPU
        best, start, t0 = tuple(st["best"]), st["epoch"], time.time() - st["t"]
        log(json.dumps({"resume_from_epoch": start, "best": best}))
    for ep in range(start, epochs):
        head.train()
        order, tot, nb = rng.permutation(N), 0.0, 0
        for i in range(0, N, batch):
            if pause is not None and pause():
                raise Paused(f"epoch {ep + 1} batch {i // batch}")
            idx = np.sort(order[i:i + batch])
            opt.zero_grad()
            b_loss = 0.0
            for part in (np.array_split(idx, accum) if accum > 1 else [idx]):  # gradient accumulation (exact mean)
                if not len(part):
                    continue
                w = len(part) / len(idx)
                f = torch.as_tensor(np.asarray(tr["feats"][part], np.float32), device=device)
                c_np = np.asarray(tr["codes"][part])
                c = torch.as_tensor(c_np.astype(np.int64), device=device)
                gates = torch.as_tensor(np.stack([g1.gate_mask(x) for x in c_np]), device=device)
                al = torch.as_tensor(labels.bits_to_bool(np.asarray(tr["bits"][part])), device=device)
                e, T = _emit(head, f, c, lnc_only), trans()
                loss = ((g1.log_partition(e, T, gates) - g1.log_partition(e, T, gates, al)) / c.shape[1]).mean()
                (loss * w).backward()
                b_loss += float(loss.detach()) * w
            torch.nn.utils.clip_grad_norm_(params, 1.0)
            opt.step()
            sched.step()
            tot, nb = tot + b_loss, nb + 1
        v = val_nll(head, trans, va, lnc_only, device)
        rec = {"arm": arm, "seed": seed, "epoch": ep + 1, "train_nll": round(tot / nb, 5), "val_nll": round(v, 5),
               "lr": opt.param_groups[0]["lr"], "t": round(time.time() - t0)}
        with open(f"{out}/metrics.jsonl", "a") as fh:
            fh.write(json.dumps(rec) + "\n")
        log(json.dumps(rec))
        if v < best[0]:
            best = (v, ep + 1)
            torch.save({"head": head.state_dict(), "trans": trans.state_dict(), "epoch": ep + 1}, f"{out}/best.pt")
        if resume:
            torch.save({"head": head.state_dict(), "trans": trans.state_dict(), "opt": opt.state_dict(),
                        "sched": sched.state_dict(), "rng": rng.bit_generator.state, "torch_rng": torch.get_rng_state(),
                        "best": list(best), "epoch": ep + 1, "t": round(time.time() - t0)}, last + ".tmp")
            os.replace(last + ".tmp", last)
    ck = torch.load(f"{out}/best.pt", map_location=device)
    head.load_state_dict(ck["head"])
    trans.load_state_dict(ck["trans"])
    res = {"arm": arm, "seed": seed, "best_epoch": best[1], "best_val_nll": best[0], "train_s": round(time.time() - t0),
           **evaluate_chromosome(head, trans, ch, lnc_only, device, with_s1=(arm == "R"))}
    json.dump(res, open(f"{out}/result.json", "w"), indent=1)
    log(json.dumps({k: (v if not isinstance(v, dict) else {kk: vv for kk, vv in v.items() if kk != "curve"})
                    for k, v in res.items()}))
    return res
