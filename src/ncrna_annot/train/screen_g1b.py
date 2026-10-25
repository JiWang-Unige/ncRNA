"""G1B arm (protocol v1.6 section 13): G1 + boundary head, otherwise identical to P1-S2 arm R (v1.4 section 11).

Boundary head: TSS / TES logits per strand from the C0 hidden state, trained with an auxiliary BCE (weight 1,
ignored bases excluded) on the P1-G0 event targets, and used as position-dependent bonuses on gene entry / exit
edges (g1.boundary_edges).  Checkpoint by validation CRF NLL per base (the BCE term is not part of the criterion).
Evaluation: S3 decoding with the same bonuses, gene confidence from G1B posteriors, metrics via screen2.summarize.
"""
import json
import math
import os
import time

import numpy as np
import torch

from ..data import labels
from ..decode import chrom
from ..eval import ruler
from ..features import generanno as G
from ..grammar import g1
from ..heads.c0 import C0BHead
from .screen2 import summarize

TOK = 6
EV_BOUNDARY = [2, 3, 6, 7]  # P1-G0 event channels (TSS+, TES+, TSS-, TES-) in g1.BOUNDARY order


def boundary_bce(logits, ev):
    """logits [B, L, 4]; ev int8 [B, L, 4] (1 pos, 0 neg, -1 ignore) -> mean BCE over kept entries."""
    keep = ev >= 0
    if not keep.any():
        return logits.sum() * 0.0
    return torch.nn.functional.binary_cross_entropy_with_logits(logits[keep], ev[keep].float())


def crf_nll(head, trans, coef, f, c, c_np, bits, device):
    e, b = head(f, c)
    T = trans()
    gates = torch.as_tensor(np.stack([g1.gate_mask(x) for x in c_np]), device=device)
    al = torch.as_tensor(labels.bits_to_bool(bits), device=device)
    edges = g1.boundary_edges(b, coef)
    nll = (g1.log_partition(e, T, gates, edges=edges) - g1.log_partition(e, T, gates, al, edges=edges)) / c.shape[1]
    return nll, b


def val_nll(head, trans, coef, va, device):
    head.eval()
    tot, n = 0.0, 0
    with torch.no_grad():
        for i in range(0, len(va["codes"]), 4):
            f = torch.as_tensor(np.asarray(va["feats"][i:i + 4], np.float32), device=device)
            c_np = np.asarray(va["codes"][i:i + 4])
            c = torch.as_tensor(c_np.astype(np.int64), device=device)
            nll, _ = crf_nll(head, trans, coef, f, c, c_np, np.asarray(va["bits"][i:i + 4]), device)
            tot += float(nll.sum())
            n += len(c_np)
    return tot / n


def evaluate_chromosome(head, trans, coef, ch, device):
    head.eval()
    feats, codes = ch["feats"], ch["codes"]
    L = len(feats) * TOK
    T = trans().detach()

    def emit_fn(s, e):
        f = torch.as_tensor(np.asarray(feats[s // TOK:e // TOK], np.float32), device=device)[None]
        c = torch.as_tensor(np.asarray(codes[s:e], np.int64), device=device)[None]
        with torch.no_grad():
            p, b = head(f, c)
        return p[0], b[0]

    t0 = time.time()
    pot = np.zeros((L, g1.S), np.float32)
    bnd = np.zeros((L, len(g1.BOUNDARY)), np.float32)
    for s in range(0, L, G.WIN * 8):
        e = min(s + G.WIN * 8, L)
        p, b = emit_fn(s, e)
        pot[s:e], bnd[s:e] = p.float().cpu().numpy(), b.float().cpu().numpy()
    post = chrom.chromosome_marginals(emit_fn, L, T, codes[:L], device, coef=coef)
    truth = ruler.read_gtf(ch["truth_gtf"], ch["chrom"])
    genes = chrom.decode_blocks(pot, codes[:L], T.double().cpu().numpy(), 20, 50, mode="s3", bnd=bnd, coef=coef)
    out = {"S3": summarize(genes, post, truth), "eval_s": round(time.time() - t0),
           "coef": {k: [float(x) for x in coef.c[k].detach().cpu()] for k in g1.TYPES}}
    if ch.get("gff_out"):
        from ..decode import genes as dec
        dec.write_gff3([g for g in genes if g["type"] != "l" or sum(b - a + 1 for a, b in g["exons"]) >= 200],
                       ch["chrom"], ch["gff_out"])
    return out


def evaluate_checkpoint(ck_path, feature_dim, ch, out, device="cuda", log=print):
    """Separate evaluation job (Kaggle input-size limits): load best.pt, evaluate on ch, write result_eval.json."""
    head, trans, coef = C0BHead(feature_dim).to(device), g1.G1Transitions().to(device), g1.BoundaryCoef().to(device)
    ck = torch.load(ck_path, map_location=device)
    head.load_state_dict(ck["head"])
    trans.load_state_dict(ck["trans"])
    coef.load_state_dict(ck["coef"])
    res = {"arm": "G1B", "best_epoch": ck["epoch"], **evaluate_chromosome(head, trans, coef, ch, device)}
    os.makedirs(out, exist_ok=True)
    json.dump(res, open(f"{out}/result.json", "w"), indent=1)
    log(json.dumps({k: (v if not isinstance(v, dict) else {kk: vv for kk, vv in v.items() if kk != "curve"})
                    for k, v in res.items()}))
    return res


class Paused(Exception):
    """Raised by run() when the pause hook asks to yield the GPU (ont: a TE or other user's job appeared)."""


def run(seed, tr, va, ch, out, epochs=12, batch=8, lr=1e-3, bce_weight=1.0, device="cuda", log=print,
        pause=None, resume=False):
    """tr / va = dict(feats, codes, bits, events [n, L, 8] int8); ch as in screen2, or None (training only;
    evaluate later with evaluate_checkpoint).
    pause: optional callable checked before every batch; returning True aborts the current epoch with Paused
    (the caller frees the GPU and restarts later).  resume: keep last.pt (model, optimiser, scheduler, RNG, best)
    after every epoch and continue from it, so an aborted run repeats at most one epoch and reproduces the
    uninterrupted run's epoch-level trajectory (same permutations, same updates)."""
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    os.makedirs(out, exist_ok=True)
    N = len(tr["codes"])
    head = C0BHead(tr["feats"].shape[-1]).to(device)
    trans = g1.G1Transitions().to(device)
    coef = g1.BoundaryCoef().to(device)
    params = list(head.parameters()) + list(trans.parameters()) + list(coef.parameters())
    opt = torch.optim.AdamW(params, lr=lr, weight_decay=1e-2)
    steps = epochs * math.ceil(N / batch)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=steps, eta_min=1e-5)
    best, t0, start = (math.inf, 0), time.time(), 0
    last = f"{out}/last.pt"
    if resume and os.path.exists(last):
        st = torch.load(last, map_location=device, weights_only=False)
        head.load_state_dict(st["head"]); trans.load_state_dict(st["trans"]); coef.load_state_dict(st["coef"])
        opt.load_state_dict(st["opt"]); sched.load_state_dict(st["sched"])
        rng.bit_generator.state = st["rng"]; torch.set_rng_state(st["torch_rng"])
        best, start, t0 = tuple(st["best"]), st["epoch"], time.time() - st["t"]
        log(json.dumps({"resume_from_epoch": start, "best": best}))
    for ep in range(start, epochs):
        head.train()
        order, tot, tot_b, nb = rng.permutation(N), 0.0, 0.0, 0
        for i in range(0, N, batch):
            if pause is not None and pause():
                raise Paused(f"epoch {ep + 1} batch {i // batch}")
            idx = np.sort(order[i:i + batch])
            f = torch.as_tensor(np.asarray(tr["feats"][idx], np.float32), device=device)
            c_np = np.asarray(tr["codes"][idx])
            c = torch.as_tensor(c_np.astype(np.int64), device=device)
            nll, b = crf_nll(head, trans, coef, f, c, c_np, np.asarray(tr["bits"][idx]), device)
            ev = torch.as_tensor(np.asarray(tr["events"][idx])[..., EV_BOUNDARY], device=device)
            bce = boundary_bce(b, ev)
            loss = nll.mean() + bce_weight * bce
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(params, 1.0)
            opt.step()
            sched.step()
            tot, tot_b, nb = tot + float(nll.mean()), tot_b + float(bce), nb + 1
        v = val_nll(head, trans, coef, va, device)
        rec = {"arm": "G1B", "seed": seed, "epoch": ep + 1, "train_nll": round(tot / nb, 5), "train_bce": round(tot_b / nb, 5),
               "val_nll": round(v, 5), "lr": opt.param_groups[0]["lr"], "t": round(time.time() - t0),
               "coef": {k: [round(float(x), 4) for x in coef.c[k].detach().cpu()] for k in g1.TYPES}}
        with open(f"{out}/metrics.jsonl", "a") as fh:
            fh.write(json.dumps(rec) + "\n")
        log(json.dumps(rec))
        if v < best[0]:
            best = (v, ep + 1)
            torch.save({"head": head.state_dict(), "trans": trans.state_dict(), "coef": coef.state_dict(),
                        "epoch": ep + 1}, f"{out}/best.pt")
        if resume:
            torch.save({"head": head.state_dict(), "trans": trans.state_dict(), "coef": coef.state_dict(),
                        "opt": opt.state_dict(), "sched": sched.state_dict(), "rng": rng.bit_generator.state,
                        "torch_rng": torch.get_rng_state(), "best": list(best), "epoch": ep + 1,
                        "t": round(time.time() - t0)}, last + ".tmp")
            os.replace(last + ".tmp", last)
    res = {"arm": "G1B", "seed": seed, "best_epoch": best[1], "best_val_nll": best[0], "train_s": round(time.time() - t0)}
    if ch is None:
        json.dump(res, open(f"{out}/train_summary.json", "w"), indent=1)
        log(json.dumps(res))
        return res
    ck = torch.load(f"{out}/best.pt", map_location=device)
    head.load_state_dict(ck["head"])
    trans.load_state_dict(ck["trans"])
    coef.load_state_dict(ck["coef"])
    res.update(evaluate_chromosome(head, trans, coef, ch, device))
    json.dump(res, open(f"{out}/result.json", "w"), indent=1)
    log(json.dumps({k: (v if not isinstance(v, dict) else {kk: vv for kk, vv in v.items() if kk != "curve"})
                    for k, v in res.items()}))
    return res
