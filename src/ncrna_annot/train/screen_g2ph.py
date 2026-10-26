"""G2-PH arm (protocol v1.7 section 14): C0BHead on the GPU, phase-type grammar forward-backward on the CPU (numba).

Training as G1B (12 epochs, cosine LR, batch 8, val-NLL checkpointing, boundary BCE weight 1) with the G2-PH CRF loss.
Windows whose observed partition is infeasible under the grammar are excluded from the loss and reported.
Evaluation: whole chromosome in 1 Mb blocks (+-100 kb): Viterbi on the expanded space (no S3 minimum lengths),
lnc >= 200 nt, gene confidence from float64 block posteriors; metrics via screen2.summarize.
"""
import json
import math
import os
import time

import numpy as np
import torch

from ..data import labels
from ..eval import ruler
from ..features import generanno as G
from ..grammar import g1, g2ph
from ..heads.c0 import C0BHead
from .screen2 import summarize
from .screen2 import Paused  # noqa: E402
from .screen_g1b import EV_BOUNDARY, boundary_bce

TOK = 6
CONFIGS = {"G2PH": (3, 2, 3), "G2PH1": (1, 1, 1)}


def _gates(c_np):
    return np.stack([g1.gate_mask(x) for x in c_np])


def feasibility(space, codes, bits, chunk=64):
    out = []
    for i in range(0, len(codes), chunk):
        c_np = np.asarray(codes[i:i + chunk])
        out.append(g2ph.feasible(space, _gates(c_np), labels.bits_to_bool(np.asarray(bits[i:i + chunk]))))
    return np.concatenate(out)


def _batch_nll(head, prm, space, feats, codes, bits, idx, feas, device):
    f = torch.as_tensor(np.asarray(feats[idx], np.float32), device=device)
    c_np = np.asarray(codes[idx])
    c = torch.as_tensor(c_np.astype(np.int64), device=device)
    x17, b = head(f, c)
    loss, (z, zo) = g2ph.nll_surrogate(space, prm, x17, b, _gates(c_np), labels.bits_to_bool(np.asarray(bits[idx])),
                                       feas[idx])
    return loss / c.shape[1], b


def val_nll(head, prm, space, va, feas, device):
    head.eval()
    tot, n = 0.0, 0
    with torch.no_grad():
        for i in range(0, len(va["codes"]), 4):
            idx = np.arange(i, min(i + 4, len(va["codes"])))
            nll, _ = _batch_nll(head, prm, space, va["feats"], va["codes"], va["bits"], idx, feas, device)
            keep = feas[idx]
            tot += float(nll[torch.as_tensor(keep, device=nll.device)].sum())
            n += int(keep.sum())
    return tot / max(n, 1)


def evaluate_chromosome(head, prm, space, ch, device, block=1_000_000, pad=100_000):
    head.eval()
    feats, codes = ch["feats"], ch["codes"]
    L = len(feats) * TOK
    t0 = time.time()
    pot = np.zeros((L, g1.S), np.float32)
    bnd = np.zeros((L, len(g1.BOUNDARY)), np.float32)
    with torch.no_grad():
        for s in range(0, L, G.WIN * 8):
            e = min(s + G.WIN * 8, L)
            f = torch.as_tensor(np.asarray(feats[s // TOK:e // TOK], np.float32), device=device)[None]
            c = torch.as_tensor(np.asarray(codes[s:e], np.int64), device=device)[None]
            p, b = head(f, c)
            pot[s:e], bnd[s:e] = p[0].float().cpu().numpy(), b[0].float().cpu().numpy()
        T = prm().detach().cpu().numpy()
    post = np.zeros((L, g1.S), np.float32)
    genes = []
    for s in range(0, L, block):
        e = min(s + block, L)
        a, b = max(s - pad, 0), min(e + pad, L)
        X = g2ph.expand(space, pot[a:b][None], g1.gate_mask(np.asarray(codes[a:b]))[None])
        with torch.no_grad():
            EV = prm.bonus(torch.as_tensor(bnd[a:b], device=prm.coef.device).double()).cpu().numpy()
        path, _ = g2ph.viterbi(space, X[0], T, EV)
        for gene in g1.path_to_genes(path, offset=a):
            if s <= gene["exons"][0][0] - 1 < e:
                genes.append(gene)
        _, p, _, _ = g2ph.forward_backward(space, X, T, EV[None])
        post[s:e] = p[0, s - a:e - a]
        del X, p
    truth = ruler.read_gtf(ch["truth_gtf"], ch["chrom"])
    out = {"PH": summarize(genes, post, truth), "eval_s": round(time.time() - t0)}
    if ch.get("gff_out"):
        from ..decode import genes as dec
        dec.write_gff3([g for g in genes if g["type"] != "l" or sum(y - x + 1 for x, y in g["exons"]) >= 200],
                       ch["chrom"], ch["gff_out"])
    return out


def _load(ck_path, arm, feature_dim, device):
    space = g2ph.Space(*CONFIGS[arm])
    head, prm = C0BHead(feature_dim).to(device), g2ph.G2PHParams(space).to(device)
    ck = torch.load(ck_path, map_location=device)
    head.load_state_dict(ck["head"])
    prm.load_state_dict(ck["grammar"])
    return head, prm, space, ck


def evaluate_checkpoint(ck_path, arm, feature_dim, ch, out, device="cuda", log=print):
    head, prm, space, ck = _load(ck_path, arm, feature_dim, device)
    res = {"arm": arm, "best_epoch": ck["epoch"], **evaluate_chromosome(head, prm, space, ch, device)}
    os.makedirs(out, exist_ok=True)
    json.dump(res, open(f"{out}/result.json", "w"), indent=1)
    log(json.dumps({k: (v if not isinstance(v, dict) else {kk: vv for kk, vv in v.items() if kk != "curve"})
                    for k, v in res.items()}))
    return res


def run(arm, seed, tr, va, out, epochs=12, batch=8, lr=1e-3, bce_weight=1.0, device="cuda", log=print,
        pause=None, resume=False):
    """Training only (evaluation is a separate job, see evaluate_checkpoint).  tr / va as in screen_g1b.
    pause: optional callable checked before every batch; True aborts the epoch with Paused.  resume: keep last.pt
    after every epoch and continue from it (see screen_g1b.run)."""
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    os.makedirs(out, exist_ok=True)
    space = g2ph.Space(*CONFIGS[arm])
    t0 = time.time()
    feas_tr = feasibility(space, tr["codes"], tr["bits"])
    feas_va = feasibility(space, va["codes"], va["bits"])
    info = {"arm": arm, "seed": seed, "states": space.n, "infeasible_train": float(1 - feas_tr.mean()),
            "infeasible_val": float(1 - feas_va.mean()), "feasibility_s": round(time.time() - t0)}
    log(json.dumps(info))
    N = len(tr["codes"])
    head = C0BHead(tr["feats"].shape[-1]).to(device)
    prm = g2ph.G2PHParams(space).to(device)
    params = list(head.parameters()) + list(prm.parameters())
    opt = torch.optim.AdamW(params, lr=lr, weight_decay=1e-2)
    steps = epochs * math.ceil(N / batch)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=steps, eta_min=1e-5)
    best, t0, start = (math.inf, 0), time.time(), 0
    last = f"{out}/last.pt"
    if resume and os.path.exists(last):
        st = torch.load(last, map_location=device, weights_only=False)
        head.load_state_dict(st["head"]); prm.load_state_dict(st["grammar"])
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
            nll, b = _batch_nll(head, prm, space, tr["feats"], tr["codes"], tr["bits"], idx, feas_tr, device)
            keep = torch.as_tensor(feas_tr[idx], device=nll.device)
            ev = torch.as_tensor(np.asarray(tr["events"][idx])[..., EV_BOUNDARY], device=device)
            bce = boundary_bce(b, ev)
            nll_mean = nll[keep].mean() if keep.any() else nll.sum() * 0.0
            loss = nll_mean.float() + bce_weight * bce
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(params, 1.0)
            opt.step()
            sched.step()
            tot, tot_b, nb = tot + float(nll_mean.detach()), tot_b + float(bce.detach()), nb + 1
        v = val_nll(head, prm, space, va, feas_va, device)
        rec = {"arm": arm, "seed": seed, "epoch": ep + 1, "train_nll": round(tot / nb, 5), "train_bce": round(tot_b / nb, 5),
               "val_nll": round(v, 5), "lr": opt.param_groups[0]["lr"], "t": round(time.time() - t0)}
        with open(f"{out}/metrics.jsonl", "a") as fh:
            fh.write(json.dumps(rec) + "\n")
        log(json.dumps(rec))
        if v < best[0]:
            best = (v, ep + 1)
            torch.save({"head": head.state_dict(), "grammar": prm.state_dict(), "epoch": ep + 1, "arm": arm},
                       f"{out}/best.pt")
        if resume:
            torch.save({"head": head.state_dict(), "grammar": prm.state_dict(), "opt": opt.state_dict(),
                        "sched": sched.state_dict(), "rng": rng.bit_generator.state, "torch_rng": torch.get_rng_state(),
                        "best": list(best), "epoch": ep + 1, "t": round(time.time() - t0)}, last + ".tmp")
            os.replace(last + ".tmp", last)
    res = {**info, "best_epoch": best[1], "best_val_nll": best[0], "train_s": round(time.time() - t0)}
    json.dump(res, open(f"{out}/train_summary.json", "w"), indent=1)
    log(json.dumps(res))
    return res
