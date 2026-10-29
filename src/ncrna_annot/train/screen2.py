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


def evaluate_chromosome(head, trans, ch, lnc_only, device, with_s1=False, mea_gammas=(), conf_grid=(), post_out=None,
                        log=None, exact_post=False):
    """ch = dict(feats [ntok, K] memmap, codes int8 [L], truth_gtf, chrom).
    mea_gammas / conf_grid / post_out / exact_post: v1.9 development, as in screen_g1b.evaluate_chromosome."""
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
    post = chrom.chromosome_marginals(emit_fn, L, T, codes[:L], device, exact=exact_post)
    if post_out:
        np.save(post_out, post)
    truth = ruler.read_gtf(ch["truth_gtf"], ch["chrom"])
    out = {"eval_s": None}
    modes = [("S3", "s3", pot, 1.0, "state")] + ([("S1", "argmax", scratch.astype(post, np.float32), 1.0, "state")] if with_s1 else [])
    modes = [m + (False,) for m in modes] + [(name, "mea", post, gm, lv, lg) for name, gm, lv, lg in mea_specs(mea_gammas)]
    for name, mode, arr, gm, lv, lg in modes:
        genes = chrom.decode_blocks(arr, codes[:L], Tn, 20, 50, mode=mode, gamma=gm, level=lv, log_gain=lg)
        out[name] = summarize(genes, post, truth, conf_grid)
        if log and mode == "mea":
            log(json.dumps({name: {k: v for k, v in out[name].items() if k != "curve"}}))
    out["eval_s"] = round(time.time() - t0)
    return out


def load_ckpt(ck_path, device):
    """One checkpoint, or (list of paths) the uniform average of their weights ("epoch" becomes the list of epochs).
    v1.9 development: averaging the last epochs of one run (same cosine-annealed basin; the head has LayerNorm only,
    no batch statistics to recompute)."""
    if isinstance(ck_path, str):
        return torch.load(ck_path, map_location=device)
    cks = [torch.load(p, map_location=device) for p in ck_path]
    avg = {k: {n: sum(c[k][n].double() for c in cks).div(len(cks)).to(cks[0][k][n].dtype) for n in cks[0][k]}
           for k, v in cks[0].items() if isinstance(v, dict)}
    avg["epoch"] = [c["epoch"] for c in cks]
    return avg


def evaluate_checkpoint(ck_path, feature_dim, ch, out, arm="R", device="cuda", log=print, lnc_only=False, with_s1=False,
                        mea_gammas=(), conf_grid=(), post_out=None, exact_post=False):
    """Eval-only job for arm R (AB0): load a checkpoint (or average a list), evaluate on ch, write result.json."""
    head, trans = C0Head(feature_dim).to(device), g1.G1Transitions().to(device)
    ck = load_ckpt(ck_path, device)
    head.load_state_dict(ck["head"])
    trans.load_state_dict(ck["trans"])
    res = {"arm": arm, "best_epoch": ck["epoch"], "lnc_only": lnc_only,
           **evaluate_chromosome(head, trans, ch, lnc_only, device, with_s1=with_s1, mea_gammas=mea_gammas,
                                 conf_grid=conf_grid, post_out=post_out, log=log, exact_post=exact_post)}
    os.makedirs(out, exist_ok=True)
    json.dump(res, open(f"{out}/result.json", "w"), indent=1)
    log(json.dumps({k: (v if not isinstance(v, dict) else {kk: vv for kk, vv in v.items() if kk != "curve"})
                    for k, v in res.items()}))
    return res


CONF_GRID = (0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)


def mea_specs(spec):
    """v1.9 development: "1,x2,pv1,pvx1" (or a sequence of gamma) -> [(name, gamma, level, log_gain)]: prefix "pv" =
    posterior-Viterbi (log gains), "x" = exon-level gains (g1.mea_duration).  Names MEA[X]_g<gamma> / PV[X]_g<gamma>."""
    out = []
    for it in ([x for x in spec.split(",") if x] if isinstance(spec, str) else list(spec)):
        it = str(it)
        lg = it.startswith("pv")
        it = it[2:] if lg else it
        lv = "exon" if it.startswith("x") else "state"
        g = float(it[1:] if lv == "exon" else it)
        out.append((("PV" if lg else "MEA") + ("X" if lv == "exon" else "") + f"_g{g:g}", g, lv, lg))
    return out


def summarize(genes, post, truth, conf_grid=()):
    """lnc >= 200 nt filter, gene confidence, gene-level PR / AUPRC and ruler L1 / L2 (protocol v1.4).
    conf_grid (v1.9 development): also ruler lnc L1 Sn / strict Pr of the lnc genes with confidence >= tau."""
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
    grid = []
    for tau in conf_grid:
        keep = [p for p, g in zip(preds, genes) if p["cls"] == "lnc" and g["conf"] >= tau]
        st = ruler.score_class("lnc", truth, keep, idx, 0)["L1_locus"]
        grid.append({"tau": tau, "n": len(keep), "sn": st["sn"], "pr": st["pr_strict"]})
    extra = {}
    if conf_grid:  # development readout: AUPRC with the gene confidence taken over exon bases only
        ex_conf = [{"strand": g["strand"], "exons": [tuple(x) for x in g["exons"]], "conf": chrom.exon_confidence(g, post)}
                   for g in genes if g["type"] == "l"]
        extra = {"conf_grid": grid, "auprc_lnc_exonconf": prcurve.pr_curve(ex_conf, truth)["auprc"]}
    return {**extra, "auprc_lnc": pr["auprc"], "best_locus_f1": pr["best_f1"], "best_point": pr["best_point"],
            "n_truth_lnc": pr["n_truth"], "n_pred_lnc": pr["n_pred"],
            "lnc_L1": {"sn": sl["L1_locus"]["sn"], "pr": sl["L1_locus"]["pr_strict"]},
            "lnc_L2": {"sn": sl["L2_intron_chain"]["sn"], "pr": sl["L2_intron_chain"]["pr_strict"]},
            "mrna_L1": {"sn": sm["L1_locus"]["sn"], "pr": sm["L1_locus"]["pr_strict"]},
            "mrna_L2": {"sn": sm["L2_intron_chain"]["sn"], "pr": sm["L2_intron_chain"]["pr_strict"]},
            "n_pred_mrna": sum(1 for g in genes if g["type"] == "m"), "curve": pr["curve"]}


def run(arm, seed, tr, va, ch, out, lnc_only=False, epochs=12, batch=8, lr=1e-3, device="cuda", log=print,
        pause=None, resume=False, accum=1, keep_epochs=False):
    """Arms R / S2 / S4 (G1 without boundary head).
    pause: optional callable checked before every batch; True aborts the epoch with Paused.  resume: keep last.pt
    after every epoch and continue from it (see screen_g1b.run).  keep_epochs: also save ep<NN>.pt every epoch."""
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
        if keep_epochs:
            torch.save({"head": head.state_dict(), "trans": trans.state_dict(), "epoch": ep + 1}, f"{out}/ep{ep + 1:02d}.pt")
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
