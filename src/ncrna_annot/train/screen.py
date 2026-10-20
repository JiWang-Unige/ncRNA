"""P1-S screening loop (protocol v1.1 section 8): proxy features, window-level validation on chr18.

arms: R (union grammar, Viterbi) also decoded as S1 (posterior argmax, no grammar constraint),
      S2 (lnc-only: mRNA states masked), S4 (F labels), S5 (v41 labels) -> selected by bits/flags.
"""
import json
import os
import time

import numpy as np
import torch

from ..data import labels
from ..decode import genes as dec
from ..eval import ruler
from ..grammar import g1
from ..heads.c0 import C0Head

M_STATES = [g1.IDX[f"m{p}{s}"] for s in g1.STRANDS for p in g1.PARTS]


def _emit(head, f, c, lnc_only):
    e = head(f, c)
    if lnc_only:
        e = e.clone()
        e[..., M_STATES] = g1.NEG
    return e


def val_decode(head, trans, feats, codes, starts, chrom, gff, lnc_only=False, mode="viterbi", device="cuda"):
    head.eval()
    T = trans().detach()
    allg = []
    with torch.no_grad():
        for i in range(len(codes)):
            f = torch.as_tensor(np.asarray(feats[i], np.float32), device=device)[None]
            c_np = np.asarray(codes[i])
            c = torch.as_tensor(c_np.astype(np.int64), device=device)[None]
            e = _emit(head, f, c, lnc_only)
            gates = g1.gate_mask(c_np)
            if mode == "viterbi":
                path, _ = g1.viterbi(e[0].cpu().numpy(), T.double().cpu().numpy(), gates)
            else:  # S1: per-base argmax of posterior marginals (no grammar constraint on the output)
                with torch.enable_grad():
                    m = g1.marginals(e, T, torch.as_tensor(gates, device=device)[None])
                path = m[0].argmax(-1).cpu().numpy().astype(np.int8)
            allg += g1.path_to_genes(path, offset=int(starts[i]))
    seen, uniq = set(), []
    for g in allg:  # overlapping windows can emit the same gene twice
        k = (g["type"], g["strand"], tuple(map(tuple, g["exons"])))
        if k not in seen:
            seen.add(k)
            uniq.append(g)
    dec.write_gff3(uniq, chrom, gff)
    return uniq


def score(gff, truth_gtf, chrom):
    truth = ruler.read_gtf(truth_gtf, chrom)
    preds = ruler.read_pred_gff3(gff, chrom)
    idx = ruler.ExonIndex(truth)
    out = {}
    for cls in ("lnc", "mrna"):
        s = ruler.score_class(cls, truth, preds, idx, 0)
        sn, pr = s["L1_locus"]["sn"] or 0.0, s["L1_locus"]["pr_strict"] or 0.0
        out[cls] = {"L1_sn": sn, "L1_pr": pr, "L1_f1": 2 * sn * pr / (sn + pr) if sn + pr else 0.0,
                    "L2_sn": s["L2_intron_chain"]["sn"], "L2_pr": s["L2_intron_chain"]["pr_strict"],
                    "n_truth": s["n_truth_genes"], "n_pred": s["n_pred_tx"]}
    return out


def run(arm, tr, va, out, lnc_only=False, epochs=6, batch=8, lr=1e-3, seed=42, device="cuda",
        extra_modes=(), log=print):
    """tr = dict(feats, codes, bits); va = dict(feats, codes, starts, truth_gtf, chrom)."""
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    os.makedirs(out, exist_ok=True)
    N = len(tr["codes"])
    head = C0Head(tr["feats"].shape[-1]).to(device)
    trans = g1.G1Transitions().to(device)
    params = list(head.parameters()) + list(trans.parameters())
    opt = torch.optim.AdamW(params, lr=lr, weight_decay=1e-2)
    best, t0 = (-1.0, 0), time.time()
    for ep in range(epochs):
        head.train()
        order, tot, nb = rng.permutation(N), 0.0, 0
        for i in range(0, N, batch):
            idx = np.sort(order[i:i + batch])
            f = torch.as_tensor(np.asarray(tr["feats"][idx], np.float32), device=device)
            c_np = np.asarray(tr["codes"][idx])
            c = torch.as_tensor(c_np.astype(np.int64), device=device)
            gates = torch.as_tensor(np.stack([g1.gate_mask(x) for x in c_np]), device=device)
            allowed = torch.as_tensor(labels.bits_to_bool(np.asarray(tr["bits"][idx])), device=device)
            e, T = _emit(head, f, c, lnc_only), trans()
            loss = ((g1.log_partition(e, T, gates) - g1.log_partition(e, T, gates, allowed)) / c.shape[1]).mean()
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(params, 1.0)
            opt.step()
            tot, nb = tot + float(loss), nb + 1
        gff = f"{out}/val_ep{ep + 1}.gff3"
        gl = val_decode(head, trans, va["feats"], va["codes"], va["starts"], va["chrom"], gff, lnc_only, "viterbi", device)
        sc = score(gff, va["truth_gtf"], va["chrom"])
        rec = {"arm": arm, "epoch": ep + 1, "loss_per_base": round(tot / max(nb, 1), 5), "t": round(time.time() - t0),
               "lnc": sc["lnc"], "mrna": sc["mrna"], "n_genes": len(gl)}
        with open(f"{out}/metrics.jsonl", "a") as fh:
            fh.write(json.dumps(rec) + "\n")
        log(json.dumps(rec))
        if sc["lnc"]["L1_f1"] > best[0]:
            best = (sc["lnc"]["L1_f1"], ep + 1)
            torch.save({"head": head.state_dict(), "trans": trans.state_dict(), "epoch": ep + 1}, f"{out}/best.pt")
    res = {"arm": arm, "best_lnc_L1_f1": best[0], "best_epoch": best[1], "seed": seed, "elapsed_s": round(time.time() - t0)}
    if extra_modes:
        ck = torch.load(f"{out}/best.pt", map_location=device)
        head.load_state_dict(ck["head"])
        trans.load_state_dict(ck["trans"])
        for mode in extra_modes:
            gff = f"{out}/val_best_{mode}.gff3"
            val_decode(head, trans, va["feats"], va["codes"], va["starts"], va["chrom"], gff, lnc_only, mode, device)
            res[f"decode_{mode}"] = score(gff, va["truth_gtf"], va["chrom"])
    json.dump(res, open(f"{out}/run.json", "w"), indent=1)
    return res
