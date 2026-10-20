"""Phase 1 training loop (protocol v1 section 2): partial-observation CRF loss, chr18 checkpoint selection.

loss = mean_b [ log Z(b) - log Z_obs(b) ] / L     (allowed-state masks from data/labels.py)
Each epoch: decode the whole validation chromosome, score with the ruler (v50, stranded) and keep
the checkpoint with the best lnc L1 F1.  One JSONL line per epoch to <out>/metrics.jsonl.
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


def lnc_l1_f1(gff, truth_gtf, chrom):
    truth = ruler.read_gtf(truth_gtf, chrom)
    preds = ruler.read_pred_gff3(gff, chrom)
    r = ruler.score_class("lnc", truth, preds, ruler.ExonIndex(truth), 0)["L1_locus"]
    m = ruler.score_class("mrna", truth, preds, ruler.ExonIndex(truth), 0)["L1_locus"]
    sn, pr = r["sn"] or 0.0, r["pr_strict"] or 0.0
    return {"lnc_L1_sn": sn, "lnc_L1_pr": pr, "lnc_L1_f1": 2 * sn * pr / (sn + pr) if sn + pr else 0.0,
            "mrna_L1_sn": m["sn"], "mrna_L1_pr": m["pr_strict"], "n_lnc_pred": sum(p["cls"] == "lnc" for p in preds)}


def train(train_feats, train_codes, train_bits, val, out, epochs=8, batch=8, lr=1e-3, seed=42,
          device="cuda", max_windows=None, log=print):
    """train_feats [N, 4096, D] (np memmap/array fp16), train_codes int8 [N, L], train_bits uint32 [N, L].
    val: dict(feats [ntok, D], codes int8 [L], truth_gtf, chrom, max_bp optional)."""
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    os.makedirs(out, exist_ok=True)
    N = len(train_codes) if max_windows is None else min(max_windows, len(train_codes))
    head = C0Head(train_feats.shape[-1]).to(device)
    trans = g1.G1Transitions().to(device)
    opt = torch.optim.AdamW(list(head.parameters()) + list(trans.parameters()), lr=lr, weight_decay=1e-2)
    best, t0 = (-1.0, None), time.time()
    for ep in range(epochs):
        head.train()
        tot, nb = 0.0, 0
        order = rng.permutation(N)
        for i in range(0, N, batch):
            idx = np.sort(order[i:i + batch])
            f = torch.as_tensor(np.asarray(train_feats[idx], np.float32), device=device)
            c_np = np.asarray(train_codes[idx])
            c = torch.as_tensor(c_np.astype(np.int64), device=device)
            gates = torch.as_tensor(np.stack([g1.gate_mask(x) for x in c_np]), device=device)
            allowed = torch.as_tensor(labels.bits_to_bool(np.asarray(train_bits[idx])), device=device)
            emit = head(f, c)
            T = trans()
            z = g1.log_partition(emit, T, gates)
            zo = g1.log_partition(emit, T, gates, allowed)
            loss = ((z - zo) / c.shape[1]).mean()
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(list(head.parameters()) + list(trans.parameters()), 1.0)
            opt.step()
            tot += float(loss)
            nb += 1
        # validation: whole chromosome (or prefix) decode + ruler
        vf, vc = val["feats"], val["codes"]
        if val.get("max_bp"):
            vf, vc = vf[:val["max_bp"] // 6], vc[:(val["max_bp"] // 6) * 6]
        pot = dec.chromosome_potentials(head, vf, vc, device)
        gl = dec.decode(pot, vc, trans().detach().double().cpu().numpy())
        gff = f"{out}/val_ep{ep + 1}.gff3"
        dec.write_gff3(gl, val["chrom"], gff)
        sc = lnc_l1_f1(gff, val["truth_gtf"], val["chrom"])
        rec = {"epoch": ep + 1, "loss_per_base": round(tot / max(nb, 1), 6), "t": round(time.time() - t0), **sc,
               "n_genes": len(gl), "gpu_mem_GB": round(torch.cuda.max_memory_allocated() / 1e9, 2) if device == "cuda" else 0}
        with open(f"{out}/metrics.jsonl", "a") as fh:
            fh.write(json.dumps(rec) + "\n")
        log(json.dumps(rec))
        if sc["lnc_L1_f1"] > best[0]:
            best = (sc["lnc_L1_f1"], ep + 1)
            torch.save({"head": head.state_dict(), "trans": trans.state_dict(), "epoch": ep + 1}, f"{out}/best.pt")
    json.dump({"best_lnc_L1_f1_val": best[0], "best_epoch": best[1], "epochs": epochs, "seed": seed,
               "n_train_windows": N, "elapsed_s": round(time.time() - t0)}, open(f"{out}/run.json", "w"), indent=1)
    return best
