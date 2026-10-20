"""0.5 heads (PROTOCOL v2), shared by the Modal and the ont runners.

run(arm, seed, data_dir, feat_dir, out_root, smoke): arm in F / Uev / Upu.
data_dir holds labels/unknown/PRIOR_pi (prep_05.py); feat_dir holds feat_{train,val,chr20}.npy.
"""
import json
import os
import time

TOK, NT, K = 6, 4096, 512


def run(arm, seed, data_dir, feat_dir, out_root, smoke=False):
    import numpy as np
    import torch
    import torch.nn as nn
    import torch.nn.functional as Fn
    torch.manual_seed(seed)
    np.random.seed(seed)
    src = data_dir
    D = feat_dir
    out = f"{out_root}/{arm}_s{seed}"
    os.makedirs(out, exist_ok=True)
    Xtr = np.load(f"{D}/feat_train.npy")  # [n, 4096, 512] fp16
    Xva = np.load(f"{D}/feat_val.npy")
    n_tr, n_va = len(Xtr), len(Xva)
    Ytr = np.load(f"{src}/labels_train.npy")[:n_tr].astype(np.int64)
    Yva = np.load(f"{src}/labels_val.npy")[:n_va].astype(np.int64)
    Utr = np.load(f"{src}/unknown_train.npy")[:n_tr]
    pi = json.load(open(f"{src}/PRIOR_pi.json"))["pi"]
    epochs = 1 if smoke else 6

    class Head(nn.Module):
        def __init__(self):
            super().__init__()
            self.conv = nn.Conv1d(K, 256, 9, padding=4)
            self.gru = nn.GRU(256, 256, batch_first=True, bidirectional=True)
            self.out = nn.Linear(512, 5)

        def forward(self, x):  # x [B, T, K]
            h = Fn.gelu(self.conv(x.transpose(1, 2))).transpose(1, 2)
            return self.out(self.gru(h)[0])

    net = Head().cuda()
    opt = torch.optim.AdamW(net.parameters(), lr=1e-3, weight_decay=1e-2)

    def loss_fn(logits, y, unk):
        lp = Fn.log_softmax(logits.float(), -1)
        if arm == "F":
            return Fn.nll_loss(lp.reshape(-1, 5), y.reshape(-1))
        if arm == "Uev":
            ce = -lp.gather(-1, y[..., None])[..., 0]
            marg = -torch.logsumexp(lp[..., :3], -1)  # background | lnc exon | lnc intron
            return torch.where(unk, marg, ce).mean()
        # Upu: nnPU on the lnc-exon dimension; unlabeled = unannotated (label 0)
        pos, unl = y == 1, y == 0
        lab_other = (~pos) & (~unl)
        p1 = lp[..., 1].exp().clamp(1e-6, 1 - 1e-6)
        terms = []
        if lab_other.any():
            terms.append(-lp[lab_other].gather(-1, y[lab_other][:, None]).mean())
        if pos.any() and unl.any():
            r_pp = -torch.log(p1[pos]).mean()
            r_pn = -torch.log1p(-p1[pos]).mean()
            r_un = -torch.log1p(-p1[unl]).mean()
            terms.append(pi * r_pp + torch.clamp(r_un - pi * r_pn, min=0))
            terms.append(-(lp[..., 0][unl] - torch.log1p(-p1[unl])).mean())  # bg among non-lnc-exon
        return sum(terms)

    def predict(X, bs=8):
        net.eval()
        P = []
        with torch.no_grad():
            for i in range(0, len(X), bs):
                x = torch.as_tensor(X[i:i + bs].astype(np.float32), device="cuda")
                P.append(Fn.softmax(net(x).float(), -1)[..., 1].cpu().numpy())
        net.train()
        return np.concatenate(P)

    def best_f1(p, y):
        best = (0.0, 0.5)
        for t in np.arange(0.05, 0.96, 0.05):
            pr, tr = p >= t, y == 1
            tp = (pr & tr).sum()
            f1 = 2 * tp / max(pr.sum() + tr.sum(), 1)
            best = max(best, (float(f1), float(t)))
        return best

    t0, best, state = time.time(), (-1.0, 0.5, 0), None
    for ep in range(epochs):
        order = np.random.permutation(n_tr)
        tot = 0.0
        for i in range(0, n_tr, 8):
            idx = np.sort(order[i:i + 8])
            x = torch.as_tensor(Xtr[idx].astype(np.float32), device="cuda")
            y = torch.as_tensor(Ytr[idx], device="cuda")
            u = torch.as_tensor(Utr[idx], device="cuda")
            logits = net(x)  # fp32: torch 2.2 has no bf16 CUDA GRU kernel
            loss = loss_fn(logits, y, u)
            opt.zero_grad()
            loss.backward()
            nn.utils.clip_grad_norm_(net.parameters(), 1.0)
            opt.step()
            tot += float(loss)
        f1, thr = best_f1(predict(Xva), Yva)
        rec = {"epoch": ep + 1, "loss": round(tot / max(n_tr // 8, 1), 5), "val_f1_v41": round(f1, 4),
               "thr": thr, "t": round(time.time() - t0)}
        open(f"{out}/metrics.jsonl", "a").write(json.dumps(rec) + "\n")
        if f1 > best[0]:
            best, state = (f1, thr, ep + 1), {k: v.clone() for k, v in net.state_dict().items()}
    net.load_state_dict(state)
    # whole-chr20 inference in 4096-token chunks with 256-token margins
    C = np.load(f"{D}/feat_chr20.npy", mmap_mode="r")
    n, m = len(C), 256
    p = np.zeros(n, np.float32)
    net.eval()
    with torch.no_grad():
        for s in range(0, n, NT - 2 * m):
            lo, hi = max(s - m, 0), min(s + NT - m, n)
            x = torch.as_tensor(np.asarray(C[lo:hi], np.float32)[None], device="cuda")
            q = Fn.softmax(net(x).float(), -1)[0, :, 1].cpu().numpy()
            p[s:min(s + NT - 2 * m, n)] = q[s - lo:s - lo + min(NT - 2 * m, n - s)]
    on = p >= best[1]
    edges = np.flatnonzero(np.diff(np.concatenate([[0], on.astype(np.int8), [0]])))
    segs = [(a * TOK + 1, b * TOK) for a, b in zip(edges[::2], edges[1::2]) if (b - a) * TOK >= 50]
    with open(f"{out}/pred.gff3", "w") as fh:
        fh.write("##gff-version 3\n")
        for i, (a, b) in enumerate(segs):
            g = f"{arm}_s{seed}_g{i}"
            fh.write(f"chr20\tp05\tgene\t{a}\t{b}\t.\t.\t.\tID={g}\n"
                     f"chr20\tp05\tlnc_RNA\t{a}\t{b}\t.\t.\t.\tID={g}.t1;Parent={g}\n"
                     f"chr20\tp05\texon\t{a}\t{b}\t.\t.\t.\tParent={g}.t1\n")
    json.dump({"arm": arm, "seed": seed, "pi": pi, "best_epoch": best[2], "val_f1_v41": best[0],
               "threshold": best[1], "n_segments": len(segs), "pred_bp": int(sum(b - a + 1 for a, b in segs)),
               "chr20_tokens_pos": int(on.sum()), "elapsed_s": round(time.time() - t0, 1), "smoke": smoke},
              open(f"{out}/run.json", "w"), indent=1)
    np.save(f"{out}/p_lnc_exon_chr20.npy", p.astype(np.float16))
    return arm, seed, best
