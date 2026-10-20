#!/usr/bin/env python3
"""Phase 0.4: per-layer linear probes of frozen DNA LMs for lncRNA structure (self-contained).

Train probes on hg38 chr21 windows, evaluate on chr20 windows; labels from GENCODE v50
(either strand). Targets per position: lnc_exon, lnc_intron, mrna_exon, mrna_intron,
lnc_splice (+-3 bp of a lnc donor/acceptor), lnc_tss / lnc_tes (+-50 bp).
Every model is read at the same positions: the centre of each 6-bp step (k-mer models: the
token centre; single-nucleotide models: every 6th base), so models see equal sample sizes.
Metric: average precision (AUPRC) per layer and target on the chr20 positions.

  probe_backbones.py --backbone ntv2 --data-dir DATA --out probe_ntv2.json [--n-train 300 --n-eval 300]
Downloads chr20/chr21 (UCSC) and GENCODE v50 into DATA when missing (Kaggle-friendly).
"""
import argparse
import gzip
import json
import os
import time
import urllib.request

import numpy as np
import torch

BACKBONES = {
    "ntv2": ("InstaDeepAI/nucleotide-transformer-v2-500m-multi-species", 6, 1),
    "generanno_cds": ("GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview", 6, 1),
    "generanno_base": ("GenerTeam/GENERanno-eukaryote-0.5b-base", 1, 1),
    "caduceus_ps": ("kuleshov-group/caduceus-ps_seqlen-131k_d_model-256_n_layer-16", 1, 0),
}  # name -> (repo, nt per token, leading special tokens)
TARGETS = ["lnc_exon", "lnc_intron", "mrna_exon", "mrna_intron", "lnc_splice", "lnc_tss", "lnc_tes"]
WIN, STEP = 6000, 6
UCSC = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/{}.fa.gz"
GTF = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gtf.gz"


def fetch(url, path):
    if not os.path.exists(path):
        urllib.request.urlretrieve(url, path + ".part")
        os.rename(path + ".part", path)


def read_fa(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        return "".join(l.strip() for l in fh if not l.startswith(">")).upper()


def labels(gtf, chrom, n):
    """uint8 bit-flags per base (index = TARGETS order)."""
    lab = np.zeros(n, np.uint8)
    tx = {}
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if not line.startswith(chrom + "\t"):
                continue
            f = line.split("\t")
            if f[2] != "exon":
                continue
            gt = f[8].split('gene_type "')[1].split('"')[0]
            if gt not in ("lncRNA", "protein_coding"):
                continue
            tid = f[8].split('transcript_id "')[1].split('"')[0]
            tx.setdefault(tid, [gt, f[6], []])[2].append((int(f[3]) - 1, int(f[4])))  # 0-based half-open
    for gt, strand, ex in tx.values():
        ex.sort()
        lnc = gt == "lncRNA"
        eb, ib = (0, 1) if lnc else (2, 3)
        for s, e in ex:
            lab[s:e] |= 1 << eb
        for (s1, e1), (s2, e2) in zip(ex, ex[1:]):
            lab[e1:s2] |= 1 << ib
            if lnc:
                for p in (e1, s2):
                    lab[max(p - 3, 0):p + 3] |= 1 << 4
        if lnc:
            tss, tes = (ex[0][0], ex[-1][1]) if strand == "+" else (ex[-1][1], ex[0][0])
            lab[max(tss - 50, 0):tss + 50] |= 1 << 5
            lab[max(tes - 50, 0):tes + 50] |= 1 << 6
    lab[(lab & 0b11) == 0b11] &= ~np.uint8(0b10)  # exon wins over intron (overlapping isoforms)
    lab[(lab & 0b1100) == 0b1100] &= ~np.uint8(0b1000)
    return lab


def pick_windows(seq, lab, n, rng):
    """half centred on random lncRNA exon bases, half uniform; drop windows with >10% N."""
    lnc_pos = np.flatnonzero(lab & 1)
    out = []
    while len(out) < n:
        if len(out) % 2 == 0 and len(lnc_pos):
            c = int(rng.choice(lnc_pos)) - WIN // 2
        else:
            c = int(rng.integers(0, len(seq) - WIN))
        c = min(max(c, 0), len(seq) - WIN)
        w = seq[c:c + WIN]
        if w.count("N") <= WIN // 10:
            out.append(c)
    return out


def load(name, device):
    from transformers import AutoModel, AutoModelForMaskedLM, AutoTokenizer
    repo = BACKBONES[name][0]
    tok = AutoTokenizer.from_pretrained(repo, trust_remote_code=True)
    model = None
    for cls in (AutoModelForMaskedLM, AutoModel):
        try:
            model = cls.from_pretrained(repo, trust_remote_code=True)
            break
        except Exception as e:  # noqa: BLE001 - try the next auto class
            last = e
    if model is None:
        raise last
    return tok, model.to(device).eval()


def hidden(name, tok, model, window, dtype):
    nt, lead = BACKBONES[name][1], BACKBONES[name][2]
    ids = torch.tensor([tok(window)["input_ids"]], device=model.device)
    with torch.no_grad(), torch.autocast("cuda", dtype=dtype):
        out = model(input_ids=ids, output_hidden_states=True)
    hs = out.hidden_states
    centres = np.arange(STEP // 2, WIN, STEP)
    idx = lead + (centres // nt)
    idx = idx[idx < hs[0].shape[1]]
    return [h[0, idx].float() for h in hs], centres[:len(idx)]


def avg_precision(score, y):
    o = np.argsort(-score)
    y = y[o]
    tp = np.cumsum(y)
    prec = tp / np.arange(1, len(y) + 1)
    return float((prec * y).sum() / max(y.sum(), 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--backbone", required=True, choices=list(BACKBONES))
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-train", type=int, default=300)
    ap.add_argument("--n-eval", type=int, default=300)
    ap.add_argument("--epochs", type=int, default=2)
    a = ap.parse_args()
    t0 = time.time()
    os.makedirs(a.data_dir, exist_ok=True)
    for c in ("chr20", "chr21"):
        fetch(UCSC.format(c), f"{a.data_dir}/{c}.fa.gz")
    fetch(GTF, f"{a.data_dir}/gencode.v50.annotation.gtf.gz")
    data = {}
    for c in ("chr21", "chr20"):
        seq = read_fa(f"{a.data_dir}/{c}.fa.gz")
        lab = labels(f"{a.data_dir}/gencode.v50.annotation.gtf.gz", c, len(seq))
        data[c] = (seq, lab)
    rng = np.random.default_rng(0)
    wins = {"chr21": pick_windows(*data["chr21"], a.n_train, rng),
            "chr20": pick_windows(*data["chr20"], a.n_eval, rng)}
    dev = "cuda"
    dtype = torch.bfloat16 if torch.cuda.is_bf16_supported() else torch.float16
    tok, model = load(a.backbone, dev)
    print(json.dumps({"setup_s": round(time.time() - t0, 1), "dtype": str(dtype)}), flush=True)

    def batch(c, start):
        seq, lab = data[c]
        feats, centres = hidden(a.backbone, tok, model, seq[start:start + WIN], dtype)
        bits = lab[start + centres]
        y = torch.tensor(((bits[:, None] >> np.arange(len(TARGETS))) & 1).astype(np.float32), device=dev)
        return feats, y

    f0, _ = batch("chr21", wins["chr21"][0])
    nl, hd = len(f0), f0[0].shape[-1]
    heads = torch.nn.ModuleList([torch.nn.Linear(hd, len(TARGETS)) for _ in range(nl)]).to(dev)
    opt = torch.optim.AdamW(heads.parameters(), lr=2e-3, weight_decay=1e-4)
    pos = torch.zeros(len(TARGETS), device=dev)
    for s in wins["chr21"][:50]:
        pos += batch("chr21", s)[1].sum(0)
    pw = ((50 * (WIN // STEP) - pos) / pos.clamp(min=1)).clamp(max=100)  # positive class weighting
    ln = lambda x: torch.nn.functional.layer_norm(x, (hd,))
    for ep in range(a.epochs):
        tot = 0.0
        for i in rng.permutation(len(wins["chr21"])):
            feats, y = batch("chr21", wins["chr21"][i])
            loss = sum(torch.nn.functional.binary_cross_entropy_with_logits(heads[k](ln(f)), y, pos_weight=pw)
                       for k, f in enumerate(feats)) / nl
            opt.zero_grad()
            loss.backward()
            opt.step()
            tot += float(loss)
        print(json.dumps({"epoch": ep + 1, "loss": round(tot / len(wins["chr21"]), 4),
                          "t": round(time.time() - t0)}), flush=True)
    scores, ys = [[] for _ in range(nl)], []
    with torch.no_grad():
        for s in wins["chr20"]:
            feats, y = batch("chr20", s)
            ys.append(y.cpu().numpy())
            for k, f in enumerate(feats):
                scores[k].append(heads[k](ln(f)).float().cpu().numpy())
    Y = np.concatenate(ys)
    res = {"backbone": a.backbone, "repo": BACKBONES[a.backbone][0], "n_layers": nl, "hidden": hd,
           "n_eval_positions": int(len(Y)), "prevalence": dict(zip(TARGETS, Y.mean(0).round(5).tolist())),
           "auprc": {}, "best": {}, "elapsed_s": round(time.time() - t0, 1),
           "gpu": torch.cuda.get_device_name(0), "dtype": str(dtype)}
    for k in range(nl):
        S = np.concatenate(scores[k])
        res["auprc"][k] = {t: round(avg_precision(S[:, j], Y[:, j]), 4) for j, t in enumerate(TARGETS)}
    for t in TARGETS:
        k = max(res["auprc"], key=lambda kk: res["auprc"][kk][t])
        res["best"][t] = {"layer": int(k), "auprc": res["auprc"][k][t], "prevalence": res["prevalence"][t]}
    json.dump(res, open(a.out, "w"), indent=1)
    print(json.dumps(res["best"]), flush=True)


if __name__ == "__main__":
    main()
