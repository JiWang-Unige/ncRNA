#!/usr/bin/env python3
"""SegmentNT (human, 30 kb) chr20 baseline with author-style threshold decoding (docs/phase1_results/PREREG_BASELINES.md,
9608a3a; card baselines-chr20, user 2026-11-03).  One read of hg38 chr20; nothing is tuned.

  segnt_chr20.py infer  [--out outputs/PHASE1/BASELINES/segnt]   -> <out>/probs.npy float16 [L, 3]
                                                                   (protein_coding_gene, lncRNA, exon), infer.json
  segnt_chr20.py decode [--out ...]                               -> <out>/pred_t<thr>_<strand rule>.gff3, SCORE.json
  segnt_chr20.py time   [--codes data/splice_own/codes_chr21.npy]  (speed / memory check on chr21 windows, no chr20)

Windows: chr20 split at N into N-free segments (upper case); segments >= 30,000 bp tiled with 30,000-bp windows,
stride 15,000, last window aligned to the segment end; each base takes the window whose centre is nearest.  Shorter
segments are not run (probabilities 0).  Forward strand only.  Probabilities: softmax over the last logit axis, [..., 1].
"""
import argparse
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1"))
MODEL = f"{ROOT}/ext/segment_nt"
REVISION = "1048ad869036bb8e50a15d82b323ebeb52e8ab92"
WIN, STRIDE = 30_000, 15_000
FEATS = ("protein_coding_gene", "lncRNA", "exon")
BASES = np.array(list("ACGTN"))


def load(device="cuda"):
    import torch
    from transformers import AutoModel, AutoTokenizer
    tok = AutoTokenizer.from_pretrained(MODEL, trust_remote_code=True)
    model = AutoModel.from_pretrained(MODEL, trust_remote_code=True).to(device).eval()
    idx = [model.config.features.index(f) for f in FEATS]
    return tok, model, idx, torch


def run_window(tok, model, idx, torch, seq, device="cuda"):
    ids = tok.batch_encode_plus([seq], return_tensors="pt", padding="max_length", max_length=WIN // 6 + 1)["input_ids"]
    assert ids.shape[1] == WIN // 6 + 1, ids.shape
    ids = ids.to(device)
    with torch.no_grad():
        out = model(ids, attention_mask=ids != tok.pad_token_id, output_hidden_states=True)
        p = torch.softmax(out.logits, dim=-1)[0, :, idx, -1]  # [WIN, 3]
    assert p.shape[0] == WIN, p.shape
    return p.float().cpu().numpy()


def segments(codes):
    """N-free [start, end) segments (0-based)."""
    isn = np.concatenate([[True], codes == 4, [True]])
    d = np.diff(isn.astype(np.int8))
    return list(zip(np.flatnonzero(d == -1), np.flatnonzero(d == 1)))


def windows(a, b):
    offs = list(range(a, b - WIN + 1, STRIDE))
    if offs[-1] != b - WIN:
        offs.append(b - WIN)
    cen = [o + WIN / 2 for o in offs]
    cuts = [a] + [int((c1 + c2) // 2) for c1, c2 in zip(cen, cen[1:])] + [b]
    return [(o, cuts[i], cuts[i + 1]) for i, o in enumerate(offs)]  # window start, own [lo, hi)


def infer(out):
    import p15_read as R
    os.makedirs(f"{ROOT}/{out}", exist_ok=True)
    if R.gate_busy():
        sys.exit(75)
    codes = np.asarray(np.load(f"{ROOT}/data/p1/chr20_codes.npy", mmap_mode="r"))
    tok, model, idx, torch = load()
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    torch.cuda.set_per_process_memory_fraction(float(os.environ.get("NCRNA_GPU_FRAC", "0.30")), 0)
    L = len(codes)
    probs = np.zeros((L, 3), np.float16)
    segs = segments(codes)
    run = [(a, b) for a, b in segs if b - a >= WIN]
    meta = {"revision": REVISION, "weights_sha256": R.sha256(f"{MODEL}/pytorch_model.bin"), "L": int(L),
            "segments": len(segs), "segments_run": len(run), "bp_not_run_nonN": int(sum(b - a for a, b in segs if b - a < WIN)),
            "started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    t0, n = time.time(), 0
    for a, b in run:
        for o, lo, hi in windows(a, b):
            seq = "".join(BASES[codes[o:o + WIN]])
            p = run_window(tok, model, idx, torch, seq)
            probs[lo:hi] = p[lo - o:hi - o]
            n += 1
            if n % 200 == 0:
                print(f"{n} windows {time.time() - t0:.0f}s", flush=True)
                if R.gate_busy():  # restart from scratch next time (the read is cheap and stateless)
                    sys.exit(75)
    meta.update(windows=n, infer_s=round(time.time() - t0), max_mem_gb=round(torch.cuda.max_memory_allocated() / 2**30, 2))
    np.save(f"{ROOT}/{out}/probs.npy", probs)
    json.dump(meta, open(f"{ROOT}/{out}/infer.json", "w"), indent=1)
    print(json.dumps(meta), flush=True)


def runs_of(mask):
    """1-based closed intervals where mask is True."""
    d = np.diff(np.concatenate([[0], mask.astype(np.int8), [0]]))
    return [(int(s) + 1, int(e)) for s, e in zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1))]


PLUS = {("GT", "AG"), ("GC", "AG")}
MINUS = {("CT", "AC"), ("CT", "GC")}


def strand_of(exons, seq):
    v = 0
    for (_, e1), (s2, _) in zip(exons, exons[1:]):
        m = (seq[e1:e1 + 2], seq[s2 - 3:s2 - 1])  # first two / last two intron bases (1-based e1 + 1 .. s2 - 1)
        v += 1 if m in PLUS else -1 if m in MINUS else 0
    return "+" if v > 0 else "-" if v < 0 else None


def decode(probs, codes, thr, strand_rule):
    seq = "".join(BASES[codes])
    exon = probs[:, 2] >= thr
    out = []
    for cls, col, min_len in (("lnc", 1, 200), ("mrna", 0, 0)):
        for a, b in runs_of(probs[:, col] >= thr):
            ex = [(s + a - 1, e + a - 1) for s, e in runs_of(exon[a - 1:b])] or [(a, b)]
            if sum(e - s + 1 for s, e in ex) < min_len:
                continue
            st = strand_of(ex, seq)
            for s in ([st] if st else (["+", "-"] if strand_rule == "both" else [])):
                out.append({"cls": cls, "strand": s, "exons": ex})
    return out


def write_gff3(preds, path):
    with open(path, "w") as fh:
        fh.write("##gff-version 3\n")
        for i, p in enumerate(preds):
            a, b, k = p["exons"][0][0], p["exons"][-1][1], {"lnc": "lnc_RNA", "mrna": "mRNA"}[p["cls"]]
            fh.write(f"chr20\tSegmentNT\tgene\t{a}\t{b}\t.\t{p['strand']}\t.\tID=g{i}\n")
            fh.write(f"chr20\tSegmentNT\t{k}\t{a}\t{b}\t.\t{p['strand']}\t.\tID=g{i}.t1;Parent=g{i}\n")
            for j, (x, y) in enumerate(p["exons"]):
                fh.write(f"chr20\tSegmentNT\texon\t{x}\t{y}\t.\t{p['strand']}\t.\tID=g{i}.t1.e{j};Parent=g{i}.t1\n")


def decode_all(out):
    import p15_score as P
    from ncrna_annot.eval import ruler
    probs = np.load(f"{ROOT}/{out}/probs.npy").astype(np.float32)
    codes = np.asarray(np.load(f"{ROOT}/data/p1/chr20_codes.npy", mmap_mode="r"))
    truth = ruler.read_gtf(f"{ROOT}/data/human/gencode.v50.annotation.gtf.gz", "chr20")
    res = {}
    for thr in (0.5, 0.3, 0.7):
        for rule in ("both", "stranded"):
            name = f"t{thr:g}_{rule}"
            preds = decode(probs, codes, thr, rule)
            write_gff3(preds, f"{ROOT}/{out}/pred_{name}.gff3")
            fin, _ = P.score_preds(truth, [dict(p, id=str(i)) for i, p in enumerate(preds)])
            res[name] = P.norm(fin)
            sl, sm = fin["lnc"], fin["mrna"]
            print(name, json.dumps({"n_lnc": fin["n_pred_lnc"], "lnc_L1": [sl["L1_locus"]["sn"], sl["L1_locus"]["pr_strict"]],
                                    "n_mrna": fin["n_pred_mrna"], "mrna_L1": [sm["L1_locus"]["sn"], sm["L1_locus"]["pr_strict"]]}),
                  flush=True)
    json.dump(res, open(f"{ROOT}/{out}/SCORE.json", "w"), indent=1)


def time_check(codes_path):
    tok, model, idx, torch = load()
    torch.cuda.set_per_process_memory_fraction(float(os.environ.get("NCRNA_GPU_FRAC", "0.30")), 0)
    c = np.asarray(np.load(f"{ROOT}/{codes_path}", mmap_mode="r"))
    starts = [o for o in range(20_000_000, 40_000_000, 1_000_000) if not (c[o:o + WIN] == 4).any()][:6]
    t = []
    for o in starts:
        seq = "".join(BASES[c[o:o + WIN]])
        t0 = time.time()
        run_window(tok, model, idx, torch, seq)
        t.append(time.time() - t0)
    print(json.dumps({"features": list(model.config.features), "s_per_window": t,
                      "max_mem_gb": torch.cuda.max_memory_allocated() / 2**30}), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("infer", "decode", "time"))
    ap.add_argument("--out", default="outputs/PHASE1/BASELINES/segnt")
    ap.add_argument("--codes", default="data/splice_own/codes_chr21.npy")
    a = ap.parse_args()
    if a.cmd == "infer":
        infer(a.out)
    elif a.cmd == "decode":
        decode_all(a.out)
    else:
        time_check(a.codes)


if __name__ == "__main__":
    main()
