#!/usr/bin/env python3
"""S11 boundary-jitter gate (docs/phase1_results/PLAN_S11.md): GENATATOR's segmentation model inside intervals around
chr18 truth lnc genes whose ends are exact, padded outward or trimmed inward.  Run with the genatator env, cwd anywhere.

  s11_gate.py [--n 200] [--out outputs/PHASE1/S11]  -> <out>/GATE.json (+ per-gene rows in <out>/rows.jsonl, resumable)
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import torch

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
PIPE = f"{ROOT}/models/genatator-pipeline"
SEG = f"{ROOT}/models/genatator-caduceus-ps-multispecies-segmentation"
sys.path.insert(0, PIPE)
sys.path.insert(0, os.path.join(ROOT, "scripts/phase0"))
import genatator_core as G  # noqa: E402
import ruler  # noqa: E402  (stdlib-only copy of the ruler)

TRUTH = f"{ROOT}/outputs/PHASE1/DIAG_L2/chr18_v50.gtf"
FASTA = f"{ROOT}/data/human/chr18.fa"
MAX_SPAN = 60_000
CONDS = [("A", 0), ("B", 0), ("O1", 1000), ("O5", 5000), ("O10", 10000), ("O20", 20000), ("I1", -1000), ("I5", -5000)]


def read_fasta(path):
    seq = []
    for line in open(path):
        if not line.startswith(">"):
            seq.append(line.strip())
    return "".join(seq).upper()


def pick_genes(truth, n):
    rows = []
    for gid in sorted(truth):
        g = truth[gid]
        if g["cls"] != "lnc":
            continue
        multi = {t: ruler._merge(ex) for t, ex in g["tx"].items() if len(ruler._merge(ex)) > 1}
        if not multi:
            continue
        gs = min(e[0] for t in g["tx"].values() for e in t)
        ge = max(e[1] for t in g["tx"].values() for e in t)
        if ge - gs + 1 > MAX_SPAN:
            continue
        tid = max(sorted(multi), key=lambda t: multi[t][-1][1] - multi[t][0][0])
        rows.append({"gene": gid, "strand": g["strand"], "gene_span": [gs, ge],
                     "tx_span": [multi[tid][0][0], multi[tid][-1][1]],
                     "chains": sorted({ruler.introns(ex) for ex in multi.values()}),
                     "introns": sorted({i for ex in multi.values() for i in ruler.introns(ex)})})
    rng = np.random.default_rng(0)
    if len(rows) > n:
        rows = [rows[i] for i in sorted(rng.choice(len(rows), n, replace=False))]
    return rows


def interval(row, cond, d, chrom_len):
    s, e = row["gene_span"] if cond == "B" else row["tx_span"]  # 1-based closed
    s, e = s - d, e + d
    s, e = max(1, s), min(chrom_len, e)
    return (s, e) if e - s + 1 >= 200 else None


def segment(seq, s, e, strand, tok, model, dev):
    """1-based closed interval -> (exons 1-based closed, reason)."""
    sub = seq[s - 1:e]
    tracks = G.infer_segmentation_tracks_with_rc(sequence=sub, tokenizer=tok, model=model, chunk_size=250000, overlap=0,
                                                 batch_size=1, device=dev, use_reverse_complement=True,
                                                 apply_sigmoid=True)
    tracks = np.asarray(tracks, dtype=np.float32)
    res = G.build_segmentation_result(interval_start=s - 1, sequence=sub, strand=strand, segmentation_tracks=tracks,
                                      splice_filter=True, is_cds_binary=False)
    lab = res.labels_argmax
    if len(lab) and (int(lab[0]) == G.INTRON_IDX or int(lab[-1]) == G.INTRON_IDX):
        return [], "intronic_boundary"
    ex = sorted((a + 1, b) for a, b in res.exon_segments)  # 0-based half-open -> 1-based closed
    return ex, None if ex else "no_exons"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--out", default="outputs/PHASE1/S11")
    a = ap.parse_args()
    out = f"{ROOT}/{a.out}"
    os.makedirs(out, exist_ok=True)
    truth = ruler.read_gtf(TRUTH, "chr18")
    genes = pick_genes(truth, a.n)
    seq = read_fasta(FASTA)
    done = {}
    rows_path = f"{out}/rows.jsonl"
    if os.path.exists(rows_path):
        for line in open(rows_path):
            r = json.loads(line)
            done[(r["gene"], r["cond"])] = r
    from transformers import AutoModelForTokenClassification, AutoTokenizer
    dev = torch.device("cuda")
    tok = AutoTokenizer.from_pretrained(SEG, trust_remote_code=True, local_files_only=True)
    model = AutoModelForTokenClassification.from_pretrained(SEG, trust_remote_code=True, local_files_only=True,
                                                            torch_dtype=torch.float32).to(dev).eval()
    t0 = time.time()
    with open(rows_path, "a") as fh:
        for k, row in enumerate(genes):
            for cond, d in CONDS:
                if (row["gene"], cond) in done:
                    continue
                iv = interval(row, cond, d, len(seq))
                r = {"gene": row["gene"], "cond": cond, "interval": iv, "exons": [], "reason": "too_short" if iv is None else None}
                if iv is not None:
                    with torch.no_grad():
                        r["exons"], r["reason"] = segment(seq, iv[0], iv[1], row["strand"], tok, model, dev)
                    torch.cuda.empty_cache()
                done[(row["gene"], cond)] = r
                fh.write(json.dumps(r) + "\n")
                fh.flush()
            if k % 10 == 0:
                print(f"gene {k + 1}/{len(genes)} {time.time() - t0:.0f}s", flush=True)

    summ = {}
    for cond, _ in CONDS:
        chain = tp = npred = ntruth = rej = single = 0
        for row in genes:
            r = done[(row["gene"], cond)]
            ex = [tuple(x) for x in r["exons"]]
            pi = set(ruler.introns(ex)) if len(ex) > 1 else set()
            ti = {tuple(i) for i in row["introns"]}
            chain += bool(pi) and tuple(ruler.introns(ex)) in {tuple(tuple(i) for i in c) for c in row["chains"]}
            tp += len(pi & ti)
            npred += len(pi)
            ntruth += len(ti)
            rej += r["reason"] is not None
            single += r["reason"] is None and len(ex) == 1
        n = len(genes)
        summ[cond] = {"n_genes": n, "chain_correct": chain / n, "junction_precision": tp / npred if npred else None,
                      "junction_recall": tp / ntruth if ntruth else None, "rejected": rej / n, "single_exon": single / n}
    ref = max(summ["A"]["chain_correct"], summ["B"]["chain_correct"])
    c = {k: v["chain_correct"] for k, v in summ.items()}
    if ref < 0.10:
        verdict = "stop: reference chain accuracy < 0.10"
    elif c["O10"] < 0.2 * ref:
        verdict = "stop: O10 < 0.2 x reference"
    elif c["O10"] >= 0.5 * ref and c["O20"] >= 0.5 * ref and c["I1"] >= 0.5 * ref:
        verdict = "pass: go to predicted loci"
    else:
        verdict = "boundary-sensitive: do not proceed"
    res = {"plan": "docs/phase1_results/PLAN_S11.md", "n_genes": len(genes), "max_span": MAX_SPAN, "summary": summ,
           "reference": ref, "verdict": verdict, "seconds": round(time.time() - t0, 1)}
    json.dump(res, open(f"{out}/GATE.json", "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
