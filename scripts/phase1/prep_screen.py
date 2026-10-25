#!/usr/bin/env python3
"""P1-S prep (protocol v1.2: canonical-transcript labels + feasibility repair), CPU on ont.

data/p1s/
  {train,val}_codes.npy      int8   [n, 24576]  real sequence codes (hg38)
  {train,val}_bits_{v50uev,v50f,v41uev,v50uevlnc}.npy   uint32 [n, 24576]
  val_windows.tsv, val_truth_v50.gtf (v50 genes fully inside a val window, chr18)
"""
import gzip
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ncrna_annot.data import evidence, labels  # noqa: E402
from ncrna_annot.grammar import g1  # noqa: E402

ROOT = os.path.expanduser("~/ncrna-annot")
P05, OUT = f"{ROOT}/data/p05", f"{ROOT}/data/p1s_v12"
GTF = {"v50": f"{ROOT}/data/human/gencode.v50.annotation.gtf.gz", "v41": f"{ROOT}/data/human/gencode.v41.annotation.gtf.gz"}
WIN = 24576


def read_genome(chroms):
    seqs, cur, buf = {}, None, []
    with gzip.open(f"{ROOT}/data/human/hg38.fa.gz", "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur in chroms:
                    seqs[cur] = "".join(buf)
                cur, buf = line[1:].split()[0], []
            elif cur in chroms:
                buf.append(line.strip())
    if cur in chroms:
        seqs[cur] = "".join(buf)
    return seqs


def main():
    global OUT
    import argparse
    ap = argparse.ArgumentParser(description="v1.2 labels for a window set (defaults: the P1-S windows -> data/p1s_v12)")
    ap.add_argument("--windows-train", default=f"{P05}/windows_train.tsv")
    ap.add_argument("--windows-val", default=f"{P05}/windows_val.tsv")
    ap.add_argument("--splits", default="train,val")
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--check-train-codes", default=None, help="assert the regenerated train codes equal this .npy")
    a = ap.parse_args()
    OUT = a.out
    os.makedirs(OUT, exist_ok=True)
    src = {"train": a.windows_train, "val": a.windows_val}
    wins = {sp: [(l.split("\t")[0], int(l.split("\t")[1])) for l in open(src[sp])] for sp in a.splits.split(",")}
    chroms = {c for w in wins.values() for c, _ in w}
    seqs = read_genome(chroms)
    tx = {v: labels.read_transcripts(p, chroms) for v, p in GTF.items()}
    spans = {v: labels.gene_spans(p, chroms) for v, p in GTF.items()}
    variants = ("v50uev", "v50f", "v41uev", "v50uevlnc")
    arr = {(sp, va): np.zeros((len(w), WIN), np.uint32) for sp, w in wins.items() for va in variants}
    codes = {sp: np.zeros((len(w), WIN), np.int8) for sp, w in wins.items()}
    for c in sorted(chroms):
        n = len(seqs[c])
        cc = g1.encode(seqs[c])
        unk = {v: evidence.unknown_mask(c, n, spans[v].get(c, []), cache=f"{ROOT}/data/p1/junctions/{c}.tsv") for v in GTF}
        b = {"v50uev": labels.allowed_bits(n, tx["v50"].get(c, []), unk["v50"]),
             "v50f": labels.allowed_bits(n, tx["v50"].get(c, []), None),
             "v41uev": labels.allowed_bits(n, tx["v41"].get(c, []), unk["v41"]),
             "v50uevlnc": labels.allowed_bits(n, [t for t in tx["v50"].get(c, []) if t[0] == "l"], unk["v50"])}
        for sp, w in wins.items():
            for i, (wc, st) in enumerate(w):
                if wc == c:
                    codes[sp][i] = cc[st:st + WIN]
                    for va in variants:
                        arr[(sp, va)][i] = b[va][st:st + WIN]
        print(c, "done", flush=True)
    if a.check_train_codes and "train" in wins:
        ref = np.load(a.check_train_codes, mmap_mode="r")
        assert ref.shape == codes["train"].shape and np.array_equal(ref, codes["train"]), "train codes differ from " + a.check_train_codes
        print("train codes identical to", a.check_train_codes, flush=True)
    stats = {}
    for sp in wins:
        np.save(f"{OUT}/{sp}_codes.npy", codes[sp])
        for va in variants:
            A = arr[(sp, va)]
            rel, viol = [], []
            for i in range(len(A)):
                A[i], r, v = labels.repair(A[i], codes[sp][i])
                rel.append(r)
                viol.append(v)
            stats[f"{sp}_{va}"] = {"windows_infeasible_before_repair": int(sum(v > 0 for v in viol)),
                                   "violations_mean": float(np.mean(viol)), "relaxed_bp_mean": float(np.mean(rel)),
                                   "relaxed_frac": float(np.sum(rel) / A.size)}
            print(sp, va, stats[f"{sp}_{va}"], flush=True)
            np.save(f"{OUT}/{sp}_bits_{va}.npy", A)
    import json
    json.dump(stats, open(f"{OUT}/repair_stats.json", "w"), indent=1)
    if "val" not in wins:
        return
    with open(f"{OUT}/val_windows.tsv", "w") as fh:
        fh.writelines(f"{c}\t{st}\n" for c, st in wins["val"])
    # truth for window-level validation: v50 genes whose full span lies inside one val window
    iv = sorted((st, st + WIN) for c, st in wins["val"])
    keep = set()
    genes = {}
    with gzip.open(GTF["v50"], "rt") as fh:
        rows = [l for l in fh if l.startswith("chr18\t")]
    for l in rows:
        f = l.split("\t")
        if f[2] == "gene":
            gid = f[8].split('gene_id "')[1].split('"')[0]
            genes[gid] = (int(f[3]) - 1, int(f[4]))
    for gid, (s, e) in genes.items():
        if any(a <= s and e <= b for a, b in iv):
            keep.add(gid)
    with open(f"{OUT}/val_truth_v50.gtf", "w") as fh:
        for l in rows:
            if l.split('gene_id "')[1].split('"')[0] in keep:
                fh.write(l)
    print("val truth genes inside windows:", len(keep), flush=True)


if __name__ == "__main__":
    main()
