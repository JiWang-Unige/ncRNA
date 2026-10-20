#!/usr/bin/env python3
"""P1.1 data prep (Phase 1 protocol v1 section 1), CPU on ont.

data/p1/
  train_windows.tsv            chrom, start (0-based), kind (lnc / pc / rand)   3,000 x 24,576 bp
  train_codes.npy              int8   [3000, 24576]  (A0 C1 G2 T3 N4)
  train_bits_{v50uev,v50f,v41uev,v41f}.npy   uint32 [3000, 24576]  allowed-state bitmasks
  chr18_codes.npy, chr18_bits_{...}.npy      whole internal-validation chromosome
  chr20_codes.npy                            test chromosome (no labels: scored with the ruler)
  prep_stats.json
"""
import gzip
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ncrna_annot.data import evidence, labels  # noqa: E402
from ncrna_annot.grammar import g1  # noqa: E402

ROOT = os.path.expanduser("~/ncrna-annot")
OUT = f"{ROOT}/data/p1"
HG38 = f"{ROOT}/data/human/hg38.fa.gz"
GTF = {"v50": f"{ROOT}/data/human/gencode.v50.annotation.gtf.gz", "v41": f"{ROOT}/data/human/gencode.v41.annotation.gtf.gz"}
TRAIN = [f"chr{i}" for i in list(range(1, 18)) + [19, 21, 22]] + ["chrX"]
VAL, TEST = "chr18", "chr20"
WIN, NTRAIN, MIX = 24576, 3000, (("lnc", 0.4), ("pc", 0.3), ("rand", 0.3))


def read_genome(chroms):
    seqs, cur, buf = {}, None, []
    with gzip.open(HG38, "rt") as fh:
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
    os.makedirs(f"{OUT}/junctions", exist_ok=True)
    chroms = set(TRAIN) | {VAL, TEST}
    seqs = read_genome(chroms)
    tx = {v: labels.read_transcripts(p, chroms) for v, p in GTF.items()}
    spans = {v: labels.gene_spans(p, chroms) for v, p in GTF.items()}
    stats = {"windows": {}, "unknown_bp": {}, "class_fraction": {}}
    codes = {c: g1.encode(seqs[c]) for c in TRAIN + [VAL]}
    # 1) windows (sampled once from v50 gene positions; identical for every label variant)
    rng = np.random.default_rng(0)
    pos = {"lnc": {}, "pc": {}}
    for c in TRAIN:
        for kind, k in (("lnc", "l"), ("pc", "m")):
            ex = [(a, b) for kk, _, exs in tx["v50"].get(c, []) if kk == k for a, b in exs]
            pos[kind][c] = np.array([(a + b) // 2 for a, b in ex], np.int64)
    lens = np.array([len(seqs[c]) for c in TRAIN], float)
    wins = []
    for kind, frac in MIX:
        want, got = int(round(NTRAIN * frac)), 0
        while got < want:
            c = TRAIN[rng.choice(len(TRAIN), p=lens / lens.sum())]
            if kind == "rand":
                st = int(rng.integers(0, len(seqs[c]) - WIN))
            else:
                if not len(pos[kind][c]):
                    continue
                st = int(rng.choice(pos[kind][c])) - WIN // 2
            st = min(max(st, 0), len(seqs[c]) - WIN)
            if (codes[c][st:st + WIN] == 4).sum() > WIN // 10:
                continue
            wins.append((c, st, kind))
            got += 1
    wins = [wins[i] for i in rng.permutation(len(wins))]
    with open(f"{OUT}/train_windows.tsv", "w") as fh:
        fh.writelines(f"{c}\t{st}\t{k}\n" for c, st, k in wins)
    np.save(f"{OUT}/train_codes.npy", np.stack([codes[c][st:st + WIN] for c, st, _ in wins]))
    # 2) labels per chromosome, sliced into windows, then discarded (memory)
    variants = [(v, arm) for v in ("v50", "v41") for arm in ("uev", "f")]
    W = {va: np.zeros((len(wins), WIN), np.uint32) for va in variants}
    for c in TRAIN + [VAL]:
        n = len(seqs[c])
        for v in ("v50", "v41"):
            unk = evidence.unknown_mask(c, n, spans[v].get(c, []), cache=f"{OUT}/junctions/{c}.tsv")
            stats["unknown_bp"][f"{c}_{v}"] = int(unk.sum())
            for arm in ("uev", "f"):
                b = labels.allowed_bits(n, tx[v].get(c, []), unk if arm == "uev" else None)
                if c == VAL:
                    np.save(f"{OUT}/chr18_bits_{v}{arm}.npy", b)
                for i, (wc, st, _) in enumerate(wins):
                    if wc == c:
                        W[(v, arm)][i] = b[st:st + WIN]
                del b
        print(c, "done", flush=True)
    for (v, arm), B in W.items():
        np.save(f"{OUT}/train_bits_{v}{arm}.npy", B)
        ig_only = float((B == labels.IG_BIT).mean())
        stats["class_fraction"][f"train_{v}{arm}"] = {
            "IG_only": ig_only, "IG_or_unknown_union": float(((B & labels.IG_BIT) > 0).mean()) - ig_only}
    np.save(f"{OUT}/chr18_codes.npy", codes[VAL])
    np.save(f"{OUT}/chr20_codes.npy", g1.encode(seqs[TEST]))
    stats["windows"] = {k: sum(1 for w in wins if w[2] == k) for k, _ in MIX}
    stats["chr18_len"], stats["chr20_len"] = len(seqs[VAL]), len(seqs[TEST])
    json.dump(stats, open(f"{OUT}/prep_stats.json", "w"), indent=1)
    print(json.dumps(stats["windows"]), json.dumps(stats["class_fraction"]), flush=True)


if __name__ == "__main__":
    main()
