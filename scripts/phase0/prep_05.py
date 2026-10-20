#!/usr/bin/env python3
"""0.5 data prep (PROTOCOL v2), CPU on ont.

Writes to data/p05/:
  windows_{train,val}.tsv          chrom<TAB>start (0-based), 24,576 bp each
  seqs_{train,val}.txt             one window sequence per line (non-ACGT -> A; grid stays 6-mer aligned)
  labels_{train,val}.npy           uint8 [n, 4096]   0 bg, 1 lnc exon, 2 lnc intron, 3 mRNA exon, 4 mRNA intron (v41)
  unknown_{train,val}.npy          bool  [n, 4096]   U-ev mask (token centre base inside an unknown interval)
  unknown_intervals.bed, evidence_stats.json, PRIOR_pi.json
"""
import gzip
import json
import os
import subprocess
import urllib.request
from collections import defaultdict

import numpy as np

ROOT = os.path.expanduser("~/ncrna-annot")
D = f"{ROOT}/data/p05"
GTF41 = f"{ROOT}/data/human/gencode.v41.annotation.gtf.gz"
HG38 = f"{ROOT}/data/human/hg38.fa.gz"
TRAIN = [f"chr{i}" for i in list(range(1, 18)) + [19, 21]] + ["chrX"]
VAL = ["chr18"]
WIN, TOK = 24576, 6
NT = WIN // TOK
SNAP = ("https://snaptron.cs.jhu.edu/srav3h/snaptron?regions={c}:1-{n}&rfilter=samples_count>:99"
        "&rfilter=coverage_sum>:499&fields=chromosome,start,end,strand,left_motif,right_motif,samples_count,coverage_sum")
CANON = {"GT-AG", "CT-AC", "GC-AG", "CT-GC", "AT-AC", "GT-AT"}


def read_genome(chroms):
    if not os.path.exists(HG38):
        urllib.request.urlretrieve("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz", HG38 + ".part")
        os.rename(HG38 + ".part", HG38)
    seqs, cur, buf = {}, None, []
    with gzip.open(HG38, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur in chroms:
                    seqs[cur] = "".join(buf).upper()
                cur, buf = line[1:].split()[0], []
            elif cur in chroms:
                buf.append(line.strip())
    if cur in chroms:
        seqs[cur] = "".join(buf).upper()
    return seqs


def read_v41(chroms):
    """per chrom: gene spans (all classes) and transcripts of lncRNA / protein_coding."""
    spans, tx = defaultdict(list), defaultdict(dict)
    gspan = {}
    with gzip.open(GTF41, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 9 or f[0] not in chroms or f[2] != "exon":
                continue
            a = f[8]
            gid = a.split('gene_id "')[1].split('"')[0]
            gt = a.split('gene_type "')[1].split('"')[0]
            s, e = int(f[3]) - 1, int(f[4])
            k = (f[0], gid)
            gspan[k] = (min(gspan.get(k, (s, e))[0], s), max(gspan.get(k, (s, e))[1], e))
            if gt in ("lncRNA", "protein_coding"):
                tid = a.split('transcript_id "')[1].split('"')[0]
                tx[f[0]].setdefault(tid, [gt, []])[1].append((s, e))
    for (c, _), se in gspan.items():
        spans[c].append(se)
    return spans, tx


def base_labels(n, txs):
    """Precedence (PROTOCOL v2): exon > intron, lnc > mRNA. Order: mRNA intron, mRNA exon, lnc intron
    (only over background / mRNA intron), lnc exon (over everything)."""
    lab = np.zeros(n, np.uint8)
    by = {"protein_coding": [], "lncRNA": []}
    for gt, ex in txs.values():
        by[gt].append(sorted(ex))
    for ex in by["protein_coding"]:
        for (_, e1), (s2, _) in zip(ex, ex[1:]):
            seg = lab[e1:s2]
            seg[seg == 0] = 4
    for ex in by["protein_coding"]:
        for s, e in ex:
            lab[s:e] = 3
    for ex in by["lncRNA"]:
        for (_, e1), (s2, _) in zip(ex, ex[1:]):
            seg = lab[e1:s2]
            seg[(seg == 0) | (seg == 4)] = 2
    for ex in by["lncRNA"]:
        for s, e in ex:
            lab[s:e] = 1
    return lab


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def subtract(iv, cut):
    cut, out = merge(cut), []
    for s, e in iv:
        cur = s
        for cs, ce in cut:
            if ce <= cur or cs >= e:
                continue
            if cs > cur:
                out.append([cur, cs])
            cur = max(cur, ce)
        if cur < e:
            out.append([cur, e])
    return out


def unknown_intervals(c, n, gene_spans, stats):
    url = SNAP.format(c=c, n=n)
    rows = urllib.request.urlopen(url, timeout=900).read().decode().splitlines()[1:]
    gs = merge(gene_spans)
    starts = np.array([s for s, _ in gs])
    keep = []
    for r in rows:
        f = r.split("\t")
        s, e = int(f[2]) - 1, int(f[3])
        motif = f"{f[5].upper()}-{f[6].upper()}"
        if motif not in CANON or e - s > 500_000:
            continue
        i = np.searchsorted(starts, e) - 1  # last gene span starting before e
        if i >= 0 and gs[i][1] > s:
            continue
        keep.append([max(s - 500, 0), min(e + 500, n)])
    stats[c] = {"junctions_returned": len(rows), "kept_outside_v41": len(keep)}
    return subtract(merge(keep), gs)


def pick(seqs, lab_by_c, chroms, n, rng):
    lnc_pos = {c: np.flatnonzero(lab_by_c[c] == 1) for c in chroms}
    lens = np.array([len(seqs[c]) for c in chroms], float)
    out = []
    while len(out) < n:
        c = chroms[rng.choice(len(chroms), p=lens / lens.sum())]
        if len(out) % 2 == 0 and len(lnc_pos[c]):
            st = int(rng.choice(lnc_pos[c])) - WIN // 2
        else:
            st = int(rng.integers(0, len(seqs[c]) - WIN))
        st = min(max(st, 0), len(seqs[c]) - WIN)
        if seqs[c][st:st + WIN].count("N") <= WIN // 10:
            out.append((c, st))
    return out


def main():
    os.makedirs(D, exist_ok=True)
    chroms = TRAIN + VAL
    seqs = read_genome(set(chroms))
    spans, tx = read_v41(set(chroms))
    labs, unk, stats = {}, {}, {}
    for c in chroms:
        labs[c] = base_labels(len(seqs[c]), tx[c])
        iv = unknown_intervals(c, len(seqs[c]), spans[c], stats)
        m = np.zeros(len(seqs[c]), bool)
        for s, e in iv:
            m[s:e] = True
        unk[c] = m
        stats[c]["unknown_bp"] = int(m.sum())
        stats[c]["unannotated_bp"] = int((labs[c] == 0).sum())
        with open(f"{D}/unknown_intervals.bed", "a") as fh:
            fh.writelines(f"{c}\t{s}\t{e}\n" for s, e in iv)
        print(c, stats[c], flush=True)
    json.dump(stats, open(f"{D}/evidence_stats.json", "w"), indent=1)
    rng = np.random.default_rng(0)
    centres = np.arange(TOK // 2, WIN, TOK)
    for split, cs, n in (("train", TRAIN, 2000), ("val", VAL, 300)):
        wins = pick(seqs, labs, cs, n, rng)
        L = np.stack([labs[c][st + centres] for c, st in wins])
        U = np.stack([unk[c][st + centres] for c, st in wins])
        np.save(f"{D}/labels_{split}.npy", L)
        np.save(f"{D}/unknown_{split}.npy", U)
        with open(f"{D}/windows_{split}.tsv", "w") as fh:
            fh.writelines(f"{c}\t{st}\n" for c, st in wins)
        with open(f"{D}/seqs_{split}.txt", "w") as fh:
            for c, st in wins:
                w = seqs[c][st:st + WIN]
                fh.write("".join(ch if ch in "ACGT" else "A" for ch in w) + "\n")
        print(split, "class fractions", np.bincount(L.ravel(), minlength=5) / L.size, "unknown frac", U.mean(), flush=True)
    L = np.load(f"{D}/labels_train.npy")
    pi = float((L == 1).sum() / max((L == 0).sum(), 1))
    json.dump({"pi": pi, "definition": "v41 lnc-exon tokens / unannotated tokens in training windows"},
              open(f"{D}/PRIOR_pi.json", "w"), indent=1)
    print("pi", pi)


if __name__ == "__main__":
    main()
