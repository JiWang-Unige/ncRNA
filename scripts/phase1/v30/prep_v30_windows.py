#!/usr/bin/env python3
"""v3.0 training windows (protocol section 18, execution detail 4), with prep_p1.py's rule: 40% centred on a random
lnc exon midpoint, 30% on a protein-coding exon midpoint, 30% uniform; chromosome drawn by length; windows with
> 10% N redrawn; order permuted at the end.

  prep_v30_windows.py human   3,000 hg38 windows from the v2.0 training chromosomes (TRAIN without chr9), v50 exons,
                              numpy default_rng(3); a window overlapping any data/p1v20 training window by > 50% is
                              redrawn  -> data/p1v30/human_new_windows.tsv
  prep_v30_windows.py mouse   3,000 mm39 windows from chr1-17 and chrX (chr18 dev and chr19 test excluded), GENCODE
                              vM39 exons, default_rng(2)  -> data/p1v30/mouse_windows.tsv
"""
import gzip
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/..")
import prep_p1 as P  # noqa: E402

OUT = f"{P.ROOT}/data/p1v30"
MM39 = f"{P.ROOT}/data/v30/mouse/mm39.fa.gz"
VM39 = f"{P.ROOT}/data/batch2/mouse/gencode.vM39.annotation.gtf.gz"
MOUSE_TRAIN = [f"chr{i}" for i in range(1, 18)] + ["chrX"]


def read_fasta(path, chroms):
    seqs, cur, buf = {}, None, []
    with gzip.open(path, "rt") as fh:
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


def sample(rng, chroms, seqs, gtf, avoid=()):
    tx = P.labels.read_transcripts(gtf, set(chroms))
    pos = {"lnc": {}, "pc": {}}
    for c in chroms:
        for kind, k in (("lnc", "l"), ("pc", "m")):
            ex = [(a, b) for kk, _, exs, *_ in tx.get(c, []) if kk == k for a, b in exs]
            pos[kind][c] = np.array([(a + b) // 2 for a, b in ex], np.int64)
    codes = {c: P.g1.encode(seqs[c]) for c in chroms}
    lens = np.array([len(seqs[c]) for c in chroms], float)
    old = {}
    for c, st in avoid:
        old.setdefault(c, []).append(st)
    old = {c: np.array(sorted(v)) for c, v in old.items()}
    wins, redrawn = [], {"n": 0, "overlap": 0}
    for kind, frac in P.MIX:
        want, got = int(round(P.NTRAIN * frac)), 0
        while got < want:
            c = chroms[rng.choice(len(chroms), p=lens / lens.sum())]
            if kind == "rand":
                st = int(rng.integers(0, len(seqs[c]) - P.WIN))
            else:
                if not len(pos[kind][c]):
                    continue
                st = int(rng.choice(pos[kind][c])) - P.WIN // 2
            st = min(max(st, 0), len(seqs[c]) - P.WIN)
            if (codes[c][st:st + P.WIN] == 4).sum() > P.WIN // 10:
                redrawn["n"] += 1
                continue
            if c in old and len(old[c]) and np.abs(old[c] - st).min() < P.WIN // 2:   # overlap > 50%
                redrawn["overlap"] += 1
                continue
            wins.append((c, st, kind))
            got += 1
    return [wins[i] for i in rng.permutation(len(wins))], redrawn


def main():
    which = sys.argv[1]
    os.makedirs(OUT, exist_ok=True)
    if which == "human":
        chroms = [c for c in P.TRAIN if c != "chr9"]
        avoid = [(l.split("\t")[0], int(l.split("\t")[1])) for l in open(f"{P.ROOT}/data/p1v20/train_windows.tsv")]
        wins, red = sample(np.random.default_rng(3), chroms, P.read_genome(set(chroms)), P.GTF["v50"], avoid)
    elif which == "mouse":
        chroms = MOUSE_TRAIN
        wins, red = sample(np.random.default_rng(2), chroms, read_fasta(MM39, set(chroms)), VM39)
    else:
        raise SystemExit("usage: prep_v30_windows.py human|mouse")
    with open(f"{OUT}/{which}{'_new' if which == 'human' else ''}_windows.tsv", "w") as fh:
        fh.writelines(f"{c}\t{st}\t{k}\n" for c, st, k in wins)
    st = {"which": which, "n": len(wins), "kinds": {k: sum(w[2] == k for w in wins) for k in ("lnc", "pc", "rand")},
          "chroms": {c: sum(w[0] == c for w in wins) for c in chroms}, "redrawn": red}
    json.dump(st, open(f"{OUT}/{which}_windows_stats.json", "w"), indent=1)
    print(json.dumps(st), flush=True)


if __name__ == "__main__":
    main()
