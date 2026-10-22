#!/usr/bin/env python3
"""P1-G0 (protocol v1.6 section 13): per-window event targets for the boundary head (G1B) and G2 candidates. CPU on ont.

data/p1g/
  {train,val}_events.npy  int8 [n, 24576, 8]   1 positive, 0 negative, -1 ignored
  events_stats.json
Channels (EV): per strand s in (+, -): donor, acceptor, TSS, TES (transcript orientation), at the base where
G1 places the corresponding state: + donor = first intron base (L), + acceptor = last intron base (R),
- donor = R position, - acceptor = L position; TSS / TES = first / last transcribed base of the transcript.
Positives: events of v50 Ensembl_canonical lncRNA / protein_coding transcripts, +-10 bp.
Ignored: events of the other (non-canonical) transcripts +-50 bp (unless positive); U-ev unknown bases (unless positive).
Windows are the P1-S windows (data/p05/windows_{train,val}.tsv), identical to data/p1s_v12.
"""
import gzip
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ncrna_annot.data import evidence, labels  # noqa: E402

ROOT = os.path.expanduser("~/ncrna-annot")
P05, OUT = f"{ROOT}/data/p05", f"{ROOT}/data/p1g"
GTF = f"{ROOT}/data/human/gencode.v50.annotation.gtf.gz"
WIN, POS_W, IGN_W = 24576, 10, 50
EV = [f"{k}{s}" for s in ("+", "-") for k in ("donor", "acceptor", "TSS", "TES")]
CH = {n: i for i, n in enumerate(EV)}


def chrom_lengths(chroms):
    n, cur, out = 0, None, {}
    with gzip.open(f"{ROOT}/data/human/hg38.fa.gz", "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur in chroms:
                    out[cur] = n
                cur, n = line[1:].split()[0], 0
            elif cur in chroms:
                n += len(line.rstrip("\n"))
    if cur in chroms:
        out[cur] = n
    return out


def tx_events(s, ex):
    """0-based positions per channel name for one transcript (exons sorted, half-open)."""
    ev = {}
    left, right = ex[0][0], ex[-1][1] - 1
    Lpos = [b for (_, b), _ in zip(ex, ex[1:])]           # first intron base
    Rpos = [a2 - 1 for (_, _), (a2, _) in zip(ex, ex[1:])]  # last intron base
    if s == "+":
        ev["TSS+"], ev["TES+"], ev["donor+"], ev["acceptor+"] = [left], [right], Lpos, Rpos
    else:
        ev["TSS-"], ev["TES-"], ev["donor-"], ev["acceptor-"] = [right], [left], Rpos, Lpos
    return ev


def paint(track, positions, w, lo, hi, value, only_if=None):
    for p in positions:
        a, b = max(p - w, lo), min(p + w + 1, hi)
        if a >= b:
            continue
        seg = track[a - lo:b - lo]
        if only_if is None:
            seg[:] = value
        else:
            seg[seg == only_if] = value


def main():
    os.makedirs(OUT, exist_ok=True)
    wins = {sp: [(l.split("\t")[0], int(l.split("\t")[1])) for l in open(f"{P05}/windows_{sp}.tsv")] for sp in ("train", "val")}
    chroms = {c for w in wins.values() for c, _ in w}
    clen = chrom_lengths(chroms)
    tx = labels.read_transcripts(GTF, chroms)
    spans = labels.gene_spans(GTF, chroms)
    stats = {}
    for sp, w in wins.items():
        out = np.zeros((len(w), WIN, len(EV)), np.int8)
        for c in sorted({c for c, _ in w}):
            unk = evidence.unknown_mask(c, clen[c], spans.get(c, []), cache=f"{ROOT}/data/p1/junctions/{c}.tsv")
            txs = tx.get(c, [])
            left = np.array([t[2][0][0] for t in txs])
            right = np.array([t[2][-1][1] for t in txs])
            for i, (wc, st) in enumerate(w):
                if wc != c:
                    continue
                lo, hi = st, st + WIN
                near = [txs[k] for k in np.flatnonzero((right + IGN_W >= lo) & (left - IGN_W < hi))]
                evs = [(t[3], tx_events(t[1], t[2])) for t in near]
                win = out[i]
                for ch, j in CH.items():
                    tr = np.zeros(WIN, np.int8)
                    for canon, ev in evs:
                        if not canon:
                            paint(tr, ev.get(ch, []), IGN_W, lo, hi, -1)
                    tr[unk[lo:hi] & (tr == 0)] = -1
                    for canon, ev in evs:
                        if canon:
                            paint(tr, ev.get(ch, []), POS_W, lo, hi, 1)
                    win[:, j] = tr
            print(sp, c, "done", flush=True)
        np.save(f"{OUT}/{sp}_events.npy", out)
        stats[sp] = {ch: {"pos_frac": float((out[..., j] == 1).mean()), "ign_frac": float((out[..., j] == -1).mean()),
                          "n_pos_events_approx": int((out[..., j] == 1).sum() // (2 * POS_W + 1))} for ch, j in CH.items()}
        print(sp, json.dumps(stats[sp]), flush=True)
    stats["channels"] = EV
    json.dump(stats, open(f"{OUT}/events_stats.json", "w"), indent=1)


if __name__ == "__main__":
    main()
