#!/usr/bin/env python3
"""P1-G0 supplement (protocol v1.6 section 13): exact canonical event / intron coordinates per window. CPU on ont.

data/p1g/{train,val}_exact.npz (window-relative 0-based coordinates; only v50 Ensembl_canonical lncRNA /
protein_coding transcripts, as for the event targets):
  ev_w, ev_pos, ev_ch, ev_k, ev_tx     events inside the window (ev_ch in prep_events.EV order; ev_k 0 = m, 1 = l)
  in_w, in_L, in_R, in_s, in_k, in_tx  introns with BOTH boundary bases inside the window (L = first intron
                                       base, R = last intron base, genomic order; in_s 0 = +, 1 = -)
  tx_w, tx_k, tx_s, tx_n_intron, tx_inside   transcripts overlapping the window (inside = whole transcript inside)
Used for G2-0 candidate thresholds / recall ceilings and to force true events into the G2 training lattice.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ncrna_annot.data import labels  # noqa: E402
from prep_events import CH, GTF, OUT, P05, WIN, tx_events  # noqa: E402


def main():
    wins = {sp: [(l.split("\t")[0], int(l.split("\t")[1])) for l in open(f"{P05}/windows_{sp}.tsv")] for sp in ("train", "val")}
    chroms = {c for w in wins.values() for c, _ in w}
    tx = labels.read_transcripts(GTF, chroms)
    for sp, w in wins.items():
        E = {k: [] for k in ("ev_w", "ev_pos", "ev_ch", "ev_k", "ev_tx")}
        I = {k: [] for k in ("in_w", "in_L", "in_R", "in_s", "in_k", "in_tx")}
        T = {k: [] for k in ("tx_w", "tx_k", "tx_s", "tx_n_intron", "tx_inside")}
        ntx = 0
        for c in sorted({c for c, _ in w}):
            canon = [t for t in tx.get(c, []) if t[3]]
            left = np.array([t[2][0][0] for t in canon])
            right = np.array([t[2][-1][1] for t in canon])
            for i, (wc, st) in enumerate(w):
                if wc != c:
                    continue
                lo, hi = st, st + WIN
                for j in np.flatnonzero((right > lo) & (left < hi)):
                    k, s, ex, _ = canon[j]
                    kk, ss = (1 if k == "l" else 0), (0 if s == "+" else 1)
                    for ch, ps in tx_events(s, ex).items():
                        for p in ps:
                            if lo <= p < hi:
                                for key, v in zip(E, (i, p - lo, CH[ch], kk, ntx)):
                                    E[key].append(v)
                    for (_, b), (a2, _) in zip(ex, ex[1:]):
                        Lp, Rp = b, a2 - 1
                        if lo <= Lp and Rp < hi:
                            for key, v in zip(I, (i, Lp - lo, Rp - lo, ss, kk, ntx)):
                                I[key].append(v)
                    for key, v in zip(T, (i, kk, ss, len(ex) - 1, int(ex[0][0] >= lo and ex[-1][1] <= hi))):
                        T[key].append(v)
                    ntx += 1
        np.savez(f"{OUT}/{sp}_exact.npz", **{k: np.asarray(v, np.int64) for d in (E, I, T) for k, v in d.items()})
        print(sp, "events", len(E["ev_w"]), "introns", len(I["in_w"]), "transcripts", ntx,
              "lnc introns", int(np.sum(np.asarray(I["in_k"]) == 1)), flush=True)


if __name__ == "__main__":
    main()
