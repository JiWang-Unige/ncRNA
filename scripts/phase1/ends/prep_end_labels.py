#!/usr/bin/env python3
"""PLAN_END.md section 2: class-specific end labels for the v2.0 training windows (CPU, ont).

data/p1ends/train_{ann,ev}.npy  int8 [3000, 24576, 8]  1 positive, 0 negative, -1 ignored
Channels (CH8): per strand s in (+, -): lncTSS, lncTES, mTSS, mTES (transcript orientation, positions as
prep_events.tx_events: TSS / TES = first / last transcribed base).
  ann  positives = Ensembl_canonical ends of that class +-10 bp (= the current event labels, split by class);
       ignored = non-canonical ends of that class +-50 bp, U-ev unknown bases.
  ev   a canonical end is supported when a same-strand FANTOM5 CAGE peak (TSS, +-200 bp) / PolyASite 2.0 cluster
       (TES, +-100 bp) is near; positive = the nearest supporting peak +-10 bp; unsupported canonical ends ignored
       +-1 kb; non-canonical ends +-50 bp and U-ev unknown bases ignored as in ann.
Paint order (later wins): 0 -> non-canonical ignore -> unknown ignore (on zeros) -> unsupported ignore (ev) -> positives.

  prep_end_labels.py [--out data/p1ends]
"""
import argparse
import json
import os
import sys

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
for d in ("scripts/phase1", "scripts/phase1/s3"):
    sys.path.insert(0, os.path.join(ROOT, d))
import prep_events as PE  # noqa: E402
import s3_features as S3  # noqa: E402
from ncrna_annot.data import evidence, labels  # noqa: E402

WIN, POS_W, IGN_W, UNSUP_W = 24576, 10, 50, 1000
CH8 = [f"{k}{e}{s}" for s in ("+", "-") for k in ("l", "m") for e in ("TSS", "TES")]
CAGE = f"{ROOT}/data/bioval/raw/fantom5_cage_peaks.bed.gz"
PAS = f"{ROOT}/data/bioval/raw/polyasite2_GRCh38.96.tsv.gz"
EV_D = {"TSS": 200, "TES": 100}


def nearest(sorted_pos, x, d):
    if not len(sorted_pos):
        return None
    i = np.searchsorted(sorted_pos, x)
    best = None
    for j in (i - 1, i):
        if 0 <= j < len(sorted_pos) and abs(int(sorted_pos[j]) - x) <= d:
            if best is None or abs(int(sorted_pos[j]) - x) < abs(best - x):
                best = int(sorted_pos[j])
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="data/p1ends")
    a = ap.parse_args()
    out_dir = f"{ROOT}/{a.out}"
    os.makedirs(out_dir, exist_ok=True)
    wins = [(l.split("\t")[0], int(l.split("\t")[1])) for l in open(f"{ROOT}/data/p1v20/train_windows.tsv")]
    codes = np.load(f"{ROOT}/data/p1v20/train_codes.npy", mmap_mode="r")
    assert len(wins) == len(codes), (len(wins), len(codes))
    chroms = {c for c, _ in wins}
    clen = PE.chrom_lengths(chroms)
    tx = labels.read_transcripts(PE.GTF, chroms)
    spans = labels.gene_spans(PE.GTF, chroms)
    outs = {arm: np.zeros((len(wins), WIN, len(CH8)), np.int8) for arm in ("ann", "ev")}
    stats = {arm: {ch: {"pos_events": 0, "unsupported_ignored": 0} for ch in CH8} for arm in ("ann", "ev")}
    for c in sorted(chroms):
        unk = evidence.unknown_mask(c, clen[c], spans.get(c, []), cache=f"{ROOT}/data/p1/junctions/{c}.tsv")
        peaks = {"TSS": S3.load_peaks(CAGE, c, 5, 6, 7),
                 "TES": S3.load_peaks(PAS, c, 5, 1, 2, skip_header=True, chrom_name=c.replace("chr", ""))}
        txs = tx.get(c, [])
        left = np.array([t[2][0][0] for t in txs])
        right = np.array([t[2][-1][1] for t in txs])
        for i, (wc, st) in enumerate(wins):
            if wc != c:
                continue
            lo, hi = st, st + WIN
            near = [txs[k] for k in np.flatnonzero((right + UNSUP_W >= lo) & (left - UNSUP_W < hi))]
            evs = [(t[0], t[1], t[3], PE.tx_events(t[1], t[2])) for t in near]
            for j, ch in enumerate(CH8):
                k, e, s = ch[0], ch[1:4], ch[4]
                for arm in ("ann", "ev"):
                    tr = np.zeros(WIN, np.int8)
                    for kk, ss, canon, ev in evs:
                        if kk == k and ss == s and not canon:
                            PE.paint(tr, ev.get(f"{e}{s}", []), IGN_W, lo, hi, -1)
                    tr[unk[lo:hi] & (tr == 0)] = -1
                    pos = []
                    for kk, ss, canon, ev in evs:
                        if kk != k or ss != s or not canon:
                            continue
                        for x in ev.get(f"{e}{s}", []):
                            if arm == "ann":
                                pos.append(x)
                                continue
                            y = nearest(peaks[e][s], x, EV_D[e])
                            if y is None:
                                PE.paint(tr, [x], UNSUP_W, lo, hi, -1)
                                if lo <= x < hi:
                                    stats[arm][ch]["unsupported_ignored"] += 1
                            else:
                                pos.append(y)
                    PE.paint(tr, pos, POS_W, lo, hi, 1)
                    stats[arm][ch]["pos_events"] += sum(lo <= x < hi for x in pos)
                    outs[arm][i, :, j] = tr
        print(c, "done", flush=True)
    for arm, o in outs.items():
        np.save(f"{out_dir}/train_{arm}.npy", o)
        for j, ch in enumerate(CH8):
            x = o[..., j]
            stats[arm][ch].update(pos_bases=int((x == 1).sum()), ign_bases=int((x == -1).sum()))
    stats["channels"] = CH8
    json.dump(stats, open(f"{out_dir}/stats.json", "w"), indent=1)
    print(json.dumps(stats, indent=1))


if __name__ == "__main__":
    main()
