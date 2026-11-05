#!/usr/bin/env python3
"""v3.0 mouse U-ev junctions (protocol section 18, execution detail 3): Snaptron srav1m (mm10) -> mm39.

Each junction's first and last intron base (Snaptron start / end, 1-based) is carried through mm10ToMm39.over.chain
with the S4 chain reader; a junction is kept when both bases map once, to the same mm39 chromosome and chain
orientation, and the mapped intron length differs from the original by <= 10%.  Kept rows are written in the cached
Snaptron row format of evidence.fetch_junctions (source, chrom, start, end, strand, left motif, right motif, samples,
coverage; motifs as in mm10) to data/v30/junctions_mm39/<chrom>.tsv; evidence.unknown_mask then applies the U-ev rule
unchanged.

  mouse_uev.py --out data/v30/junctions_mm39
"""
import argparse
import glob
import gzip
import json
import os
import sys

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1/s4"))
import s4_synteny as S  # noqa: E402

FLIP = {"+": "-", "-": "+"}


def lift_junction(ch, s, e, strand, max_rel=0.10):
    """1-based inclusive intron [s, e] on mm10 -> (qname index, strand, s', e') on mm39 (1-based), or a reason string."""
    a, b = S.map_interval(ch, s - 1, s), S.map_interval(ch, e - 1, e)
    if len(a) != 1 or len(b) != 1:
        return "unmapped"
    (qa, sa, pa, _), (qb, sb, pb, _) = a[0], b[0]
    if qa != qb or sa != sb:
        return "split"
    lo, hi = min(pa, pb), max(pa, pb)
    if abs((hi - lo) - (e - s)) > max_rel * (e - s):
        return "length"
    st = strand if sa > 0 or strand not in FLIP else FLIP[strand]
    return qa, st, lo + 1, hi + 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--srav1m", default="data/v30/srav1m")
    ap.add_argument("--chain", default="data/v30/chains/mm10ToMm39.over.chain.gz")
    ap.add_argument("--out", default="data/v30/junctions_mm39")
    a = ap.parse_args()
    os.makedirs(f"{ROOT}/{a.out}", exist_ok=True)
    stats, rows_by = {}, {}
    for path in sorted(glob.glob(f"{ROOT}/{a.srav1m}/*.tsv")):
        c10 = os.path.basename(path)[:-4]
        with gzip.open(f"{ROOT}/{a.chain}", "rt") as fh:
            ch = S.read_chains(fh, c10)
        st = {"rows": 0, "kept": 0, "unmapped": 0, "split": 0, "length": 0, "other_chrom": 0}
        for line in open(path):
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or not f[2].isdigit():
                continue
            st["rows"] += 1
            r = lift_junction(ch, int(f[2]), int(f[3]), f[4])
            if isinstance(r, str):
                st[r] += 1
                continue
            qi, strand, s2, e2 = r
            c39 = ch["qnames"][qi]
            if c39 != c10:
                st["other_chrom"] += 1
            st["kept"] += 1
            rows_by.setdefault(c39, []).append("\t".join([f[0], c39, str(s2), str(e2), strand] + f[5:9]))
        st["kept_frac"] = st["kept"] / max(st["rows"], 1)
        stats[c10] = st
        print(c10, st, flush=True)
    for c39, rows in rows_by.items():
        rows.sort(key=lambda r: (int(r.split("\t")[2]), int(r.split("\t")[3])))
        open(f"{ROOT}/{a.out}/{c39}.tsv", "w").write("\n".join(rows))
    json.dump(stats, open(f"{ROOT}/{a.out}/LIFT_STATS.json", "w"), indent=1)


if __name__ == "__main__":
    main()
