#!/usr/bin/env python3
"""M1 experiment 3 (docs/m1/PLAN_M1_AG.md): start table, AlphaGenome fold assignment and Colab job windows (CPU).

Starts: exactly the DEDEND rows (dedend_metrics.rows_of: GENCODE v50 Ensembl_canonical lncRNA / protein_coding
transcripts, 0-based TSS, FANTOM5 CAGE +-200 bp, promoter CpG o/e) on chr8, chr9 (main), chr18, chr20 (secondary).
Folds (alphagenome/data/fold_intervals.py): model fold-k holds out two Borzoi blocks (validation, test).  A start is
read by model k when it lies inside an interval of k's held-out blocks AND is >= 430 kb from every interval of k's
training blocks (training windows are 1 Mb around an interval's centre plus shift augmentation); otherwise it is
excluded.  test_only marks starts inside an interval of k's test block (sensitivity analysis).
Windows as in S2-AG (1,048,576 bp, window k starts at k*524288 - 262144); each start goes to the window whose centre
is nearest, so start +- 20 kb lies >= 242 kb inside the window.
  ag_fold_prep.py -> outputs/M1/AG/{tss.json, job_windows.json, PREP.json}
"""
import gzip
import hashlib
import json
import os
import sys
from collections import defaultdict

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
for d in ("src", "scripts/phase1", "scripts/phase1/s3", "scripts/phase1/ends", "scripts/phase1/dedend"):
    sys.path.insert(0, os.path.join(ROOT, d))
import dedend_metrics as DM  # noqa: E402
import s3_features as S3  # noqa: E402
from ncrna_annot.grammar import g1  # noqa: E402

OUT = f"{ROOT}/outputs/M1/AG"
BED = f"{ROOT}/data/alphagenome/borzoi_sequences_human.bed.gz"
VALID = {"fold-0": "fold0", "fold-1": "fold3", "fold-2": "fold2", "fold-3": "fold6"}
TEST = {"fold-0": "fold1", "fold-1": "fold4", "fold-2": "fold5", "fold-3": "fold7"}
MARGIN = 430_000
CHROMS = ("chr8", "chr9", "chr18", "chr20")
W, STEP, H, HALF = 1 << 20, 1 << 19, 1 << 18, 20_000


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def read_bed():
    iv = defaultdict(lambda: defaultdict(list))
    with gzip.open(BED, "rt") as fh:
        for line in fh:
            c, s, e, f = line.split()[:4]
            iv[c][f].append((int(s), int(e)))
    return {c: {f: np.array(sorted(v), np.int64) for f, v in d.items()} for c, d in iv.items()}


def inside(arr, x):
    return bool(len(arr)) and bool(((arr[:, 0] <= x) & (x < arr[:, 1])).any())


def dist(arr, x):
    if not len(arr):
        return 10 ** 12
    d = np.maximum(np.maximum(arr[:, 0] - x, x - (arr[:, 1] - 1)), 0)
    return int(d.min())


def assign(iv, x):
    for m in TEST:
        held = (VALID[m], TEST[m])
        if not any(inside(iv.get(f, np.zeros((0, 2), np.int64)), x) for f in held):
            continue
        train = [f for f in iv if f not in held]
        dmin = min(dist(iv[f], x) for f in train) if train else 10 ** 12
        if dmin >= MARGIN:
            return m, inside(iv.get(TEST[m], np.zeros((0, 2), np.int64)), x), dmin
        return None, False, dmin
    return None, False, None


def main():
    os.makedirs(OUT, exist_ok=True)
    bed = read_bed()
    rows, windows, prep = [], defaultdict(lambda: defaultdict(lambda: defaultdict(list))), {"codes_sha256": {}, "chroms": {}}
    for c in CHROMS:
        seq = S3.read_fasta(f"{ROOT}/data/human/{c}.fa")
        prep["codes_sha256"][c] = sha256_bytes(g1.encode(seq).astype(np.uint8).tobytes())
        rr = DM.rows_of(c)
        L = len(seq)
        n = defaultdict(int)
        for i, r in enumerate(rr):
            m, test_only, dmin = assign(bed.get(c, {}), r["TSS"])
            r.update(id=f"{c}:{i}", chrom=c, model=m, test_only=bool(test_only), train_dist=dmin)
            n[(r["k"], m)] += 1
            if m is not None and HALF <= r["TSS"] < L - HALF:
                k = int(round((r["TSS"] - H) / STEP))
                off = r["TSS"] - (k * STEP - H)
                assert HALF <= off < W - HALF
                windows[m][c][str(k)].append([r["id"], off, r["s"]])
            rows.append(r)
        prep["chroms"][c] = {"L": L, "n_rows": len(rr), "by_class_model": {f"{k}|{m}": v for (k, m), v in sorted(n.items(), key=str)}}
        print(c, json.dumps(prep["chroms"][c]), flush=True)
    jw = {m: {c: dict(sorted(d.items(), key=lambda t: int(t[0]))) for c, d in v.items()} for m, v in sorted(windows.items())}
    prep["windows"] = {m: sum(len(d) for d in v.values()) for m, v in jw.items()}
    prep["starts"] = {m: sum(len(x) for d in v.values() for x in d.values()) for m, v in jw.items()}
    prep["bed_sha256"] = hashlib.sha256(open(BED, "rb").read()).hexdigest()
    prep["margin"] = MARGIN
    json.dump(rows, open(f"{OUT}/tss.json", "w"))
    json.dump(jw, open(f"{OUT}/job_windows.json", "w"))
    json.dump(prep, open(f"{OUT}/PREP.json", "w"), indent=1)
    print(json.dumps({k: prep[k] for k in ("windows", "starts", "bed_sha256")}), flush=True)


if __name__ == "__main__":
    main()
