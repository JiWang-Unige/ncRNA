#!/usr/bin/env python3
"""AUTH-E1 (docs/auth/e1/PLAN_E1.md, author approved 2026-11-24): molecule-linked endpoint labels from PacBio rep1
(ENCFF373TKM: ENCODE long-read pipeline, minimap2 -ax splice -uf --secondary=no, CCS/FLNC reads, primary only), and the
section 6 step 2 feasibility gate on the training chromosomes only.

Read ends: transcript strand = alignment strand (-uf); 5' end = first aligned base on that strand, 3' end = last.
Candidates: AUTH pairs cache (outputs/AUTH/pairs/<c>.npz, K20 ends, length >= 200), v41 same-strand coding pairs removed.
Per TSS group: cov = reads whose 5' end is within W of the TSS; per candidate end: n_link = of those, reads whose 3' end is
within W of the end.  Supported end: n_link >= 2.  Competing end: n_link == 0.  Group enters the main analysis if it has
>= 1 supported end, >= 1 competing end and cov >= 10.  W = 100 (main), 50 / 200 reported.
Gate: training groups < 162 (x 0.37 = test < 60) -> stop (inconclusive, archive), the test is not read.
  e1_labels.py --mode feasibility   -> outputs/AUTH_E1/FEASIBILITY.json, outputs/AUTH_E1/labels/<train chrom>.npz
"""
import argparse
import json
import os
import re
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import auth_lib as A                       # noqa: E402

BAM = f"{A.ROOT}/data/wtc11/encode/ENCFF373TKM.bam"
SAMTOOLS = os.path.expanduser("~/envs/lrtools/bin/samtools")
OUT = f"{A.ROOT}/outputs/AUTH_E1"
WINDOWS = (50, 100, 200)
MIN_LINK, MIN_COV, GATE_TRAIN = 2, 10, 162
CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def read_ends(c):
    """{strand: (five, three)} arrays of read 5' / 3' ends (1-based) on chromosome c."""
    A.guard(c)
    out = {"+": ([], []), "-": ([], [])}
    p = subprocess.Popen([SAMTOOLS, "view", "-F", "2308", BAM, c], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for line in p.stdout:
        f = line.split("\t", 6)
        flag, pos = int(f[1]), int(f[3])
        ref = sum(int(n) for n, op in CIG.findall(f[5]) if op in "MDN=X")
        end = pos + ref - 1
        if flag & 16:
            out["-"][0].append(end)
            out["-"][1].append(pos)
        else:
            out["+"][0].append(pos)
            out["+"][1].append(end)
    assert p.wait() == 0
    res = {}
    for st, (a, b) in out.items():
        a, b = np.array(a, np.int64), np.array(b, np.int64)
        o = np.argsort(a, kind="stable")
        res[st] = (a[o], b[o])
    return res


def label_chrom(c):
    z = np.load(f"{A.OUT}/pairs/{c}.npz")
    keep = z["X"][:, A.FI["pc_same"]] == 0
    ends = read_ends(c)
    n = len(z["s"])
    cov = {w: np.zeros(n, np.int64) for w in WINDOWS}
    link = {w: np.zeros(n, np.int64) for w in WINDOWS}
    for k, st in enumerate(A.STR):
        five, three = ends[st]
        idx = np.nonzero((z["strand"] == k) & keep)[0]
        for i in idx:
            t, e = int(z["tss"][i]), int(z["endpos"][i])
            for w in WINDOWS:
                lo, hi = np.searchsorted(five, t - w, "left"), np.searchsorted(five, t + w, "right")
                cov[w][i] = hi - lo
                link[w][i] = int(np.sum(np.abs(three[lo:hi] - e) <= w))
    return {"keep": keep, "strand": z["strand"], "row": z["row"], "tss": z["tss"], "endpos": z["endpos"],
            **{f"cov{w}": cov[w] for w in WINDOWS}, **{f"link{w}": link[w] for w in WINDOWS}}


def groups(L, w=100):
    """Qualifying TSS groups: list of index arrays (kept candidates of one TSS)."""
    out = []
    key = L["strand"].astype(np.int64) * 10 ** 9 + L["row"]
    for k in np.unique(key[L["keep"]]):
        idx = np.nonzero((key == k) & L["keep"])[0]
        lk, cv = L[f"link{w}"][idx], L[f"cov{w}"][idx]
        if (lk >= MIN_LINK).any() and (lk == 0).any() and cv.max() >= MIN_COV:
            out.append(idx)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True, choices=["feasibility"])
    ap.parse_args()
    os.makedirs(f"{OUT}/labels", exist_ok=True)
    per, tot = {}, {w: 0 for w in WINDOWS}
    for c in A.TRAIN:
        L = label_chrom(c)
        np.savez(f"{OUT}/labels/{c}.npz", **L)
        per[c] = {w: len(groups(L, w)) for w in WINDOWS}
        for w in WINDOWS:
            tot[w] += per[c][w]
        print(c, per[c], flush=True)
    res = {"rule": __doc__.split("Gate:")[0].strip(), "gate_train_min": GATE_TRAIN, "train_groups": tot, "per_chrom": per,
           "expected_test_groups_main": round(tot[100] * 312 / 842, 1), "gate_pass": tot[100] >= GATE_TRAIN}
    json.dump(res, open(f"{OUT}/FEASIBILITY.json", "w"), indent=1)
    print(json.dumps({k: res[k] for k in ("train_groups", "expected_test_groups_main", "gate_pass")}))


if __name__ == "__main__":
    main()
