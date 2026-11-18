#!/usr/bin/env python3
"""M1 weeks 2-4: integrity check of the test-read queue outputs after ont's unplanned reboot (2026-11-18, ~20:11 UTC;
the journal ends without a shutdown, so writes still in the page cache may have been lost).
  - data/p1/w24_*_codes.npy: md5 against outputs/M1/W24/PARTS.json
  - every finished read under outputs/M1/W24/read/w24_*/<run>/: read_meta.json and genes.json parse; post.npy is finite
    and every row sums to 1 (+-0.02) -- lost writes would show up as zero rows.  One chr9 v2.0 read (written weeks ago)
    is checked the same way as a control for the row-sum rule.
  - Tiberius GTFs of the test chromosomes: every non-comment line has 9 tab-separated fields.
  integrity_check.py -> JSON on stdout, exit 1 on any failure
"""
import glob
import hashlib
import json
import os
import sys

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
W = f"{ROOT}/outputs/M1/W24"


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def check_post(path):
    post = np.load(path, mmap_mode="r")
    bad = nonfin = 0
    first_bad = None
    for a in range(0, len(post), 4_000_000):
        s = np.asarray(post[a:a + 4_000_000], np.float32).sum(1)
        nonfin += int((~np.isfinite(s)).sum())
        idx = np.flatnonzero(~(np.abs(s - 1) <= 0.02))
        bad += len(idx)
        if len(idx) and first_bad is None:
            first_bad = int(a + idx[0])
    return {"rows": int(len(post)), "bad_rows": bad, "nonfinite": nonfin, "first_bad": first_bad}


def main():
    res = {"codes": {}, "reads": {}, "control": {}, "tiberius": {}}
    for name, p in json.load(open(f"{W}/PARTS.json"))["parts"].items():
        res["codes"][name] = md5(f"{ROOT}/data/p1/{name}_codes.npy") == p["codes_md5"]
    for d in sorted(glob.glob(f"{W}/read/w24_*/*/")):
        json.load(open(d + "read_meta.json"))
        r = check_post(d + "post.npy")
        r["n_genes"] = len(json.load(open(d + "genes.json")))
        r["ok"] = r["bad_rows"] == 0 and r["nonfinite"] == 0
        res["reads"][os.path.relpath(d, W)] = r
    ctl = sorted(glob.glob(f"{ROOT}/outputs/PHASE1/V20/read/chr9/*/post.npy"))[:1]
    for f in ctl:
        res["control"][os.path.relpath(f, ROOT)] = check_post(f)
    for f in sorted(glob.glob(f"{W}/tiberius/*.gtf")):
        n, ok = 0, True
        for line in open(f):
            if line.startswith("#"):
                continue
            n += 1
            ok &= len(line.rstrip("\n").split("\t")) == 9
        res["tiberius"][os.path.basename(f)] = {"lines": n, "ok": ok}
    res["all_ok"] = (all(res["codes"].values()) and all(r["ok"] for r in res["reads"].values())
                     and all(t["ok"] for t in res["tiberius"].values()))
    print(json.dumps(res, indent=1))
    sys.exit(0 if res["all_ok"] else 1)


if __name__ == "__main__":
    main()
