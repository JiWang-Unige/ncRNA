#!/usr/bin/env python3
"""M1 weeks 2-4: configs/m1_w24.json from configs/v20.json (same runs, weights, decoding, calibration) for the DNA-only
reads of the test parts (outputs/M1/W24/PARTS.json).

chr18 reproduction reference = the v2.0 chr18 reads (outputs/PHASE1/V20/read/chr18/<run>: post sha256 from
read_meta.json, genes.json sha256).  Code files: every file also listed in configs/v20.json must be unchanged.
  w24_config.py --commit <sha>
"""
import argparse
import glob
import hashlib
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
from ncrna_annot.features import generanno as G  # noqa: E402



def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--commit", required=True)
    a = ap.parse_args()
    pp = f"{ROOT}/configs/v20.json"
    prev = json.load(open(pp))
    runs = {}
    for name in prev["main"]["runs"].values():
        r = prev["runs"][name]
        assert all(sha256(f"{ROOT}/{c}") == h for c, h in zip(r["ckpts"], r["ckpt_sha256"])), name
        runs[name] = r
    ref = {}
    for name in runs:
        rd = f"{ROOT}/{prev['out_dir']}/chr18/{name}"
        m = json.load(open(f"{rd}/read_meta.json"))
        assert m["config_sha256"] == sha256(pp), name
        ref[name] = {"read": os.path.relpath(rd, ROOT), "post_sha256": m["post_sha256"], "genes_sha256": sha256(f"{rd}/genes.json")}
    targets = {"chr18": dict(prev["targets"]["chr18"], runs=sorted(runs))}
    pj = json.load(open(f"{ROOT}/outputs/M1/W24/PARTS.json"))
    for name, p in pj["parts"].items():
        codes = f"data/p1/{name}_codes.npy"
        L = len(np.load(f"{ROOT}/{codes}", mmap_mode="r"))
        targets[name] = {"feats": f"data/p1_feat/{name}.npy", "cover": f"data/p1_feat/{name}.cover.json",
                         "tiles": len(G.chromosome_tiles(L)), "codes": codes,
                         "codes_meta": json.load(open(f"{ROOT}/{codes[:-4]}.json")), "offset0": p["offset0"],
                         "chrom": p["chrom"], "seqid": name, "runs": sorted(runs)}
    code = sorted(glob.glob(f"{ROOT}/src/ncrna_annot/**/*.py", recursive=True)
                  + [f"{ROOT}/scripts/phase1/p15_read.py", f"{ROOT}/scripts/phase1/p16_read.py"])
    files = {os.path.relpath(p, ROOT): sha256(p) for p in code}
    changed = [p for p, h in files.items() if p in prev["code_files"] and prev["code_files"][p] != h]
    assert not changed, f"code differs from configs/v20.json: {changed}"
    cfg = {"protocol": "docs/m1/PLAN_M1_W24.md (M1 weeks 2-4, DNA-only reads of the test parts); decoding and calibration of configs/v20.json",
           "previous": {"config": "configs/v20.json", "config_sha256": sha256(pp)}, "out_dir": "outputs/M1/W24/read",
           "code_commit": a.commit, "created": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "decode": prev["decode"], "calibration": prev["calibration"], "main": prev["main"], "runs": runs,
           "chr18_reference": ref, "targets": targets, "code_files": files}
    json.dump(cfg, open(f"{ROOT}/configs/m1_w24.json", "w"), indent=1)
    print(json.dumps({"runs": sorted(runs), "targets": {k: v["tiles"] for k, v in targets.items()}}))


if __name__ == "__main__":
    main()
