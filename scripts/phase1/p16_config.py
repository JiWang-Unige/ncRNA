#!/usr/bin/env python3
"""P1.6 zero-shot panel (protocol v1 section 5, descriptive, no gate): write the read configuration on ont.

  p16_config.py --commit <git sha of the code> --previous configs/p15_v19b.json --out configs/p16_v19.json

The decoding (section 16, v1.9: exon-level MEA gamma* = 3, tau* = 0.6, exact posteriors, epoch 10-12 weight average)
and the three main runs are taken from the P1.5 batch-2 config; their weights must still have the recorded sha256.
Targets: chr18 (reproduction: the posteriors and MEAX genes of every main run must equal the P1.5 batch-2 chr18 read
bit for bit) and the four panel chromosomes (Ensembl 116 / 63 truth, as Phase 0 0.1a / 0.1b).  Reading code: sha256
of src/ncrna_annot/**/*.py, scripts/phase1/p15_read.py and p16_read.py (p16_read.py refuses other code); the scoring
script records its own sha256 in its output.  Outputs go to outputs/PHASE1/P16.
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

PANEL = {  # name: (truth GTF, seqid)
    "chicken_5": ("data/chicken/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.116.gtf.gz", "5"),
    "zebrafish_22": ("data/zebrafish/Danio_rerio.GRCz11.116.gtf.gz", "22"),
    "fly_2L": ("data/fly/Drosophila_melanogaster.BDGP6.54.63.gtf.gz", "2L"),
    "worm_I": ("data/worm/Caenorhabditis_elegans.WBcel235.63.gtf.gz", "I"),
    "mouse_chr19": ("data/batch2/mouse/gencode.vM39.annotation.gtf.gz", "chr19")}   # S5 (PLAN_BATCH2.md B), --only mouse_chr19


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--commit", required=True)
    ap.add_argument("--previous", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--only", default="chicken_5,zebrafish_22,fly_2L,worm_I", help="panel targets besides the chr18 reproduction")
    ap.add_argument("--out-dir", default="outputs/PHASE1/P16")
    a = ap.parse_args()
    pp = a.previous if os.path.isabs(a.previous) else f"{ROOT}/{a.previous}"
    prev = json.load(open(pp))
    main_runs = prev["main"]["runs"]
    runs = {}
    for name in main_runs.values():
        r = prev["runs"][name]
        bad = [c for c, h in zip(r["ckpts"], r["ckpt_sha256"]) if sha256(f"{ROOT}/{c}") != h]
        assert not bad, f"{name}: weights changed {bad}"
        runs[name] = r
    ref = {}
    for name in runs:
        rd = f"{ROOT}/{prev['out_dir']}/chr18/{name}"
        m = json.load(open(f"{rd}/read_meta.json"))
        assert m["config_sha256"] == sha256(pp), f"{name}: chr18 read of another config"
        ref[name] = {"read": os.path.relpath(rd, ROOT), "post_sha256": m["post_sha256"], "genes_sha256": sha256(f"{rd}/genes.json")}
    targets = {"chr18": dict(prev["targets"]["chr18"], runs=sorted(runs))}
    for name, (gtf, seqid) in PANEL.items():
        if name not in a.only.split(","):
            continue
        codes = f"data/p1/{name}_codes.npy"
        L = len(np.load(f"{ROOT}/{codes}", mmap_mode="r"))
        targets[name] = {"feats": f"data/p1_feat/{name}.npy", "cover": f"data/p1_feat/{name}.cover.json",
                         "tiles": len(G.chromosome_tiles(L)), "codes": codes, "codes_meta": json.load(open(f"{ROOT}/{codes[:-4]}.json")),
                         "truth": gtf, "truth_chrom": seqid, "seqid": seqid, "runs": sorted(runs)}
    code = sorted(glob.glob(f"{ROOT}/src/ncrna_annot/**/*.py", recursive=True)
                  + [f"{ROOT}/scripts/phase1/p15_read.py", f"{ROOT}/scripts/phase1/p16_read.py"])
    cfg = {"protocol": "PHASE1_PROTOCOL.md section 5 P1.6 (zero-shot panel, descriptive); decoding of section 16 (v1.9)",
           "previous": {"config": a.previous, "config_sha256": sha256(pp)}, "out_dir": a.out_dir,
           "code_commit": a.commit, "created": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "decode": prev["decode"], "main": prev["main"], "runs": runs, "chr18_reference": ref, "targets": targets,
           "code_files": {os.path.relpath(p, ROOT): sha256(p) for p in code}}
    out = a.out if os.path.isabs(a.out) else f"{ROOT}/{a.out}"
    json.dump(cfg, open(out, "w"), indent=1)
    print(json.dumps({"out": out, "runs": sorted(runs), "targets": {k: v["tiles"] for k, v in targets.items()}}, indent=1))


if __name__ == "__main__":
    main()
