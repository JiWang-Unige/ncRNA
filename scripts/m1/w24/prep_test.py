#!/usr/bin/env python3
"""M1 weeks 2-4 (docs/m1/PLAN_M1_W24.md section 6): sequences and parts of the test chromosomes for the DNA-only v2.0
reads.  chr9 already has v2.0 reads (outputs/PHASE1/V20/read/chr9).  chr12 and chr13 come from data/human/hg38.fa.gz
(md5 checked against data/human/MD5SUMS) and are cut in two equal halves (disk); chr20 is read whole from
data/human/chr20.fa.  Every part gets its own name (w24_*) so that no earlier cache or shard record is touched.
  prep_test.py -> data/human/{chr12,chr13}.fa, data/p1/w24_<part>_codes.npy + .json, outputs/M1/W24/PARTS.json
"""
import gzip
import hashlib
import json
import os
import sys

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
from ncrna_annot.grammar import g1  # noqa: E402

HG = f"{ROOT}/data/human/hg38.fa.gz"
HALVES = ("chr12", "chr13")


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def main():
    want = [l.split() for l in open(f"{ROOT}/data/human/MD5SUMS") if l.strip().endswith("hg38.fa.gz")]
    assert want and want[0][0] == md5(HG), "hg38.fa.gz md5 differs from MD5SUMS"
    seqs, cur = {}, None
    with gzip.open(HG, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                cur = line[1:].split()[0]
                cur = cur if cur in HALVES else None
                if cur:
                    seqs[cur] = []
            elif cur:
                seqs[cur].append(line.strip())
    parts = {}
    for c in HALVES:
        s = "".join(seqs[c])
        with open(f"{ROOT}/data/human/{c}.fa", "w") as out:
            out.write(f">{c}\n")
            for i in range(0, len(s), 60):
                out.write(s[i:i + 60] + "\n")
        L = len(s)
        cut = L // 2
        for x, (a, b) in zip("ab", ((0, cut), (cut, L))):
            parts[f"w24_{c}{x}"] = {"chrom": c, "offset0": a, "end0": b}
        seqs[c] = s
    s20 = "".join(l.strip() for l in open(f"{ROOT}/data/human/chr20.fa") if not l.startswith(">"))
    seqs["chr20"] = s20
    parts["w24_chr20"] = {"chrom": "chr20", "offset0": 0, "end0": len(s20)}
    for name, p in parts.items():
        codes = g1.encode(seqs[p["chrom"]][p["offset0"]:p["end0"]])
        out = f"{ROOT}/data/p1/{name}_codes.npy"
        np.save(out, codes)
        meta = dict(p, name=name, length=len(codes), n_frac=float((codes == 4).mean()))
        json.dump(meta, open(out[:-4] + ".json", "w"), indent=1)
        p["length"], p["codes_md5"] = len(codes), md5(out)
        print(json.dumps(meta), flush=True)
    os.makedirs(f"{ROOT}/outputs/M1/W24", exist_ok=True)
    json.dump({"hg38_md5": want[0][0], "parts": parts, "chr9": "outputs/PHASE1/V20/read/chr9 (existing v2.0 reads)"},
              open(f"{ROOT}/outputs/M1/W24/PARTS.json", "w"), indent=1)


if __name__ == "__main__":
    main()
