#!/usr/bin/env python3
"""M1 audit option A (docs/m1/PLAN_M1_UNSEEN.md): giraffe chr14 (NC_137288.1) -> three equal parts a / b / c.

Writes data/p1/giraffe_14{a,b,c}_codes.npy + .json (g1.encode, as prep_panel.py), data/giraffe/chr14.gff3 (the
NC_137288.1 rows of the RefSeq GFF) and data/giraffe/PARTS.json (0-based half-open part coordinates, file md5s).
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

SEQID = "NC_137288.1"
FA = f"{ROOT}/data/giraffe/{SEQID}.fa"
GFF = f"{ROOT}/data/giraffe/genomic.gff.gz"
GFF_MD5 = "8876963915e310b064fe1007f21c2948"   # NCBI FTP md5checksums.txt
LENGTH = 72_936_098                            # NCBI sequence report


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def main():
    assert md5(GFF) == GFF_MD5
    with open(FA) as fh:
        header = fh.readline().strip()
        seq = "".join(line.strip() for line in fh)
    assert header.split()[0] == ">" + SEQID and ">" not in seq and len(seq) == LENGTH, (header, len(seq))
    cuts = [0, LENGTH // 3, 2 * LENGTH // 3, LENGTH]
    parts = {}
    for i, x in enumerate("abc"):
        name = f"giraffe_14{x}"
        s0, e0 = cuts[i], cuts[i + 1]
        codes = g1.encode(seq[s0:e0])
        out = f"{ROOT}/data/p1/{name}_codes.npy"
        np.save(out, codes)
        meta = {"name": name, "seqid": SEQID, "offset0": s0, "end0": e0, "length": len(codes),
                "n_frac": float((codes == 4).mean()), "fasta": os.path.relpath(FA, ROOT), "fasta_md5": md5(FA)}
        json.dump(meta, open(out[:-4] + ".json", "w"), indent=1)
        parts[name] = [s0, e0]
        print(json.dumps(meta), flush=True)
    n = 0
    with gzip.open(GFF, "rt") as fh, open(f"{ROOT}/data/giraffe/chr14.gff3", "w") as out:
        for line in fh:
            if line.startswith(SEQID + "\t"):
                out.write(line)
                n += 1
    json.dump({"seqid": SEQID, "length": LENGTH, "parts": parts, "fasta_md5": md5(FA), "gff_md5": GFF_MD5,
               "chr14_gff3_rows": n, "chr14_gff3_md5": md5(f"{ROOT}/data/giraffe/chr14.gff3")},
              open(f"{ROOT}/data/giraffe/PARTS.json", "w"), indent=1)
    print("gff rows", n)


if __name__ == "__main__":
    main()
