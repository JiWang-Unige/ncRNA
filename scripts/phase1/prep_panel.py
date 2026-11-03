#!/usr/bin/env python3
"""P1.6 zero-shot panel (protocol v1 section 1: chicken 5, zebrafish 22, fly 2L, worm I; Ensembl 116 / 63): codes.

  prep_panel.py
For each panel chromosome: check the FASTA's md5 against data/<species>/MD5SUMS (written when Phase 0 fetched it),
then write data/p1/<name>_codes.npy (int8, A0 C1 G2 T3 N4, g1.encode as prep_p1.py) and <name>_codes.json
(length, N fraction, md5).  Afterwards `NAMES="..." CHROM_ONLY=1 bash scripts/phase1/p1_feat_ont.sh` extracts features.
"""
import hashlib
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ncrna_annot.grammar import g1  # noqa: E402

ROOT = os.path.expanduser("~/ncrna-annot")
PANEL = {"chicken_5": ("chicken", "5.fa", "5"), "zebrafish_22": ("zebrafish", "22.fa", "22"),
         "fly_2L": ("fly", "2L.fa", "2L"), "worm_I": ("worm", "I.fa", "I")}


def main():
    for name, (sp, fa, seqid) in PANEL.items():
        path = f"{ROOT}/data/{sp}/{fa}"
        md5 = hashlib.md5(open(path, "rb").read()).hexdigest()
        want = {ln.split()[1].lstrip("./"): ln.split()[0] for ln in open(f"{ROOT}/data/{sp}/MD5SUMS") if ln.strip()}
        assert want[fa] == md5, (name, md5, want)
        seq, header = [], None
        with open(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    assert header is None, "expected exactly one record"
                    header = line.strip()
                else:
                    seq.append(line.strip())
        assert header.split()[0] == ">" + seqid, header
        codes = g1.encode("".join(seq))
        out = f"{ROOT}/data/p1/{name}_codes.npy"
        np.save(out, codes)
        meta = {"name": name, "seqid": seqid, "header": header, "length": len(codes), "n_frac": float((codes == 4).mean()),
                "fasta": os.path.relpath(path, ROOT), "fasta_md5": md5}
        json.dump(meta, open(out[:-4] + ".json", "w"), indent=1)
        print(json.dumps(meta), flush=True)


if __name__ == "__main__":
    main()
