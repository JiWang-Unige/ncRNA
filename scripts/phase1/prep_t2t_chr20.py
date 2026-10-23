#!/usr/bin/env python3
"""T2T-CHM13 chr20 (NC_060944.1) codes for the GENATATOR-comparable test (protocol v1 section 1).

  prep_t2t_chr20.py [--fasta data/human/t2t_chr20.fa]
Downloads the sequence from NCBI E-utilities (efetch, ~66 MB) if --fasta is missing, then writes
data/p1/t2t_chr20_codes.npy (int8, A0 C1 G2 T3 N4, same encoding as prep_p1.py) and records length + md5.
Run on ont; afterwards `bash scripts/phase1/p1_feat_ont.sh` picks the chromosome up (its tiles were skipped before).
"""
import argparse
import hashlib
import json
import os
import sys
import urllib.request

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ncrna_annot.grammar import g1  # noqa: E402

ROOT = os.path.expanduser("~/ncrna-annot")
ACC = "NC_060944.1"
URL = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={ACC}&rettype=fasta&retmode=text"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", default=f"{ROOT}/data/human/t2t_chr20.fa")
    ap.add_argument("--out", default=f"{ROOT}/data/p1/t2t_chr20_codes.npy")
    a = ap.parse_args()
    if not os.path.exists(a.fasta):
        os.makedirs(os.path.dirname(a.fasta), exist_ok=True)
        print("downloading", URL, flush=True)
        urllib.request.urlretrieve(URL, a.fasta + ".part")
        os.rename(a.fasta + ".part", a.fasta)
    seq, header = [], None
    with open(a.fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                assert header is None, "expected exactly one record"
                header = line.strip()
            else:
                seq.append(line.strip())
    seq = "".join(seq).upper()
    assert ACC in header, header
    codes = g1.encode(seq)
    np.save(a.out, codes)
    meta = {"accession": ACC, "header": header, "length": len(seq), "n_frac": float((codes == 4).mean()),
            "fasta_md5": hashlib.md5(open(a.fasta, "rb").read()).hexdigest(), "codes": a.out}
    json.dump(meta, open(a.out[:-4] + ".json", "w"), indent=1)
    print(json.dumps(meta))


if __name__ == "__main__":
    main()
