#!/usr/bin/env python3
"""v2.0 (protocol section 17 detail 6): chr9 sequence codes for the single read -> data/p1/chr9_codes.npy (+ .json).
The uppercase md5 of chr9 must equal the lock record configs/v20_testdata.json.  Sequence only: no annotation, no
model output."""
import hashlib
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/..")
import prep_p1 as P  # noqa: E402


def main():
    lock = json.load(open(f"{P.ROOT}/configs/v20_testdata.json"))
    seq = P.read_genome({"chr9"})["chr9"]
    md5 = hashlib.md5(seq.upper().encode()).hexdigest()
    assert len(seq) == lock["genome"]["chr9_len"] and md5 == lock["genome"]["chr9_md5_uppercase"], (len(seq), md5)
    codes = P.g1.encode(seq)
    with open(f"{P.ROOT}/data/human/chr9.fa", "w") as fh:   # soft-masked sequence for s3_features.py (repeats)
        fh.write(">chr9\n")
        fh.writelines(seq[i:i + 60] + "\n" for i in range(0, len(seq), 60))
    np.save(f"{P.ROOT}/data/p1/chr9_codes.npy", codes)
    meta = {"name": "chr9", "length": len(codes), "n_frac": float((codes == 4).mean()), "md5_uppercase": md5,
            "genome_sha256": lock["genome"]["sha256"]}
    json.dump(meta, open(f"{P.ROOT}/data/p1/chr9_codes.json", "w"), indent=1)
    print(json.dumps(meta), flush=True)


if __name__ == "__main__":
    main()
