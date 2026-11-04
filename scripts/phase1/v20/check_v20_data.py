#!/usr/bin/env python3
"""v2.0 (protocol section 17, execution detail 2): the regenerated training arrays in data/p1v20 must equal
data/p1v12 bit for bit on the 2,846 rows that were not replaced; writes data/p1v20/train_new_codes.npy (the codes of
the 154 new windows, for feature extraction) and data/p1v20/CHECK.json.  Exit 1 on any difference."""
import json
import os
import sys

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
A, B = f"{ROOT}/data/p1v12", f"{ROOT}/data/p1v20"


def main():
    rows = [r["row"] for r in json.load(open(f"{B}/replaced.json"))["replaced"]]
    keep = np.setdiff1d(np.arange(3000), rows)
    res, ok = {"n_replaced": len(rows), "n_kept": int(len(keep))}, True
    for f in ["train_codes.npy"] + [f"train_bits_{v}.npy" for v in ("v50uev", "v50f", "v41uev", "v50uevlnc")] + ["train_events.npy"]:
        x, y = np.load(f"{A}/{f}", mmap_mode="r"), np.load(f"{B}/{f}", mmap_mode="r")
        same = x.shape == y.shape and all(np.array_equal(x[i:i + 100][np.isin(np.arange(i, min(i + 100, 3000)), keep)],
                                                         y[i:i + 100][np.isin(np.arange(i, min(i + 100, 3000)), keep)])
                                          for i in range(0, 3000, 100))
        res[f] = bool(same)
        ok = ok and same
    codes = np.load(f"{B}/train_codes.npy", mmap_mode="r")
    new = np.asarray(codes[rows])
    res["new_n_frac_max"] = float((new == 4).mean(1).max())
    np.save(f"{B}/train_new_codes.npy", new)
    res["ok"] = ok
    json.dump(res, open(f"{B}/CHECK.json", "w"), indent=1)
    print(json.dumps(res), flush=True)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
