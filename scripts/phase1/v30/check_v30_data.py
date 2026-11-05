#!/usr/bin/env python3
"""v3.0 data checks (protocol section 18, execution details 1 and 5).

  check_v30_data.py regress   the label / event scripts with their new options left at the defaults must reproduce
                              data/p1v20 bit for bit (codes, four bit variants, events) -> data/p1v30/REGRESS.json
  check_v30_data.py mm39      chr19 of the downloaded mm39 must equal data/mouse/chr19.fa (S5 test sequence)
  check_v30_data.py summary   window sets, N fractions, label repair and U-ev statistics -> data/p1v30/PREP_SUMMARY.json
Exit 1 on any failed check.
"""
import gzip
import hashlib
import json
import os
import sys

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
D = f"{ROOT}/data/p1v30"


def seq_md5(path, name=None):
    h, on = hashlib.md5(), name is None
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if on and name is not None:
                    break
                on = name is None or line[1:].split()[0] == name
                continue
            if on:
                h.update(line.strip().upper().encode())
    return h.hexdigest()


def main():
    what = sys.argv[1]
    if what == "regress":
        res, ok = {}, True
        for f in ["train_codes.npy"] + [f"train_bits_{v}.npy" for v in ("v50uev", "v50f", "v41uev", "v50uevlnc")] + ["train_events.npy"]:
            x, y = np.load(f"{ROOT}/data/p1v20/{f}", mmap_mode="r"), np.load(f"{D}/regress/{f}", mmap_mode="r")
            same = x.shape == y.shape and all(np.array_equal(x[i:i + 100], y[i:i + 100]) for i in range(0, len(x), 100))
            res[f] = bool(same)
            ok = ok and same
        res["ok"] = ok
        json.dump(res, open(f"{D}/REGRESS.json", "w"), indent=1)
    elif what == "mm39":
        a = seq_md5(f"{ROOT}/data/v30/mouse/mm39.fa.gz", "chr19")
        b = seq_md5(f"{ROOT}/data/mouse/chr19.fa")
        res = {"mm39_chr19_md5": a, "s5_chr19_md5": b, "ok": a == b}
        ok = res["ok"]
        json.dump(res, open(f"{D}/MM39_CHECK.json", "w"), indent=1)
    elif what == "summary":
        res = {}
        for part in ("human_new", "mouse"):
            codes = np.load(f"{D}/{part}/train_codes.npy", mmap_mode="r")
            res[part] = {"n": int(len(codes)), "n_frac_max": float(max((codes[i:i + 100] == 4).mean(1).max() for i in range(0, len(codes), 100))),
                         "windows": json.load(open(f"{D}/{part.split('_')[0]}_windows_stats.json")),
                         "repair": json.load(open(f"{D}/{part}/repair_stats.json")).get("train_v50uev"),
                         "events": {k: v for k, v in json.load(open(f"{D}/{part}/events_stats.json")).get("train", {}).items()
                                    if k in ("donor+", "acceptor+", "TSS+")}}
            ub = json.load(open(f"{D}/{part}/repair_stats.json")).get("unknown_bp", {})
            res[part]["unknown_bp_total"] = int(sum(ub.values()))
        lift = f"{ROOT}/data/v30/junctions_mm39/LIFT_STATS.json"
        if os.path.exists(lift):
            st = json.load(open(lift))
            res["mouse_uev_lift"] = {"rows": sum(v["rows"] for v in st.values()), "kept": sum(v["kept"] for v in st.values()),
                                     "kept_frac": sum(v["kept"] for v in st.values()) / max(sum(v["rows"] for v in st.values()), 1)}
        ok = all(res[p]["n"] == 3000 and res[p]["n_frac_max"] <= 0.1 for p in ("human_new", "mouse"))
        res["ok"] = ok
        json.dump(res, open(f"{D}/PREP_SUMMARY.json", "w"), indent=1)
    else:
        raise SystemExit("usage: check_v30_data.py regress|mm39|summary")
    print(json.dumps(res)[:2000], flush=True)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
