#!/usr/bin/env python3
"""POST HOC, descriptive (RESULT_BASELINES.md): what the pre-registered SegmentNT decoding produced on chr20.
Per class and threshold: spans (genes before strand duplication), span length quantiles, share < 200 / < 1,000 bp,
share multi-exon, share stranded; base-level fraction of chr20 above threshold.  Reads pred_t*_both.gff3 + probs.npy."""
import json, os
from collections import defaultdict
import numpy as np
D = os.path.expanduser("~/ncrna-annot/outputs/PHASE1/BASELINES/segnt")
probs = np.load(f"{D}/probs.npy", mmap_mode="r")
out = {}
for thr in ("0.3", "0.5", "0.7"):
    genes = {}
    cur = None
    for line in open(f"{D}/pred_t{thr}_both.gff3"):
        f = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(f) < 9:
            continue
        if f[2] in ("lnc_RNA", "mRNA"):
            cur = (f[2], int(f[3]), int(f[4]))
            genes.setdefault(cur, {"strands": set(), "n_ex": 0})
            genes[cur]["strands"].add(f[6])
            genes[cur]["n_ex_" + f[6]] = 0
        elif f[2] == "exon":
            genes[cur]["n_ex_" + f[6]] += 1
    res = {}
    for cls in ("lnc_RNA", "mRNA"):
        g = [(b - a + 1, v) for (c, a, b), v in genes.items() if c == cls]
        L = np.array([x for x, _ in g])
        multi = np.array([max(v[k] for k in v if k.startswith("n_ex_")) > 1 for _, v in g])
        stranded = np.array([len(v["strands"]) == 1 for _, v in g])
        res[cls] = {"spans": len(g), "len_q10_50_90": [int(np.quantile(L, q)) for q in (0.1, 0.5, 0.9)] if len(L) else None,
                    "frac_lt_200": float((L < 200).mean()) if len(L) else None,
                    "frac_lt_1000": float((L < 1000).mean()) if len(L) else None,
                    "frac_multi_exon": float(multi.mean()) if len(L) else None,
                    "frac_stranded": float(stranded.mean()) if len(L) else None}
    t = float(thr)
    res["base_frac_ge_thr"] = {n: float((np.asarray(probs[:, i], np.float32) >= t).mean()) for i, n in enumerate(("protein_coding_gene", "lncRNA", "exon"))}
    out[thr] = res
    print(thr, json.dumps(res), flush=True)
json.dump(out, open(f"{D}/DESC.json", "w"), indent=1)
