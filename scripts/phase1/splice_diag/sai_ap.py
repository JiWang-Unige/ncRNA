#!/usr/bin/env python3
"""POST HOC (2026-11-03 evening, RESULT_SPLICE_OWN.md section 2): chr18 lnc / mRNA site AP of SpliceAI and of the selected
own net, both with own_splice.site_ap.  -> outputs/PHASE1/SPLICE_OWN/AP_COMPARE.json"""
import json, os, sys
import numpy as np
ROOT = os.path.expanduser("~/ncrna-annot")
sys.path.insert(0, f"{ROOT}/src"); sys.path.insert(0, f"{ROOT}/scripts/phase1"); sys.path.insert(0, f"{ROOT}/scripts/phase1/splice_diag")
import own_splice, sai_score
from ncrna_annot.eval import ruler
codes = np.asarray(np.load(f"{ROOT}/data/p1/chr18_codes.npy", mmap_mode="r"))
truth = ruler.read_gtf(f"{ROOT}/outputs/PHASE1/DIAG_L2/chr18_v50.gtf", "chr18")
ts = {c: sai_score.truth_sites(truth, c) for c in ("mrna", "lnc")}
out = {}
for name, f in (("spliceai", "outputs/PHASE1/SPLICE_DIAG/sai_chr18.npy"), ("own_s10k", "outputs/PHASE1/SPLICE_OWN/own_chr18.npy")):
    sc = np.asarray(np.load(f"{ROOT}/{f}", mmap_mode="r"), np.float32)
    out[name] = own_splice.site_ap(sc, codes, lambda cls, col: ts[cls][col])
    print(name, round(out[name]["lnc_ap"], 4), round(out[name]["mrna_ap"], 4), flush=True)
json.dump(out, open(f"{ROOT}/outputs/PHASE1/SPLICE_OWN/AP_COMPARE.json", "w"), indent=1)
