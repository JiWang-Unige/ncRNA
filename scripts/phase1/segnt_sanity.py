#!/usr/bin/env python3
"""POST HOC sanity check (RESULT_BASELINES.md): per-base agreement of the SegmentNT chr20 probabilities with v50 masks,
to rule out a channel / orientation mix-up in segnt_chr20.py.  AUROC on every 37th base (N excluded) of
  protein_coding_gene vs v50 mRNA-class gene bodies, lncRNA vs v50 lnc gene bodies, exon vs exons of any v50 gene;
plus the share of 1-bp runs among protein_coding_gene >= 0.5 runs and the mean |p - 0.5| at their bases."""
import json, os, sys
import numpy as np
ROOT = os.path.expanduser("~/ncrna-annot")
sys.path.insert(0, f"{ROOT}/src")
from ncrna_annot.eval import ruler
D = f"{ROOT}/outputs/PHASE1/BASELINES/segnt"
probs = np.load(f"{D}/probs.npy", mmap_mode="r")
codes = np.load(f"{ROOT}/data/p1/chr20_codes.npy", mmap_mode="r")
L = len(codes)
truth = ruler.read_gtf(f"{ROOT}/data/human/gencode.v50.annotation.gtf.gz", "chr20")
body = {"mrna": np.zeros(L, bool), "lnc": np.zeros(L, bool)}
exon = np.zeros(L, bool)
for g in truth.values():
    for ex in g["tx"].values():
        for s, e in ex:
            exon[s - 1:e] = True
        if g["cls"] in body:
            lo, hi = min(s for s, _ in ex), max(e for _, e in ex)
            body[g["cls"]][lo - 1:hi] = True
idx = np.arange(0, L, 37)
idx = idx[np.asarray(codes[idx]) != 4]
def auroc(score, lab):
    o = np.argsort(score, kind="stable"); r = np.empty(len(o)); r[o] = np.arange(1, len(o) + 1)
    n1 = lab.sum(); n0 = len(lab) - n1
    return float((r[lab].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))
P = np.asarray(probs[idx], np.float32)
out = {"n_bases": int(len(idx)),
       "truth_frac": {"mrna_body": float(body["mrna"][idx].mean()), "lnc_body": float(body["lnc"][idx].mean()), "exon_any": float(exon[idx].mean())},
       "auroc": {"protein_coding_gene_vs_mrna_body": auroc(P[:, 0], body["mrna"][idx]),
                 "lncRNA_vs_lnc_body": auroc(P[:, 1], body["lnc"][idx]),
                 "exon_vs_exon_any": auroc(P[:, 2], exon[idx])}}
pc = np.asarray(probs[:, 0], np.float32)
m = pc >= 0.5
d = np.diff(np.concatenate([[0], m.astype(np.int8), [0]]))
st, en = np.flatnonzero(d == 1), np.flatnonzero(d == -1)
ln = en - st
out["pcg_runs"] = {"n": int(len(ln)), "n_len1": int((ln == 1).sum()), "n_ge_1000": int((ln >= 1000).sum()),
                   "bp_in_runs_ge_1000_frac": float(ln[ln >= 1000].sum() / ln.sum())}
one = st[ln == 1]
out["pcg_runs"]["mean_abs_p_minus_half_at_len1"] = float(np.abs(pc[one] - 0.5).mean()) if len(one) else None
print(json.dumps(out, indent=1))
json.dump(out, open(f"{D}/SANITY.json", "w"), indent=1)
