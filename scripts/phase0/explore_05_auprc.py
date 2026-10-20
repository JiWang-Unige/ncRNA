#!/usr/bin/env python3
"""0.5 EXPLORATORY (post hoc, not part of the pre-registered gate): threshold-free ranking quality.

Per run, on chr20 tokens (6-bp, centre base), using saved p_lnc_exon_chr20.npy:
  A) delta-vs-intergenic: positives = exon tokens of delta genes (v50-new lncRNA loci),
     negatives = tokens outside every v50 gene (any class, any strand)
  B) all-v50-lnc-exon vs rest (whole chromosome)
Reports average precision; arm means over seeds and paired seed-mean differences vs F.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ruler  # noqa: E402

ROOT = os.path.expanduser("~/ncrna-annot")
OUT = f"{ROOT}/outputs/PHASE0/0.5_unknown"
new = ruler.read_gtf(f"{ROOT}/data/human/gencode.v50.annotation.gtf.gz", "chr20")
old = ruler.read_gtf(f"{ROOT}/data/human/gencode.v41.annotation.gtf.gz", "chr20")
delta = ruler.delta_genes(new, old)
ntok = len(np.load(f"{OUT}/F_s42/p_lnc_exon_chr20.npy", mmap_mode="r"))
centre = np.arange(ntok) * 6 + 3 + 1  # 1-based centre base


def mask(genes, exons_only=True):
    m = np.zeros(ntok, bool)
    for g in genes.values():
        for ex in g["tx"].values():
            ivs = ex if exons_only else [(min(s for s, _ in ex), max(e for _, e in ex))]
            for s, e in ivs:
                m[(s - 4) // 6 + 1 if s > 4 else 0:(e - 4) // 6 + 1] = True
    return m


def ap(score, y):
    o = np.argsort(-score, kind="stable")
    y = y[o]
    tp = np.cumsum(y)
    return float((tp / np.arange(1, len(y) + 1))[y].sum() / max(y.sum(), 1))


pos_delta = mask({g: new[g] for g in delta})
in_any_gene = mask(new, exons_only=False)
lnc_ex = mask({g: v for g, v in new.items() if v["cls"] == "lnc"})
selA = pos_delta | ~in_any_gene
res = {}
for arm in ("F", "Uev", "Upu"):
    for s in (42, 123, 456):
        p = np.load(f"{OUT}/{arm}_s{s}/p_lnc_exon_chr20.npy").astype(np.float32)
        res[f"{arm}_s{s}"] = {"A_delta_vs_intergenic_ap": ap(p[selA], pos_delta[selA]),
                              "B_v50_lnc_exon_ap": ap(p, lnc_ex)}
out = {"note": "EXPLORATORY, post hoc; not used for GATE_0.5",
       "prevalence": {"A": float(pos_delta[selA].mean()), "B": float(lnc_ex.mean())}, "runs": res, "arms": {}}
for arm in ("F", "Uev", "Upu"):
    for k in ("A_delta_vs_intergenic_ap", "B_v50_lnc_exon_ap"):
        v = [res[f"{arm}_s{s}"][k] for s in (42, 123, 456)]
        out["arms"].setdefault(arm, {})[k] = [round(float(np.mean(v)), 4), round(float(np.std(v)), 4)]
for arm in ("Uev", "Upu"):
    for k in ("A_delta_vs_intergenic_ap", "B_v50_lnc_exon_ap"):
        d = [res[f"{arm}_s{s}"][k] - res[f"F_s{s}"][k] for s in (42, 123, 456)]
        out.setdefault("diff_vs_F", {})[f"{arm}|{k}"] = {"mean": round(float(np.mean(d)), 4),
                                                         "per_seed": [round(x, 4) for x in d]}
json.dump(out, open(f"{OUT}/EXPLORATORY_auprc.json", "w"), indent=1)
print(json.dumps({k: out[k] for k in ("prevalence", "arms", "diff_vs_F")}, indent=1))
