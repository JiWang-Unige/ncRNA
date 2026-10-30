"""v1.9 development (chr18 only): §3 guardrails at one exon-level MEA operating point -- mRNA L1 Sn and the number of
truth lnc genes hit only by predicted mRNAs.  guard.py root name post gamma tau out.json"""
import json
import sys

import numpy as np

root, name, path, gamma, tau, out = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), float(sys.argv[5]), sys.argv[6]
sys.path.insert(0, root + "/src")
from ncrna_annot.decode import chrom
from ncrna_annot.eval import ruler
from ncrna_annot.grammar import g1

post = np.load(path, mmap_mode="r")
codes = np.load(f"{root}/data/p1/chr18_codes.npy", mmap_mode="r")[:len(post)]
T = g1.G1Transitions()().detach().double().numpy()
truth = ruler.read_gtf(f"{root}/outputs/PHASE1/DIAG_L2/chr18_v50.gtf", "chr18")
genes = chrom.decode_blocks(post, codes, T, 20, 50, mode="mea", gamma=gamma, level="exon", log_gain=False)
genes = [g for g in genes if g["type"] != "l" or sum(b - a + 1 for a, b in g["exons"]) >= 200]
for g in genes:
    g["conf"] = chrom.gene_confidence(g, post)
genes = [g for g in genes if g["type"] != "l" or g["conf"] >= tau]
preds = [{"id": str(i), "strand": g["strand"], "cls": "lnc" if g["type"] == "l" else "mrna",
          "exons": [tuple(x) for x in g["exons"]]} for i, g in enumerate(genes)]
idx = ruler.ExonIndex(truth)
sm = ruler.score_class("mrna", truth, preds, idx, 0)
sl = ruler.score_class("lnc", truth, preds, idx, 0)
lnc_ids = {k for k, v in truth.items() if v["cls"] == "lnc"}
mrna_hit_lnc = set()  # truth lnc genes overlapped (exon, same strand) by a predicted mRNA that hits no truth mRNA
for p in preds:
    if p["cls"] != "mrna":
        continue
    h = idx.hits(p["strand"], p["exons"])
    if h and not any(truth[x]["cls"] == "mrna" for x in h):
        mrna_hit_lnc |= h & lnc_ids
res = {"name": name, "gamma": gamma, "tau": tau, "lnc_L1": sl["L1_locus"], "mrna_L1": sm["L1_locus"],
       "n_lnc_genes_reported_as_mrna": len(mrna_hit_lnc), "n_pred_lnc": sum(p["cls"] == "lnc" for p in preds),
       "n_pred_mrna": sum(p["cls"] == "mrna" for p in preds)}
print(json.dumps({k: v for k, v in res.items() if k not in ("lnc_L1", "mrna_L1")} | {"lnc": [sl["L1_locus"]["sn"], sl["L1_locus"]["pr_strict"]], "mrna": [sm["L1_locus"]["sn"], sm["L1_locus"]["pr_strict"]], "mrna_other": sm["L1_locus"]["pred_hits_other_class_only"]}), flush=True)
json.dump(res, open(out, "w"), indent=1)
