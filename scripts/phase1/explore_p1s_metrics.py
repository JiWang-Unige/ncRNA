#!/usr/bin/env python3
"""EXPLORATORY (post hoc, not a P1-S decision input): fragment-robust view of P1-S val GFFs.

For a GFF on chr18 val windows:
  - raw L1 (as in the pre-registered screen)
  - locus-merged L1: predicted lnc transcripts overlapping on the same strand are merged into one predicted
    locus before computing precision (fragments inside one true gene count once)
  - length filter (lncRNA definition): drop predicted lnc transcripts with mature length < 200 nt
usage: explore_p1s_metrics.py <truth_gtf> <gff> [<gff> ...]
"""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ncrna_annot.eval import ruler  # noqa: E402


def merged_loci(preds):
    by = {}
    for p in preds:
        by.setdefault(p["strand"], []).append((p["exons"][0][0], p["exons"][-1][1], p))
    loci = []
    for s, items in by.items():
        items.sort(key=lambda x: x[0])
        cur = None
        for a, b, p in items:
            if cur and a <= cur["end"]:
                cur["end"] = max(cur["end"], b)
                cur["exons"] += p["exons"]
            else:
                if cur:
                    loci.append(cur)
                cur = {"strand": s, "end": b, "exons": list(p["exons"])}
        if cur:
            loci.append(cur)
    return loci


def evaluate(truth, preds):
    lnc_truth = {g: v for g, v in truth.items() if v["cls"] == "lnc"}
    idx_all, idx_lnc = ruler.ExonIndex(truth), ruler.ExonIndex(lnc_truth)
    lp = [p for p in preds if p["cls"] == "lnc"]
    hit_genes = set()
    for p in lp:
        hit_genes |= idx_lnc.hits(p["strand"], p["exons"])
    sn = len(hit_genes) / max(len(lnc_truth), 1)
    raw_pr = sum(bool(idx_lnc.hits(p["strand"], p["exons"])) for p in lp) / max(len(lp), 1)
    loci = merged_loci(lp)
    loc_pr = sum(bool(idx_lnc.hits(l["strand"], l["exons"])) for l in loci) / max(len(loci), 1)
    f = lambda a, b: 2 * a * b / (a + b) if a + b else 0.0
    return {"sn": round(sn, 3), "raw_pr": round(raw_pr, 3), "raw_f1": round(f(sn, raw_pr), 3),
            "n_pred_tx": len(lp), "n_loci": len(loci), "locus_pr": round(loc_pr, 3), "locus_f1": round(f(sn, loc_pr), 3)}


def main():
    truth = ruler.read_gtf(sys.argv[1], "chr18")
    out = {}
    for gff in sys.argv[2:]:
        preds = ruler.read_pred_gff3(gff, "chr18")
        keep = [p for p in preds if p["cls"] != "lnc" or sum(b - a + 1 for a, b in p["exons"]) >= 200]
        out[os.path.basename(os.path.dirname(gff)) + "/" + os.path.basename(gff)] = {
            "all": evaluate(truth, preds), "len>=200": evaluate(truth, keep)}
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
