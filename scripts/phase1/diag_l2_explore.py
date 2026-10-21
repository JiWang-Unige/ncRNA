#!/usr/bin/env python3
"""v1.5 diagnostic, EXPLORATORY follow-up (not part of the pre-registered A/B/C decision).

usage: diag_l2_explore.py <truth.gtf> <chrom> <pred.gff3> <out.json>
For L1-hit predictions (pairing as in diag_l2.py): how many truth genes each prediction spans (merging),
junction deviation for predicted introns inside / outside the paired truth span, exact recovery of truth
introns (recall side), length ratios; and how many predictions hit each truth gene (fragmentation).
"""
import json
import os
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ncrna_annot.eval import ruler  # noqa: E402
from diag_l2 import BINS, introns, ov  # noqa: E402


def span(ex):
    return ex[0][0], ex[-1][1]


def binof(d):
    return next(n for n, lo, hi in BINS if lo <= d <= hi)


def analyse(truth, preds, cls):
    genes = {g: v for g, v in truth.items() if v["cls"] == cls}
    idx, allidx = ruler.ExonIndex(genes), ruler.ExonIndex(truth)
    tspan = {g: (min(e[0] for ex in v["tx"].values() for e in ex), max(e[1] for ex in v["tx"].values() for e in ex))
             for g, v in truth.items()}
    n_genes_spanned, inside, outside, recov, len_ratio, intr_per_pred = [], Counter(), Counter(), [], [], []
    hits_per_gene = Counter()
    for p in (p for p in preds if p["cls"] == cls):
        hits = idx.hits(p["strand"], p["exons"])
        if not hits:
            continue
        for g in hits:
            hits_per_gene[g] += 1
        best, bo = None, -1
        for gid in hits:
            for ex in genes[gid]["tx"].values():
                ex = sorted(ex)
                o = ov(p["exons"], ex)
                if o > bo:
                    best, bo = ex, o
        q = p["exons"]
        a, b = span(q)
        # truth genes (any class, same strand) whose span overlaps the prediction span
        n_genes_spanned.append(sum(1 for g, v in truth.items() if v["strand"] == p["strand"]
                                   and tspan[g][0] <= b and tspan[g][1] >= a))
        ta, tb = span(best)
        ti, qi = introns(best), introns(q)
        intr_per_pred.append(len(qi))
        for x, y in qi:
            d = min(max(abs(x - u), abs(y - v)) for u, v in ti) if ti else None
            tgt = inside if (x >= ta and y <= tb) else outside
            tgt[binof(d) if d is not None else "no_truth_intron"] += 1
        if ti:
            recov.append(sum(1 for t in ti if t in set(qi)) / len(ti))
        len_ratio.append((b - a + 1) / (tb - ta + 1))
    pct = lambda v, q: float(np.percentile(v, q)) if v else None
    frag = Counter(min(c, 5) for c in hits_per_gene.values())
    norm = lambda c: {k: {"n": v, "frac": v / max(sum(c.values()), 1)} for k, v in sorted(c.items())}
    return {
        "n_paired": len(len_ratio),
        "truth_genes_spanned_per_pred": {"median": pct(n_genes_spanned, 50), "frac_ge2": float(np.mean([n >= 2 for n in n_genes_spanned])) if n_genes_spanned else None},
        "pred_introns_inside_truth_span": norm(inside),
        "pred_introns_outside_truth_span": norm(outside),
        "truth_intron_exact_recovery": {"median": pct(recov, 50), "mean": float(np.mean(recov)) if recov else None},
        "pred_span_over_truth_span": {"median": pct(len_ratio, 50), "q25": pct(len_ratio, 25), "q75": pct(len_ratio, 75)},
        "pred_introns_per_pred_median": pct(intr_per_pred, 50),
        "preds_per_hit_truth_gene(5=5+)": {str(k): v for k, v in sorted(frag.items())},
    }


def main():
    tgtf, chrom, gff, out = sys.argv[1:5]
    truth, preds = ruler.read_gtf(tgtf, chrom), ruler.read_pred_gff3(gff, chrom)
    res = {"exploratory": True, "lnc": analyse(truth, preds, "lnc"), "mrna": analyse(truth, preds, "mrna")}
    json.dump(res, open(out, "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
