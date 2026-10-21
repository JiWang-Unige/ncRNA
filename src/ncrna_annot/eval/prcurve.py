"""Gene-level PR curve for lncRNA predictions with locus-merged precision (protocol v1.4 section 11).

preds: list of dict(strand, exons [(a, b) 1-based], conf); truth: ruler.read_gtf genes.
Thresholds = every distinct confidence (descending).  At each:
  Sn = fraction of truth lnc genes hit (same-strand exon overlap) by kept predictions
  Pr = kept predictions merged into same-strand overlapping loci; fraction of loci hitting a truth lnc gene
AUPRC = sum_i (Sn_i - Sn_{i-1}) * Pr_i  (step-wise, Sn_0 = 0).
"""
import numpy as np

from . import ruler


def _merge(preds):
    loci = []
    for s in ("+", "-", "."):
        items = sorted((p["exons"][0][0], p["exons"][-1][1], p) for p in preds if p["strand"] == s)
        cur = None
        for a, b, p in items:
            if cur and a <= cur["end"]:
                cur["end"] = max(cur["end"], b)
                cur["members"].append(p)
            else:
                if cur:
                    loci.append(cur)
                cur = {"end": b, "members": [p]}
        if cur:
            loci.append(cur)
    return loci


def pr_curve(preds, truth, max_points=200):
    lnc = {g: v for g, v in truth.items() if v["cls"] == "lnc"}
    idx = ruler.ExonIndex(lnc)
    for p in preds:
        p["_hits"] = idx.hits(p["strand"], p["exons"])
    order = sorted(preds, key=lambda p: -p["conf"])
    confs = sorted({p["conf"] for p in order}, reverse=True)
    if len(confs) > max_points:
        confs = list(np.quantile(confs, np.linspace(1, 0, max_points)))
    pts, prev_sn, auprc = [], 0.0, 0.0
    for c in confs:
        kept = [p for p in order if p["conf"] >= c]
        hit = set().union(*(p["_hits"] for p in kept)) if kept else set()
        sn = len(hit) / max(len(lnc), 1)
        loci = _merge(kept)
        pr = sum(any(m["_hits"] for m in l["members"]) for l in loci) / max(len(loci), 1)
        auprc += max(sn - prev_sn, 0.0) * pr
        prev_sn = max(prev_sn, sn)
        pts.append({"conf": float(c), "sn": sn, "pr": pr, "n_kept": len(kept), "n_loci": len(loci)})
    best = max(pts, key=lambda q: 2 * q["sn"] * q["pr"] / (q["sn"] + q["pr"]) if q["sn"] + q["pr"] else 0) if pts else None
    f1 = 2 * best["sn"] * best["pr"] / (best["sn"] + best["pr"]) if best and best["sn"] + best["pr"] else 0.0
    for p in preds:
        p.pop("_hits", None)
    return {"auprc": auprc, "best_f1": f1, "best_point": best, "n_truth": len(lnc), "n_pred": len(preds), "curve": pts}
