#!/usr/bin/env python3
"""v1.5 diagnostic step 2 (protocol section 12): why are intron chains wrong?

usage: diag_l2.py <truth.gtf(.gz)> <chrom> <pred.gff3> <out.json>
Per class (lnc, mrna), for predictions that hit a same-class same-strand truth gene (L1), each prediction is
paired with the truth isoform of largest exon overlap.  Reports items 1-6 and the A/B/C decision for lnc.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ncrna_annot.eval import ruler  # noqa: E402


def ov(a, b):
    return sum(max(0, min(y, d) - max(x, c) + 1) for x, y in a for c, d in b)


def introns(ex):
    return [(ex[i][1] + 1, ex[i + 1][0] - 1) for i in range(len(ex) - 1)]


BINS = [("exact", 0, 0), ("1-5", 1, 5), ("6-20", 6, 20), ("21-100", 21, 100), (">100", 101, 10 ** 12)]


def analyse(truth, preds, cls):
    genes = {g: v for g, v in truth.items() if v["cls"] == cls}
    idx = ruler.ExonIndex(genes)
    cp = [p for p in preds if p["cls"] == cls]
    single_vs_multi = paired = paired_multi = 0
    jbins = {b[0]: 0 for b in BINS}
    jbins["no_truth_intron"] = 0
    cover, extra, tss, tes, spancov, n_pred_introns, truth_intr_len, pred_intr_len = [], [], [], [], [], 0, [], []
    for p in cp:
        hits = idx.hits(p["strand"], p["exons"])
        if not hits:
            continue
        best, bo = None, -1
        for gid in hits:
            for ex in genes[gid]["tx"].values():
                ex = sorted(ex)
                o = ov(p["exons"], ex)
                if o > bo:
                    best, bo = ex, o
        paired += 1
        t, q = best, p["exons"]
        if len(t) > 1:
            paired_multi += 1
            if len(q) == 1:
                single_vs_multi += 1
        ti = introns(t)
        truth_intr_len += [b - a + 1 for a, b in ti]
        for a, b in introns(q):
            n_pred_introns += 1
            pred_intr_len.append(b - a + 1)
            if not ti:
                jbins["no_truth_intron"] += 1
                continue
            dev = min(max(abs(a - x), abs(b - y)) for x, y in ti)
            for name, lo, hi in BINS:
                if lo <= dev <= hi:
                    jbins[name] += 1
                    break
        cover.append(sum(1 for e in t if ov([e], q) > 0) / len(t))
        extra.append(sum(1 for e in q if ov([e], t) == 0))
        plus = p["strand"] == "+"
        tss.append(abs((q[0][0] - t[0][0]) if plus else (q[-1][1] - t[-1][1])))
        tes.append(abs((q[-1][1] - t[-1][1]) if plus else (q[0][0] - t[0][0])))
        spancov.append(max(0, min(q[-1][1], t[-1][1]) - max(q[0][0], t[0][0]) + 1) / (t[-1][1] - t[0][0] + 1))
    # composition of truth
    n_multi = sum(1 for v in genes.values() if any(len(ex) > 1 for ex in v["tx"].values()))
    pct = lambda v, q: float(np.percentile(v, q)) if v else None
    res = {
        "n_pred": len(cp), "n_paired": paired,
        "1_single_exon_pred_for_multi_exon_truth": {"n": single_vs_multi, "n_paired_multi_exon_truth": paired_multi,
                                                    "frac": single_vs_multi / max(paired_multi, 1)},
        "2_junction_deviation": {"n_pred_introns": n_pred_introns,
                                 "bins": {k: {"n": v, "frac": v / max(n_pred_introns, 1)} for k, v in jbins.items()}},
        "3_exon_coverage": {"median": pct(cover, 50), "q25": pct(cover, 25), "q75": pct(cover, 75),
                            "extra_pred_exons_mean": float(np.mean(extra)) if extra else None},
        "4_ends": {"tss_dist_median": pct(tss, 50), "tes_dist_median": pct(tes, 50),
                   "span_coverage_median": pct(spancov, 50)},
        "6_truth_composition": {"n_genes": len(genes), "frac_multi_exon": n_multi / max(len(genes), 1),
                                "truth_intron_len_median": pct(truth_intr_len, 50), "pred_intron_len_median": pct(pred_intr_len, 50)},
    }
    return res


def strand_type_errors(truth, preds):
    idx = ruler.ExonIndex(truth)
    flip = {"+": "-", "-": "+"}
    out = {}
    for cls in ("lnc", "mrna"):
        cp = [p for p in preds if p["cls"] == cls]
        same = wrong_type = wrong_strand = 0
        for p in cp:
            h = idx.hits(p["strand"], p["exons"])
            if any(truth[g]["cls"] == cls for g in h):
                same += 1
            elif h:
                wrong_type += 1
            elif idx.hits(flip[p["strand"]], p["exons"]):
                wrong_strand += 1
        out[cls] = {"n": len(cp), "correct_class_and_strand": same, "other_class_same_strand": wrong_type,
                    "opposite_strand_only": wrong_strand}
    return out


def decide(r):
    j = r["2_junction_deviation"]["bins"]
    near = j["1-5"]["frac"] + j["6-20"]["frac"]
    A = near >= 0.5
    B = r["1_single_exon_pred_for_multi_exon_truth"]["frac"] >= 0.5 or (r["3_exon_coverage"]["median"] or 0) < 0.5
    C = (r["4_ends"]["span_coverage_median"] or 0) < 0.5 or max(r["4_ends"]["tss_dist_median"] or 0, r["4_ends"]["tes_dist_median"] or 0) > 1000
    order = [k for k, v in (("B", B), ("C", C), ("A", A)) if v]
    return {"A_boundary_precision": A, "near_miss_frac": near, "B_structure": B, "C_ends": C,
            "priority": order, "G2_indicated": B or C}


def main():
    tgtf, chrom, gff, out = sys.argv[1:5]
    truth = ruler.read_gtf(tgtf, chrom)
    preds = ruler.read_pred_gff3(gff, chrom)
    res = {"lnc": analyse(truth, preds, "lnc"), "mrna": analyse(truth, preds, "mrna"),
           "5_strand_type": strand_type_errors(truth, preds)}
    res["decision_lnc"] = decide(res["lnc"])
    json.dump(res, open(out, "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
