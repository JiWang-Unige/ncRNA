#!/usr/bin/env python3
"""M1 weeks 2-4, chr18 development: score the reference-free long-read baselines (StringTie -L, IsoQuant de novo on PacBio
rep1) against the WTC11 truth (rep2 & rep3).  Loci = the tool's genes (transcript spans merged per gene id); lnc
candidates = stranded loci not overlapping a Tiberius coding gene on the same strand (as the pairing baseline).
Also stratified by truth expression (rep2 + rep3 reads: <= median vs > median).
  dev_longread_score.py [--chrom chr18] -> outputs/M1/W24/dev/LONGREAD_<chrom>.json
"""
import argparse
import json
import re
from collections import defaultdict

import numpy as np

import w24_lib as L
from dev_baselines import tib_hit


def loci_of(gtf, chrom):
    g = defaultdict(lambda: [None, 10 ** 12, 0])
    for line in open(gtf):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[0] != chrom or f[2] != "transcript" or f[6] not in "+-":
            continue
        gid = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
        x = g[(gid, f[6])]
        x[0], x[1], x[2] = f[6], min(x[1], int(f[3])), max(x[2], int(f[4]))
    return [{"strand": st, "exons": [(s, e)]} for st, s, e in g.values() if e - s + 1 >= 200]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", default="chr18")
    a = ap.parse_args()
    assert a.chrom == "chr18"
    import simple_det as D
    truth = L.load_truth(a.chrom)
    tib = D.read_tiberius(a.chrom)
    O = f"{L.W}/longread/{a.chrom}"
    expr = {i: x["count2"] + x["count3"] for i, x in truth["lnc"].items()}
    med = float(np.median(list(expr.values())))
    res = {"chrom": a.chrom, "expr_median_reads": med}
    for tool, gtf in (("stringtie", f"{O}/stringtie.gtf"), ("isoquant", f"{O}/isoquant/rep1/rep1.transcript_models.gtf")):
        loci = loci_of(gtf, a.chrom)
        lnc = [p for p in loci if not tib_hit(tib, p["strand"], p["exons"][0][0], p["exons"][-1][1])]
        r = L.score(truth, lnc)
        hit = set(r.pop("S50_ids"))
        lo = [i for i, v in expr.items() if v <= med]
        hi = [i for i, v in expr.items() if v > med]
        r["S50_sn_low_expr"] = sum(i in hit for i in lo) / len(lo)
        r["S50_sn_high_expr"] = sum(i in hit for i in hi) / len(hi)
        r["n_loci_all"] = len(loci)
        res[tool] = r
    json.dump(res, open(f"{L.W}/dev/LONGREAD_{a.chrom}.json", "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
