#!/usr/bin/env python3
"""S3 post hoc checks (written after S3_chr18 / S3_chr20 were read; descriptive, labelled post hoc in RESULT_S3.md).
1. u1_pas = u1_kb - pas_kb made the frozen R1 design collinear: refit without u1_pas.
2. Promoter CpG may stand for lncRNAs that share a promoter with a protein-coding gene.  Add
   divergent  = a v50 mRNA gene on the opposite strand has a TSS within 1 kb of the lncRNA's outermost TSS (head to head)
   antisense  = the lncRNA body overlaps an opposite-strand v50 mRNA gene body
   log_d_mrna = log10(1 + distance from the lncRNA TSS to the nearest v50 mRNA TSS on either strand)
and refit; detection rate by divergent / antisense / neither.

  s3_posthoc.py --chrom chr18 --out outputs/PHASE1/S3/S3_posthoc_chr18.json
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s3_features as F  # noqa: E402
from s3_features import ROOT, RUNS, ruler, P  # noqa: E402


def tss_of(g):
    txs = [ruler._merge(sorted(ex)) for ex in g["tx"].values()]
    return [t[0][0] - 1 if g["strand"] == "+" else t[-1][1] - 1 for t in txs]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True, choices=("chr18", "chr20"))
    ap.add_argument("--config", default="configs/p15_v19b.json")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    cfg = json.load(open(f"{ROOT}/{a.config}"))
    tg = cfg["targets"][a.chrom]
    base = json.load(open(f"{ROOT}/outputs/PHASE1/S3/S3_{a.chrom}.json"))
    truth = ruler.read_gtf(f"{ROOT}/{tg['truth']}", tg["truth_chrom"])
    v41 = ruler.read_gtf(f"{ROOT}/data/human/gencode.v41.annotation.gtf.gz", a.chrom)
    ids41 = {gid.split(".")[0] for gid in v41}
    seq = F.read_fasta(f"{ROOT}/data/human/{a.chrom}.fa")
    cage = F.load_peaks(f"{ROOT}/data/bioval/raw/fantom5_cage_peaks.bed.gz", a.chrom, 5, 6, 7)
    polya = F.load_peaks(f"{ROOT}/data/bioval/raw/polyasite2_GRCh38.96.tsv.gz", a.chrom, 5, 1, 2, skip_header=True,
                         chrom_name=a.chrom.replace("chr", ""))
    lnc = {gid: g for gid, g in truth.items() if g["cls"] == "lnc"}
    mrna = [g for g in truth.values() if g["cls"] == "mrna"]
    gids = sorted(lnc)
    m_tss = {s: np.sort(np.array([x for g in mrna if g["strand"] == s for x in tss_of(g)], np.int64)) for s in "+-"}
    all_m_tss = np.sort(np.concatenate([m_tss["+"], m_tss["-"]]))
    m_body = {s: [(min(a_ for ex in g["tx"].values() for a_, _ in ex), max(b for ex in g["tx"].values() for _, b in ex))
                  for g in mrna if g["strand"] == s] for s in "+-"}
    rows, extra = [], []
    for gid in gids:
        g = lnc[gid]
        r = F.gene_features(g, seq, cage, polya)
        r["era"] = float(gid.split(".")[0] in ids41)
        rows.append(r)
        opp = "-" if g["strand"] == "+" else "+"
        t = min(tss_of(g)) if g["strand"] == "+" else max(tss_of(g))
        div = F.near(m_tss[opp], t, 1000)
        s0 = min(a_ for ex in g["tx"].values() for a_, _ in ex)
        e0 = max(b for ex in g["tx"].values() for _, b in ex)
        anti = any(a_ <= e0 and s0 <= b for a_, b in m_body[opp])
        i = np.searchsorted(all_m_tss, t)
        d = min([abs(all_m_tss[j] - t) for j in (i - 1, i) if 0 <= j < len(all_m_tss)] or [10 ** 7])
        extra.append({"divergent": float(div), "antisense": float(anti), "log_d_mrna": float(np.log10(1 + d))})
    idx = ruler.ExonIndex(lnc)
    Y = {}
    for run in RUNS:
        genes = P.load_genes(base["pred_files"][run])
        hit = set()
        for p in P.final_preds(genes, base["tau"]):
            if p["cls"] == "lnc":
                hit |= set(idx.hits(p["strand"], p["exons"]))
        Y[run] = np.array([gid in hit for gid in gids], float)
    assert all(abs(Y[r].mean() - base["detected_frac"][r]) < 1e-12 for r in RUNS)   # same labels as the frozen run
    rng = np.random.default_rng(0)
    boots = [rng.integers(0, len(gids), len(gids)) for _ in range(F.N_BOOT)]

    def fit(names, table):
        X = np.array([[t[k] for k in names] for t in table], float)
        sd = X.std(0)
        sd[sd == 0] = 1
        Z = (X - X.mean(0)) / sd
        coefs = {r: F.fit_logit(Z, Y[r])[1:] for r in RUNS}
        bm = np.array([np.mean([F.fit_logit(Z[b], Y[r][b])[1:] for r in RUNS], 0) for b in boots])
        lo, hi = np.percentile(bm, [2.5, 97.5], 0)
        return {k: {"per_run": [float(coefs[r][i]) for r in RUNS], "mean": float(np.mean([coefs[r][i] for r in RUNS])),
                    "ci95": [float(lo[i]), float(hi[i])],
                    "associated": bool((lo[i] > 0 or hi[i] < 0) and len({np.sign(coefs[r][i]) for r in RUNS}) == 1)}
                for i, k in enumerate(names)}
    no_upas = [k for k in F.FEATS if k != "u1_pas"]
    table = [{**r, **e} for r, e in zip(rows, extra)]
    ymean = np.mean([Y[r] for r in RUNS], 0)
    cat = np.array(["divergent" if e["divergent"] else "antisense" if e["antisense"] else "neither" for e in extra])
    res = {"chrom": a.chrom, "post_hoc": True,
           "R1_without_u1_pas": fit(no_upas, table),
           "R1_with_promoter_sharing": fit(no_upas + ["divergent", "antisense", "log_d_mrna"], table),
           "detected_by_category": {c: {"n": int((cat == c).sum()), "detected": float(ymean[cat == c].mean()) if (cat == c).any() else None,
                                        "cpg_oe_mean": float(np.mean([table[i]["cpg_oe"] for i in np.flatnonzero(cat == c)])) if (cat == c).any() else None}
                                    for c in ("divergent", "antisense", "neither")}}
    json.dump(res, open(f"{ROOT}/{a.out}", "w"), indent=1)
    show = lambda d: {k: (round(v["mean"], 3), [round(x, 3) for x in v["ci95"]], v["associated"]) for k, v in d.items()}
    print(json.dumps({"chrom": a.chrom, "without_u1_pas": show(res["R1_without_u1_pas"]),
                      "with_sharing": show(res["R1_with_promoter_sharing"]), "by_category": res["detected_by_category"]}), flush=True)


if __name__ == "__main__":
    main()
