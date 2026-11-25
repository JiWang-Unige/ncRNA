"""DEC1: Decima lncRNA blind-spot diagnostic (frozen per PROTOCOL_DEC1.md).

Usage: python dec1.py <Supplementary_file_1.h5ad> <out_dir>
"""
import json
import sys

import anndata as ad
import numpy as np
import pandas as pd
import patsy

H5, OUT = sys.argv[1], sys.argv[2]
NB = 2000
rng = np.random.default_rng(0)

a = ad.read_h5ad(H5)
var = a.var.copy()
keep = var.gene_type.isin(["protein_coding", "lncRNA"]).values
a = a[:, keep].copy()
var = var.loc[keep].copy()

X = np.asarray(a.X, dtype=np.float64)
P = np.asarray(a.layers["preds"], dtype=np.float64)


def colpearson(X, P):
    out = np.full(X.shape[1], np.nan)
    for j in range(X.shape[1]):
        x, p = X[:, j], P[:, j]
        m = np.isfinite(x) & np.isfinite(p)
        if m.sum() > 10 and x[m].std() > 0 and p[m].std() > 0:
            out[j] = np.corrcoef(x[m], p[m])[0, 1]
    return out


var["r"] = colpearson(X, P)
var["sd_obs"] = np.nanstd(X, axis=0)
var["log_mean"] = np.log10(var["mean_counts"].astype(float) + 1e-3)
var["log_len"] = np.log10(var["gene_length"].astype(float))
var["is_lnc"] = (var.gene_type == "lncRNA").astype(int)
res = {"check_recomputed_vs_var_pearson_median_absdiff":
       float(np.nanmedian(np.abs(var["r"] - var["pearson"])))}

FORM = "r ~ is_lnc + cr(log_mean, df=4) + cr(sd_obs, df=4) + log_len + frac_nan"


def ols_coef(d):
    y, Xd = patsy.dmatrices(FORM, d, return_type="dataframe")
    b, *_ = np.linalg.lstsq(Xd.values, y.values[:, 0], rcond=None)
    return b[list(Xd.columns).index("is_lnc")]


def analyse(split):
    d = var[(var.dataset == split) & var.r.notna()].copy()
    lnc, pc = d[d.is_lnc == 1], d[d.is_lnc == 0]
    o = {"n_lnc": int(len(lnc)), "n_pc": int(len(pc)),
         "median_r_lnc": float(lnc.r.median()), "median_r_pc": float(pc.r.median()),
         "median_mean_counts_lnc": float(lnc.mean_counts.median()),
         "median_mean_counts_pc": float(pc.mean_counts.median())}
    gap = o["median_r_pc"] - o["median_r_lnc"]
    bl, bp = lnc.r.values, pc.r.values
    gb = [np.median(rng.choice(bp, len(bp))) - np.median(rng.choice(bl, len(bl)))
          for _ in range(NB)]
    o["H1_gap_pc_minus_lnc"] = float(gap)
    o["H1_ci"] = [float(np.percentile(gb, 2.5)), float(np.percentile(gb, 97.5))]
    c = ols_coef(d)
    idx = np.arange(len(d))
    cb = []
    for _ in range(NB):
        cb.append(ols_coef(d.iloc[rng.choice(idx, len(idx))].reset_index(drop=True)))
    o["H2_coef_is_lnc"] = float(c)
    o["H2_ci"] = [float(np.percentile(cb, 2.5)), float(np.percentile(cb, 97.5))]
    # decile-matched
    d["dec"] = pd.qcut(d.log_mean, 10, labels=False, duplicates="drop")
    rows = []
    for k, g in d.groupby("dec"):
        gl, gp = g[g.is_lnc == 1], g[g.is_lnc == 0]
        rows.append({"decile": int(k), "n_lnc": int(len(gl)), "n_pc": int(len(gp)),
                     "med_r_lnc": float(gl.r.median()) if len(gl) else None,
                     "med_r_pc": float(gp.r.median()) if len(gp) else None})
    o["decile_matched"] = rows
    return o


for s in ["test", "val"]:
    res[s] = analyse(s)

t = res["test"]
c, (lo, hi) = t["H2_coef_is_lnc"], t["H2_ci"]
if c <= -0.05 and hi < 0:
    verdict = "pass"
elif -0.05 < c <= -0.02 and hi < 0:
    verdict = "grey"
else:
    verdict = "fail"
res["H2_verdict"] = verdict

# H3 coverage
res["H3"] = {"decima_lnc": int((var.gene_type == "lncRNA").sum()),
             "decima_pc": int((var.gene_type == "protein_coding").sum()),
             "gencode_lnc_approx": 18039, "gencode_pc_approx": 20017}
res["H3"]["frac_lnc"] = res["H3"]["decima_lnc"] / 18039
res["H3"]["frac_pc"] = res["H3"]["decima_pc"] / 20017

# exploratory: lncRNA relation to nearest protein-coding gene
pc = var[var.is_lnc == 0]
cls = {}
for g, row in var[var.is_lnc == 1].iterrows():
    q = pc[pc.chrom == row.chrom]
    s, e = row.gene_start, row.gene_end
    ov = q[(q.gene_start < e) & (q.gene_end > s)]
    if len(ov) and (ov.strand != row.strand).any():
        cls[g] = "antisense_overlap"
    elif len(ov):
        cls[g] = "sense_overlap"
    else:
        dist = np.minimum(np.abs(q.gene_start.values - e), np.abs(s - q.gene_end.values))
        cls[g] = "near_10kb" if len(dist) and dist.min() <= 10000 else "intergenic"
var["lnc_class"] = pd.Series(cls)
ex = {}
for s in ["test", "val", "train"]:
    d = var[(var.dataset == s) & (var.is_lnc == 1)]
    ex[s] = {k: {"n": int(len(g)), "median_r": float(g.r.median()),
                 "median_mean_counts": float(g.mean_counts.median())}
             for k, g in d.groupby("lnc_class")}
res["explore_lnc_class"] = ex

json.dump(res, open(f"{OUT}/DEC1.json", "w"), indent=1)
var[["gene_id", "gene_type", "dataset", "mean_counts", "sd_obs", "gene_length",
     "frac_nan", "pearson", "r", "lnc_class"]].to_csv(f"{OUT}/dec1_genes.tsv.gz", sep="\t")
print("done", verdict)
