"""DEC3: ENTEx total RNA-seq vs GTEx poly(A) vs Decima, organ level (PROTOCOL_DEC3.md).

python dec3.py <OUT>   (ctx18 env)
"""
import gzip
import json
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.expanduser("~/ncrna-annot")
EX = f"{ROOT}/data/ctx18/entex"
DEC = f"{ROOT}/data/ctx18/decima/Supplementary_file_1.h5ad"
GTF = f"{ROOT}/data/human/gencode.v50.annotation.gtf.gz"
OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
NB = 2000
rng = np.random.default_rng(0)

MAP = {  # organ: (ENTEx tissues, GTEx tissues)
    "adipose": (["subcutaneous adipose tissue", "omental fat pad"],
                ["Adipose - Subcutaneous", "Adipose - Visceral (Omentum)"]),
    "adrenal": (["adrenal gland"], ["Adrenal Gland"]),
    "breast": (["breast epithelium"], ["Breast - Mammary Tissue"]),
    "gut": (["transverse colon", "sigmoid colon", "stomach", "Peyer's patch",
             "esophagus muscularis mucosa", "esophagus squamous epithelium",
             "gastroesophageal sphincter"],
            ["Colon - Transverse", "Colon - Sigmoid", "Stomach", "Small Intestine - Terminal Ileum",
             "Esophagus - Muscularis", "Esophagus - Mucosa", "Esophagus - Gastroesophageal Junction"]),
    "heart": (["heart left ventricle", "right atrium auricular region"],
              ["Heart - Left Ventricle", "Heart - Atrial Appendage"]),
    "liver": (["right lobe of liver"], ["Liver"]),
    "lung": (["upper lobe of left lung"], ["Lung"]),
    "muscle": (["gastrocnemius medialis"], ["Muscle - Skeletal"]),
    "ovary": (["ovary"], ["Ovary"]),
    "pancreas": (["body of pancreas"], ["Pancreas"]),
    "prostate": (["prostate gland"], ["Prostate"]),
    "skin": (["lower leg skin", "suprapubic skin"],
             ["Skin - Sun Exposed (Lower leg)", "Skin - Not Sun Exposed (Suprapubic)"]),
    "spleen": (["spleen"], ["Spleen"]),
    "uterus": (["uterus"], ["Uterus"]),
    "vasculature": (["thoracic aorta", "ascending aorta"], ["Artery - Aorta"]),
}


def strip(g):
    return str(g).split(".")[0]


def gtf_info():
    res = {}
    with gzip.open(GTF, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            p = line.split("\t", 9)
            if p[2] != "gene":
                continue
            a = p[8]
            g = strip(a.split('gene_id "')[1].split('"')[0])
            bt = a.split('gene_type "')[1].split('"')[0]
            res[g] = (bt, int(p[4]) - int(p[3]) + 1)
    return res


def main():
    import anndata as ad
    import patsy
    organs = list(MAP)
    # ENTEx
    man = pd.read_csv(f"{EX}/manifest.tsv", sep="\t", header=None,
                      names=["exp", "tissue", "donor", "acc", "md5", "href"])
    tis = {}
    for t, g in man.groupby("tissue"):
        mats = []
        for acc in g.acc:
            q = pd.read_csv(f"{EX}/gq/{acc}.tsv", sep="\t", usecols=["gene_id", "TPM"])
            q = q[q.gene_id.str.startswith("ENSG")]
            q.index = q.gene_id.map(strip)
            mats.append(q.TPM[~q.index.duplicated()])
        tis[t] = pd.concat(mats, axis=1).mean(1)
    tot = pd.DataFrame({o: pd.concat([tis[t] for t in MAP[o][0] if t in tis], axis=1).mean(1)
                        for o in organs})
    # GTEx
    gt = pd.read_csv(f"{EX}/gtex_v8_median_tpm.gct.gz", sep="\t", skiprows=2)
    gt.index = gt.Name.map(strip)
    gt = gt[~gt.index.duplicated()]
    pa = pd.DataFrame({o: gt[MAP[o][1]].mean(1) for o in organs})
    # Decima
    a = ad.read_h5ad(DEC)
    obs, var = a.obs, a.var
    P = np.asarray(a.layers["preds"], dtype=np.float32)
    dec = {}
    for o in organs:
        idx = np.where((obs.organ == o).values)[0]
        dec[o] = P[idx].mean(0) if len(idx) else None
    organs = [o for o in organs if dec[o] is not None]
    res = {"organs": organs, "n_organs": len(organs),
           "tracks_per_organ": {o: int((obs.organ == o).sum()) for o in organs}}
    if len(organs) < 10:
        res["stop"] = "organs < 10"
        json.dump(res, open(f"{OUT}/DEC3.json", "w"), indent=1)
        return
    D = pd.DataFrame(np.vstack([dec[o] for o in organs]).T, columns=organs,
                     index=var.gene_id.map(strip).values)
    D = D[~D.index.duplicated()]
    dset = pd.Series(var.dataset.astype(str).values, index=var.gene_id.map(strip).values)
    dset = dset[~dset.index.duplicated()]
    info = gtf_info()
    genes = [g for g in D.index if g in tot.index and g in pa.index and g in info
             and info[g][0] in ("lncRNA", "protein_coding")]
    T, A, Dm = np.log1p(tot.loc[genes, organs].values), np.log1p(pa.loc[genes, organs].values), \
        D.loc[genes, organs].values
    ok = ((np.expm1(T) >= 1).sum(1) >= 5) & ((np.expm1(A) >= 1).sum(1) >= 5) & \
        (T.std(1) > 0) & (A.std(1) > 0) & (Dm.std(1) > 0)

    def pr(x, y):
        x = x - x.mean(1, keepdims=True)
        y = y - y.mean(1, keepdims=True)
        return (x * y).sum(1) / np.sqrt((x * x).sum(1) * (y * y).sum(1))

    G = pd.DataFrame({
        "gene": np.array(genes)[ok],
        "is_lnc": [int(info[g][0] == "lncRNA") for g in np.array(genes)[ok]],
        "dataset": [dset[g] for g in np.array(genes)[ok]],
        "log_len": [np.log10(info[g][1]) for g in np.array(genes)[ok]],
        "r_tot": pr(Dm[ok], T[ok]), "r_pA": pr(Dm[ok], A[ok]),
        "lm_tot": np.log10(np.expm1(T[ok]).mean(1) + 1e-3), "sd_tot": T[ok].std(1),
        "lm_pA": np.log10(np.expm1(A[ok]).mean(1) + 1e-3), "sd_pA": A[ok].std(1)})
    G["d"] = G.r_tot - G.r_pA
    G.to_csv(f"{OUT}/dec3_genes.tsv.gz", sep="\t", index=False)

    def coef(d, form, term="is_lnc"):
        y, X = patsy.dmatrices(form, d, return_type="dataframe")
        b, *_ = np.linalg.lstsq(X.values, y.values[:, 0], rcond=None)
        return b[list(X.columns).index(term)]

    def boot(d, form):
        c = coef(d, form)
        ix = np.arange(len(d))
        cb = [coef(d.iloc[rng.choice(ix, len(ix))].reset_index(drop=True), form) for _ in range(NB)]
        return float(c), [float(np.percentile(cb, 2.5)), float(np.percentile(cb, 97.5))]

    F1 = "r_tot ~ is_lnc + cr(lm_tot, df=4) + cr(sd_tot, df=4) + log_len"
    F1p = "r_pA ~ is_lnc + cr(lm_pA, df=4) + cr(sd_pA, df=4) + log_len"
    F2 = ("d ~ is_lnc + cr(lm_tot, df=3) + cr(sd_tot, df=3) + cr(lm_pA, df=3) + cr(sd_pA, df=3)"
          " + log_len")
    H = G[G.dataset.isin(["test", "val"])].reset_index(drop=True)
    Tr = G[G.dataset == "train"].reset_index(drop=True)

    def summ(d):
        return {str(k): {"n": int(len(g)), "median_r_tot": float(g.r_tot.median()),
                         "median_r_pA": float(g.r_pA.median())} for k, g in d.groupby("is_lnc")}

    res["groups_heldout"], res["groups_train"] = summ(H), summ(Tr)
    c, ci = boot(H, F1)
    res["R1"] = {"coef": c, "ci": ci, "pass": bool(c <= -0.05 and ci[1] < 0)}
    c, ci = boot(H, F1p)
    res["R1_polyA_report"] = {"coef": c, "ci": ci}
    c, ci = boot(H, F2)
    res["R2"] = {"coef": c, "ci": ci, "pass": bool(c <= -0.03 and ci[1] < 0)}
    res["train_report"] = {"R1": float(coef(Tr, F1)), "R1_polyA": float(coef(Tr, F1p)),
                           "R2": float(coef(Tr, F2))}
    if not res["R1"]["pass"]:
        res["decision"] = "not_replicated"
    elif res["R2"]["pass"]:
        res["decision"] = "replicated_assay_specific"
    else:
        res["decision"] = "replicated_external_generalisation"
    json.dump(res, open(f"{OUT}/DEC3.json", "w"), indent=1)
    print("done", res["decision"])


main()
