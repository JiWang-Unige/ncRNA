#!/usr/bin/env python
"""D0: gene x cell interaction audit for lncRNA nuclear/cytoplasmic localisation (PROTOCOL_D0.md v1, frozen).
ENCODE Gingeras polyA+ nucleus/cytosol, GRCh38 GENCODE V29 RSEM gene quantifications, 13 cell lines.
Run with ~/envs/rnafate/bin/python (sklearn 1.3.2)."""
import gzip, json, os, sys
import numpy as np, pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

ROOT = os.path.expanduser("~/ncrna-annot")
DATA = f"{ROOT}/data/context_pilot"
OUT = f"{ROOT}/outputs/context_pilot/D0"
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(0)

# ---------------- labels
man = pd.read_csv(f"{DATA}/encode/manifest.tsv", sep="\t")
tpm = {}
for r in man.itertuples():
    q = pd.read_csv(f"{DATA}/encode/{r.file}.tsv", sep="\t", usecols=["gene_id", "TPM"])
    q["gid"] = q.gene_id.str.split(".").str[0]
    q = q[~q.gene_id.str.contains("PAR_Y")].groupby("gid").TPM.sum()
    tpm.setdefault((r.cell, r.fraction), []).append(q)
GENES = pd.read_csv(f"{DATA}/genes_v41.tsv", sep="\t")  # built from gencode.v41 GTF (gene_id, gene_type, gene_name)
gtype = GENES.drop_duplicates("gene_id").set_index("gene_id").gene_type

cells = sorted({c for c, _ in tpm})
Y, NOISE, EXPR = {}, {}, {}
for c in cells:
    nu, cy = tpm[(c, "nucleus")], tpm[(c, "cytosol")]
    nmean = pd.concat(nu, axis=1).mean(1); cmean = pd.concat(cy, axis=1).mean(1)
    idx = nmean.index.intersection(cmean.index)
    ok = ((nmean[idx] + cmean[idx]) / 2) >= 1
    reps = [np.log2((cy[min(i, len(cy) - 1)][idx] + 0.1) / (nu[min(i, len(nu) - 1)][idx] + 0.1)) for i in range(2)]
    Y[c] = ((reps[0] + reps[1]) / 2)[ok]
    if len(nu) == 2 and len(cy) == 2:
        NOISE[c] = (reps[0] - reps[1])[ok]
    EXPR[c] = np.log2((nmean[idx] + cmean[idx]) / 2 + 1)
Ymat = pd.DataFrame(Y)
lnc = Ymat.index[Ymat.index.isin(gtype[gtype == "lncRNA"].index)]
pc = Ymat.index[Ymat.index.isin(gtype[gtype == "protein_coding"].index)]


def two_way(M, iters=50):
    mu = M.mean(1); al = pd.Series(0.0, index=M.columns)
    for _ in range(iters):
        al = (M.sub(mu, axis=0)).median(0)
        mu = (M.sub(al, axis=1)).mean(1)
    return mu, al


def analyse(genes, label):
    M = Ymat.loc[genes]
    M = M[M.notna().sum(1) >= 6]
    res = {"set": label, "n_genes": int(len(M)), "n_obs": int(M.notna().sum().sum()),
           "cells_per_gene_median": float(M.notna().sum(1).median())}
    mu, al = two_way(M)
    delta = M.sub(mu, axis=0).sub(al, axis=1)
    resid_mu = M.sub(mu, axis=0)
    # ---- G1: variance decomposition with replicate noise
    nz = []
    for c, d in NOISE.items():
        dd = d.reindex(M.index).dropna()
        if len(dd) > 50:
            e = EXPR[c].reindex(dd.index)
            q = pd.qcut(e, 5, labels=False, duplicates="drop")
            v = dd.groupby(q).var() / 4
            w = q.value_counts(normalize=True).reindex(v.index)
            nz.append(float((v * w).sum()))
    noise = float(np.mean(nz))
    dvals = delta.stack().dropna().values
    rvals = resid_mu.stack().dropna().values
    gene_idx = delta.stack().dropna().index.get_level_values(0).values
    ug = np.unique(gene_idx)
    pos = {g: np.nonzero(gene_idx == g)[0] for g in ug}
    boots = []
    for _ in range(2000):
        pick = rng.choice(ug, len(ug))
        ix = np.concatenate([pos[g] for g in pick])
        boots.append(dvals[ix].var() - noise)
    true_int = float(dvals.var() - noise)
    res["G1"] = {"var_y_minus_mu": float(rvals.var()), "var_delta_hat": float(dvals.var()), "noise_var_mean": noise,
                 "true_interaction_var": true_int, "ci": [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))],
                 "fraction_beyond_mu": true_int / float(rvals.var()), "var_mu": float(mu.var()), "var_alpha": float(al.var()),
                 "noise_cells": sorted(NOISE)}
    res["G1"]["pass"] = bool(res["G1"]["ci"][0] > 0 and res["G1"]["fraction_beyond_mu"] >= 0.20)
    return res, M


def predictions(M, X):
    """Leave-one-cell-line-out prediction of delta for M0..M3."""
    allrbp = [l.strip() for l in open(f"{DATA}/eclip_rbps.txt") if l.strip()]
    rbpexp = rbp_matrix(allrbp)
    out = {m: [] for m in ("M0", "M1", "M2", "M3")}
    per_cell = []
    for c in M.columns:
        tr = [x for x in M.columns if x != c]
        Mt = M[tr]
        mu, al = two_way(Mt)
        dtr = Mt.sub(mu, axis=0).sub(al, axis=1)
        test = M[c].dropna()
        test = test[test.index.isin(mu.index) & test.index.isin(X.index)]
        a_c = float((test - mu[test.index]).median())
        dte = test - mu[test.index] - a_c
        pca = PCA(n_components=8, random_state=0).fit(StandardScaler().fit(rbpexp.loc[tr]).transform(rbpexp.loc[tr]))
        sc = StandardScaler().fit(rbpexp.loc[tr])
        Z = pd.DataFrame(pca.transform(sc.transform(rbpexp)), index=rbpexp.index)
        # M1: nearest 3 training cells by Z
        dist = np.linalg.norm(Z.loc[tr].values - Z.loc[c].values, axis=1)
        nn = np.argsort(dist)[:3]; w = 1 / (dist[nn] + 1e-6)
        m1 = sum(w[k] * dtr[tr[nn[k]]].reindex(dte.index).fillna(0) for k in range(3)) / w.sum()
        # training pairs
        st = dtr.stack().dropna(); st = st[st.index.get_level_values(0).isin(X.index)]
        g_tr, c_tr = st.index.get_level_values(0), st.index.get_level_values(1)
        Xg_tr, Zc_tr = X.loc[g_tr].values, Z.loc[c_tr].values
        F_tr = np.einsum("ni,nj->nij", Xg_tr, Zc_tr).reshape(len(st), -1)
        F_te = np.einsum("ni,nj->nij", X.loc[dte.index].values, np.repeat(Z.loc[[c]].values, len(dte), 0)).reshape(len(dte), -1)
        alpha = choose_alpha(F_tr, st.values, c_tr)
        m2 = Ridge(alpha=alpha).fit(F_tr, st.values).predict(F_te)
        m3 = kernel_pred(Xg_tr, Zc_tr, st.values, X.loc[dte.index].values, Z.loc[c].values, c_tr)
        for name, p in (("M0", np.zeros(len(dte))), ("M1", m1.values), ("M2", m2), ("M3", m3)):
            out[name].append(pd.DataFrame({"cell": c, "gene": dte.index, "y": dte.values, "p": p}))
        per_cell.append(c)
    return {k: pd.concat(v) for k, v in out.items()}


def rbp_matrix(rbps):
    # whole-cell log expression of eCLIP RBPs, symbols mapped via GENCODE v41 gene names
    m = dict(zip(GENES.gene_name, GENES.gene_id))
    ids = [m[r] for r in rbps if r in m]
    return pd.DataFrame({c: EXPR[c].reindex(ids) for c in cells}).T.fillna(0)


def choose_alpha(F, y, groups, cands=(1e-2, 1e0, 1e2, 1e4)):
    best, ba = -np.inf, cands[0]
    ug = np.unique(groups)
    for a in cands:
        s = []
        for g in ug:
            tr, te = groups != g, groups == g
            p = Ridge(alpha=a).fit(F[tr], y[tr]).predict(F[te])
            s.append(1 - ((y[te] - p) ** 2).sum() / (y[te] ** 2).sum())
        if np.mean(s) > best:
            best, ba = np.mean(s), a
    return ba


def kernel_pred(Xg, Zc, y, Xte, zte, groups, lam_cands=(1e-1, 1e0, 1e1, 1e2), n=6000):
    idx = rng.choice(len(y), min(n, len(y)), replace=False)
    Xg, Zc, y, groups = Xg[idx], Zc[idx], y[idx], groups[idx]
    d2 = ((Xg[:2000, None, :] - Xg[None, :2000, :]) ** 2).sum(-1)
    bw = np.median(d2[d2 > 0])
    def K(A, B, za, zb):
        dd = (A ** 2).sum(1)[:, None] + (B ** 2).sum(1)[None, :] - 2 * A @ B.T
        return np.exp(-dd / bw) * (za @ zb.T)
    Ktr = K(Xg, Xg, Zc, Zc)
    lam = lam_cands[1]
    best = -np.inf
    ug = np.unique(groups)
    for l in lam_cands:
        s = []
        for g in ug[:4]:
            tr, te = groups != g, groups == g
            a = np.linalg.solve(Ktr[np.ix_(tr, tr)] + l * np.eye(tr.sum()), y[tr])
            p = Ktr[np.ix_(te, tr)] @ a
            s.append(1 - ((y[te] - p) ** 2).sum() / (y[te] ** 2).sum())
        if np.mean(s) > best:
            best, lam = np.mean(s), l
    a = np.linalg.solve(Ktr + lam * np.eye(len(y)), y)
    Kte = K(Xte, Xg, np.repeat(zte[None, :], len(Xte), 0), Zc)
    return Kte @ a


def score(pred):
    r = {}
    for m, d in pred.items():
        pc_ = d.groupby("cell").apply(lambda x: np.corrcoef(x.y, x.p)[0, 1] if x.p.std() > 0 else 0.0)
        r[m] = {"pearson_median": float(pc_.median()), "pearson_by_cell": pc_.round(3).to_dict(),
                "R2_interaction": float(1 - ((d.y - d.p) ** 2).sum() / (d.y ** 2).sum())}
    return r


if __name__ == "__main__":
    feats = pd.read_csv(f"{ROOT}/data/rna_fate/e1_features.tsv.gz", sep="\t")
    feats["gid"] = feats.gene.str.split(".").str[0]
    feats = feats.drop_duplicates("gid").set_index("gid")
    fcols = [c for c in feats.columns if c.startswith(("b_", "s_"))]
    report = {}
    for label, genes in (("lncRNA", lnc), ("protein_coding_control", pc)):
        res, M = analyse(genes, label)
        if label == "lncRNA":
            F = feats.reindex(M.index)[fcols].dropna()
            X = pd.DataFrame(PCA(n_components=20, random_state=0).fit_transform(StandardScaler().fit_transform(F)), index=F.index)
            pred = predictions(M.loc[M.index.isin(X.index)], X)
            res["G2"] = score(pred)
            best = max(("M2", "M3"), key=lambda m: res["G2"][m]["pearson_median"])
            res["G2"]["best"] = best
            res["G2"]["pass"] = bool(res["G2"][best]["pearson_median"] >= 0.30 and
                                    res["G2"][best]["R2_interaction"] - max(0.0, res["G2"]["M1"]["R2_interaction"]) >= 0.05)
            for m, d in pred.items():
                d.to_csv(f"{OUT}/pred_{m}.tsv.gz", sep="\t", index=False)
        report[label] = res
    json.dump(report, open(f"{OUT}/D0.json", "w"), indent=1, default=str)
    print(json.dumps(report, indent=1, default=str)[:6000])
