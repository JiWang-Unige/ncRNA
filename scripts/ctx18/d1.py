#!/usr/bin/env python
"""D1 (PROTOCOL_D1.md): lncRNA nuclear-enrichment gene x context interaction on Tabula Sapiens total-RNA.
Steps:  build   -> pseudobulk + bootstrap noise, saves npz (no statistics printed)
        analyse -> G0..G3, writes D1.json (run ONCE after freezing)
Run with ~/envs/ctx18/bin/python on ont."""
import os, sys, json, gzip
import numpy as np, pandas as pd, h5py, anndata as ad
from anndata.io import read_elem
from scipy import sparse
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge

ROOT = os.path.expanduser("~/ncrna-annot")
H5 = f"{ROOT}/data/ctx18/ts_total/TSP33_total_RNA.h5ad"
OUT = f"{ROOT}/outputs/ctx18/D1"
os.makedirs(OUT, exist_ok=True)
NBOOT, MINC, MIN_CTX, MIN_CPM = 50, 50, 15, 2.0
rng = np.random.default_rng(0)


# ---------------------------------------------------------------- build
def build():
    with h5py.File(H5, "r") as f:
        obs = read_elem(f["obs"]); var = read_elem(f["var"])
        X = read_elem(f["layers/raw_counts"])
    X = sparse.csr_matrix(X)
    keep = var.biotype.isin(["lncRNA", "protein_coding"]).values
    X = X[:, keep].tocsc().tocsr(); var = var[keep]
    t = obs.groupby(["tissue", "cell_type", "method"], observed=True).size().unstack(fill_value=0)
    ctx = t[(t.get("Cells", 0) >= MINC) & (t.get("Nuclei", 0) >= MINC)].index.tolist()
    G = X.shape[1]
    cpm = {m: np.zeros((len(ctx), G), np.float32) for m in ("Cells", "Nuclei")}
    boot = {m: np.zeros((len(ctx), NBOOT, G), np.float32) for m in ("Cells", "Nuclei")}
    for k, (ti, ct) in enumerate(ctx):
        for m in ("Cells", "Nuclei"):
            rows = np.nonzero(((obs.tissue == ti) & (obs.cell_type == ct) & (obs.method == m)).values)[0]
            Xr = X[rows]
            s = np.asarray(Xr.sum(0)).ravel()
            cpm[m][k] = s / s.sum() * 1e6
            W = rng.multinomial(len(rows), np.full(len(rows), 1 / len(rows)), size=NBOOT).astype(np.float32)
            B = np.asarray(Xr.T @ W.T).T                         # NBOOT x G
            boot[m][k] = B / B.sum(1, keepdims=True) * 1e6
        print("ctx", k, ti, ct, flush=True)
    np.savez_compressed(f"{OUT}/pseudobulk.npz", ctx=np.array(["|".join(c) for c in ctx]),
                        gene=var.gene_ids.str.split(".").str[0].values, sym=var.index.values,
                        biotype=var.biotype.values, cpm_cell=cpm["Cells"], cpm_nuc=cpm["Nuclei"],
                        boot_cell=boot["Cells"], boot_nuc=boot["Nuclei"])
    print("saved", len(ctx), "contexts", G, "genes")


# ---------------------------------------------------------------- analyse helpers
def two_way(M, iters=50):
    mu = M.mean(1); al = pd.Series(0.0, index=M.columns)
    for _ in range(iters):
        al = (M.sub(mu, axis=0)).median(0)
        mu = (M.sub(al, axis=1)).mean(1)
    return mu, al


def labels(z, biotype):
    cc, cn = z["cpm_cell"], z["cpm_nuc"]
    y = np.log2((cn + 1) / (cc + 1)); obsd = (cn + cc) / 2 >= MIN_CPM
    yb = np.log2((z["boot_nuc"] + 1) / (z["boot_cell"] + 1))      # ctx x boot x gene
    nv = yb.var(1)                                                # ctx x gene
    sel = (z["biotype"] == biotype) & (obsd.sum(0) >= MIN_CTX)
    Y = pd.DataFrame(np.where(obsd, y, np.nan)[:, sel].T, index=z["gene"][sel], columns=z["ctx"])
    NV = pd.DataFrame(np.where(obsd, nv, np.nan)[:, sel].T, index=z["gene"][sel], columns=z["ctx"])
    Y = Y[~Y.index.duplicated()]; NV = NV.loc[Y.index]
    return Y, NV


def g1(Y, NV):
    mu, al = two_way(Y)
    d = Y.sub(mu, axis=0).sub(al, axis=1)
    st = d.stack().dropna(); r = Y.sub(mu, axis=0).stack().dropna()
    noise = float(NV.stack().dropna().mean())
    gi = st.index.get_level_values(0).values; ug = np.unique(gi)
    pos = {g: np.nonzero(gi == g)[0] for g in ug}; v = st.values
    boots = [v[np.concatenate([pos[g] for g in rng.choice(ug, len(ug))])].var() - noise for _ in range(2000)]
    ti = float(v.var() - noise)
    res = {"n_genes": int(len(Y)), "n_obs": int(Y.notna().sum().sum()), "var_mu": float(mu.var()), "var_alpha": float(al.var()),
           "var_y_minus_mu": float(r.var()), "var_delta_hat": float(v.var()), "noise_var": noise, "true_interaction_var": ti,
           "ci": [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))], "fraction_beyond_mu": ti / float(r.var())}
    res["pass"] = bool(res["ci"][0] > 0 and res["fraction_beyond_mu"] >= 0.20)
    return res, mu


def encode_mu():
    D = f"{ROOT}/data/context_pilot"
    man = pd.read_csv(f"{D}/encode/manifest.tsv", sep="\t"); tpm = {}
    for r in man.itertuples():
        q = pd.read_csv(f"{D}/encode/{r.file}.tsv", sep="\t", usecols=["gene_id", "TPM"])
        q = q[~q.gene_id.str.contains("PAR_Y")]; q["gid"] = q.gene_id.str.split(".").str[0]
        tpm.setdefault((r.cell, r.fraction), []).append(q.groupby("gid").TPM.sum())
    cols = {}
    for c in sorted({c for c, _ in tpm}):
        n = pd.concat(tpm[(c, "nucleus")], axis=1).mean(1); y = pd.concat(tpm[(c, "cytosol")], axis=1).mean(1)
        i = n.index.intersection(y.index); ok = (n[i] + y[i]) / 2 >= 1
        cols[c] = np.log2((n[i] + 0.1) / (y[i] + 0.1))[ok]
    M = pd.DataFrame(cols)
    return M[M.notna().sum(1) >= 6].mean(1)


def ctx_features(z, train_idx):
    pc = z["biotype"] == "protein_coding"
    L = np.log2(z["cpm_cell"][:, pc] + 1)
    v = L[train_idx].var(0) > 0
    sc = StandardScaler().fit(L[train_idx][:, v]); p = PCA(n_components=10, random_state=0).fit(sc.transform(L[train_idx][:, v]))
    return p.transform(sc.transform(L[:, v]))


def mlp_pred(Gtr, Ctr, ytr, Gte, Cte, seed=0):
    import torch, torch.nn as nn
    torch.manual_seed(seed)
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    class T(nn.Module):
        def __init__(s, dg, dc):
            super().__init__()
            s.g = nn.Sequential(nn.Linear(dg, 64), nn.ReLU(), nn.Linear(64, 32)); s.c = nn.Sequential(nn.Linear(dc, 64), nn.ReLU(), nn.Linear(64, 32))
            s.bg = nn.Linear(dg, 1); s.bc = nn.Linear(dc, 1)
        def forward(s, g, c):
            return (s.g(g) * s.c(c)).sum(1) / 32 + s.bg(g).squeeze(1) + s.bc(c).squeeze(1)
    m = T(Gtr.shape[1], Ctr.shape[1]).to(dev); opt = torch.optim.Adam(m.parameters(), 1e-3, weight_decay=1e-4)
    g, c, y = (torch.tensor(a, dtype=torch.float32, device=dev) for a in (Gtr, Ctr, ytr))
    for _ in range(300):
        perm = torch.randperm(len(y), device=dev)
        for b in range(0, len(y), 4096):
            ix = perm[b:b + 4096]; opt.zero_grad(); ((m(g[ix], c[ix]) - y[ix]) ** 2).mean().backward(); opt.step()
    m.eval()
    with torch.no_grad():
        return m(torch.tensor(Gte, dtype=torch.float32, device=dev), torch.tensor(Cte, dtype=torch.float32, device=dev)).cpu().numpy()


def choose_alpha(F, y, groups, cands=(1e-1, 1e1, 1e3, 1e5)):
    best, ba = -np.inf, cands[0]
    for a in cands:
        s = []
        for g in np.unique(groups):
            tr, te = groups != g, groups == g
            p = Ridge(alpha=a).fit(F[tr], y[tr]).predict(F[te]); s.append(1 - ((y[te] - p) ** 2).sum() / (y[te] ** 2).sum())
        if np.mean(s) > best: best, ba = np.mean(s), a
    return ba


def run_split(Y, z, GF, tissue_of, gene_fold=None, shuffle_ctx=False, models=("M0", "M1", "M2", "M3", "M4")):
    """leave-one-tissue-out; if gene_fold given, test genes are also held out (G3)."""
    ctx = list(Y.columns); ctx_index = {c: i for i, c in enumerate(z["ctx"])}
    out = {m: [] for m in models}
    for T in sorted(set(tissue_of.values())):
        te_c = [c for c in ctx if tissue_of[c] == T]; tr_c = [c for c in ctx if tissue_of[c] != T]
        CF = ctx_features(z, [ctx_index[c] for c in tr_c])
        CF = pd.DataFrame(CF, index=z["ctx"])
        if shuffle_ctx:
            perm = rng.permutation(len(tr_c)); CF.loc[tr_c] = CF.loc[tr_c].values[perm]
        folds = [None] if gene_fold is None else sorted(set(gene_fold.values))
        for fo in folds:
            g_te = Y.index if fo is None else gene_fold.index[gene_fold == fo]
            g_tr = Y.index if fo is None else gene_fold.index[gene_fold != fo]
            Mt = Y.loc[g_tr, tr_c]; mu, al = two_way(Mt)
            dtr = Mt.sub(mu, axis=0).sub(al, axis=1)
            if fo is None:
                mu_te = mu
            else:   # unseen genes: their mean comes from training contexts of the same gene (labels of other tissues)
                mu_te, _ = two_way(Y.loc[g_te, tr_c])
            test = Y.loc[g_te, te_c].stack().dropna()
            test = test[test.index.get_level_values(0).isin(mu_te.dropna().index)]
            gi, ci = test.index.get_level_values(0), test.index.get_level_values(1)
            resid = test.values - mu_te.loc[gi].values
            a_c = pd.Series(resid, index=ci).groupby(level=0).median()
            yte = resid - a_c.loc[ci].values
            st = dtr.stack().dropna(); gtr, ctr = st.index.get_level_values(0), st.index.get_level_values(1)
            gfe = lambda g, m: np.column_stack([GF.loc[g].values, m.loc[g].values])
            Gtr, Ctr = gfe(gtr, mu), CF.loc[ctr].values
            Gte, Cte = gfe(gi, mu_te), CF.loc[ci].values
            sg = StandardScaler().fit(Gtr); Gtr, Gte = sg.transform(Gtr), sg.transform(Gte)
            preds = {"M0": np.zeros(len(yte))}
            if "M1" in models:
                d = np.linalg.norm(CF.loc[tr_c].values[None] - CF.loc[te_c].values[:, None], axis=2)
                m1 = {}
                for j, c in enumerate(te_c):
                    nn_ = np.argsort(d[j])[:3]; w = 1 / (d[j][nn_] + 1e-6)
                    m1[c] = sum(w[k] * dtr[tr_c[nn_[k]]] for k in range(3)) / w.sum()
                preds["M1"] = np.array([m1[c].get(g, 0.0) if not np.isnan(m1[c].get(g, 0.0)) else 0.0 for g, c in zip(gi, ci)])
            if "M2" in models:
                F = lambda Gm, Cm: np.einsum("ni,nj->nij", Gm, Cm).reshape(len(Gm), -1)
                Ftr = F(Gtr, Ctr); grp = np.array([tissue_of[c] for c in ctr])
                a = choose_alpha(Ftr, st.values, grp)
                preds["M2"] = Ridge(alpha=a).fit(Ftr, st.values).predict(F(Gte, Cte))
            if "M3" in models:
                preds["M3"] = mlp_pred(Gtr, Ctr, st.values, Gte, Cte)
            if "M4" in models:
                p4 = np.zeros(len(yte))
                for g, ix in pd.Series(np.arange(len(gi)), index=gi).groupby(level=0):
                    if g not in dtr.index: continue
                    row = dtr.loc[g, tr_c].dropna()
                    if len(row) < 5: continue
                    r = Ridge(alpha=10).fit(CF.loc[row.index].values, row.values)
                    p4[ix.values] = r.predict(CF.loc[[ci[i] for i in ix.values]].values)
                preds["M4"] = p4
            for m in models:
                out[m].append(pd.DataFrame({"tissue": T, "gene": gi, "ctx": ci, "y": yte, "p": preds[m]}))
    return {m: pd.concat(v, ignore_index=True) for m, v in out.items()}


def score(pred, ref="M1"):
    r = {}
    for m, d in pred.items():
        pt = d.groupby("tissue").apply(lambda x: np.corrcoef(x.y, x.p)[0, 1] if x.p.std() > 0 else 0.0)
        r[m] = {"pearson_median_tissue": float(pt.median()), "pearson_by_tissue": pt.round(3).to_dict(),
                "R2": float(1 - ((d.y - d.p) ** 2).sum() / (d.y ** 2).sum())}
    return r


def gain_ci(pred, best, ref="M1"):
    a, b = pred[best], pred[ref] if ref in pred else None
    ug = a.gene.unique(); idx = {g: np.nonzero(a.gene.values == g)[0] for g in ug}
    ya, pa = a.y.values, a.p.values; pb = b.p.values if b is not None else np.zeros(len(ya))
    def stat(ix):
        ss = (ya[ix] ** 2).sum()
        return (1 - ((ya[ix] - pa[ix]) ** 2).sum() / ss) - max(0.0, 1 - ((ya[ix] - pb[ix]) ** 2).sum() / ss)
    bs = [stat(np.concatenate([idx[g] for g in rng.choice(ug, len(ug))])) for _ in range(1000)]
    return float(stat(np.arange(len(ya)))), [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]


def analyse():
    z = dict(np.load(f"{OUT}/pseudobulk.npz", allow_pickle=True))
    tissue_of = {c: c.split("|")[0] for c in z["ctx"]}
    rep = {}
    Yl, NVl = labels(z, "lncRNA"); Yp, NVp = labels(z, "protein_coding")
    # G0
    enc = encode_mu()
    g1l, mul = g1(Yl, NVl); g1p, mup = g1(Yp, NVp)
    common = mul.index.intersection(enc.index)
    from scipy.stats import spearmanr
    rho = spearmanr(mul.loc[common], enc.loc[common]).correlation
    rep["G0"] = {"n_common": int(len(common)), "spearman": float(rho), "pass": bool(rho >= 0.30)}
    rep["G1_lncRNA"] = g1l; rep["G1_protein_coding_control"] = g1p
    # gene features
    feats = pd.read_csv(f"{ROOT}/data/rna_fate/e1_features.tsv.gz", sep="\t")
    feats["gid"] = feats.gene.str.split(".").str[0]; feats = feats.drop_duplicates("gid").set_index("gid")
    fcols = [c for c in feats.columns if c.startswith(("b_", "s_"))]
    F = feats.reindex(Yl.index)[fcols].dropna()
    GF = pd.DataFrame(PCA(n_components=20, random_state=0).fit_transform(StandardScaler().fit_transform(F)), index=F.index)
    Y = Yl.loc[GF.index]
    rep["G2_n_genes_with_features"] = int(len(Y))
    pred = run_split(Y, z, GF, tissue_of)
    sc = score(pred)
    best = max(("M2", "M3", "M4"), key=lambda m: sc[m]["pearson_median_tissue"])
    gain, ci = gain_ci(pred, best)
    sc.update({"best": best, "R2_gain_vs_M1": gain, "R2_gain_ci": ci})
    sc["pass"] = bool(sc[best]["pearson_median_tissue"] >= 0.20 and gain >= 0.05 and ci[0] > 0)
    sc["kill"] = bool(sc[best]["pearson_median_tissue"] < 0.10)
    shuf = run_split(Y, z, GF, tissue_of, shuffle_ctx=True, models=("M0", best))
    sc["diag_shuffled_ctx_R2"] = score(shuf)[best]["R2"]
    rep["G2"] = sc
    for m, d in pred.items(): d.to_csv(f"{OUT}/pred_G2_{m}.tsv.gz", sep="\t", index=False)
    # G3: unseen genes
    gf = pd.Series(rng.integers(0, 5, len(Y)), index=Y.index)
    p3 = run_split(Y, z, GF, tissue_of, gene_fold=gf, models=("M0", "M2", "M3"))
    s3 = score(p3)
    b3 = max(("M2", "M3"), key=lambda m: s3[m]["R2"])
    g3, c3 = gain_ci(p3, b3, ref="none")
    s3.update({"best": b3, "R2": g3, "R2_ci": c3, "pass": bool(g3 > 0 and c3[0] > 0)})
    rep["G3"] = s3
    json.dump(rep, open(f"{OUT}/D1.json", "w"), indent=1, default=str)
    print(json.dumps(rep, indent=1, default=str)[:8000])


if __name__ == "__main__":
    {"build": build, "analyse": analyse}[sys.argv[1]]()
