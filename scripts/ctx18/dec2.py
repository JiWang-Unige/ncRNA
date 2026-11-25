"""DEC2: Decima vs Tabula Sapiens total-RNA, lncRNA blind spot (PROTOCOL_DEC2.md).

python dec2.py prep    <OUT>   (ctx18 env)  -> matching, gene sets, measured matrix
python dec2.py infer   <OUT>   (decima env) -> rep0 predictions for GPU gene set
python dec2.py analyse <OUT>   (ctx18 env)  -> DEC2.json
"""
import gzip
import json
import os
import sys

import numpy as np

ROOT = os.path.expanduser("~/ncrna-annot")
PB = f"{ROOT}/outputs/ctx18/D1/pseudobulk.npz"
DEC = f"{ROOT}/data/ctx18/decima/Supplementary_file_1.h5ad"
GTF = f"{ROOT}/data/human/gencode.v50.annotation.gtf.gz"
FA = f"{ROOT}/data/rna_fate/hg38.fa"
STEP, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
NB = 2000
CTX, TSS_POS = 524_288, 163_840

ORG = {"Bladder": ["urinary"], "Large_intestine": ["gut"], "Small_intestine": ["gut"],
       "Lung": ["lung", "airway"], "Mucosa": ["mouth"], "Salivary_gland": ["mouth"],
       "Tongue": ["mouth"], "Muscle": ["muscle"], "Prostate": ["prostate"],
       "Vasculature": ["vasculature"]}
CT = {"Alveolar type I cell": ["type I pneumocyte"],
      "Alveolar type II cell": ["type II pneumocyte"],
      "Basal epithelial cell": ["basal cell", "respiratory basal cell"],
      "Ciliated epithelial cell": ["ciliated cell"],
      "Club cell": ["club cell"],
      "Colonocyte": ["enterocyte"], "Enterocyte": ["enterocyte"],
      "BEST4+ enterocyte": ["enterocyte"],
      "Endothelial cell": ["blood vessel endothelial cell", "capillary endothelial cell",
                           "endothelial cell of vascular tree"],
      "Lymphatic endothelial cell": ["endothelial cell of lymphatic vessel"],
      "Fibroblast": ["fibroblast"], "Myofibroblast": ["myofibroblast cell"],
      "Goblet cell": ["goblet cell"], "Macrophage": ["macrophage"],
      "Mast cell": ["mast cell"], "Monocyte": ["classical monocyte"],
      "Plasma cell": ["plasma cell"], "Smooth muscle cell": ["smooth muscle cell"],
      "T cell": ["CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell"],
      "Luminal epithelial cell": ["luminal cell of prostate epithelium"],
      "Neuron": ["neuron"]}


def strip(g):
    return str(g).split(".")[0]


def prep():
    import anndata as ad
    z = np.load(PB, allow_pickle=True)
    ctx, gid = z["ctx"], np.array([strip(g) for g in z["gene"]])
    bt, cell, nuc = z["biotype"], z["cpm_cell"], z["cpm_nuc"]
    a = ad.read_h5ad(DEC, backed="r")
    o = a.obs
    tracks, keep = [], []
    for i, c in enumerate(ctx):
        t, ct = c.split("|")
        if t not in ORG or ct not in CT:
            continue
        idx = np.where((o.disease == "healthy").values & o.organ.isin(ORG[t]).values
                       & o.cell_type.isin(CT[ct]).values)[0]
        if len(idx):
            keep.append(i)
            tracks.append(idx.tolist())
    print("matched contexts", len(keep))
    if len(keep) < 20:
        json.dump({"stop": "matched contexts < 20", "n": len(keep)}, open(f"{OUT}/DEC2.json", "w"))
        sys.exit(0)
    C = cell[keep]
    ok = np.isin(bt, ["lncRNA", "protein_coding"]) & ((C >= 1).sum(0) >= 10)
    M = np.log1p(C)
    ok &= M.std(0) > 0
    # nuclear enrichment over all 73 contexts
    y = np.log2((nuc + 1) / (cell + 1))
    obs = (nuc + cell) / 2 >= 2
    nobs = obs.sum(0)
    mu = np.where(nobs >= 10, np.where(obs, y, 0).sum(0) / np.maximum(nobs, 1), np.nan)
    var = a.var
    dec_ids = {strip(g): j for j, g in enumerate(var.gene_id)}
    dset = var.dataset.astype(str).values
    in_set = np.array([g in dec_ids for g in gid])
    G = np.where(ok)[0]
    rng = np.random.default_rng(0)
    out_lnc = G[(bt[G] == "lncRNA") & ~in_set[G]]
    held = lambda j: dset[dec_ids[gid[j]]] in ("test", "val")
    in_lnc_held = np.array([j for j in G if in_set[j] and bt[j] == "lncRNA" and held(j)])
    in_pc_held = np.array([j for j in G if in_set[j] and bt[j] == "protein_coding" and held(j)])
    gpu = np.concatenate([rng.choice(out_lnc, min(3000, len(out_lnc)), replace=False),
                          in_lnc_held,
                          rng.choice(in_pc_held, min(500, len(in_pc_held)), replace=False)])
    np.savez_compressed(f"{OUT}/prep.npz", keep=np.array(keep), ctx=ctx[keep],
                        tracks=np.array(json.dumps(tracks)), M=M.astype(np.float32),
                        ok=ok, mu=mu, gid=gid, bt=bt, in_set=in_set, gpu=gpu)
    print("genes ok", int(ok.sum()), "out_lnc", len(out_lnc), "gpu set", len(gpu))


def read_gtf(ids):
    want, res = set(ids), {}
    with gzip.open(GTF, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            p = line.split("\t", 9)
            if p[2] != "gene":
                continue
            g = strip(p[8].split('gene_id "')[1].split('"')[0])
            if g in want:
                res[g] = (p[0], int(p[3]) - 1, int(p[4]), p[6])
    return res


class Fasta:
    def __init__(self, fa):
        self.f = open(fa, "rb")
        self.idx = {}
        for line in open(fa + ".fai"):
            n, ln, off, lb, lby = line.split("\t")[:5]
            self.idx[n] = (int(ln), int(off), int(lb), int(lby))

    def fetch(self, chrom, s, e):
        ln, off, lb, lby = self.idx[chrom]
        if s < 0 or e > ln:
            return None
        a = off + (s // lb) * lby + s % lb
        b = off + (e // lb) * lby + e % lb
        self.f.seek(a)
        return self.f.read(b - a).replace(b"\n", b"").decode().upper()


COMP = str.maketrans("ACGTN", "TGCAN")


def infer():
    import torch
    from decima.data.dataset import SeqDataset
    from decima.hub import load_decima_model
    p = np.load(f"{OUT}/prep.npz", allow_pickle=True)
    gid, gpu = p["gid"], p["gpu"]
    coords = read_gtf(gid[gpu])
    fa = Fasta(FA)
    seqs, st, en, names = [], [], [], []
    for j in gpu:
        g = gid[j]
        if g not in coords:
            continue
        ch, s, e, strand = coords[g]
        L = e - s
        if strand == "+":
            ws = s - TSS_POS
            seq = fa.fetch(ch, ws, ws + CTX)
        else:
            we = e + TSS_POS
            seq = fa.fetch(ch, we - CTX, we)
            seq = seq.translate(COMP)[::-1] if seq else None
        if seq is None or len(seq) != CTX:
            continue
        seqs.append(seq)
        st.append(TSS_POS)
        en.append(min(TSS_POS + L, CTX))
        names.append(int(j))
    print("windows", len(seqs), flush=True)
    model = load_decima_model("v1_rep0", device="cuda")
    preds = []
    B = 100
    for k in range(0, len(seqs), B):
        ds = SeqDataset(seqs[k:k + B], st[k:k + B], en[k:k + B])
        r = model.predict_on_dataset(ds, device="cuda", batch_size=2, num_workers=2,
                                     float_precision="bf16-mixed")
        r = r["expression"] if isinstance(r, dict) else r
        r = np.asarray(r)
        if r.shape[0] == 8856 and r.shape[1] != 8856:
            r = r.T
        preds.append(r.astype(np.float32))
        print("done", k + len(r), flush=True)
    np.savez_compressed(f"{OUT}/rep0_preds.npz", idx=np.array(names), preds=np.vstack(preds))


def analyse():
    import anndata as ad
    import pandas as pd
    import patsy
    rng = np.random.default_rng(1)
    p = np.load(f"{OUT}/prep.npz", allow_pickle=True)
    tracks = json.loads(str(p["tracks"]))
    M, ok, mu, gid, bt, in_set = p["M"], p["ok"], p["mu"], p["gid"], p["bt"], p["in_set"]
    a = ad.read_h5ad(DEC)
    var = a.var
    dec_ids = {strip(g): j for j, g in enumerate(var.gene_id)}
    dset = var.dataset.astype(str).values
    P = np.asarray(a.layers["preds"], dtype=np.float32)  # tracks x genes
    coords = read_gtf(gid[ok])

    def ctxpred(trackpred):  # trackpred: (8856,) -> (n_ctx,)
        return np.array([trackpred[t].mean() for t in tracks])

    def pear(x, y):
        return float(np.corrcoef(x, y)[0, 1]) if x.std() > 0 and y.std() > 0 else np.nan

    def row(j, pr):
        m = M[:, j]
        ch = coords.get(gid[j])
        return {"j": int(j), "r": pear(pr, m), "log_mean": float(np.log10(np.expm1(m).mean() + 1e-3)),
                "sd": float(m.std()), "log_len": float(np.log10(ch[2] - ch[1])) if ch else np.nan,
                "is_lnc": int(bt[j] == "lncRNA"), "mu": float(mu[j])}

    res = {"n_ctx": int(len(tracks)), "n_genes_ok": int(ok.sum())}
    rows = []
    for j in np.where(ok & in_set)[0]:
        k = dec_ids[gid[j]]
        d = row(j, ctxpred(P[:, k]))
        d["dataset"] = dset[k]
        rows.append(d)
    D = pd.DataFrame(rows).dropna(subset=["r", "log_len"])
    COV = "cr(log_mean, df=4) + cr(sd, df=4) + log_len"

    def coef(d, form, term):
        y, X = patsy.dmatrices(form, d, return_type="dataframe")
        b, *_ = np.linalg.lstsq(X.values, y.values[:, 0], rcond=None)
        return b[list(X.columns).index(term)]

    def boot(d, form, term, scale=1.0):
        c = coef(d, form, term) * scale
        idx = np.arange(len(d))
        cb = [coef(d.iloc[rng.choice(idx, len(idx))].reset_index(drop=True), form, term) * scale
              for _ in range(NB)]
        return float(c), [float(np.percentile(cb, 2.5)), float(np.percentile(cb, 97.5))]

    def verdict(c, ci):
        if c <= -0.05 and ci[1] < 0:
            return "pass"
        if -0.05 < c <= -0.02 and ci[1] < 0:
            return "grey"
        return "fail"

    def summ(d, by):
        return {str(k): {"n": int(len(g)), "median_r": float(g.r.median())} for k, g in d.groupby(by)}

    # H1
    H = D[D.dataset.isin(["test", "val"])].reset_index(drop=True)
    c, ci = boot(H, "r ~ is_lnc + " + COV, "is_lnc")
    res["H1"] = {"groups_heldout": summ(H, "is_lnc"), "groups_train": summ(D[D.dataset == "train"], "is_lnc"),
                 "coef_is_lnc": c, "ci": ci, "verdict": verdict(c, ci)}
    Tr = D[D.dataset == "train"].reset_index(drop=True)
    res["H1"]["train_coef_is_lnc"] = float(coef(Tr, "r ~ is_lnc + " + COV, "is_lnc"))
    # H2
    res["H2"] = {}
    for name, sub in [("lncRNA", D[D.is_lnc == 1]), ("protein_coding", D[D.is_lnc == 0])]:
        s = sub.dropna(subset=["mu"]).reset_index(drop=True)
        iqr = float(np.subtract(*np.percentile(s.mu, [75, 25])))
        c, ci = boot(s, "r ~ mu + C(dataset) + " + COV, "mu", scale=iqr)
        res["H2"][name] = {"n": int(len(s)), "iqr_mu": iqr, "effect_per_iqr": c, "ci": ci}
    res["H2"]["verdict"] = verdict(res["H2"]["lncRNA"]["effect_per_iqr"], res["H2"]["lncRNA"]["ci"])
    # S0 + H3
    q = np.load(f"{OUT}/rep0_preds.npz")
    idx, R = q["idx"], q["preds"]
    # deviation 1: compare own rep0 with the published single-model layer v1_rep0
    import h5py
    meta = f"{ROOT}/data/ctx18/decima/metadata.h5ad"
    mv = ad.read_h5ad(meta, backed="r").var
    mids = {strip(g): j for j, g in enumerate(mv.gene_id)}
    sel = [(i, mids[gid[j]]) for i, j in enumerate(idx) if in_set[j]]
    cols = sorted(set(c for _, c in sel))
    pos = {c: k for k, c in enumerate(cols)}
    r0 = h5py.File(meta)["layers/v1_rep0"][:, cols]
    s0 = [pear(R[i], r0[:, pos[c]]) for i, c in sel]
    s0e = [pear(R[i], P[:, dec_ids[gid[j]]]) for i, j in enumerate(idx) if in_set[j]]
    res["S0"] = {"n": len(s0), "median_pearson_rep0_vs_published_rep0": float(np.nanmedian(s0)),
                 "median_pearson_rep0_vs_published_ensemble": float(np.nanmedian(s0e))}
    res["S0"]["pass"] = res["S0"]["median_pearson_rep0_vs_published_rep0"] >= 0.95
    if res["S0"]["pass"]:
        rows = []
        for i, j in enumerate(idx):
            d = row(j, ctxpred(R[i]))
            d["is_out"] = int(not in_set[j])
            rows.append(d)
        E = pd.DataFrame(rows).dropna(subset=["r", "log_len"])
        L = E[E.is_lnc == 1].reset_index(drop=True)
        c, ci = boot(L, "r ~ is_out + " + COV, "is_out")
        res["H3"] = {"groups_lnc": summ(L, "is_out"), "pc_in_set_rep0": summ(E[E.is_lnc == 0], "is_lnc"),
                     "coef_is_out": c, "ci": ci, "verdict": verdict(c, ci)}
        E.to_csv(f"{OUT}/dec2_gpu_genes.tsv.gz", sep="\t", index=False)
    D.to_csv(f"{OUT}/dec2_inset_genes.tsv.gz", sep="\t", index=False)
    v = [res["H1"]["verdict"], res["H2"]["verdict"], res.get("H3", {}).get("verdict", "not_run")]
    res["decision"] = "blind_spot" if "pass" in v else "close"
    json.dump(res, open(f"{OUT}/DEC2.json", "w"), indent=1)
    print("done", v, res["decision"])


{"prep": prep, "infer": infer, "analyse": analyse}[STEP]()
