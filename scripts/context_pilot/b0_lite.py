#!/usr/bin/env python
"""B0-lite: predict compartment change in icSHAPE reactivity (PROTOCOL_B0.md v1 section 5, frozen).
Step 'build': per-nucleotide table (needs RNAplfold from ~/envs/xaudit). Step 'model': baselines + HistGBR (~/envs/rnafate).
CNN is in b0_cnn.py (~/envs/te_benchmark)."""
import gzip, json, os, subprocess, sys, tempfile
import numpy as np, pandas as pd

ROOT = os.path.expanduser("~/ncrna-annot")
DATA = f"{ROOT}/data/context_pilot"
OUT = f"{ROOT}/outputs/context_pilot/B0"
os.makedirs(OUT, exist_ok=True)
TASKS = {"main_np_to_cy": ("293_np_vivo", "293_cy_vivo", "human"),
         "ch_to_np": ("293_ch_vivo", "293_np_vivo", "human"),
         "vitro_to_vivo_cy": ("293_cy_vitro", "293_cy_vivo", "human"),
         "mes_np_to_cy": ("mes_np_vivo", "mes_cy_vivo", "mouse")}
FILES = {f.split("_", 1)[1].replace(".out.txt.gz", ""): f for f in os.listdir(f"{DATA}/icshape")}
TRIM, W_SEQ, W_GC = 30, 12, 25
RNAPLFOLD = os.path.expanduser("~/envs/xaudit/bin/RNAplfold")


def read_out(name):
    d = {}
    with gzip.open(f"{DATA}/icshape/{FILES[name]}", "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            tid = f[0].split(".")[0]
            vals = np.array([np.nan if v == "NULL" else float(v) for v in f[3:]], dtype=np.float32)
            d[tid] = (int(f[1]), float(f[2]), vals)
    return d


def fetch_seqs(ex, fa):
    """samtools faidx in batches; returns {tid: seq}."""
    sam = os.path.expanduser("~/envs/lrtools/bin/samtools")
    comp = str.maketrans("ACGTN", "TGCAN")
    out = {}
    tids = list(ex)
    for i in range(0, len(tids), 2000):
        batch = tids[i:i + 2000]
        regs = [f"{c}:{s + 1}-{e}" for t in batch for c, s, e, _ in ex[t]]
        res = subprocess.run([sam, "faidx", fa] + regs, capture_output=True, text=True).stdout
        seqs, cur = {}, None
        for line in res.split("\n"):
            if line.startswith(">"):
                cur = line[1:]; seqs[cur] = []
            elif cur:
                seqs[cur].append(line.strip())
        for t in batch:
            parts = sorted(ex[t], key=lambda x: x[1])
            s = "".join("".join(seqs.get(f"{c}:{a + 1}-{b}", [])) for c, a, b, _ in parts).upper()
            if parts[0][3] == "-":
                s = s.translate(comp)[::-1]
            out[t] = s
    return out


def plfold_unpaired(seq_by_tid):
    res = {}
    with tempfile.TemporaryDirectory() as td:
        fa = os.path.join(td, "in.fa")
        with open(fa, "w") as fh:
            for t, s in seq_by_tid.items():
                fh.write(f">{t}\n{s.replace('T', 'U')}\n")
        subprocess.run(f"cd {td} && nice -n 19 {RNAPLFOLD} -W 150 -L 100 -u 1 < in.fa", shell=True, check=True, capture_output=True)
        for t in seq_by_tid:
            p = os.path.join(td, f"{t}_lunp")
            if not os.path.exists(p):
                continue
            v = np.loadtxt(p, comments="#", usecols=1, ndmin=1)
            res[t] = v.astype(np.float32)
    return res


def build(task):
    ref, tgt, species = TASKS[task]
    A, B = read_out(ref), read_out(tgt)
    tx = load_tx(species)
    rows = []
    common = [t for t in A if t in B and t in tx and len(tx[t][0]) == A[t][0] == B[t][0]]
    seqs = {t: tx[t][0] for t in common}
    unp = {}
    items = list(seqs.items())
    for i in range(0, len(items), 500):
        unp.update(plfold_unpaired(dict(items[i:i + 500])))
    for t in common:
        L = A[t][0]; ra, rb = A[t][2], B[t][2]
        n = min(len(ra), len(rb), L)
        idx = np.arange(TRIM, n - TRIM)
        ok = idx[~np.isnan(ra[idx]) & ~np.isnan(rb[idx])]
        if len(ok) == 0 or t not in unp:
            continue
        s = seqs[t]
        for i in ok:
            rows.append((t, tx[t][2], tx[t][1], i, float(ra[i]), float(rb[i]), float(A[t][1]), L, s[max(0, i - W_SEQ):i + W_SEQ + 1].rjust(W_SEQ + 1 + min(i, W_SEQ)).ljust(2 * W_SEQ + 1),
                         s[max(0, i - W_GC):i + W_GC + 1].count("G") + s[max(0, i - W_GC):i + W_GC + 1].count("C"), float(unp[t][i]) if i < len(unp[t]) else np.nan))
    df = pd.DataFrame(rows, columns=["tid", "gid", "ttype", "pos", "r_ref", "r_tgt", "rpkm_ref", "tlen", "win", "gc_cnt", "unpaired"])
    df["delta"] = df.r_tgt - df.r_ref
    df.to_pickle(f"{OUT}/{task}.pkl")
    return df


def load_tx(species):
    cache = f"{DATA}/tx_{species}.pkl"
    if os.path.exists(cache):
        return pd.read_pickle(cache)
    gtf = f"{ROOT}/data/human/gencode.v41.annotation.gtf.gz" if species == "human" else f"{ROOT}/data/batch2/mouse/gencode.vM39.annotation.gtf.gz"
    fa = f"{ROOT}/data/rna_fate/hg38.fa" if species == "human" else f"{ROOT}/data/context_pilot/mm39.fa"
    ex, info = {}, {}
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[2] != "exon":
                continue
            a = f[8]
            tid = a.split('transcript_id "')[1].split('"')[0].split(".")[0]
            gid = a.split('gene_id "')[1].split('"')[0].split(".")[0]
            tt = a.split('transcript_type "')[1].split('"')[0]
            ex.setdefault(tid, []).append((f[0], int(f[3]) - 1, int(f[4]), f[6]))
            info[tid] = (tt, gid)
    seqs = fetch_seqs(ex, fa)
    tx = {t: (seqs[t], info[t][0], info[t][1]) for t in seqs}
    pd.to_pickle(tx, cache)
    return tx


# ------------------------------------------------------------------ models
def featurize(df, with_ref):
    onehot = np.zeros((len(df), (2 * W_SEQ + 1) * 4), np.float32)
    lut = {"A": 0, "C": 1, "G": 2, "T": 3}
    for k, w in enumerate(df.win.values):
        for j, ch in enumerate(w):
            if ch in lut:
                onehot[k, j * 4 + lut[ch]] = 1
    tt = pd.get_dummies(df.ttype.where(df.ttype.isin(["protein_coding", "lncRNA"]), "other")).reindex(columns=["protein_coding", "lncRNA", "other"], fill_value=0).values
    cols = [onehot, (df.gc_cnt.values / (2 * W_GC + 1))[:, None], df.unpaired.fillna(0.5).values[:, None], tt, np.log(df.tlen.values)[:, None]]
    if with_ref:
        cols += [df.r_ref.values[:, None], np.log1p(df.rpkm_ref.values)[:, None]]
    return np.hstack(cols).astype(np.float32)


def model(task):
    from sklearn.ensemble import HistGradientBoostingRegressor
    from sklearn.linear_model import LinearRegression
    from sklearn.model_selection import GroupKFold
    df = pd.read_pickle(f"{OUT}/{task}.pkl")
    y = df.delta.values
    folds = list(GroupKFold(5).split(df, groups=df.gid.values))
    preds = {k: np.zeros(len(df)) for k in ("zero", "mean_shift", "ref_only", "HGB_E", "HGB_S")}
    XE, XS = featurize(df, True), featurize(df, False)
    for tr, te in folds:
        preds["mean_shift"][te] = y[tr].mean()
        preds["ref_only"][te] = LinearRegression().fit(df.r_ref.values[tr, None], y[tr]).predict(df.r_ref.values[te, None])
        for name, X in (("HGB_E", XE), ("HGB_S", XS)):
            m = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.05, random_state=0).fit(X[tr], y[tr])
            preds[name][te] = m.predict(X[te])
    pd.DataFrame({"tid": df.tid, "gid": df.gid, "ttype": df.ttype, "y": y, **preds}).to_pickle(f"{OUT}/{task}_preds.pkl")
    return summarize(task)


def summarize(task, extra=None):
    d = pd.read_pickle(f"{OUT}/{task}_preds.pkl")
    if extra is not None:
        for k, v in extra.items():
            d[k] = v
    models = [c for c in d.columns if c not in ("tid", "gid", "ttype", "y")]
    rng = np.random.default_rng(0)
    def stats(sub):
        o = {}
        for m in models:
            p = sub[m].values
            o[m] = {"pearson": float(np.corrcoef(sub.y, p)[0, 1]) if p.std() > 0 else 0.0,
                    "R2": float(1 - ((sub.y - p) ** 2).sum() / ((sub.y - sub.y.mean()) ** 2).sum())}
        return o
    res = {"task": task, "n_nt": int(len(d)), "n_genes": int(d.gid.nunique()),
           "n_lnc_genes": int(d[d.ttype == "lncRNA"].gid.nunique()), "n_lnc_nt": int((d.ttype == "lncRNA").sum()),
           "all": stats(d), "lncRNA": stats(d[d.ttype == "lncRNA"]) if (d.ttype == "lncRNA").sum() > 100 else None}
    simple = ["zero", "mean_shift", "ref_only"]
    learned = [m for m in models if m not in simple]
    bs = max(simple, key=lambda m: res["all"][m]["pearson"]); bl = max(learned, key=lambda m: res["all"][m]["pearson"])
    genes = d.gid.values; ug = np.unique(genes); pos = pd.Series(np.arange(len(d))).groupby(genes).apply(np.array)
    diffs = []
    for _ in range(2000):
        ix = np.concatenate(pos.loc[rng.choice(ug, len(ug))].values)
        a = d.y.values[ix]
        pl, ps = d[bl].values[ix], d[bs].values[ix]
        diffs.append(np.corrcoef(a, pl)[0, 1] - (np.corrcoef(a, ps)[0, 1] if ps.std() > 0 else 0.0))
    res["best_simple"], res["best_learned"] = bs, bl
    res["gain_pearson"] = res["all"][bl]["pearson"] - res["all"][bs]["pearson"]
    res["gain_ci"] = [float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))]
    res["G0"] = bool(res["n_lnc_genes"] >= 300 and res["n_lnc_nt"] >= 30000)
    res["G2"] = bool(res["all"][bl]["pearson"] >= 0.30 and res["all"][bl]["R2"] > 0)
    res["G3"] = bool(res["gain_pearson"] >= 0.05 and res["gain_ci"][0] > 0)
    json.dump(res, open(f"{OUT}/{task}_summary.json", "w"), indent=1)
    return res


if __name__ == "__main__":
    step, task = sys.argv[1], sys.argv[2]
    if step == "build":
        df = build(task); print(task, len(df), df.gid.nunique())
    elif step == "model":
        print(json.dumps(model(task), indent=1))
    elif step == "cnnexport":
        df = pd.read_pickle(f"{OUT}/{task}.pkl"); ref, tgt, sp = TASKS[task]
        tx = load_tx(sp); A = read_out(ref); ts = df.tid.unique()
        np.savez(f"{OUT}/{task}_cnn_in.npz", tid=df.tid.values.astype(str), pos=df.pos.values, gid=df.gid.values.astype(str),
                 delta=df.delta.values, tx_seq=np.array({t: tx[t][0] for t in ts}, dtype=object),
                 ref_tracks=np.array({t: A[t][2] for t in ts}, dtype=object))
    elif step == "cnnsum":
        print(json.dumps(summarize(task, {"CNN_E": np.load(f"{OUT}/{task}_cnn_E.npy"), "CNN_S": np.load(f"{OUT}/{task}_cnn_S.npy")}), indent=1))
