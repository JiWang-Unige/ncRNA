#!/usr/bin/env python3
"""S3 first half (PLAN_S3.md, card subjournal-portfolio 2026-11-04): which sequence / annotation features go with
v1.9 lncRNA locus detection.  Reads the frozen P1.5 predictions only (no decoding, no training).

Unit = v50 lnc gene of the chromosome; y = strict L1 hit (same-strand exon overlap) by a final lnc prediction
(confidence >= tau*) of each main run.  Features (sequence and public annotation only):
  cpg_oe   CpG observed / expected in outermost TSS +- 500 bp      cage   FANTOM5 CAGE peak within 200 bp of any TSS
  u1_kb    GT[AG]AGT per kb in the first 1 kb downstream of the TSS (sense)   pas_kb  AATAAA / ATTAAA per kb, same window
  u1_pas   u1_kb - pas_kb                                           polya  PolyASite 2.0 cluster within 100 bp of any TES
  repeat   soft-masked fraction of the merged exons                 log_span, log_len (merged exons), n_exon, multi
  era      gene id (no version) present in GENCODE v41 (1) or not (0)
R1 per run: logistic regression (standardised, L2 C = 100) on all features, gene bootstrap 2,000 (seed 0) CIs;
"associated" = CI of the 3-run mean coefficient excludes 0 and the 3 run coefficients share its sign.
R2 detection rate (mean over runs) by feature quartile; R3 = R1 with era x feature interactions.

  s3_features.py --chrom chr18 --out outputs/PHASE1/S3/S3_chr18.json
"""
import argparse
import gzip
import json
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1"))
from ncrna_annot.eval import ruler  # noqa: E402
import p15_score as P  # noqa: E402

FEATS = ("cpg_oe", "cage", "u1_kb", "pas_kb", "u1_pas", "polya", "repeat", "log_span", "log_len", "n_exon", "multi", "era")
RUNS = ("g1b-42-ke", "g1b-123-ke", "g1b-456-ke")
RC = str.maketrans("ACGTacgtN", "TGCAtgcaN")
N_BOOT = 2000


def read_fasta(path):
    seq = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq)


def count_motifs(s, motifs):
    n = 0
    for m in motifs:
        i = s.find(m)
        while i >= 0:
            n += 1
            i = s.find(m, i + 1)
    return n


def cpg_oe(s):
    s = s.upper()
    c, g, cg = s.count("C"), s.count("G"), s.count("CG")
    return cg * len(s) / (c * g) if c and g else 0.0


def load_peaks(path, chrom, strand_col, start_col, end_col, chrom_col=0, skip_header=False, chrom_name=None):
    out = {"+": [], "-": []}
    with gzip.open(path, "rt") as fh:
        if skip_header:
            fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[chrom_col] != (chrom_name or chrom):
                continue
            if f[strand_col] in out:
                out[f[strand_col]].append((int(f[start_col]) + int(f[end_col])) // 2)   # 0-based mid point
    return {k: np.sort(np.array(v, np.int64)) for k, v in out.items()}


def near(sorted_pos, x, d):
    if not len(sorted_pos):
        return False
    i = np.searchsorted(sorted_pos, x - d)
    return i < len(sorted_pos) and sorted_pos[i] <= x + d


def gene_features(g, seq, cage, polya):
    st = g["strand"]
    txs = [ruler._merge(sorted(ex)) for ex in g["tx"].values()]
    s0 = min(t[0][0] for t in txs)
    e0 = max(t[-1][1] for t in txs)
    tss_all = [t[0][0] - 1 if st == "+" else t[-1][1] - 1 for t in txs]          # 0-based
    tes_all = [t[-1][1] - 1 if st == "+" else t[0][0] - 1 for t in txs]
    tss = s0 - 1 if st == "+" else e0 - 1
    prom = seq[max(tss - 500, 0):tss + 501]
    down = seq[tss:tss + 1000] if st == "+" else seq[max(tss - 999, 0):tss + 1].translate(RC)[::-1]
    down = down.upper()
    u1 = count_motifs(down, ("GTAAGT", "GTGAGT")) / max(len(down), 1) * 1000
    pas = count_motifs(down, ("AATAAA", "ATTAAA")) / max(len(down), 1) * 1000
    merged = ruler._merge(sorted(x for t in txs for x in t))
    exseq = "".join(seq[a - 1:b] for a, b in merged)
    longest = max(txs, key=lambda t: sum(b - a + 1 for a, b in t))
    return {"cpg_oe": cpg_oe(prom), "cage": float(any(near(cage[st], x, 200) for x in tss_all)),
            "u1_kb": u1, "pas_kb": pas, "u1_pas": u1 - pas, "polya": float(any(near(polya[st], x, 100) for x in tes_all)),
            "repeat": sum(ch.islower() for ch in exseq) / max(len(exseq), 1),
            "log_span": float(np.log10(e0 - s0 + 1)), "log_len": float(np.log10(sum(b - a + 1 for a, b in longest))),
            "n_exon": float(len(longest)), "multi": float(len(longest) > 1)}


def fit_logit(X, y, C=100.0):
    n, p = X.shape
    Xa = np.hstack([np.ones((n, 1)), X])

    def f(w):
        z = Xa @ w
        ll = np.logaddexp(0, z) - y * z
        g = Xa.T @ (1 / (1 + np.exp(-z)) - y)
        reg = 0.5 / C * (w[1:] @ w[1:])
        g[1:] += w[1:] / C
        return ll.sum() + reg, g
    return minimize(f, np.zeros(p + 1), jac=True, method="L-BFGS-B").x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True, choices=("chr18", "chr20"))
    ap.add_argument("--config", default="configs/p15_v19b.json")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    t0 = time.time()
    cfg = json.load(open(f"{ROOT}/{a.config}"))
    tau = cfg["decode"]["tau"]
    tg = cfg["targets"][a.chrom]
    truth = ruler.read_gtf(f"{ROOT}/{tg['truth']}", tg["truth_chrom"])
    v41 = ruler.read_gtf(f"{ROOT}/data/human/gencode.v41.annotation.gtf.gz", a.chrom)
    ids41 = {gid.split(".")[0] for gid in v41}
    seq = read_fasta(f"{ROOT}/data/human/{a.chrom}.fa")
    cage = load_peaks(f"{ROOT}/data/bioval/raw/fantom5_cage_peaks.bed.gz", a.chrom, 5, 6, 7)
    polya = load_peaks(f"{ROOT}/data/bioval/raw/polyasite2_GRCh38.96.tsv.gz", a.chrom, 5, 1, 2, skip_header=True,
                       chrom_name=a.chrom.replace("chr", ""))
    lnc = {gid: g for gid, g in truth.items() if g["cls"] == "lnc"}
    gids = sorted(lnc)
    rows = []
    for gid in gids:
        r = gene_features(lnc[gid], seq, cage, polya)
        r["era"] = float(gid.split(".")[0] in ids41)
        rows.append(r)
    X = np.array([[r[k] for k in FEATS] for r in rows], float)
    idx = ruler.ExonIndex(lnc)
    Y, res_paths = {}, {}
    for run in RUNS:
        path = f"{ROOT}/{cfg['out_dir']}/{a.chrom}/{run}/genes.json"   # P15b holds the chr18 reproduction read
        if not os.path.exists(path):
            path = f"{ROOT}/outputs/PHASE1/P15/{a.chrom}/{run}/genes.json"  # main runs were read on chr20 in batch 1
        res_paths[run] = path
        genes = P.load_genes(path)
        hit = set()
        for p in P.final_preds(genes, tau):
            if p["cls"] == "lnc":
                hit |= set(idx.hits(p["strand"], p["exons"]))
        Y[run] = np.array([gid in hit for gid in gids], float)
    mu, sd = X.mean(0), X.std(0)
    sd[sd == 0] = 1
    Z = (X - mu) / sd
    rng = np.random.default_rng(0)
    boots = [rng.integers(0, len(gids), len(gids)) for _ in range(N_BOOT)]
    res = {"chrom": a.chrom, "n_genes": len(gids), "tau": tau, "features": FEATS,
           "feature_mean_sd": {k: [float(mu[i]), float(sd[i])] for i, k in enumerate(FEATS)},
           "detected_frac": {r: float(Y[r].mean()) for r in RUNS}, "pred_files": res_paths}

    def r1(Zm, names):
        coefs = {r: fit_logit(Zm, Y[r])[1:] for r in RUNS}
        mean = np.mean([coefs[r] for r in RUNS], 0)
        bm = np.array([np.mean([fit_logit(Zm[b], Y[r][b])[1:] for r in RUNS], 0) for b in boots])
        lo, hi = np.percentile(bm, [2.5, 97.5], 0)
        out = {}
        for i, k in enumerate(names):
            signs = {np.sign(coefs[r][i]) for r in RUNS}
            assoc = (lo[i] > 0 or hi[i] < 0) and len(signs) == 1
            out[k] = {"per_run": [float(coefs[r][i]) for r in RUNS], "mean": float(mean[i]), "ci95": [float(lo[i]), float(hi[i])],
                      "associated": bool(assoc)}
        return out
    res["R1"] = r1(Z, FEATS)
    ymean = np.mean([Y[r] for r in RUNS], 0)
    r2 = {}
    for i, k in enumerate(FEATS):
        v = X[:, i]
        if len(np.unique(v)) <= 2:
            r2[k] = {str(u): {"n": int((v == u).sum()), "detected": float(ymean[v == u].mean())} for u in np.unique(v)}
        else:
            q = np.quantile(v, [0.25, 0.5, 0.75])
            b = np.searchsorted(q, v, "right")
            r2[k] = {f"Q{j + 1}": {"n": int((b == j).sum()), "detected": float(ymean[b == j].mean()) if (b == j).any() else None,
                                   "range": [float(v[b == j].min()), float(v[b == j].max())] if (b == j).any() else None} for j in range(4)}
    res["R2"] = r2
    era = X[:, FEATS.index("era")]
    other = [k for k in FEATS if k != "era"]
    inter = np.hstack([Z, Z[:, [FEATS.index(k) for k in other]] * ((era - era.mean()) / max(era.std(), 1e-9))[:, None]])
    res["R3"] = {k: v for k, v in r1(inter, list(FEATS) + [f"{k}_x_era" for k in other]).items() if k.endswith("_x_era")}
    res["seconds"] = round(time.time() - t0)
    os.makedirs(os.path.dirname(f"{ROOT}/{a.out}"), exist_ok=True)
    json.dump(res, open(f"{ROOT}/{a.out}", "w"), indent=1)
    print(json.dumps({"chrom": a.chrom, "detected_frac": res["detected_frac"],
                      "R1": {k: (round(v["mean"], 3), [round(x, 3) for x in v["ci95"]], v["associated"]) for k, v in res["R1"].items()},
                      "R3_associated": [k for k, v in res["R3"].items() if v["associated"]], "seconds": res["seconds"]}), flush=True)


if __name__ == "__main__":
    main()
