#!/usr/bin/env python3
"""v2.0 chr18 development (protocol section 17 / DRAFT_V20 sections 3.3-3.4): coding-potential features, gene-level
calibration and tau, on the chr18 reads of the three retrained main runs (outputs/PHASE1/V20/dev/chr18/<run>/).

  v20_calib.py hexamer      -> data/p1v20/hexamer_table.json  (training chromosomes only: chr1-8, 10-17, 19, 21, 22, X)
  v20_calib.py calibrate    -> outputs/PHASE1/V20/dev/CALIB.json

Hexamer table (CPAT-like): in-frame hexamer counts (step 3 from the first base) of the CDS of v50 Ensembl_canonical
protein-coding transcripts and of the spliced exons of v50 Ensembl_canonical lncRNA transcripts.  Hexamer score of a
predicted lnc transcript = mean over the in-frame hexamers of its longest ORF of log((c+1)/(Nc+4096)) -
log((n+1)/(Nn+4096)) (0 when the ORF has no hexamer).
Features per predicted lnc gene (before tau): conf (mean lnc posterior), mrna_mass (mean same-strand mRNA posterior
on the exon bases), log_len, n_exon, log_orf, orf_frac, fickett, hexamer.  y = strict L1 hit (same-strand exon
overlap with a v50 lnc gene of chr18).
Calibration: logistic, L2 C = 1, standardised, pooled over the three runs (vip_diag.fit_predict).  Cross-fit on the
chr18 halves (gene start < / >= 40 Mb): out-of-fold ECE (10 equal-width bins) and AUPRC of the calibrated
probability vs the raw conf.  Keep rule (DRAFT 3.4): keep iff oof AUPRC_cal >= oof AUPRC_raw and oof ECE_cal <
oof ECE_raw.  Final coefficients from all of chr18.
tau (section 16 rule): grid 0.2-0.9 step 0.1; among the tau with three-run mean strict Pr >= 0.43 the one with the
largest mean lnc L1 Sn (P1.5 ruler, mRNA genes unfiltered); if no tau reaches 0.43, the tau with the largest mean
strict Pr.  Chosen separately for the baseline (raw conf) and for v2.0 (calibrated probability; raw conf if the
calibration is not kept).
"""
import gzip
import json
import math
import os
import sys
import time

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1/vipsania"))
import vip_diag as V  # noqa: E402  (fickett, longest_orf, transcript, fit_predict; it also sets up src/ paths)
import p15_score as P  # noqa: E402
from ncrna_annot.eval import ruler  # noqa: E402
from ncrna_annot.grammar import g1  # noqa: E402

DEV = f"{ROOT}/outputs/PHASE1/V20/dev"
HEX = f"{ROOT}/data/p1v20/hexamer_table.json"
RUNS = ("v20-g1b-42-ke", "v20-g1b-123-ke", "v20-g1b-456-ke")
TRAIN = [f"chr{i}" for i in list(range(1, 9)) + list(range(10, 18)) + [19, 21, 22]] + ["chrX"]
FEATS = ["conf", "mrna_mass", "log_len", "n_exon", "log_orf", "orf_frac", "fickett", "hexamer"]
SPLIT, TAU_GRID, PR_MIN = 40_000_000, [round(0.2 + 0.1 * i, 1) for i in range(8)], 0.43
COMP = str.maketrans("ACGTN", "TGCAN")


def hexamer_table():
    import prep_p1
    seqs = prep_p1.read_genome(set(TRAIN))
    tx = {}
    with gzip.open(f"{ROOT}/data/human/gencode.v50.annotation.gtf.gz", "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 9 or f[0] not in seqs or 'tag "Ensembl_canonical"' not in f[8]:
                continue
            gt = f[8].split('gene_type "')[1].split('"')[0]
            if (gt, f[2]) not in (("protein_coding", "CDS"), ("lncRNA", "exon")):
                continue
            tid = f[8].split('transcript_id "')[1].split('"')[0]
            tx.setdefault(tid, ["c" if gt == "protein_coding" else "n", f[0], f[6], []])[3].append((int(f[3]), int(f[4])))
    counts = {"c": {}, "n": {}}
    ntx = {"c": 0, "n": 0}
    for kind, c, s, iv in tx.values():
        seq = "".join(seqs[c][a - 1:b] for a, b in sorted(iv)).upper()
        if s == "-":
            seq = seq.translate(COMP)[::-1]
        ntx[kind] += 1
        d = counts[kind]
        for i in range(0, len(seq) - 5, 3):
            h = seq[i:i + 6]
            if "N" not in h:
                d[h] = d.get(h, 0) + 1
    hexs = sorted(set(counts["c"]) | set(counts["n"]))
    tab = {"rule": __doc__.split("Hexamer table")[1].split("Features")[0].strip(), "chroms": TRAIN, "n_transcripts": ntx,
           "total": {k: int(sum(v.values())) for k, v in counts.items()},
           "counts": {h: [counts["c"].get(h, 0), counts["n"].get(h, 0)] for h in hexs}}
    json.dump(tab, open(HEX, "w"))
    print(json.dumps({"n_transcripts": ntx, "total": tab["total"], "n_hexamers": len(hexs)}), flush=True)


def hexamer_score(orf, tab):
    if len(orf) < 6:
        return 0.0
    Nc, Nn, cnt = tab["total"]["c"], tab["total"]["n"], tab["counts"]
    v = []
    for i in range(0, len(orf) - 5, 3):
        c, n = cnt.get(orf[i:i + 6], (0, 0))
        v.append(math.log((c + 1) / (Nc + 4096)) - math.log((n + 1) / (Nn + 4096)))
    return float(np.mean(v)) if v else 0.0


def orf_seq(s):
    best, seq = 0, ""
    for f in range(3):
        start = None
        for i in range(f, len(s) - 2, 3):
            c = s[i:i + 3]
            if start is None and c == "ATG":
                start = i
            elif start is not None and c in V.STOPS:
                if i + 3 - start > best:
                    best, seq = i + 3 - start, s[start:i + 3]
                start = None
    return seq


def gene_rows(genes, post, codes, tab, lidx=None):
    rows = []
    for k, g in enumerate(genes):
        if g["type"] != "l":
            continue
        s = g["strand"]
        ms = [g1.IDX[f"m{p}{s}"] for p in g1.PARTS]
        exb = np.concatenate([np.arange(a - 1, b) for a, b in g["exons"]])
        tx = V.transcript(codes, g)
        orf = orf_seq(tx)
        r = {"k": k, "start": g["exons"][0][0], "conf": g["conf"],
             "mrna_mass": float(np.asarray(post[exb][:, ms], np.float64).sum(1).mean()),
             "log_len": math.log(len(tx)), "n_exon": float(len(g["exons"])), "log_orf": math.log1p(len(orf)),
             "orf_frac": len(orf) / len(tx), "fickett": V.fickett(tx), "hexamer": hexamer_score(orf, tab)}
        if lidx is not None:
            r["y"] = 1.0 if lidx.hits(s, g["exons"]) else 0.0
        rows.append(r)
    return rows


def ece(p, y, bins=10):
    b = np.minimum((np.asarray(p) * bins).astype(int), bins - 1)
    return float(sum(abs(p[b == i].mean() - y[b == i].mean()) * (b == i).sum() for i in range(bins) if (b == i).any()) / len(p))


def fit(X, y, C=1.0):
    """coefficients of vip_diag.fit_predict's model, so they can be written into the config."""
    from scipy.optimize import minimize
    mu, sd = X.mean(0), X.std(0)
    sd[sd == 0] = 1.0
    A1 = np.hstack([(X - mu) / sd, np.ones((len(X), 1))])

    def f(w):
        z = A1 @ w
        p = 1 / (1 + np.exp(-z))
        g = A1.T @ (p - y)
        g[:-1] += w[:-1] / C
        return np.sum(np.logaddexp(0, z) - y * z) + 0.5 / C * np.dot(w[:-1], w[:-1]), g
    w = minimize(f, np.zeros(A1.shape[1]), jac=True, method="L-BFGS-B", options={"maxiter": 2000}).x
    return {"features": FEATS, "mean": mu.tolist(), "sd": sd.tolist(), "coef": w[:-1].tolist(), "intercept": float(w[-1])}


def apply(model, X):
    Z = (X - np.array(model["mean"])) / np.array(model["sd"])
    return 1 / (1 + np.exp(-(Z @ np.array(model["coef"]) + model["intercept"])))


def pick_tau(genes_by_run, scores_by_run, truth):
    grid = []
    for t in TAU_GRID:
        sn, pr = [], []
        for run, genes in genes_by_run.items():
            gs = [dict(g, conf=scores_by_run[run].get(k, g["conf"])) for k, g in enumerate(genes)]
            fin, _ = P.score_preds(truth, P.final_preds(gs, t), boot=0)
            sn.append(fin["lnc"]["L1_locus"]["sn"])
            pr.append(fin["lnc"]["L1_locus"]["pr_strict"] or 0.0)
        grid.append({"tau": t, "sn": float(np.mean(sn)), "pr_strict": float(np.mean(pr)), "sn_runs": sn, "pr_runs": pr})
    ok = [x for x in grid if x["pr_strict"] >= PR_MIN]
    best = max(ok, key=lambda x: x["sn"]) if ok else max(grid, key=lambda x: x["pr_strict"])
    return best["tau"], grid, bool(ok)


def calibrate():
    t0 = time.time()
    tab = json.load(open(HEX))
    truth = ruler.read_gtf(f"{ROOT}/outputs/PHASE1/DIAG_L2/chr18_v50.gtf", "chr18")
    lidx = ruler.ExonIndex({k: v for k, v in truth.items() if v["cls"] == "lnc"})
    codes = np.load(f"{ROOT}/data/p1/chr18_codes.npy", mmap_mode="r")
    genes_by_run, rows = {}, []
    for run in RUNS:
        rd = f"{DEV}/chr18/{run}"
        genes = P.load_genes(f"{rd}/genes.json")
        genes_by_run[run] = genes
        for r in gene_rows(genes, np.load(f"{rd}/post.npy", mmap_mode="r"), codes, tab, lidx):
            r["run"] = run
            rows.append(r)
    X = np.array([[r[f] for f in FEATS] for r in rows], float)
    y = np.array([r["y"] for r in rows])
    raw = X[:, FEATS.index("conf")]
    half = np.array([r["start"] >= SPLIT for r in rows])
    oof = np.zeros(len(y))
    for h in (False, True):
        oof[half == h] = V.fit_predict(X[half != h], y[half != h], X[half == h])
    ap = lambda s: V.own_splice.average_precision(s, y > 0)
    res = {"n_genes": len(rows), "n_pos": int(y.sum()), "by_run": {r: sum(1 for x in rows if x["run"] == r) for r in RUNS},
           "oof": {"auprc_raw": ap(raw), "auprc_cal": ap(oof), "ece_raw": ece(raw, y), "ece_cal": ece(oof, y)},
           "oof_by_half": {}}
    for h, name in ((False, "first"), (True, "second")):
        m = half == h
        res["oof_by_half"][name] = {"n": int(m.sum()), "auprc_raw": V.own_splice.average_precision(raw[m], y[m] > 0),
                                    "auprc_cal": V.own_splice.average_precision(oof[m], y[m] > 0), "ece_raw": ece(raw[m], y[m]),
                                    "ece_cal": ece(oof[m], y[m])}
    o = res["oof"]
    res["keep"] = bool(o["auprc_cal"] >= o["auprc_raw"] and o["ece_cal"] < o["ece_raw"])
    model = fit(X, y)
    res["model"] = model
    cal = apply(model, X)
    res["in_sample"] = {"auprc_cal": ap(cal), "ece_cal": ece(cal, y)}
    scores = {run: {} for run in RUNS}
    for r, p in zip(rows, cal):
        scores[r["run"]][r["k"]] = float(p)
    tb, grid_b, ok_b = pick_tau(genes_by_run, {run: {} for run in RUNS}, truth)
    res["tau_baseline"], res["tau_grid_baseline"], res["tau_baseline_reached_pr"] = tb, grid_b, ok_b
    if res["keep"]:
        tv, grid_v, ok_v = pick_tau(genes_by_run, scores, truth)
    else:
        tv, grid_v, ok_v = tb, grid_b, ok_b
    res["tau_v20"], res["tau_grid_v20"], res["tau_v20_reached_pr"] = tv, grid_v, ok_v
    res["hexamer_table_sha256"] = P.sha256(HEX)
    res["seconds"] = round(time.time() - t0)
    json.dump(res, open(f"{DEV}/CALIB.json", "w"), indent=1)
    print(json.dumps({k: res[k] for k in ("n_genes", "n_pos", "oof", "keep", "tau_baseline", "tau_v20", "seconds")}), flush=True)


if __name__ == "__main__":
    {"hexamer": hexamer_table, "calibrate": calibrate}[sys.argv[1]]()
