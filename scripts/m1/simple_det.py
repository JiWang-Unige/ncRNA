#!/usr/bin/env python3
"""M1 experiment 2 (docs/m1/PLAN_M1_SIMPLE.md): a simple sequence-composition lncRNA locus detector as a control for
the 1.2B backbone + grammar.  CPU only, no model, no annotation input except Tiberius (optional, sequence-only).

Detector
  * every 50-bp grid boundary x on both strands gets composition features (CpG o/e, GC, CpG density, repeats in
    +-500 bp; GC in +-5 kb, CpG o/e in +-2 kb; sense GC / AT skew, CpG density, U1 / PAS motif density in 1 kb down-
    and upstream; 16 sense dinucleotide frequencies in +-500 bp);
  * variant "lr": logistic regression (chr18 only) scores x as a lncRNA TSS (positives: grid points within 100 bp of
    any v50 lncRNA transcript start on that strand); variant "cgi": CpG-island rule (GC >= 0.5, o/e >= 0.6 in +-500 bp,
    score = CpG count), both strands;
  * peaks = local maxima, greedy non-maximum suppression within R per strand; optional Tiberius filter (drop starts
    inside a same-strand Tiberius gene, stop spans at the next same-strand Tiberius gene);
  * span = TSS + fixed length L downstream (length prior), stopped before the next selected same-strand start.
  The top-N starts (N = our own lncRNA call count on that chromosome, 3-run mean) are the predictions.

Evaluation (same rulers as B1): S50 = one-to-one reciprocal span overlap >= 0.5 (b1_strict); L1s = same-strand overlap
of the predicted span with any truth exon (span-level L1, applied identically to our predictions converted to spans).

  simple_det.py dev                         chr18: cross-fitted (halves at 40 Mb) grid search -> outputs/M1/SIMPLE/DEV.json
                                            + RULES.json (final model fitted on all chr18)
  simple_det.py read --chrom chr20|chr9     once per chromosome, needs the frozen plan -> outputs/M1/SIMPLE/<chrom>/
"""
import argparse
import hashlib
import json
import os
import re
import sys
import time
from collections import defaultdict

import numpy as np
from scipy.ndimage import maximum_filter1d
from scipy.optimize import minimize

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
for p in ("src", "scripts/phase1", "scripts/phase1/v20", "scripts/phase1/b1", "scripts/phase1/s3"):
    sys.path.insert(0, os.path.join(ROOT, p))
import p15_score as P  # noqa: E402
from ncrna_annot.eval import ruler  # noqa: E402

OUT = f"{ROOT}/outputs/M1/SIMPLE"
TIB = f"{ROOT}/outputs/PHASE1/BASELINES/tiberius"
BIN = 50
SPLIT = 40_000_000
N_BOOT = 2000
GRID_R = (1000, 2000, 5000)
GRID_L = (2000, 5000, 10000, 20000, 40000, 80000)
DINUC = [a + b for a in "ACGT" for b in "ACGT"]
COMP = str.maketrans("ACGT", "TGCA")
# chromosome -> truth, FASTA, our reads (model, dir, runs, tau source)
TGT = {
    "chr18": dict(truth="outputs/PHASE1/DIAG_L2/chr18_v50.gtf", ours=[
        ("v1.9", "outputs/PHASE1/P15b/chr18", ("g1b-42-ke", "g1b-123-ke", "g1b-456-ke")),
        ("v2.0", "outputs/PHASE1/V20/read/chr18", ("v20-g1b-42-ke", "v20-g1b-123-ke", "v20-g1b-456-ke"))]),
    "chr20": dict(truth="data/human/gencode.v50.annotation.gtf.gz", ours=[
        ("v1.9", "outputs/PHASE1/P15/chr20", ("g1b-42-ke", "g1b-123-ke", "g1b-456-ke"))]),
    "chr9": dict(truth="data/human/gencode.v50.annotation.gtf.gz", ours=[
        ("v2.0", "outputs/PHASE1/V20/read/chr9", ("v20-g1b-42-ke", "v20-g1b-123-ke", "v20-g1b-456-ke"))]),
}
FEATS = (["gc", "oe", "cpgd", "rep", "gc5k", "oe2k", "gc2", "oe2", "cpgd2",
          "gcskew_dn", "gcskew_up", "atskew_dn", "atskew_up", "cpg_dn_up", "u1_dn", "pas_dn"]
         + [f"di_{d}" for d in DINUC])


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


# ------------------------------------------------------------------ sequence -> per-bin counts
def load_seq(chrom):
    with open(f"{ROOT}/data/human/{chrom}.fa", "rb") as fh:
        raw = b"".join(line.strip() for line in fh if not line.startswith(b">"))
    arr = np.frombuffer(raw, np.uint8)
    lower = (arr >= 97) & (arr <= 122)
    lut = np.full(256, 4, np.uint8)
    for i, c in enumerate(b"ACGT"):
        lut[c] = i
        lut[c + 32] = i
    return raw.upper(), lut[arr], lower


def bin_counts(chrom):
    up, codes, lower = load_seq(chrom)
    L = len(codes)
    nb = -(-L // BIN)
    pad = nb * BIN - L
    c = np.concatenate([codes, np.full(pad, 4, np.uint8)]).reshape(nb, BIN)
    lw = np.concatenate([lower, np.zeros(pad, bool)]).reshape(nb, BIN)
    cnt = {b: (c == i).sum(1).astype(np.int32) for i, b in enumerate("ACGTN")}
    cnt["low"] = lw.sum(1).astype(np.int32)
    d = codes[:-1].astype(np.int16) * 5 + codes[1:]
    d = np.concatenate([d, np.full(pad + 1, 24, np.int16)]).reshape(nb, BIN)
    for x in DINUC:
        cnt[x] = (d == "ACGT".index(x[0]) * 5 + "ACGT".index(x[1])).sum(1).astype(np.int32)
    for name, motifs in (("u1p", ("GTAAGT", "GTGAGT")), ("u1m", ("ACTTAC", "ACTCAC")),
                         ("pasp", ("AATAAA", "ATTAAA")), ("pasm", ("TTTATT", "TTTAAT"))):
        pos = [m.start() for mo in motifs for m in re.finditer(f"(?={mo})".encode(), up)]
        cnt[name] = np.bincount(np.array(pos, np.int64) // BIN, minlength=nb).astype(np.int32)[:nb]
    return L, nb, cnt


def features(chrom):
    """Feature matrix for grid boundaries k = 0..nb (x = k * BIN), strands '+' (rows 0..K-1) and '-' (K..2K-1)."""
    L, nb, cnt = bin_counts(chrom)
    cs = {k: np.concatenate([[0], np.cumsum(v, dtype=np.int64)]) for k, v in cnt.items()}
    K = nb + 1
    k = np.arange(K)

    def S(name, a, b):  # sum of bins [k+a, k+b), clipped to the chromosome
        lo = np.clip(k + a, 0, nb)
        hi = np.clip(k + b, 0, nb)
        return (cs[name][hi] - cs[name][lo]).astype(np.float64)

    def comp(a, b):
        A, C, G, T, N = (S(x, a, b) for x in "ACGTN")
        n = A + C + G + T
        return A, C, G, T, N, n

    def div(x, y):
        return np.divide(x, y, out=np.zeros_like(x), where=y > 0)

    A, C, G, T, N, n = comp(-10, 10)
    CG = S("CG", -10, 10)
    base = {"gc": div(C + G, n), "oe": div(CG * n, C * G), "cpgd": div(CG, n), "rep": div(S("low", -10, 10), n + N),
            "nfrac": div(N, n + N)}
    A5, C5, G5, T5, _, n5 = comp(-100, 100)
    base["gc5k"] = div(C5 + G5, n5)
    _, C2, G2, _, _, n2 = comp(-40, 40)
    base["oe2k"] = div(S("CG", -40, 40) * n2, C2 * G2)
    base["cpg_count"] = CG
    di = {x: S(x, -10, 10) for x in DINUC}
    ndi = sum(di.values())
    rows = []
    for strand in "+-":
        dn, upw = ((0, 20), (-20, 0)) if strand == "+" else ((-20, 0), (0, 20))
        Ad, Cd, Gd, Td, _, nd = comp(*dn)
        Au, Cu, Gu, Tu, _, nu = comp(*upw)
        if strand == "-":   # sense strand = complement
            Ad, Cd, Gd, Td = Td, Gd, Cd, Ad
            Au, Cu, Gu, Tu = Tu, Gu, Cu, Au
        f = dict(base)
        f["gc2"], f["oe2"], f["cpgd2"] = base["gc"] ** 2, np.minimum(base["oe"], 3) ** 2, base["cpgd"] ** 2
        f["gcskew_dn"], f["gcskew_up"] = div(Gd - Cd, Gd + Cd), div(Gu - Cu, Gu + Cu)
        f["atskew_dn"], f["atskew_up"] = div(Ad - Td, Ad + Td), div(Au - Tu, Au + Tu)
        f["cpg_dn_up"] = div(S("CG", *dn), nd) - div(S("CG", *upw), nu)
        f["u1_dn"] = div(S("u1p" if strand == "+" else "u1m", *dn), nd) * 1000
        f["pas_dn"] = div(S("pasp" if strand == "+" else "pasm", *dn), nd) * 1000
        for x in DINUC:
            sx = x if strand == "+" else x[::-1].translate(COMP)
            f[f"di_{x}"] = div(di[sx], ndi)
        rows.append(f)
    X = np.vstack([np.column_stack([r[name] for name in FEATS]) for r in rows]).astype(np.float32)
    aux = {"nfrac": np.concatenate([base["nfrac"]] * 2), "gc": np.concatenate([base["gc"]] * 2),
           "oe": np.concatenate([base["oe"]] * 2), "cpg_count": np.concatenate([base["cpg_count"]] * 2)}
    return {"L": L, "K": K, "X": X, "aux": aux}


# ------------------------------------------------------------------ labels / truth
def lnc_tss(truth):
    out = {"+": set(), "-": set()}
    for g in truth.values():
        if g["cls"] != "lnc":
            continue
        for ex in g["tx"].values():
            ex = sorted(ex)
            out[g["strand"]].add(ex[0][0] - 1 if g["strand"] == "+" else ex[-1][1])   # 0-based boundary
    return out


def labels(K, tss, pos_bins=2, neg_excl=20):
    """y: 1 within pos_bins of a same-strand lnc TSS, 0 if farther than neg_excl bins, -1 in between."""
    y = np.zeros(2 * K, np.int8)
    for si, s in enumerate("+-"):
        near = np.zeros(K, np.int16)
        for x in tss[s]:
            c = int(round(x / BIN))
            near[max(c - neg_excl, 0):c + neg_excl + 1] = np.maximum(near[max(c - neg_excl, 0):c + neg_excl + 1], 1)
            near[max(c - pos_bins, 0):c + pos_bins + 1] = 2
        y[si * K:(si + 1) * K] = np.where(near == 2, 1, np.where(near == 1, -1, 0))
    return y


def fit_logit(X, y, C=1.0):
    n, p = X.shape
    Xa = np.hstack([np.ones((n, 1)), X])

    def f(w):
        z = Xa @ w
        ll = np.logaddexp(0, z) - y * z
        g = Xa.T @ (1 / (1 + np.exp(-z)) - y)
        g[1:] += w[1:] / C
        return ll.sum() + 0.5 / C * (w[1:] @ w[1:]), g
    return minimize(f, np.zeros(p + 1), jac=True, method="L-BFGS-B", options={"maxiter": 2000}).x


def train_lr(F, y, mask, seed=0, n_neg=300_000):
    rng = np.random.default_rng(seed)
    valid = (F["aux"]["nfrac"] < 0.5) & mask
    pos = np.flatnonzero(valid & (y == 1))
    neg = np.flatnonzero(valid & (y == 0))
    neg = rng.choice(neg, min(n_neg, len(neg)), replace=False)
    idx = np.concatenate([pos, neg])
    X = F["X"][idx].astype(np.float64)
    mu, sd = X.mean(0), X.std(0)
    sd[sd == 0] = 1
    w = fit_logit((X - mu) / sd, (y[idx] == 1).astype(float))
    return {"features": list(FEATS), "mu": mu.tolist(), "sd": sd.tolist(), "w": w.tolist(),
            "n_pos": int(len(pos)), "n_neg": int(len(neg))}


def score_lr(F, model):
    X = (F["X"].astype(np.float64) - np.array(model["mu"])) / np.array(model["sd"])
    w = np.array(model["w"])
    s = w[0] + X @ w[1:]
    s[F["aux"]["nfrac"] >= 0.5] = -np.inf
    return s


def score_cgi(F):
    a = F["aux"]
    s = np.where((a["gc"] >= 0.5) & (a["oe"] >= 0.6) & (a["nfrac"] < 0.5), a["cpg_count"], -np.inf)
    return s.astype(np.float64)


# ------------------------------------------------------------------ peaks -> predictions
def read_tiberius(chrom):
    spans = {}
    for line in open(f"{TIB}/{chrom}.gtf"):
        f = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(f) < 9 or f[0] != chrom or f[2] not in ("CDS", "exon"):
            continue
        gid = f[8].split('gene_id "', 1)[1].split('"', 1)[0]
        s, e = int(f[3]), int(f[4])
        if gid in spans:
            spans[gid] = (f[6], min(spans[gid][1], s), max(spans[gid][2], e))
        else:
            spans[gid] = (f[6], s, e)
    out = {"+": [], "-": []}
    for st, s, e in spans.values():
        out[st].append((s, e))
    return {k: np.array(sorted(v), np.int64).reshape(-1, 2) for k, v in out.items()}


def peaks(score, K, R, tib=None, max_keep=40_000):
    """Greedy NMS per strand over local maxima -> list of (score, strand, x) sorted by score desc."""
    rb = max(1, R // BIN)
    out = []
    for si, st in enumerate("+-"):
        s = score[si * K:(si + 1) * K].copy()
        if tib is not None and len(tib[st]):          # drop starts inside a same-strand Tiberius gene
            ins = np.zeros(K, bool)
            for a, b in tib[st]:
                ins[max(a // BIN, 0):min(b // BIN + 1, K)] = True
            s[ins] = -np.inf
        fin = np.isfinite(s)
        if not fin.any():
            continue
        mf = maximum_filter1d(np.where(fin, s, -1e300), size=2 * rb + 1, mode="constant", cval=-1e300)
        cand = np.flatnonzero(fin & (s >= mf))
        cand = cand[np.argsort(-s[cand], kind="stable")]
        taken = np.zeros(K, bool)
        kept = 0
        for c in cand:
            if taken[max(c - rb, 0):c + rb + 1].any():
                continue
            taken[c] = True
            out.append((float(s[c]), st, int(c) * BIN))
            kept += 1
            if kept >= max_keep:
                break
    out.sort(key=lambda t: -t[0])
    return out


def spans_for(sel, L_prior, chrom_len, tib=None):
    """sel: [(score, strand, x)] -> predictions with 1-based inclusive spans."""
    by = {"+": sorted(x for _, s, x in sel if s == "+"), "-": sorted(x for _, s, x in sel if s == "-")}
    preds = []
    for i, (sc, st, x) in enumerate(sel):
        if st == "+":
            s, e = x + 1, min(x + L_prior, chrom_len)
            j = np.searchsorted(by["+"], x, side="right")
            if j < len(by["+"]):
                e = min(e, by["+"][j])
            if tib is not None and len(tib["+"]):
                j = np.searchsorted(tib["+"][:, 0], x + 1, side="right")
                if j < len(tib["+"]):
                    e = min(e, tib["+"][j, 0] - 1)
        else:
            s, e = max(x - L_prior + 1, 1), x
            j = np.searchsorted(by["-"], x, side="left") - 1
            if j >= 0:
                s = max(s, by["-"][j] + 1)
            if tib is not None and len(tib["-"]):
                ends = np.sort(tib["-"][:, 1])
                j = np.searchsorted(ends, x, side="left") - 1
                if j >= 0:
                    s = max(s, int(ends[j]) + 1)
        if e - s + 1 >= 200:
            preds.append({"id": f"sd{i}", "cls": "lnc", "strand": st, "exons": [(int(s), int(e))], "score": sc})
    return preds


def top_n(pk, n, L_prior, chrom_len, tib=None):
    """Exactly n predictions (spans < 200 bp are dropped by spans_for, so take more starts until n survive)."""
    m = n
    while True:
        preds = spans_for(pk[:m], L_prior, chrom_len, tib)
        if len(preds) >= n or m >= len(pk):
            return sorted(preds, key=lambda p: -p["score"])[:n]
        m += n - len(preds)


def predict(score, K, chrom_len, cfg, n=None, thr=None, tib=None):
    t = tib if cfg["tib"] else None
    pk = peaks(score, K, cfg["R"], t)
    if thr is not None:
        return spans_for([p for p in pk if p[0] >= thr], cfg["L"], chrom_len, t)
    return top_n(pk, n, cfg["L"], chrom_len, t)


# ------------------------------------------------------------------ rulers (per-gene hit vectors)
def truth_spans(truth):
    out = []
    for gid, g in truth.items():
        if g["cls"] == "lnc":
            s = min(e[0] for t in g["tx"].values() for e in t)
            e = max(e[1] for t in g["tx"].values() for e in t)
            out.append((gid, g["strand"], s, e))
    return sorted(out, key=lambda t: t[2])


def s50(ts, preds, th=0.5):
    import b1_strict as B
    starts = [t[2] for t in ts]
    maxlen = max(t[3] - t[2] for t in ts)
    adj = []
    for p in preds:
        st, s, e = p["strand"], p["exons"][0][0], p["exons"][-1][1]
        lo, hi = np.searchsorted(starts, s - maxlen - 1), np.searchsorted(starts, e, side="right")
        a = []
        for i in range(lo, hi):
            _, gst, gs, ge = ts[i]
            if gst != st:
                continue
            ov = min(e, ge) - max(s, gs) + 1
            if ov > 0 and ov / (ge - gs + 1) >= th and ov / (e - s + 1) >= th:
                a.append(i)
        adj.append(a)
    ml = B.matching(adj, len(ts))
    return {ts[i][0] for i in ml if i != -1}


def l1s(idx, truth, preds):
    hit, nhit = set(), 0
    for p in preds:
        h = {x for x in idx.hits(p["strand"], [(p["exons"][0][0], p["exons"][-1][1])]) if truth[x]["cls"] == "lnc"}
        hit |= h
        nhit += bool(h)
    return hit, nhit


def metrics(truth, ts, idx, preds):
    m50 = s50(ts, preds)
    h1, n1 = l1s(idx, truth, preds)
    nt, n = len(ts), len(preds)
    return {"n_pred": n, "S50_sn": len(m50) / nt, "S50_pr": len(m50) / n if n else None,
            "L1s_sn": len(h1) / nt, "L1s_pr": n1 / n if n else None}, m50, h1


# ------------------------------------------------------------------ our predictions
def load_ours(chrom):
    """{model: {run: [genes with type 'l' and (calibrated) conf]}}, tau per model."""
    out = {}
    for model, rd, runs in TGT[chrom]["ours"]:
        if model == "v1.9":
            tau, cal = json.load(open(f"{ROOT}/configs/p15_v19.json"))["decode"]["tau"], None
        else:
            import v20_calib as Cc
            c20 = json.load(open(f"{ROOT}/configs/v20.json"))
            cal, tau = c20["calibration"], c20["calibration"]["tau_v20"]
            assert cal["keep"] and P.sha256(f"{ROOT}/{cal['hexamer_table']}") == cal["hexamer_table_sha256"]
            tab = json.load(open(f"{ROOT}/{cal['hexamer_table']}"))
            codes = np.load(f"{ROOT}/data/p1/{chrom}_codes.npy", mmap_mode="r")
        rr = {}
        for run in runs:
            genes = P.load_genes(f"{ROOT}/{rd}/{run}/genes.json")
            if cal:
                post = np.load(f"{ROOT}/{rd}/{run}/post.npy", mmap_mode="r")
                rows = Cc.gene_rows(genes, post, codes, tab)
                X = np.array([[r[f] for f in cal["model"]["features"]] for r in rows], float)
                sc = dict(zip((r["k"] for r in rows), Cc.apply(cal["model"], X)))
                for k, g in enumerate(genes):
                    if k in sc:
                        g["conf"] = float(sc[k])
            rr[run] = [{"id": f"o{k}", "cls": "lnc", "strand": g["strand"], "conf": g["conf"],
                        "exons": [tuple(x) for x in g["exons"]]} for k, g in enumerate(genes) if g["type"] == "l"]
        out[model] = {"tau": tau, "runs": rr}
    return out


def as_span(p):
    return dict(p, exons=[(p["exons"][0][0], p["exons"][-1][1])])


def exon_l1(idx, truth, preds):
    hit, nhit = set(), 0
    for p in preds:
        h = {x for x in idx.hits(p["strand"], p["exons"]) if truth[x]["cls"] == "lnc"}
        hit |= h
        nhit += bool(h)
    return hit, nhit


# ------------------------------------------------------------------ dev (chr18)
def dev(a):
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)
    chrom = "chr18"
    truth = ruler.read_gtf(f"{ROOT}/{TGT[chrom]['truth']}", chrom)
    ts, idx = truth_spans(truth), ruler.ExonIndex(truth)
    F = features(chrom)
    K, Lc = F["K"], F["L"]
    print("features", F["X"].shape, round(time.time() - t0), "s", flush=True)
    y = labels(K, lnc_tss(truth))
    xpos = np.concatenate([np.arange(K) * BIN] * 2)
    half = xpos < SPLIT
    # cross-fitted LR scores: model trained on one half scores the other half
    s_cf = np.full(2 * K, -np.inf)
    for train_mask in (half, ~half):
        m = train_lr(F, y, train_mask)
        s = score_lr(F, m)
        s_cf[~train_mask] = s[~train_mask]
    full = train_lr(F, y, np.ones(2 * K, bool))
    s_full = score_lr(F, full)
    s_cgi = score_cgi(F)
    tib = read_tiberius(chrom)
    ours = load_ours(chrom)
    ref = {}
    for model, o in ours.items():
        per = {}
        for run, genes in o["runs"].items():
            fin = [g for g in genes if g["conf"] >= o["tau"]]
            mt, _, _ = metrics(truth, ts, idx, [as_span(g) for g in fin])
            he, ne = exon_l1(idx, truth, fin)
            mt["L1_exon_sn"], mt["L1_exon_pr"] = len(he) / len(ts), ne / len(fin) if fin else None
            per[run] = mt
        ref[model] = {"tau": o["tau"], "runs": per,
                      "mean": {k: float(np.mean([p[k] for p in per.values()])) for k in per[next(iter(per))]}}
    n18 = int(round(ref["v1.9"]["mean"]["n_pred"]))
    print("ours chr18", json.dumps({m: r["mean"] for m, r in ref.items()}), flush=True)
    grid = []
    for variant, sc in (("lr", s_cf), ("cgi", s_cgi)):
        for R in GRID_R:
            for tb in (0, 1):
                pk = peaks(sc, K, R, tib if tb else None)
                for Lp in GRID_L:
                    cfg = {"variant": variant, "R": R, "L": Lp, "tib": tb}
                    preds = top_n(pk, n18, Lp, Lc, tib if tb else None)
                    mt, _, _ = metrics(truth, ts, idx, preds)
                    grid.append(dict(cfg, **mt))
                    print(json.dumps(grid[-1]), flush=True)
    sel, rules_sel = {}, {}
    for metric, other_m in (("S50", "L1s"), ("L1s", "S50")):
        best = max(grid, key=lambda g: (g[f"{metric}_sn"], g[f"{other_m}_sn"]))
        cfg = {k: best[k] for k in ("variant", "R", "L", "tib")}
        s_sel = s_full if cfg["variant"] == "lr" else s_cgi
        t = tib if cfg["tib"] else None
        pk = peaks(s_sel, K, cfg["R"], t)
        preds = top_n(pk, n18, cfg["L"], Lc, t)
        thr = min(p["score"] for p in preds)
        insample, _, _ = metrics(truth, ts, idx, preds)
        other = max((g for g in grid if g["variant"] != cfg["variant"]), key=lambda g: (g[f"{metric}_sn"], g[f"{other_m}_sn"]))
        sel[metric] = {"config": cfg, "cv": best, "fullfit_insample": insample, "best_other_variant": other,
                       "ratio_cv_vs_ours": {m: best[f"{metric}_sn"] / ref[m]["mean"][f"{metric}_sn"] for m in ref}}
        rules_sel[metric] = {"config": cfg, "threshold_matching_chr18_n": thr}
        print("selected", metric, json.dumps(best), "other", json.dumps(other), flush=True)
    rules = {"plan": "docs/m1/PLAN_M1_SIMPLE.md", "bin": BIN, "features": list(FEATS), "selected": rules_sel,
             "lr_model": full, "n_chr18": n18, "code_sha256": sha256(os.path.abspath(__file__)),
             "tiberius_chr18_sha256": sha256(f"{TIB}/chr18.gtf")}
    json.dump(rules, open(f"{OUT}/RULES.json", "w"), indent=1)
    res = {"chrom": chrom, "n_truth_lnc": len(ts), "ours": ref, "n_matched": n18, "grid": grid, "selected": sel,
           "lr_cv_models_note": "cross-fitted halves split at 40 Mb", "seconds": round(time.time() - t0, 1)}
    json.dump(P.norm(res), open(f"{OUT}/DEV.json", "w"), indent=1)
    print("seconds", res["seconds"], flush=True)


# ------------------------------------------------------------------ read (chr20 / chr9), once
def cpg_promoter(truth, chrom):
    import s3_features as S3
    seq = S3.read_fasta(f"{ROOT}/data/human/{chrom}.fa")
    out = {}
    for gid, g in truth.items():
        if g["cls"] != "lnc":
            continue
        s0 = min(e[0] for t in g["tx"].values() for e in t)
        e0 = max(e[1] for t in g["tx"].values() for e in t)
        t = s0 - 1 if g["strand"] == "+" else e0 - 1
        out[gid] = S3.cpg_oe(seq[max(t - 500, 0):t + 501])
    return out


def boot(vec_det, vecs_ours, rng_seed=0):
    n = len(vec_det)
    rng = np.random.default_rng(rng_seed)
    od = np.mean(vecs_ours, 0)
    d, r = [], []
    for _ in range(N_BOOT):
        b = rng.integers(0, n, n)
        sd, so = vec_det[b].mean(), od[b].mean()
        d.append(so - sd)
        r.append(sd / so if so > 0 else np.nan)
    return {"ours_minus_det": float(od.mean() - vec_det.mean()), "diff_ci95": np.nanpercentile(d, [2.5, 97.5]).tolist(),
            "ratio": float(vec_det.mean() / od.mean()) if od.mean() > 0 else None,
            "ratio_ci95": np.nanpercentile(r, [2.5, 97.5]).tolist()}


def read(a):
    chrom = a.chrom
    t0 = time.time()
    plan = f"{ROOT}/docs/m1/PLAN_M1_SIMPLE.md"
    frozen = f"{OUT}/FROZEN"   # written by hand after the author approves the plan: '<rules sha256>'
    rules_path = f"{OUT}/RULES.json"
    if a.dry:   # code check on the development chromosome only (v1.9 as 'ours'); never a test chromosome
        assert chrom == "chr18", "--dry is for chr18 only"
        od = f"{OUT}/dryrun_chr18"
    else:
        assert chrom in ("chr20", "chr9")
        if not (os.path.exists(plan) and os.path.exists(frozen)):
            sys.exit("plan not frozen: refusing to read a test chromosome")
        if open(frozen).read().split()[0] != sha256(rules_path):
            sys.exit("RULES.json changed after freezing")
        od = f"{OUT}/{chrom}"
        if os.path.exists(f"{od}/SCORE.json"):
            sys.exit(f"{chrom} already read")
    os.makedirs(od, exist_ok=True)
    rules = json.load(open(rules_path))
    truth = ruler.read_gtf(f"{ROOT}/{TGT[chrom]['truth']}", chrom)
    ts, idx = truth_spans(truth), ruler.ExonIndex(truth)
    gids = [t[0] for t in ts]
    F = features(chrom)
    tib = read_tiberius(chrom)
    ours = load_ours(chrom)
    if a.dry:
        ours = {"v1.9": ours["v1.9"]}
    (model, o), = ours.items()
    per, vec = {}, {"S50": [], "L1s": [], "L1e": []}
    for run, genes in o["runs"].items():
        fin = [g for g in genes if g["conf"] >= o["tau"]]
        mt, m50, h1 = metrics(truth, ts, idx, [as_span(g) for g in fin])
        he, ne = exon_l1(idx, truth, fin)
        mt["L1_exon_sn"], mt["L1_exon_pr"] = len(he) / len(ts), ne / len(fin) if fin else None
        per[run] = mt
        vec["S50"].append(np.array([g in m50 for g in gids], float))
        vec["L1s"].append(np.array([g in h1 for g in gids], float))
        vec["L1e"].append(np.array([g in he for g in gids], float))
    n_match = int(round(np.mean([p["n_pred"] for p in per.values()])))
    cpg = cpg_promoter(truth, chrom)
    cv = np.array([cpg[g] for g in gids])
    qs = np.quantile(cv, [0.25, 0.5, 0.75])
    q = np.searchsorted(qs, cv, side="right")
    om = {k: np.mean(vec[k], 0) for k in vec}
    det_res = {}
    for metric in ("S50", "L1s"):
        cfg = rules["selected"][metric]["config"]
        sc = score_lr(F, rules["lr_model"]) if cfg["variant"] == "lr" else score_cgi(F)
        det = predict(sc, F["K"], F["L"], cfg, n=n_match, tib=tib)
        dm, d50, d1 = metrics(truth, ts, idx, det)
        dmt, _, _ = metrics(truth, ts, idx, predict(sc, F["K"], F["L"], cfg,
                                                    thr=rules["selected"][metric]["threshold_matching_chr18_n"], tib=tib))
        vd = {"S50": np.array([g in d50 for g in gids], float), "L1s": np.array([g in d1 for g in gids], float)}
        t = tib if cfg["tib"] else None
        pk = peaks(sc, F["K"], cfg["R"], t)
        curve = [metrics(truth, ts, idx, top_n(pk, n, cfg["L"], F["L"], t))[0] for n in range(100, 3001, 100)]
        union = np.max(vec[metric], 0) > 0
        det_res[metric] = {
            "config": cfg, "matched": dm, "threshold_from_chr18": dmt,
            "compare": boot(vd[metric], vec[metric]),
            "compare_other_metric": boot(vd["L1s" if metric == "S50" else "S50"], vec["L1s" if metric == "S50" else "S50"]),
            **({"compare_vs_ours_exon_L1": boot(vd["L1s"], vec["L1e"])} if metric == "L1s" else {}),
            "overlap": {"ours_any_run_hits": int(union.sum()), "det_hits": int(vd[metric].sum()),
                        "both": int((union & (vd[metric] > 0)).sum()),
                        "frac_ours_hits_also_det": float((union & (vd[metric] > 0)).sum() / max(union.sum(), 1))},
            "compare_lowcpg_q1q2": boot(vd[metric][q <= 1], [v[q <= 1] for v in vec[metric]]),
            "cpg_strata": {f"Q{i + 1}": {"n": int((q == i).sum()), "det": float(vd[metric][q == i].mean()),
                                         "ours": float(om[metric][q == i].mean())} for i in range(4)},
            "curve": curve}
        json.dump([dict(p, exons=[list(x) for x in p["exons"]]) for p in det], open(f"{od}/pred_{metric}.json", "w"))
    curve_ours = {}
    for run, genes in o["runs"].items():
        cv_ = []
        for t in [round(0.05 * i, 2) for i in range(1, 20)]:
            fin = [g for g in genes if g["conf"] >= t]
            if fin:
                cv_.append(dict(metrics(truth, ts, idx, [as_span(g) for g in fin])[0], tau=t))
        curve_ours[run] = cv_
    ratios = {m: det_res[m]["compare"]["ratio"] for m in ("S50", "L1s")}
    res = {"chrom": chrom, "plan": "docs/m1/PLAN_M1_SIMPLE.md", "rules_sha256": sha256(rules_path),
           "n_truth_lnc": len(ts), "ours_model": model, "ours_tau": o["tau"], "ours": per,
           "ours_mean": {k: float(np.mean([p[k] for p in per.values()])) for k in per[next(iter(per))]},
           "n_matched": n_match, "det": det_res, "ratios": ratios, "cpg_quartile_cuts": qs.tolist(),
           "curve_ours": curve_ours, "seconds": round(time.time() - t0, 1)}
    json.dump(P.norm(res), open(f"{od}/SCORE.json", "w"), indent=1)
    print(json.dumps({"chrom": chrom, "n": n_match, "ours_mean": res["ours_mean"], "ratios": ratios,
                      "S50": det_res["S50"]["compare"], "L1s": det_res["L1s"]["compare"]}), flush=True)


def verdict(a):
    """Pre-registered rule (PLAN_M1_SIMPLE.md section 4) over the two test chromosomes."""
    sc = {c: json.load(open(f"{OUT}/{c}/SCORE.json")) for c in ("chr20", "chr9")}
    r = {c: s["ratios"] for c, s in sc.items()}
    retract = all(r[c]["S50"] >= 0.8 and r[c]["L1s"] >= 0.8 for c in sc)
    keep = all(r[c]["L1s"] < 0.5 and sc[c]["det"]["S50"]["compare"]["diff_ci95"][0] > 0 for c in sc)
    v = "RETRACT" if retract else ("KEEP" if keep else "GREY")
    low = {m: all(sc[c]["det"][m]["compare_lowcpg_q1q2"]["diff_ci95"][0] > 0 for c in sc) for m in ("S50", "L1s")}
    out = {"verdict": v, "ratios": r, "S50_diff_ci95": {c: sc[c]["det"]["S50"]["compare"]["diff_ci95"] for c in sc},
           "secondary_lowcpg_ours_better_both_chroms": low,
           "lowcpg": {c: {m: sc[c]["det"][m]["compare_lowcpg_q1q2"] for m in ("S50", "L1s")} for c in sc}}
    json.dump(out, open(f"{OUT}/VERDICT.json", "w"), indent=1)
    print(json.dumps(out))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("dev", "read", "verdict"))
    ap.add_argument("--chrom", choices=("chr18", "chr20", "chr9"))
    ap.add_argument("--dry", action="store_true", help="read path check on chr18 (development chromosome)")
    a = ap.parse_args()
    {"dev": dev, "read": read, "verdict": verdict}[a.cmd](a)


if __name__ == "__main__":
    main()
