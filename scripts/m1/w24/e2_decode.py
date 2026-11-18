#!/usr/bin/env python3
"""M1 weeks 2-4, E2 (evidence in the model; docs/m1/PLAN_M1_W24.md section 3): exact posteriors from the cached head
outputs with an evidence term added to the G1B gene entry / exit edges, then v2.0's MEA decoding, gene confidence,
calibration and tau.  CPU only (tiles in parallel).

Evidence -> edges (g1.boundary_edges; row t-1 of vals = transition into base t, 0-based within the tile):
  ssCAGE TSS on + at 1-based X : IG -> kE+ into base X-1        (gene starts at X)
  ssCAGE TSS on - at X         : kE- -> IG into base X           (gene's right end is X)
  QuantSeq site on + at X      : kE+ -> IG into base X           (gene's last base is X)
  QuantSeq site on - at X      : IG -> kE- into base X-1         (gene's left end is X)
for both types k (m, l); +beta on every base within +-tol of X.  Condition 'dna' adds nothing and must reproduce v2.0.
  e2_decode.py --target chr18 --run <run> --cond dna|oracle|ev --frac 1.0 --seed 0 --beta 4 --tol 20 [--check]
    -> outputs/M1/W24/e2/<target>/<run>/<tag>.json (scores) + genes
"""
import argparse
import hashlib
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
for d in ("src", "scripts/phase1", "scripts/m1", "scripts/m1/w24"):
    sys.path.insert(0, os.path.join(ROOT, d))
os.environ.setdefault("NUMBA_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

_G = {}


def _init(cache, codes_path, coef, T, add):
    import torch
    from ncrna_annot.grammar import g1
    torch.set_num_threads(1)
    _G.update(em=np.load(f"{cache}/emit.npy", mmap_mode="r"), bd=np.load(f"{cache}/bnd.npy", mmap_mode="r"),
              codes=np.load(codes_path, mmap_mode="r"), T=T, add=add, g1=g1, torch=torch)
    c = g1.BoundaryCoef()
    with torch.no_grad():
        for k, v in coef.items():
            c.c[k].copy_(torch.tensor(v))
    _G["coef"] = c.double()


def _tile(args):
    k, s, lo, hi = args
    from ncrna_annot.features import generanno as G
    g1, torch = _G["g1"], _G["torch"]
    e = np.asarray(_G["em"][k], np.float64)
    b = torch.as_tensor(np.asarray(_G["bd"][k]), dtype=torch.float32)
    gates = g1.gate_mask(np.asarray(_G["codes"][s:s + G.WIN]))
    with torch.no_grad():
        frm, to, vals = g1.boundary_edges(b.double(), _G["coef"])
    vals = vals.numpy().copy()
    for col, rows, beta in _G["add"].get(k, []):
        vals[rows, col] += beta
    m = g1.marginals_np(e, _G["T"], gates, edges=(frm, to, vals))[0]
    return k, m[lo:hi].astype(np.float16)


def edge_cols():
    from ncrna_annot.grammar import g1
    cols = {}
    i = 0
    for k in g1.TYPES:
        for name in ("in+", "out+", "in-", "out-"):
            cols.setdefault(name, []).append(i)
            i += 1
    return cols


def evidence_points(target, cond, frac, seed):
    """-> list of (kind, strand, X 1-based) with kind 'tss' / 'pas'."""
    import w24_lib as L
    pts = []
    if cond == "ev":
        ev = L.evidence(target, frac, seed)
        for st in "+-":
            pts += [("tss", st, int(x)) for x in ev["tss"][st][:, 0]] if len(ev["tss"][st]) else []
            pts += [("pas", st, int(x)) for x in ev["pas"][st]]
    elif cond == "oracle":
        d = json.load(open(f"{L.W}/truth/{target}.json"))
        for loc in d["loci"]:
            for t in loc["tx"]:
                s, e = t["exons"][0][0], t["exons"][-1][1]
                if loc["strand"] == "+":
                    pts += [("tss", "+", s), ("pas", "+", e)]
                else:
                    pts += [("tss", "-", e), ("pas", "-", s)]
        pts = sorted(set(pts))
    return pts


def build_add(pts, starts, win, beta, tol):
    """Per tile: list of (edge column, rows array, beta)."""
    cols = edge_cols()
    which = {("tss", "+"): ("in+", -1), ("tss", "-"): ("out-", 0), ("pas", "+"): ("out+", 0), ("pas", "-"): ("in-", -1)}
    starts = np.asarray(starts)
    add = {}
    for kind, st, X in pts:
        name, shift = which[(kind, st)]
        x0 = X - 1                                   # 0-based base of the evidence
        for k in np.flatnonzero((starts <= x0 + tol) & (x0 - tol < starts + win)):
            s = int(starts[k])
            rows = np.arange(x0 - tol, x0 + tol + 1) - s + shift   # vals row = (base into which we transition) - 1
            rows = rows[(rows >= 0) & (rows < win - 1)]
            if len(rows):
                for c in cols[name]:
                    add.setdefault(int(k), []).append((c, rows, beta))
    return add


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", default="chr18")
    ap.add_argument("--run", required=True)
    ap.add_argument("--cond", default="dna", choices=["dna", "oracle", "ev"])
    ap.add_argument("--frac", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--beta", type=float, default=0.0)
    ap.add_argument("--tol", type=int, default=0)
    ap.add_argument("--procs", type=int, default=24)
    ap.add_argument("--check", action="store_true", help="dna: compare post / genes with the v2.0 read")
    a = ap.parse_args()
    assert a.target == "chr18", "development is chr18 only"
    from ncrna_annot.decode import chrom
    from ncrna_annot.features import generanno as G
    from ncrna_annot.grammar import g1
    import w24_lib as L
    import simple_det as D
    cfg = json.load(open(f"{ROOT}/configs/v20.json"))
    d = cfg["decode"]
    cache = f"{L.W}/cache/{a.target}/{a.run}"
    meta = json.load(open(f"{cache}/meta.json"))
    T = np.load(f"{cache}/T.npy")
    starts, Ln = meta["starts"], meta["L"]
    codes_path = f"{ROOT}/data/p1/{a.target}_codes.npy"
    codes = np.load(codes_path, mmap_mode="r")
    tag = "dna" if a.cond == "dna" else f"{a.cond}_f{a.frac}_s{a.seed}_b{a.beta}_t{a.tol}"
    pts = evidence_points(a.target, a.cond, a.frac, a.seed) if a.cond != "dna" else []
    add = build_add(pts, starts, G.WIN, a.beta, a.tol) if pts else {}
    t0 = time.time()
    jobs = []
    for k, s in enumerate(starts):
        lo = 0 if k == 0 else G.MARGIN * G.TOK
        hi = G.WIN if k == len(starts) - 1 else G.WIN - G.MARGIN * G.TOK
        jobs.append((k, s, lo, hi))
    post = np.zeros((Ln, g1.S), np.float16)
    with Pool(a.procs, initializer=_init, initargs=(cache, codes_path, meta["coef"], T, add)) as pool:
        for k, m in pool.imap_unordered(_tile, jobs, chunksize=8):
            s, lo, hi = jobs[k][1], jobs[k][2], jobs[k][3]
            post[s + lo:s + hi] = m
    t_post = time.time() - t0
    out = f"{L.W}/e2/{a.target}/{a.run}"
    os.makedirs(out, exist_ok=True)
    res = {"run": a.run, "cond": a.cond, "frac": a.frac, "seed": a.seed, "beta": a.beta, "tol": a.tol, "n_points": len(pts),
           "post_s": round(t_post)}
    if a.check:
        ref = np.load(f"{ROOT}/outputs/PHASE1/V20/read/{a.target}/{a.run}/post.npy", mmap_mode="r")
        res["post_equal"] = bool(ref.shape == post.shape and np.array_equal(np.asarray(ref), post))
        res["post_maxdiff"] = float(np.abs(np.asarray(ref, np.float32) - post.astype(np.float32)).max())
    t1 = time.time()
    Tg = g1.G1Transitions()().detach().double().numpy()
    genes = chrom.decode_blocks(post, codes[:Ln], Tg, d["min_exon"], d["min_intron"], mode="mea", gamma=d["gamma"],
                                level=d["level"], log_gain=False)
    genes = [g for g in genes if g["type"] != "l" or sum(b - a_ + 1 for a_, b in g["exons"]) >= d["lnc_min_len"]]
    for g in genes:
        g["conf"] = chrom.gene_confidence(g, post)
    res["decode_s"] = round(time.time() - t1)
    gj = [{"type": g["type"], "strand": g["strand"], "exons": [list(x) for x in g["exons"]], "conf": g["conf"]} for g in genes]
    if a.check:
        ref_g = open(f"{ROOT}/outputs/PHASE1/V20/read/{a.target}/{a.run}/genes.json", "rb").read()
        res["genes_equal"] = hashlib.sha256(json.dumps(gj).encode()).hexdigest() == hashlib.sha256(ref_g).hexdigest()
    # calibration (v2.0) and tau
    import v20_calib as Cc
    cal = cfg["calibration"]
    tab = json.load(open(f"{ROOT}/{cal['hexamer_table']}"))
    rows = Cc.gene_rows(gj, post, codes, tab)
    X = np.array([[r[f] for f in cal["model"]["features"]] for r in rows], float)
    sc = dict(zip((r["k"] for r in rows), Cc.apply(cal["model"], X)))
    lnc = [{"strand": g["strand"], "exons": [tuple(x) for x in g["exons"]], "conf": float(sc.get(k, g["conf"]))}
           for k, g in enumerate(gj) if g["type"] == "l"]
    truth = L.load_truth(a.target)
    res["n_lnc_all"] = len(lnc)
    res["by_tau"] = {}
    for tau in (cal["tau_v20"], 0.2, 0.1, 0.0):
        p = [x for x in lnc if x["conf"] >= tau]
        r = L.score_unknown(truth, p, a.target)
        r.pop("S50_ids")
        res["by_tau"][str(tau)] = r
    json.dump(res, open(f"{out}/{tag}.json", "w"), indent=1)
    json.dump(lnc, open(f"{out}/{tag}.lnc.json", "w"))
    print(json.dumps({k: v for k, v in res.items() if k != "by_tau"}), flush=True)
    print(json.dumps(res["by_tau"][str(cal["tau_v20"])]), flush=True)


if __name__ == "__main__":
    main()
