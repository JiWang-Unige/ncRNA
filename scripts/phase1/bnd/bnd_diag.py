#!/usr/bin/env python3
"""BND (docs/phase1_results/PLAN_BND.md): why lnc gene ends are mislocated.  chr18, v1.9 main runs, CPU, existing
posteriors and predictions only.

  bnd_diag.py [--out outputs/PHASE1/BND] [--skip-s1]   -> <out>/BND.json
"""
import argparse
import json
import os
import statistics
import sys
import time

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
for d in ("scripts/phase1", "scripts/phase1/b1", "scripts/phase1/s3"):
    sys.path.insert(0, os.path.join(ROOT, d))
import b1_strict as B  # noqa: E402
import p15_score as P  # noqa: E402
import s3_features as S3  # noqa: E402
from ncrna_annot.eval import ruler  # noqa: E402
from ncrna_annot.features import generanno as G  # noqa: E402
from ncrna_annot.grammar import g1  # noqa: E402

SEAM_TOL, FLANK, SEARCH = 100, 1000, 50_000
KIND = {"lnc": "l", "mrna": "m"}


def span(g):
    return min(e[0] for t in g["tx"].values() for e in t), max(e[1] for t in g["tx"].values() for e in t)


def seam_frac(pos0):
    """fraction of 0-based positions within SEAM_TOL of a tile seam (k * STEP + MARGIN * TOK, k >= 1)."""
    if not len(pos0):
        return None
    off = (np.asarray(pos0) - G.MARGIN * G.TOK) % G.STEP
    d = np.minimum(off, G.STEP - off)
    ok = np.asarray(pos0) >= G.STEP  # no seam before the second tile
    return float(((d <= SEAM_TOL) & ok).mean())


def d0(preds, truth):
    out = {"expected_uniform": (2 * SEAM_TOL + 1) / G.STEP}
    for c in ("lnc", "mrna"):
        ps = [p for p in preds if p["cls"] == c]
        out[f"pred_{c}_start"] = seam_frac([p["exons"][0][0] - 1 for p in ps])
        out[f"pred_{c}_end"] = seam_frac([p["exons"][-1][1] - 1 for p in ps])
        ts = [span(g) for g in truth.values() if g["cls"] == c]
        out[f"truth_{c}_start"] = seam_frac([s - 1 for s, _ in ts])
        out[f"truth_{c}_end"] = seam_frac([e - 1 for _, e in ts])
        out[f"n_pred_{c}"] = len(ps)
    return out


def d1(preds, truth):
    idx = ruler.ExonIndex(truth)
    lp = [p for p in preds if p["cls"] == "lnc"]
    hit = set()
    for p in lp:
        hit |= {x for x in idx.hits(p["strand"], p["exons"]) if truth[x]["cls"] == "lnc"}
    cnt = {"ok": 0, "pred_too_long": 0, "pred_too_short": 0, "shifted": 0}
    for gid in hit:
        g = truth[gid]
        gs, ge = span(g)
        best = None
        for p in lp:
            if p["strand"] != g["strand"]:
                continue
            s, e = p["exons"][0][0], p["exons"][-1][1]
            ov = min(e, ge) - max(s, gs) + 1
            if ov > 0 and (best is None or ov > best[0]):
                best = (ov, ov / (ge - gs + 1), ov / (e - s + 1))
        _, rg, rp = best
        k = "ok" if rg >= 0.5 and rp >= 0.5 else "pred_too_long" if rg >= 0.5 else "pred_too_short" if rp >= 0.5 else "shifted"
        cnt[k] += 1
    n = len(hit)
    return {"n_L1_hit_truth": n, **{k: v / n for k, v in cnt.items()}}


def edge_dist(ind, center, rising):
    """distance from index `center` to the nearest rising (ind[i] and not ind[i-1]) / falling (ind[i] and not ind[i+1])
    position of the boolean vector, None if none."""
    if rising:
        pos = np.flatnonzero(ind[1:] & ~ind[:-1]) + 1
    else:
        pos = np.flatnonzero(ind[:-1] & ~ind[1:])
    return int(np.abs(pos - center).min()) if len(pos) else None


def d3(post, truth, cls):
    L = len(post)
    k = KIND[cls]
    rows = []
    for gid, g in truth.items():
        if g["cls"] != cls:
            continue
        st = g["strand"]
        cols = [g1.IDX[f"{k}{p}{st}"] for p in g1.PARTS]
        s, e = span(g)  # 1-based closed
        if s - 1 - SEARCH < 0 or e + SEARCH > L:
            continue

        def mass(a, b):
            return post[a:b][:, cols].astype(np.float32).sum(1)

        inside = float(mass(s - 1, e).mean())
        row = {"detected": inside >= 0.5}
        for name, at_start in (("tss", st == "+"), ("tes", st == "-")):
            c = s - 1 if at_start else e - 1  # 0-based end position on the genome
            m = mass(c - SEARCH, c + SEARCH)
            ind = m >= 0.5
            if at_start:
                step = float(m[SEARCH:SEARCH + FLANK].mean() - m[SEARCH - FLANK:SEARCH].mean())
                dist = edge_dist(ind, SEARCH, True)
            else:
                step = float(m[SEARCH - FLANK + 1:SEARCH + 1].mean() - m[SEARCH + 1:SEARCH + 1 + FLANK].mean())
                dist = edge_dist(ind, SEARCH, False)
            row[f"{name}_step"], row[f"{name}_dist"] = step, dist
        rows.append(row)
    det = [r for r in rows if r["detected"]]
    out = {"n": len(rows), "frac_detected": len(det) / len(rows) if rows else None}
    for name in ("tss", "tes"):
        out[f"{name}_step_all"] = float(np.mean([r[f"{name}_step"] for r in rows]))
        if det:
            dd = [r[f"{name}_dist"] for r in det]
            fin = [x for x in dd if x is not None]
            out[f"{name}_step_detected"] = float(np.mean([r[f"{name}_step"] for r in det]))
            out[f"{name}_dist_median_detected"] = statistics.median(fin) if fin else None
            for lim in (200, 1000, 5000):
                out[f"{name}_within_{lim}_detected"] = sum(x is not None and x <= lim for x in dd) / len(dd)
    return out


def matched_ids(truth, preds, cls, th=0.5):
    ts = B.truth_spans(truth, cls)
    ps = [(p["strand"], p["exons"][0][0], p["exons"][-1][1]) for p in preds if p["cls"] == cls]
    adj = [[] for _ in ps]
    for j, (st, s, e) in enumerate(ps):
        for i, (_, gst, gs, ge) in enumerate(ts):
            if st != gst:
                continue
            ov = min(e, ge) - max(s, gs) + 1
            if ov > 0 and ov / (ge - gs + 1) >= th and ov / (e - s + 1) >= th:
                adj[j].append(i)
    return {ts[i][0] for i in B.matching(adj, len(ts)) if i != -1}


def d4(preds, truth, cage, polya):
    def tss_tes(st, s, e):
        return (s - 1, e - 1) if st == "+" else (e - 1, s - 1)

    sup = {}
    for gid, g in truth.items():
        if g["cls"] == "lnc":
            t5, t3 = tss_tes(g["strand"], *span(g))
            sup[gid] = (bool(S3.near(cage[g["strand"]], t5, 200)), bool(S3.near(polya[g["strand"]], t3, 100)))
    m = matched_ids(truth, preds, "lnc")
    out = {"n_truth_lnc": len(sup), "frac_tss_cage": float(np.mean([a for a, _ in sup.values()])),
           "frac_tes_polya": float(np.mean([b for _, b in sup.values()]))}
    for name, f in (("both", lambda a, b: a and b), ("one", lambda a, b: a != b), ("neither", lambda a, b: not a and not b)):
        ids = [gid for gid, (a, b) in sup.items() if f(a, b)]
        out[f"S50_sn_{name}"] = {"n": len(ids), "sn": sum(gid in m for gid in ids) / len(ids) if ids else None}
    idx = ruler.ExonIndex(truth)
    pc, tc, rc = [], [], []
    for p in preds:
        if p["cls"] != "lnc":
            continue
        hit = [x for x in idx.hits(p["strand"], p["exons"]) if truth[x]["cls"] == "lnc"]
        if not hit:
            continue
        st = p["strand"]
        p5, _ = tss_tes(st, p["exons"][0][0], p["exons"][-1][1])
        pc.append(bool(S3.near(cage[st], p5, 200)))
        rc.append(bool(S3.near(cage[st], p5 + 50_000, 200)))
        tc.append(any(sup[x][0] for x in hit))
    out["L1_hit_pred_start_cage"] = float(np.mean(pc)) if pc else None
    out["their_truth_tss_cage"] = float(np.mean(tc)) if tc else None
    out["pred_start_shifted_50kb_cage"] = float(np.mean(rc)) if rc else None
    return out


def mean_of(runs, path):
    v = []
    for r in runs:
        x = r
        for k in path:
            x = x.get(k) if isinstance(x, dict) else None
            if x is None:
                break
        if isinstance(x, (int, float)) and not isinstance(x, bool):
            v.append(x)
    return sum(v) / len(v) if v else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="outputs/PHASE1/BND")
    ap.add_argument("--skip-s1", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    os.makedirs(f"{ROOT}/{a.out}", exist_ok=True)
    cfg = json.load(open(f"{ROOT}/configs/p15_v19.json"))
    d, tg = cfg["decode"], cfg["targets"]["chr18"]
    truth = ruler.read_gtf(f"{ROOT}/{tg['truth']}", tg["truth_chrom"])
    codes = np.load(f"{ROOT}/{tg['codes']}", mmap_mode="r")
    cage = S3.load_peaks(f"{ROOT}/data/bioval/raw/fantom5_cage_peaks.bed.gz", "chr18", 5, 6, 7)
    polya = S3.load_peaks(f"{ROOT}/data/bioval/raw/polyasite2_GRCh38.96.tsv.gz", "chr18", 5, 1, 2, skip_header=True,
                          chrom_name="18")
    res = {"plan": "docs/phase1_results/PLAN_BND.md", "runs": {}}
    for run in cfg["main"]["runs"].values():
        rd = f"{ROOT}/{cfg['out_dir']}/chr18/{run}"
        genes = P.load_genes(f"{rd}/genes.json")
        post = np.load(f"{rd}/post.npy", mmap_mode="r")
        preds = P.final_preds(genes, d["tau"])
        r = {"D0": d0(preds, truth), "D1": d1(preds, truth),
             "D3": {c: d3(post, truth, c) for c in ("lnc", "mrna")}, "D4": d4(preds, truth, cage, polya),
             "D2": {"meax": B.score(truth, preds)}}
        if not a.skip_s1:
            s1 = P.s1_genes(post, codes, d)
            r["D2"]["s1"] = B.score(truth, P.final_preds(s1, d["tau"]))
            r["D0_s1"] = d0(P.final_preds(s1, d["tau"]), truth)
        res["runs"][run] = r
        print(run, round(time.time() - t0), "s", json.dumps({"D0": r["D0"], "D1": r["D1"]}), flush=True)
        json.dump(P.norm(res), open(f"{ROOT}/{a.out}/BND.json", "w"), indent=1)
    runs = list(res["runs"].values())
    m = lambda *p: mean_of(runs, p)  # noqa: E731
    seam = np.mean([m("D0", k) for k in ("pred_lnc_start", "pred_lnc_end")])
    s1_sn, mx_sn = m("D2", "s1", "lnc", "S50", "sn"), m("D2", "meax", "lnc", "S50", "sn")
    l_tss, m_tss = m("D3", "lnc", "tss_within_1000_detected"), m("D3", "mrna", "tss_within_1000_detected")
    both, neither = m("D4", "S50_sn_both", "sn"), m("D4", "S50_sn_neither", "sn")
    res["verdict"] = {
        "seam": {"frac_pred_lnc_ends_at_seam": float(seam), "holds": bool(seam >= 0.20)},
        "decoder": {"s1_S50_sn": s1_sn, "meax_S50_sn": mx_sn,
                    "holds": bool(s1_sn is not None and s1_sn >= 2 * mx_sn and s1_sn >= 0.10)},
        "no_posterior_signal": {"lnc_tss_within_1kb": l_tss, "mrna_tss_within_1kb": m_tss,
                                "holds": bool(l_tss < 0.30 and m_tss >= 0.60)},
        "truth_ends": {"S50_sn_both": both, "S50_sn_neither": neither,
                       "holds": bool(both is not None and neither is not None and both >= 2 * neither)}}
    res["seconds"] = round(time.time() - t0, 1)
    json.dump(P.norm(res), open(f"{ROOT}/{a.out}/BND.json", "w"), indent=1)
    print(json.dumps(res["verdict"], indent=1))


if __name__ == "__main__":
    main()
