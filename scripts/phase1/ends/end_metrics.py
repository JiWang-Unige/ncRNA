#!/usr/bin/env python3
"""PLAN_END.md sections 3-4: end localisation on chr18 (CPU).

Score sources: stage 0 = raw boundary logits of the v2.0 models (<ENDS>/stage0/<run>_bnd.npy, [L, 4] TSS+, TES+, TSS-, TES-,
class-agnostic: the same channel scores lnc and mRNA ends); stage 1 = probes (<ENDS>/end-<arm>-42/chr18_end.npy, [L, 8]
prep_end_labels.CH8, class-specific).
Truth: chr18 GENCODE v50 Ensembl_canonical transcripts (one per gene) of lncRNA / protein_coding genes.
Loc1k / Loc200: share of ends whose score argmax in [end - 20 kb, end + 20 kb] (same strand, own channel) lies within
1 kb / 200 bp.  Strata: CAGE (TSS, same strand +-200 bp) / polyA (TES, +-100 bp) support; lnc TSS by promoter CpG o/e
quartile (TSS +-500 bp, sense).  Secondary: peak-level AP (local maxima, +-250 bp non-maximum suppression, top 20,000
per channel, a peak matches an unmatched end of its class within 500 bp, greedy by score).

  end_metrics.py [--ends outputs/PHASE1/ENDS] -> <ends>/ENDS.json
"""
import argparse
import json
import os
import sys
import time

import numpy as np
from scipy.ndimage import maximum_filter1d

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
for d in ("scripts/phase1", "scripts/phase1/s3", "scripts/phase1/ends"):
    sys.path.insert(0, os.path.join(ROOT, d))
import prep_end_labels as PL  # noqa: E402
import prep_events as PE  # noqa: E402
import s3_features as S3  # noqa: E402
from ncrna_annot.data import labels  # noqa: E402

TRUTH = f"{ROOT}/outputs/PHASE1/DIAG_L2/chr18_v50.gtf"
HALF, TOPK, MATCH = 20_000, 20_000, 500
STAGE0_RUNS = ("v20-g1b-42-ke", "v20-g1b-123-ke", "v20-g1b-456-ke")
BND = {"TSS+": 0, "TES+": 1, "TSS-": 2, "TES-": 3}


def ends_table(seq):
    tx = labels.read_transcripts(TRUTH, {"chr18"})["chr18"]
    cage = S3.load_peaks(PL.CAGE, "chr18", 5, 6, 7)
    pas = S3.load_peaks(PL.PAS, "chr18", 5, 1, 2, skip_header=True, chrom_name="18")
    rows = []
    for k, s, ex, canon in tx:
        if not canon:
            continue
        ev = PE.tx_events(s, ex)
        t5, t3 = ev[f"TSS{s}"][0], ev[f"TES{s}"][0]
        prom = seq[max(0, t5 - 500):t5 + 501]
        rows.append({"k": k, "s": s, "TSS": t5, "TES": t3,
                     "TSS_ev": bool(S3.near(cage[s], t5, 200)), "TES_ev": bool(S3.near(pas[s], t3, 100)),
                     "cpg": S3.cpg_oe(prom)})
    return rows


def channel(src, k, e, s):
    if src == "stage0":
        return BND[f"{e}{s}"]
    return PL.CH8.index(f"{k}{e}{s}")


def loc(track, rows, e):
    L = len(track)
    d = []
    for r in rows:
        x = r[e]
        a, b = max(0, x - HALF), min(L, x + HALF + 1)
        w = np.asarray(track[a:b], np.float32)
        d.append(abs(a + int(np.argmax(w)) - x))
    d = np.array(d)
    return {"n": len(d), "loc1k": float((d <= 1000).mean()) if len(d) else None,
            "loc200": float((d <= 200).mean()) if len(d) else None,
            "median_dist": float(np.median(d)) if len(d) else None}


def peak_ap(track, pos):
    x = np.asarray(track, np.float32)
    mx = maximum_filter1d(x, size=501, mode="nearest")
    cand = np.flatnonzero(x == mx)
    if len(cand) > TOPK:
        cand = cand[np.argpartition(-x[cand], TOPK)[:TOPK]]
    cand = cand[np.argsort(-x[cand])]
    pos = np.sort(np.asarray(pos))
    used = np.zeros(len(pos), bool)
    tp = []
    for p in cand:
        i = np.searchsorted(pos, p - MATCH)
        hit = False
        while i < len(pos) and pos[i] <= p + MATCH:
            if not used[i]:
                used[i] = hit = True
                break
            i += 1
        tp.append(hit)
    tp = np.array(tp, float)
    if not len(pos) or not len(tp):
        return None
    prec = np.cumsum(tp) / np.arange(1, len(tp) + 1)
    return {"ap": float((prec * tp).sum() / len(pos)), "recall_at_topk": float(tp.sum() / len(pos)), "n_true": len(pos)}


def score_source(src, track, rows):
    out = {}
    for k in ("l", "m"):
        cls = "lnc" if k == "l" else "mrna"
        for e in ("TSS", "TES"):
            res = {}
            # localisation (both strands pooled, each end on its own strand channel)
            def run(sel):
                parts = []
                for s in ("+", "-"):
                    rr = [r for r in rows if r["k"] == k and r["s"] == s and sel(r)]
                    parts.append((loc(track[:, channel(src, k, e, s)], rr, e), len(rr)))
                n = sum(p["n"] for p, _ in parts)
                if not n:
                    return {"n": 0}
                comb = {"n": n}
                for m in ("loc1k", "loc200"):
                    comb[m] = sum((p[m] or 0) * p["n"] for p, _ in parts) / n
                return comb
            res["all"] = run(lambda r: True)
            res["with_evidence"] = run(lambda r: r[f"{e}_ev"])
            res["without_evidence"] = run(lambda r: not r[f"{e}_ev"])
            if e == "TSS" and k == "l":
                q = np.quantile([r["cpg"] for r in rows if r["k"] == "l"], [0.25, 0.5, 0.75])
                for qi in range(4):
                    lo_ = -np.inf if qi == 0 else q[qi - 1]
                    hi_ = np.inf if qi == 3 else q[qi]
                    res[f"cpg_q{qi + 1}"] = run(lambda r, lo_=lo_, hi_=hi_: lo_ < r["cpg"] <= hi_)
            aps = []
            for s in ("+", "-"):
                pos = [r[e] for r in rows if r["k"] == k and r["s"] == s]
                aps.append(peak_ap(track[:, channel(src, k, e, s)], pos))
            res["peak_ap_mean_strands"] = float(np.mean([x["ap"] for x in aps if x])) if any(aps) else None
            out[f"{cls}_{e}"] = res
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ends", default="outputs/PHASE1/ENDS")
    a = ap.parse_args()
    t0 = time.time()
    E = f"{ROOT}/{a.ends}"
    seq = S3.read_fasta(f"{ROOT}/data/human/chr18.fa").upper()
    rows = ends_table(seq)
    res = {"plan": "docs/phase1_results/PLAN_END.md", "n_ends": {
        f"{c}_{e}{'_ev' if ev else ''}": sum(r["k"] == c and (not ev or r[f'{e}_ev']) for r in rows)
        for c in ("l", "m") for e in ("TSS", "TES") for ev in (False, True)}, "sources": {}}
    for run in STAGE0_RUNS:
        p = f"{E}/stage0/{run}_bnd.npy"
        if os.path.exists(p + ".done"):
            res["sources"][f"stage0/{run}"] = score_source("stage0", np.load(p, mmap_mode="r"), rows)
            print("scored", run, round(time.time() - t0), flush=True)
    for arm in ("ann", "ev"):
        d = f"{E}/end-{arm}-42"
        if os.path.exists(f"{d}/CHR18_DONE"):
            res["sources"][f"P-{arm}"] = score_source("probe", np.load(f"{d}/chr18_end.npy", mmap_mode="r"), rows)
            print("scored", arm, round(time.time() - t0), flush=True)

    S = res["sources"]
    g = lambda src, key, sub="with_evidence", m="loc1k": S[src][key][sub].get(m) if src in S else None  # noqa: E731
    s0 = [g(f"stage0/{r}", "lnc_TSS") for r in STAGE0_RUNS if f"stage0/{r}" in S]
    st0 = float(np.mean(s0)) if s0 else None
    pa, pe = g("P-ann", "lnc_TSS"), g("P-ev", "lnc_TSS")
    ma, me = g("P-ann", "mrna_TSS"), g("P-ev", "mrna_TSS")
    if None not in (pa, pe, ma, me, st0):
        A = pe >= 0.50 and pe >= pa + 0.15
        B = pa >= 0.50 and st0 <= pa - 0.15
        C = pa < 0.30 and pe < 0.30 and ma >= 0.60 and me >= 0.60
        res["verdict"] = {"lnc_TSS_loc1k_with_CAGE": {"stage0_mean": st0, "P-ann": pa, "P-ev": pe},
                          "mrna_TSS_loc1k_with_CAGE": {"P-ann": ma, "P-ev": me},
                          "A_labels": A, "B_structure": B, "C_not_locatable": C,
                          "outcome": "A" if A else "B" if B else "C" if C else "grey"}
        print(json.dumps(res["verdict"], indent=1))
    res["seconds"] = round(time.time() - t0, 1)
    json.dump(res, open(f"{E}/ENDS.json", "w"), indent=1)


if __name__ == "__main__":
    main()
