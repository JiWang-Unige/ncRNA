#!/usr/bin/env python3
"""M1 weeks 2-4: the one-time test scoring of docs/m1/PLAN_M1_W24.md (sections 3-7), pooled over chr9 + chr12 + chr20 + chr13.

Methods (all with the Tiberius same-strand coding filter, as in development):
  ours        : candidates = each ssCAGE TSS x the K = 3 nearest downstream same-strand QuantSeq sites within 50 kb; end = the
                farthest of them; units accepted greedily by ssCAGE tags (no same-strand overlap); kept if the run's DNA-only
                span score (mean lnc posterior - mean mRNA posterior on that strand) >= 0.  One result per v2.0 run.
  ev_farthest : the same without the model filter (evidence only).
  pair        : nearest site within 50 kb (evidence only; the development pairing baseline).
  snap        : the run's DNA-only lnc calls (calibrated, tau_v20) with each end snapped to the nearest evidence within 1 kb.
  dna         : the run's DNA-only lnc calls.
  stringtie / isoquant : reference-free, PacBio rep1 (outputs/M1/W24/longread/<chrom>).
Evidence: 100 % and the gradient 30 % / 10 % (seeds 0-2).  Truth: WTC11 rep2 & rep3 loci outside the exclusion zone; loci
crossing the chr12 / chr13 half cut are dropped, and predictions crossing it are dropped for every method.
  w24_score.py --mode test   (needs outputs/M1/W24/FROZEN = {file: sha256} of this file and its three modules; runs once)
                              -> outputs/M1/W24/RESULT.json
  w24_score.py --mode dev    (chr18 only, development check)                       -> outputs/M1/W24/dev/SCORE_DEV_chr18.json
"""
import argparse
import hashlib
import json
import os
import re
from collections import defaultdict

import numpy as np

import w24_lib as L
from dev_baselines import pair, snap, tib_hit
from e3_anchor import candidates, select_evidence, select_hybrid, span_mean

TEST = ("chr9", "chr12", "chr20", "chr13")
FRACS = ((1.0, (0,)), (0.3, (0, 1, 2)), (0.1, (0, 1, 2)))
DMAX, K, THETA, SNAP_W = 50_000, 3, 0.0, 1_000
N_BOOT = 2000


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def read_tib(path, chrom):
    spans = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(f) < 9 or f[0] != chrom or f[2] not in ("CDS", "exon"):
            continue
        gid = f[8].split('gene_id "', 1)[1].split('"', 1)[0]
        s, e = int(f[3]), int(f[4])
        spans[gid] = (f[6], min(spans.get(gid, (0, s, e))[1], s), max(spans.get(gid, (0, s, e))[2], e))
    out = {"+": [], "-": []}
    for st, s, e in spans.values():
        out[st].append((s, e))
    return {k: np.array(sorted(v), np.int64).reshape(-1, 2) for k, v in out.items()}


def sources(chrom, cfg):
    """[(part dir template, offset0, length, codes path)] for the run reads covering the chromosome."""
    if chrom in ("chr9", "chr18"):
        n = len(np.load(f"{L.ROOT}/data/p1/{chrom}_codes.npy", mmap_mode="r"))
        return [(f"{L.ROOT}/outputs/PHASE1/V20/read/{chrom}/{{run}}", 0, n, f"{L.ROOT}/data/p1/{chrom}_codes.npy")]
    pj = json.load(open(f"{L.W}/PARTS.json"))["parts"]
    out = []
    for name, p in sorted(pj.items(), key=lambda kv: kv[1]["offset0"]):
        if p["chrom"] == chrom:
            out.append((f"{L.W}/read/{name}/{{run}}", p["offset0"], p["length"], f"{L.ROOT}/data/p1/{name}_codes.npy"))
    return out


def cuts_of(chrom):
    return [o for _, o, _, _ in sources(chrom, None)[1:]] if chrom not in ("chr9", "chr18") else []


def run_inputs(chrom, run, cfg, tab):
    """Chromosome-wide cumulative posteriors (lnc / mRNA, per strand) and calibrated DNA-only lnc calls of one run."""
    import v20_calib as Cc
    from ncrna_annot.grammar import g1
    srcs = sources(chrom, cfg)
    Ltot = srcs[-1][1] + srcs[-1][2]
    cs = {(k, st): np.zeros(Ltot + 1, np.float64) for k in "lm" for st in "+-"}
    calls = []
    cal = cfg["calibration"]
    for tmpl, off, n, codes_p in srcs:
        d = tmpl.format(run=run)
        post = np.load(f"{d}/post.npy", mmap_mode="r")
        m = len(post)                      # reads cover the 6-mer token grid: up to 5 trailing bases have no posterior
        assert 0 <= n - m < 6, (d, m, n)
        for k in "lm":
            for st in "+-":
                idx = [g1.IDX[f"{k}{p}{st}"] for p in ("E", "L", "I", "R")]
                for a in range(0, m, 5_000_000):
                    b = min(a + 5_000_000, m)
                    cs[(k, st)][off + a + 1:off + b + 1] = np.asarray(post[a:b][:, idx], np.float32).sum(1)
        genes = json.load(open(f"{d}/genes.json"))
        codes = np.load(codes_p, mmap_mode="r")
        rows = Cc.gene_rows(genes, post, codes, tab)
        X = np.array([[r[f] for f in cal["model"]["features"]] for r in rows], float)
        sc = dict(zip((r["k"] for r in rows), Cc.apply(cal["model"], X)))
        for j, g in enumerate(genes):
            conf = float(sc.get(j, g["conf"]))   # as simple_det.load_ours: genes without a calibration row keep their raw conf
            if g["type"] == "l" and conf >= cal["tau_v20"]:
                calls.append({"strand": g["strand"], "exons": [(a + off, b + off) for a, b in g["exons"]], "conf": conf})
    for key in cs:
        cs[key] = np.cumsum(cs[key])
    return cs, calls


def longread(chrom, tool):
    O = f"{L.W}/longread/{chrom}"
    gtf = f"{O}/stringtie.gtf" if tool == "stringtie" else f"{O}/isoquant/rep1/rep1.transcript_models.gtf"
    g = defaultdict(lambda: [None, 10 ** 12, 0])
    for line in open(gtf):
        f = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(f) < 9 or f[0] != chrom or f[2] != "transcript" or f[6] not in "+-":
            continue
        x = g[(re.search(r'gene_id "([^"]+)"', f[8]).group(1), f[6])]
        x[0], x[1], x[2] = f[6], min(x[1], int(f[3])), max(x[2], int(f[4]))
    return [{"strand": st, "exons": [(s, e)]} for st, s, e in g.values() if e - s + 1 >= 200]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True, choices=["test", "dev"])
    a = ap.parse_args()
    if a.mode == "test":
        here = os.path.dirname(os.path.abspath(__file__))
        frozen = json.load(open(f"{L.W}/FROZEN"))
        for f in ("w24_score.py", "w24_lib.py", "e3_anchor.py", "dev_baselines.py"):
            assert frozen[f] == sha(f"{here}/{f}"), f"{f} differs from the frozen version"
        assert not os.path.exists(f"{L.W}/RESULT.json"), "the test has already been read once"
        chroms = TEST
    else:
        chroms = ("chr18",)
    cfg = json.load(open(f"{L.ROOT}/configs/v20.json"))
    tab = json.load(open(f"{L.ROOT}/{cfg['calibration']['hexamer_table']}"))
    runs = cfg["targets"]["chr18"]["runs"]
    tibdir = {c: (f"{L.W}/tiberius/{c}.gtf" if c in ("chr12", "chr13") else f"{L.ROOT}/outputs/PHASE1/BASELINES/tiberius/{c}.gtf") for c in chroms}
    per_locus = defaultdict(dict)     # method -> {locus id: matched (0/1, or mean over runs)}
    counts = defaultdict(lambda: defaultdict(float))
    info = {"chroms": {}}
    for c in chroms:
        truth = L.load_truth(c)
        cuts = cuts_of(c)
        truth["ts"] = [t for t in truth["ts"] if not any(t[2] - 1 < x < t[3] for x in cuts)]   # 0-based cut inside the span
        truth["lnc"] = {i: v for i, v in truth["lnc"].items() if i in {t[0] for t in truth["ts"]}}
        truth["zone"] = list(truth["zone"]) + [(x - 1, x + 1) for x in cuts]
        tib = read_tib(tibdir[c], c)
        sig = L.quant_signal(c)
        info["chroms"][c] = {"n_truth_lnc": len(truth["ts"]), "cuts": cuts}
        res_c = {}
        lr = {t: [p for p in longread(c, t) if not tib_hit(tib, p["strand"], p["exons"][0][0], p["exons"][-1][1])]
              for t in ("stringtie", "isoquant")} if os.path.exists(f"{L.W}/longread/{c}/stringtie.gtf") else {}
        ri = {r: run_inputs(c, r, cfg, tab) for r in runs}
        for frac, seeds in FRACS:
            for sd in seeds:
                ev = L.evidence(c, frac, sd)
                tag = f"f{frac}_s{sd}"
                cands = candidates(ev, DMAX, K)
                cands1 = candidates(ev, DMAX, 1)
                info["chroms"][c][f"coverage_{tag}"] = L.coverage(truth, ev, 100, 100)
                sets = {"ev_farthest": [select_evidence(cands, sig, "farthest", tib)],
                        "pair": [select_evidence(cands1, sig, "nearest", tib)],
                        "ours": [select_hybrid(cands, ri[r][0], sig, "farthest", THETA, tib) for r in runs],
                        "snap": [snap(ri[r][1], ev, SNAP_W) for r in runs]}
                if tag == "f1.0_s0":
                    sets["dna"] = [ri[r][1] for r in runs]
                    sets.update({t: [v] for t, v in lr.items()})
                for m, plist in sets.items():
                    rs = [L.score_unknown(truth, p, c) for p in plist]
                    key = f"{m}|{tag}"
                    ids = [set(r["S50_ids"]) for r in rs]
                    for t in truth["ts"]:
                        per_locus[key][t[0]] = float(np.mean([t[0] in s for s in ids]))
                    for r in rs:
                        n_tp = len(r["S50_ids"])
                        counts[key]["tp"] += n_tp / len(rs)
                        counts[key]["n_pred"] += r["n_pred"] / len(rs)
                        counts[key]["n_unknown"] += r["n_unknown"] / len(rs)
                    res_c[key] = {"S50_sn": float(np.mean([r["S50_sn"] for r in rs])),
                                  "S50_pr_known": float(np.mean([r["S50_pr_known"] for r in rs if r["S50_pr_known"] is not None])) if any(r["S50_pr_known"] is not None for r in rs) else None,
                                  "n_pred": float(np.mean([r["n_pred"] for r in rs]))}
        info["chroms"][c]["by_method"] = res_c
    # pooled
    loci = sorted(per_locus["ours|f1.0_s0"])
    nt = len(loci)
    pooled = {}
    for key, cnt in counts.items():
        tp, n, unk = cnt["tp"], cnt["n_pred"], cnt["n_unknown"]
        pooled[key] = {"S50_sn": tp / nt if nt else None, "S50_pr": tp / n if n else None,
                       "S50_pr_known": tp / (n - unk) if n - unk > 0 else None, "n_pred": n}
    rng = np.random.default_rng(0)
    B = [rng.integers(0, nt, nt) for _ in range(N_BOOT)]

    def diff(m1, m2):
        x = np.array([per_locus[m1][i] - per_locus[m2][i] for i in loci])
        return {"diff": float(x.mean()), "ci95": np.percentile([x[b].mean() for b in B], [2.5, 97.5]).tolist()}

    def ci(m):
        x = np.array([per_locus[m][i] for i in loci])
        return np.percentile([x[b].mean() for b in B], [2.5, 97.5]).tolist()

    out = {"plan": "docs/m1/PLAN_M1_W24.md", "mode": a.mode, "n_truth_lnc": nt, "info": info, "pooled": pooled,
           "ci95_S50_sn": {k: ci(k) for k in per_locus}, "paired": {}}
    for tag in [f"f{f}_s{s}" for f, ss in FRACS for s in ss]:
        out["paired"][f"ours-ev_farthest|{tag}"] = diff(f"ours|{tag}", f"ev_farthest|{tag}")
        out["paired"][f"ours-snap|{tag}"] = diff(f"ours|{tag}", f"snap|{tag}")
        out["paired"][f"ours-pair|{tag}"] = diff(f"ours|{tag}", f"pair|{tag}")
    if "isoquant|f1.0_s0" in per_locus:
        only = [i for i in loci if per_locus["ours|f1.0_s0"][i] >= 0.5 and per_locus["isoquant|f1.0_s0"][i] == 0]
        out["ours_not_isoquant"] = {"n": len(only), "frac_of_truth": len(only) / nt}
    o = pooled["ours|f1.0_s0"]
    gate = {"S50_sn": o["S50_sn"], "S50_pr_known": o["S50_pr_known"], "vs_snap": out["paired"]["ours-snap|f1.0_s0"]}
    gate["met"] = bool(o["S50_sn"] >= 0.30 and (o["S50_pr_known"] or 0) >= 0.40 and gate["vs_snap"]["ci95"][0] > 0)
    e = out["paired"]["ours-ev_farthest|f1.0_s0"]
    dp = (o["S50_pr_known"] or 0) - (pooled["ev_farthest|f1.0_s0"]["S50_pr_known"] or 0)
    seq = bool(e["ci95"][0] > 0 or (dp >= 0.05 and e["diff"] >= -0.02))
    out["verdict"] = {"gate": gate, "sequence_prior_contributes": seq, "precision_known_minus_ev": dp,
                      "decision_hint": "B" if gate["met"] and seq else ("D-leaning (evidence works, our prior adds nothing)" if not seq else "gate not met")}
    path = f"{L.W}/RESULT.json" if a.mode == "test" else f"{L.W}/dev/SCORE_DEV_chr18.json"
    json.dump(out, open(path, "w"), indent=1)
    print(json.dumps({"n_truth": nt, "verdict": out["verdict"], "pooled_f1": {k: v for k, v in pooled.items() if k.endswith("f1.0_s0")}}, indent=1))


if __name__ == "__main__":
    main()
