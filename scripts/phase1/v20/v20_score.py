#!/usr/bin/env python3
"""v2.0 scoring (protocol section 17 / DRAFT_V20 section 4.1) of the reads of configs/v20.json on one target.

  v20_score.py --config configs/v20.json --target chr9    -> <out_dir>/chr9/SCORE_V20.json
  (--target chr18 is allowed as a dry run on the development chromosome)

Per run, from <out_dir>/<target>/<run>/{genes.json, post.npy}:
  baseline = v1.9 decoding of the retrained run, lnc genes kept when raw confidence >= decode.tau (tau_baseline);
  v2.0     = same genes, lnc genes kept when the calibrated probability (calibration.model on the v20_calib.py
             features, hexamer table checked by sha256) >= calibration.tau_v20; raw confidence when keep is false.
mRNA genes are never filtered, so the two differ only in which lnc genes are reported.  Ruler: p15_score.score_preds
(P1.5).  Paired comparisons, 2,000 resamples, numpy default_rng(0), statistic = three-run mean:
  Sn (L1 locus sensitivity): resample the truth lnc genes; per resample Sn = mean hit indicator over the sample;
  strict Pr: resample the decoded lnc genes (the union both arms choose from); per resample strict Pr of each arm =
             hits / reported among the sampled genes.
ECE: lnc genes before tau, 10 equal-width bins, raw confidence (baseline) vs calibrated probability (v2.0).
Section 4.1 verdicts are written with their inputs; nothing is tuned here.
"""
import argparse
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1/v20"))
import v20_calib as C  # noqa: E402
import p15_score as P  # noqa: E402
from ncrna_annot.eval import ruler  # noqa: E402

N_BOOT = 2000


def hits_of(preds, idx):
    return [bool(idx.hits(p["strand"], p["exons"])) for p in preds]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/v20.json")
    ap.add_argument("--target", default="chr9")
    a = ap.parse_args()
    t0 = time.time()
    cfg_path = f"{ROOT}/{a.config}"
    cfg = json.load(open(cfg_path))
    cal, d = cfg["calibration"], cfg["decode"]
    assert P.sha256(f"{ROOT}/{cal['hexamer_table']}") == cal["hexamer_table_sha256"], "hexamer table changed"
    tab = json.load(open(f"{ROOT}/{cal['hexamer_table']}"))
    tg = cfg["targets"][a.target]
    truth = ruler.read_gtf(f"{ROOT}/{tg['truth']}", tg["truth_chrom"])
    lnc_truth = {k: v for k, v in truth.items() if v["cls"] == "lnc"}
    lidx = ruler.ExonIndex(lnc_truth)
    codes = np.load(f"{ROOT}/{tg['codes']}", mmap_mode="r")
    out_dir = f"{ROOT}/{cfg['out_dir']}/{a.target}"
    tau_b, tau_v = d["tau"], cal["tau_v20"]
    res = {"config": a.config, "config_sha256": P.sha256(cfg_path), "target": a.target, "tau_baseline": tau_b,
           "tau_v20": tau_v, "calibration_kept": cal["keep"], "runs": {}}
    tids = sorted(lnc_truth)
    hit_b, hit_v, pairs = {}, {}, {}
    for run in tg["runs"]:
        rd = f"{out_dir}/{run}"
        meta = json.load(open(f"{rd}/read_meta.json"))
        assert meta["config_sha256"] == res["config_sha256"], f"{run}: read with another config"
        genes = P.load_genes(f"{rd}/genes.json")
        post = np.load(f"{rd}/post.npy", mmap_mode="r")
        rows = C.gene_rows(genes, post, codes, tab)
        X = np.array([[r[f] for f in cal["model"]["features"]] for r in rows], float)
        p = C.apply(cal["model"], X) if cal["keep"] else X[:, cal["model"]["features"].index("conf")]
        score = {r["k"]: float(s) for r, s in zip(rows, p)}
        gv = [dict(g, conf=score.get(k, g["conf"])) for k, g in enumerate(genes)]
        fb, _ = P.score_preds(truth, P.final_preds(genes, tau_b), boot=N_BOOT)
        fv, _ = P.score_preds(truth, P.final_preds(gv, tau_v), boot=N_BOOT)
        lnc_rows = [(k, g) for k, g in enumerate(genes) if g["type"] == "l"]
        y = np.array([1.0 if lidx.hits(g["strand"], g["exons"]) else 0.0 for _, g in lnc_rows])
        raw = np.array([g["conf"] for _, g in lnc_rows])
        calp = np.array([score[k] for k, _ in lnc_rows])
        res["runs"][run] = {"baseline": P.norm(fb), "v20": P.norm(fv),
                            "ece_baseline": C.ece(raw, y), "ece_v20": C.ece(calp, y),
                            "n_lnc_decoded": len(lnc_rows), "n_lnc_reported": {"baseline": int((raw >= tau_b).sum()), "v20": int((calp >= tau_v).sum())}}
        # per-truth-gene hits (Sn) and per-decoded-gene (reported?, hit) for the paired bootstrap
        for arm, gs, tau, store in (("b", genes, tau_b, hit_b), ("v", gv, tau_v, hit_v)):
            rep = [g for g in gs if g["type"] == "l" and g["conf"] >= tau]
            hit_genes = set()
            for g in rep:
                hit_genes |= set(lidx.hits(g["strand"], g["exons"]))
            store[run] = np.array([t in hit_genes for t in tids], float)
        pairs[run] = (raw >= tau_b, calp >= tau_v, y)
    rng = np.random.default_rng(0)
    runs = tg["runs"]
    sn_b = np.mean([hit_b[r] for r in runs], 0)
    sn_v = np.mean([hit_v[r] for r in runs], 0)
    d_sn = []
    for _ in range(N_BOOT):
        i = rng.integers(0, len(tids), len(tids))
        d_sn.append(sn_v[i].mean() - sn_b[i].mean())
    d_pr = []
    for _ in range(N_BOOT):
        vals = []
        for r in runs:
            rb, rv, y = pairs[r]
            i = rng.integers(0, len(y), len(y))
            pb = y[i][rb[i]].mean() if rb[i].any() else 0.0
            pv = y[i][rv[i]].mean() if rv[i].any() else 0.0
            vals.append(pv - pb)
        d_pr.append(np.mean(vals))
    mean = lambda arm, key: float(np.mean([res["runs"][r][arm]["lnc"]["L1_locus"][key] or 0.0 for r in runs]))
    m = {"lnc_sn": {"baseline": mean("baseline", "sn"), "v20": mean("v20", "sn")},
         "lnc_pr_strict": {"baseline": mean("baseline", "pr_strict"), "v20": mean("v20", "pr_strict")},
         "lnc_L2_sn": {arm: float(np.mean([res["runs"][r][arm]["lnc"]["L2_intron_chain"]["sn"] for r in runs])) for arm in ("baseline", "v20")},
         "mrna_sn": {arm: float(np.mean([res["runs"][r][arm]["mrna"]["L1_locus"]["sn"] for r in runs])) for arm in ("baseline", "v20")},
         "guard_lnc_as_mrna": {arm: float(np.mean([res["runs"][r][arm]["guards"]["lnc_genes_reported_as_mrna"] for r in runs])) for arm in ("baseline", "v20")},
         "ece": {"baseline": float(np.mean([res["runs"][r]["ece_baseline"] for r in runs])), "v20": float(np.mean([res["runs"][r]["ece_v20"] for r in runs]))}}
    m["delta_sn"] = {"point": m["lnc_sn"]["v20"] - m["lnc_sn"]["baseline"], "ci95": np.percentile(d_sn, [2.5, 97.5]).tolist()}
    m["delta_pr_strict"] = {"point": m["lnc_pr_strict"]["v20"] - m["lnc_pr_strict"]["baseline"], "ci95": np.percentile(d_pr, [2.5, 97.5]).tolist()}
    res["mean"] = m
    res["verdicts"] = {
        "go_replication": bool(m["lnc_sn"]["v20"] >= 0.20 and m["lnc_pr_strict"]["v20"] >= 0.40),
        "structure_L2": bool(m["lnc_L2_sn"]["v20"] >= 0.05),
        "precision": bool(m["delta_pr_strict"]["ci95"][0] > 0 and m["lnc_sn"]["baseline"] - m["lnc_sn"]["v20"] <= 0.03),
        "calibration": bool(m["ece"]["v20"] <= 0.15),
        "guardrails": bool(m["mrna_sn"]["v20"] >= 0.80 and m["mrna_sn"]["baseline"] - m["mrna_sn"]["v20"] <= 0.02
                           and m["guard_lnc_as_mrna"]["v20"] <= m["guard_lnc_as_mrna"]["baseline"])}
    res["seconds"] = round(time.time() - t0)
    json.dump(res, open(f"{out_dir}/SCORE_V20.json", "w"), indent=1, default=lambda x: sorted(x) if isinstance(x, set) else str(x))
    print(json.dumps({"mean": m, "verdicts": res["verdicts"], "seconds": res["seconds"]}), flush=True)


if __name__ == "__main__":
    main()
