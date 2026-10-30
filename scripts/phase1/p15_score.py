#!/usr/bin/env python3
"""Protocol v1.9 section 16, P1.5 scoring (CPU) of the p15_read.py outputs (outputs/PHASE1/P15/<target>/<run>/).

  p15_score.py --config configs/p15_v19.json --target chr18       reproduction of the section 16 chr18 numbers
  p15_score.py --config configs/p15_v19.json --target t2t_chr20   main runs: ruler v1.3 + leaderboard metric (k = 250)
  p15_score.py --config configs/p15_v19.json --target chr20       section 3 endpoints, Go / Strong / guardrails
(score t2t_chr20 before chr20: the chr20 verdict takes the T2T part of Strong from P15/t2t_chr20/SCORE.json)

chr18: for every run of the config's "chr18_reference" the read's post.npy must have the sha256 of the development
dump, and the MEAX gamma* summary (screen2.summarize, conf_grid at tau*) must equal the development result in every
number; only then P15/chr18/REPRO_OK (config sha256) is written, which p15_read.py requires for chr20 / t2t_chr20.
One run on chr20 / t2t_chr20:
  core    screen2.summarize on all MEAX genes: AUPRC_lnc without threshold (as in development), conf_grid at tau*
  final   ruler v1.3 L1-L6 for lnc and mRNA of the final predictions (lnc confidence >= tau*, every mRNA), with
          2,000-resample CIs
  guards  truth lnc genes hit by a predicted mRNA that hits no truth mRNA (Phase 0 T2T-panel definition) and
          predicted mRNAs whose truth hits are lnc only (the Phase 0 hg38 chr20 count, 46 for GENATATOR)
  ece     predicted lnc genes (all, before tau*): confidence vs L1 hit, 10 equal-width bins, raw posteriors
  strata  (chr20) final-prediction L1 Sn by truth lnc length (longest merged transcript < 1 kb, 1-3 kb, >= 3 kb),
          exon count (1, 2, 3-4, >= 5), best transcript_support_level (1, 2-3, 4-5, NA), the v50-new subset
          (ruler.delta_genes against v41) and the non-homologous subsets (P15/homology_chr20.json; "masked" primary)
Verdict (chr20): 3-seed means of the main runs (config "main"); paired bootstrap over truth lnc genes (2,000
resamples, numpy rng seed 0) of mean_seed(hit) - GENATATOR(hit) for the L1 Sn gain.  Ablations: every arm's runs
as above with Delta AUPRC_lnc against the same-seed main run; AB1 = S1 (posterior-argmax segments) of the main
runs' posteriors with tau*; AB6 = its own evaluation (result.json).
"""
import argparse
import gzip
import hashlib
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
P15 = f"{ROOT}/outputs/PHASE1/P15"
BOOT = 2000
G_CHR20 = f"{ROOT}/outputs/PHASE0/0.1_genatator/human_chr20/pred.gff3"
G_CHR20_SCORE = f"{ROOT}/outputs/PHASE0/0.1_genatator/human_chr20/score.json"
V41 = f"{ROOT}/data/human/gencode.v41.annotation.gtf.gz"

from ncrna_annot.decode import chrom  # noqa: E402
from ncrna_annot.eval import ruler  # noqa: E402
from ncrna_annot.grammar import g1  # noqa: E402
from ncrna_annot.train import screen2 as s2  # noqa: E402


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def norm(x):
    return json.loads(json.dumps(x))


def load_genes(path):
    return [dict(type=g["type"], strand=g["strand"], exons=[tuple(x) for x in g["exons"]], conf=g["conf"])
            for g in json.load(open(path))]


def final_preds(genes, tau):
    return [{"id": str(i), "strand": g["strand"], "cls": "lnc" if g["type"] == "l" else "mrna", "exons": g["exons"]}
            for i, g in enumerate(genes) if g["type"] != "l" or g["conf"] >= tau]


def core(genes, post, truth, tau):
    r = s2.summarize([dict(g) for g in genes], post, truth, (tau,))
    r.pop("curve")
    return r


def hit_set(cls, truth, preds, idx):
    out = set()
    for p in preds:
        if p["cls"] == cls:
            out |= {x for x in idx.hits(p["strand"], p["exons"]) if truth[x]["cls"] == cls}
    return out


def guards(truth, preds, idx, sm):
    lnc_ids = {k for k, v in truth.items() if v["cls"] == "lnc"}
    as_mrna = set()
    for p in preds:
        if p["cls"] == "mrna":
            h = idx.hits(p["strand"], p["exons"])
            if h and not any(truth[x]["cls"] == "mrna" for x in h):
                as_mrna |= h & lnc_ids
    return {"lnc_genes_reported_as_mrna": len(as_mrna),
            "mrna_preds_hitting_lnc_only": sm["L1_locus"]["pred_hits_other_class_only"].get("lnc", 0)}


def ece(genes, truth, bins=10):
    idx = ruler.ExonIndex({k: v for k, v in truth.items() if v["cls"] == "lnc"})
    xs = np.array([(g["conf"], 1.0 if idx.hits(g["strand"], g["exons"]) else 0.0) for g in genes if g["type"] == "l"])
    if not len(xs):
        return None
    b = np.minimum((xs[:, 0] * bins).astype(int), bins - 1)
    rel = [{"bin": [i / bins, (i + 1) / bins], "n": int((b == i).sum()), "conf": float(xs[b == i, 0].mean()),
            "frac_hit": float(xs[b == i, 1].mean())} for i in range(bins) if (b == i).any()]
    return {"ece": float(sum(r["n"] * abs(r["conf"] - r["frac_hit"]) for r in rel) / len(xs)), "n": len(xs),
            "reliability": rel}


def score_preds(truth, preds, boot=BOOT):
    idx = ruler.ExonIndex(truth)
    sl, sm = (ruler.score_class(c, truth, preds, idx, boot) for c in ("lnc", "mrna"))
    return {"lnc": sl, "mrna": sm, "guards": guards(truth, preds, idx, sm),
            "n_pred_lnc": sum(p["cls"] == "lnc" for p in preds), "n_pred_mrna": sum(p["cls"] == "mrna" for p in preds)}, \
        hit_set("lnc", truth, preds, idx)


def score_run(run_dir, truth, tau, post=None):
    genes = load_genes(f"{run_dir}/genes.json")
    post = np.load(f"{run_dir}/post.npy", mmap_mode="r") if post is None else post
    fin, hits = score_preds(truth, final_preds(genes, tau))
    return {"core": core(genes, post, truth, tau), "final": fin, "ece": ece(genes, truth)}, hits


def s1_genes(post, codes, d):
    T = g1.G1Transitions()().detach().double().numpy()
    gs = chrom.decode_blocks(np.asarray(post, np.float32), codes[:len(post)], T, d["min_exon"], d["min_intron"], mode="argmax")
    gs = [g for g in gs if g["type"] != "l" or sum(b - a + 1 for a, b in g["exons"]) >= d["lnc_min_len"]]
    for g in gs:
        g["conf"] = chrom.gene_confidence(g, post)
    return gs


def tsl_by_gene(gtf, chrom_name):
    best = {}
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if not line.startswith(chrom_name + "\t"):
                continue
            f = line.split("\t", 9)
            if f[2] != "transcript":
                continue
            a = ruler._attrs_gtf(f[8])
            t = a.get("transcript_support_level", "NA").split()[0]
            v = int(t) if t.isdigit() else 6
            best[a["gene_id"]] = min(best.get(a["gene_id"], 9), v)
    return best


def strata(truth, gtf, chrom_name):
    lnc = {k: v for k, v in truth.items() if v["cls"] == "lnc"}
    tsl = tsl_by_gene(gtf, chrom_name)
    out = {}
    for gid, g in lnc.items():
        merged = [ruler._merge(t) for t in g["tx"].values()]
        ln, ne = max(sum(e - s + 1 for s, e in m) for m in merged), max(len(m) for m in merged)
        keys = ["len<1kb" if ln < 1000 else "len1-3kb" if ln < 3000 else "len>=3kb",
                "exons1" if ne == 1 else "exons2" if ne == 2 else "exons3-4" if ne <= 4 else "exons>=5"]
        t = tsl.get(gid, 6)
        keys.append("tsl1" if t == 1 else "tsl2-3" if t <= 3 else "tsl4-5" if t <= 5 else "tslNA")
        for k in keys:
            out.setdefault(k, set()).add(gid)
    old = ruler.read_gtf(V41, chrom_name)
    out["delta_v50_new"] = ruler.delta_genes(truth, old)
    hp = f"{P15}/homology_{chrom_name}.json"
    if os.path.exists(hp):
        h = json.load(open(hp))
        for var, v in h["variants"].items():
            out[f"nonhomologous_{var}"] = set(lnc) - set(v["homologous"])
    return out


def sn_in(hits, genes):
    return len(hits & genes) / len(genes) if genes else None


def mean(xs):
    xs = [x for x in xs if x is not None]
    return float(np.mean(xs)) if xs else None


def repro(cfg, cfg_sha):
    d, rep, ok = cfg["decode"], {}, True
    tgt = cfg["targets"]["chr18"]
    truth = ruler.read_gtf(f"{ROOT}/{tgt['truth']}", tgt["truth_chrom"])
    for name, ref in cfg["chr18_reference"].items():
        rd = f"{P15}/chr18/{name}"
        assert json.load(open(f"{rd}/read_meta.json"))["config_sha256"] == cfg_sha, f"{name}: read with another config"
        same_post = sha256(f"{rd}/post.npy") == ref["post_sha256"]
        got = core(load_genes(f"{rd}/genes.json"), np.load(f"{rd}/post.npy", mmap_mode="r"), truth, d["tau"])
        want = dict(ref["summary"])
        want["conf_grid"] = [x for x in want["conf_grid"] if x["tau"] == d["tau"]]
        diff = {k: (norm(got.get(k)), want[k]) for k in want if norm(got.get(k)) != want[k]}
        rep[name] = {"post_sha256_equal": same_post, "diff": diff}
        ok &= same_post and not diff
        print(json.dumps({name: rep[name]}), flush=True)
    out = {"config_sha256": cfg_sha, "ok": ok, "runs": rep, "time": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    json.dump(out, open(f"{P15}/chr18/REPRO.json", "w"), indent=1)
    if ok:
        json.dump(out, open(f"{P15}/chr18/REPRO_OK", "w"), indent=1)
    return out


def t2t(cfg):
    import p15_t2t_official as off
    d, tgt = cfg["decode"], cfg["targets"]["t2t_chr20"]
    truth = ruler.read_gff3_truth(f"{ROOT}/{tgt['truth']}", tgt["truth_chrom"])
    res = {"runs": {}, "GENATATOR": {}}
    g_fin, _ = score_preds(truth, ruler.read_pred_gff3(off.GENATATOR, tgt["truth_chrom"]))
    res["GENATATOR"] = {"final": g_fin, "official": off.official(off.GENATATOR)}
    for name in tgt["runs"]:
        rd = f"{P15}/t2t_chr20/{name}"
        r, _ = score_run(rd, truth, d["tau"])
        r["official"] = off.official(f"{rd}/pred.gff3")
        res["runs"][name] = r
    main = [res["runs"][n] for n in cfg["main"]["runs"].values() if n in res["runs"]]
    seg = lambda r, m: r["official"]["lnc_RNA|k250"]["segmentation-level"][m]
    res["mean"] = {"lnc_seg_recall_k250": mean([seg(r, "recall") for r in main]),
                   "lnc_seg_precision_k250": mean([seg(r, "precision") for r in main]),
                   "lnc_L1_sn": mean([r["final"]["lnc"]["L1_locus"]["sn"] for r in main]),
                   "lnc_L1_pr": mean([r["final"]["lnc"]["L1_locus"]["pr_strict"] for r in main])}
    res["genatator_seg_k250"] = {m: seg(res["GENATATOR"], m) for m in ("recall", "precision")}
    return res


def chr20(cfg):
    d, tgt = cfg["decode"], cfg["targets"]["chr20"]
    truth = ruler.read_gtf(f"{ROOT}/{tgt['truth']}", tgt["truth_chrom"])
    lnc_ids = sorted(k for k, v in truth.items() if v["cls"] == "lnc")
    st = strata(truth, f"{ROOT}/{tgt['truth']}", tgt["truth_chrom"])
    g_fin, g_hits = score_preds(truth, ruler.read_pred_gff3(G_CHR20, "chr20"))
    g_ph0 = json.load(open(G_CHR20_SCORE))["scores"]["new"]["lnc"]["L1_locus"]
    res = {"truth": {"n_lnc": len(lnc_ids), "strata_n": {k: len(v) for k, v in st.items()}},
           "GENATATOR": {"final": g_fin, "strata_sn": {k: sn_in(g_hits, v) for k, v in st.items()},
                         "phase0_lnc_L1": {"sn": g_ph0["sn"], "pr_strict": g_ph0["pr_strict"]}},
           "runs": {}, "ablations": {}}
    hits = {}
    for name in tgt["runs"]:
        rd = f"{P15}/chr20/{name}"
        if cfg["runs"][name]["arm"] == "G2PH":
            r = json.load(open(f"{rd}/result.json"))
            res["runs"][name] = {"own_evaluation": {k: (v if not isinstance(v, dict) else {kk: vv for kk, vv in v.items()
                                                                                          if kk != "curve"}) for k, v in r.items()}}
            continue
        r, h = score_run(rd, truth, d["tau"])
        r["strata_sn"] = {k: sn_in(h, v) for k, v in st.items()}
        res["runs"][name], hits[name] = r, h
        print(json.dumps({name: {"auprc": r["core"]["auprc_lnc"], "lnc_L1": [r["final"]["lnc"]["L1_locus"][k] for k in ("sn", "pr_strict")],
                                 "n_lnc": r["final"]["n_pred_lnc"]}}), flush=True)
    codes = np.load(f"{ROOT}/{tgt['codes']}", mmap_mode="r")
    for seed, name in cfg["main"]["runs"].items():  # AB1: S1 on the main runs' posteriors, tau*
        post = np.load(f"{P15}/chr20/{name}/post.npy", mmap_mode="r")
        gs = s1_genes(post, codes, d)
        fin, _ = score_preds(truth, final_preds(gs, d["tau"]))
        res["ablations"].setdefault("AB1", {})[seed] = {"run": name, "core": core(gs, post, truth, d["tau"]), "final": fin}
    main = {s: res["runs"][n] for s, n in cfg["main"]["runs"].items()}
    for arm, by_seed in cfg["arms"].items():
        for seed, name in by_seed.items():
            r = res["runs"].get(name)
            if r and "core" in r:
                res["ablations"].setdefault(arm, {})[seed] = {
                    "run": name, "auprc_lnc": r["core"]["auprc_lnc"], "delta_auprc_vs_main": r["core"]["auprc_lnc"] - main[seed]["core"]["auprc_lnc"],
                    "lnc_L1": [r["final"]["lnc"]["L1_locus"][k] for k in ("sn", "pr_strict")]}
    # verdict (main runs, 3-seed means)
    L1 = lambda r, k: r["final"]["lnc"]["L1_locus"][k]
    mn = {"lnc_L1_sn": mean([L1(r, "sn") for r in main.values()]), "lnc_L1_pr_strict": mean([L1(r, "pr_strict") for r in main.values()]),
          "lnc_L2_sn": mean([r["final"]["lnc"]["L2_intron_chain"]["sn"] for r in main.values()]),
          "mrna_L1_sn": mean([r["final"]["mrna"]["L1_locus"]["sn"] for r in main.values()]),
          "lnc_genes_reported_as_mrna": mean([r["final"]["guards"]["lnc_genes_reported_as_mrna"] for r in main.values()]),
          "mrna_preds_hitting_lnc_only": mean([r["final"]["guards"]["mrna_preds_hitting_lnc_only"] for r in main.values()]),
          "auprc_lnc": mean([r["core"]["auprc_lnc"] for r in main.values()])}
    for k in st:
        mn[f"sn_{k}"] = mean([r["strata_sn"][k] for r in main.values()])
    H = np.array([[g in hits[n] for g in lnc_ids] for n in cfg["main"]["runs"].values()], float)
    dvec = H.mean(0) - np.array([g in g_hits for g in lnc_ids], float)
    bi = np.random.default_rng(0).integers(0, len(lnc_ids), (BOOT, len(lnc_ids)))
    boots = dvec[bi].mean(1)
    gain = {"mean": float(dvec.mean()), "ci95": [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))]}
    gg = g_fin["guards"]
    tp = f"{P15}/t2t_chr20/SCORE.json"
    t2 = json.load(open(tp)) if os.path.exists(tp) else None
    go = {"lnc_L1_sn>=0.20": mn["lnc_L1_sn"] >= 0.20, "lnc_L1_pr_strict>=0.40": mn["lnc_L1_pr_strict"] >= 0.40,
          "sn_gain_vs_GENATATOR_ci_low>0": gain["ci95"][0] > 0}
    strong = {"lnc_L2_sn>=0.05": mn["lnc_L2_sn"] >= 0.05,
              "nonhomologous_masked_sn>=0.15": (mn.get("sn_nonhomologous_masked") or 0) >= 0.15,
              "t2t_seg_recall>GENATATOR": None if t2 is None else t2["mean"]["lnc_seg_recall_k250"] > t2["genatator_seg_k250"]["recall"],
              "t2t_seg_precision>=GENATATOR": None if t2 is None else t2["mean"]["lnc_seg_precision_k250"] >= t2["genatator_seg_k250"]["precision"]}
    guard = {"mrna_L1_sn>=0.80": mn["mrna_L1_sn"] >= 0.80,
             "lnc_genes_reported_as_mrna<GENATATOR": mn["lnc_genes_reported_as_mrna"] < gg["lnc_genes_reported_as_mrna"],
             "mrna_preds_hitting_lnc_only<GENATATOR": mn["mrna_preds_hitting_lnc_only"] < gg["mrna_preds_hitting_lnc_only"]}
    res["verdict"] = {"main": cfg["main"], "mean": mn, "sn_gain_vs_GENATATOR": gain, "go": go, "go_pass": all(go.values()),
                      "strong": strong, "strong_pass": all(go.values()) and all(v is True for v in strong.values()),
                      "guardrails": guard, "guardrails_pass": all(guard.values()),
                      "genatator": {"lnc_L1": [g_fin["lnc"]["L1_locus"][k] for k in ("sn", "pr_strict")], "guards": gg}}
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--target", required=True, choices=["chr18", "chr20", "t2t_chr20"])
    a = ap.parse_args()
    cfg_path = a.config if os.path.isabs(a.config) else f"{ROOT}/{a.config}"
    cfg, cfg_sha = json.load(open(cfg_path)), sha256(cfg_path)
    t0 = time.time()
    if a.target == "chr18":
        res = repro(cfg, cfg_sha)
        print(json.dumps({"repro_ok": res["ok"]}))
        sys.exit(0 if res["ok"] else 1)
    for name in cfg["targets"][a.target]["runs"]:
        m = json.load(open(f"{P15}/{a.target}/{name}/read_meta.json"))
        assert m["config_sha256"] == cfg_sha, f"{name}: read with another config"
    res = t2t(cfg) if a.target == "t2t_chr20" else chr20(cfg)
    res.update(config=a.config, config_sha256=cfg_sha, score_s=round(time.time() - t0))
    json.dump(res, open(f"{P15}/{a.target}/SCORE.json", "w"), indent=1, default=lambda x: sorted(x) if isinstance(x, set) else str(x))
    if a.target == "chr20":
        print(json.dumps({k: res["verdict"][k] for k in ("mean", "sn_gain_vs_GENATATOR", "go", "go_pass", "strong", "guardrails")}))
    else:
        print(json.dumps(res["mean"]))


if __name__ == "__main__":
    main()
