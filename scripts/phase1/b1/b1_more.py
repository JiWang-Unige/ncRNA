#!/usr/bin/env python3
"""B1 extension (manuscript section L, 2026-11-14): strict one-to-one span matching (b1_strict.score, S50 / S25) for the
other already-read sets - T2T chr20 (P1.5), mouse chr19 (S5), the P1.6 zero-shot panel - and GENATATOR on the same
chromosomes (Phase 0 predictions).  Existing predictions only (final lnc conf >= tau 0.6, all mRNA); CPU.

  b1_more.py [--out outputs/PHASE1/B1]   -> <out>/B1_MORE.json
"""
import argparse
import json
import os
import sys
import time

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
for d in ("scripts/phase1", "scripts/phase1/b1"):
    sys.path.insert(0, os.path.join(ROOT, d))
import b1_strict as B  # noqa: E402
import p15_score as P  # noqa: E402
import p15_t2t_official as OFF  # noqa: E402
from ncrna_annot.eval import ruler  # noqa: E402

G0 = f"{ROOT}/outputs/PHASE0/0.1_genatator"


def runs_of(cfg_path, target, reader):
    cfg = json.load(open(f"{ROOT}/{cfg_path}"))
    tg, tau = cfg["targets"][target], cfg["decode"]["tau"]
    truth = reader(f"{ROOT}/{tg['truth']}", tg["truth_chrom"])
    runs = {}
    for run in cfg["main"]["runs"].values():
        genes = P.load_genes(f"{ROOT}/{cfg['out_dir']}/{target}/{run}/genes.json")
        runs[run] = B.score(truth, P.final_preds(genes, tau))
    return truth, tg, {"tau": tau, "runs": runs, "mean": {c: B.mean_runs(list(runs.values()), c) for c in ("lnc", "mrna")}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="outputs/PHASE1/B1")
    a = ap.parse_args()
    t0 = time.time()
    res = {"plan": "docs/phase1_results/PLAN_B1.md (definitions); manuscript section L", "sets": {}}
    truth, tg, r = runs_of("configs/p15_v19.json", "t2t_chr20", ruler.read_gff3_truth)
    res["sets"]["v19_t2t_chr20"] = r
    res["sets"]["genatator_t2t_chr20"] = B.score(truth, B.gff_preds(OFF.GENATATOR, tg["truth_chrom"]))
    truth, tg, r = runs_of("configs/p17_mouse.json", "mouse_chr19", ruler.read_gtf)
    res["sets"]["v19_mouse_chr19"] = r
    for tgt in ("chicken_5", "zebrafish_22", "fly_2L", "worm_I"):
        truth, tg, r = runs_of("configs/p16_v19.json", tgt, ruler.read_gtf)
        res["sets"][f"v19_{tgt}"] = r
        gp = f"{G0}/{tgt}/pred.gff3"
        if os.path.exists(gp):
            res["sets"][f"genatator_{tgt}"] = B.score(truth, B.gff_preds(gp, tg["seqid"]))
    for k, v in res["sets"].items():
        m = v.get("mean") or {c: {f"{key}_{f}": v[c][key].get(f) for key in ("S50", "S25", "L1") for f in ("sn", "pr", "pr_strict") if v[c][key].get(f) is not None} for c in ("lnc", "mrna")}
        print(k, json.dumps({c: {x: round(y, 3) for x, y in m[c].items()} for c in m}), flush=True)
    res["seconds"] = round(time.time() - t0, 1)
    json.dump(P.norm(res), open(f"{ROOT}/{a.out}/B1_MORE.json", "w"), indent=1)


if __name__ == "__main__":
    main()
