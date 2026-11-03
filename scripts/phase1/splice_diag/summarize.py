#!/usr/bin/env python3
"""chr18 splice diagnostic, step 3: one table of every decoding variant (splice_diag.py result.json) next to the
baseline (w = 0: the P1.5 batch-2 chr18 reads of the same runs, scored by the same functions without bootstrap).

  summarize.py [--dir outputs/PHASE1/SPLICE_DIAG] -> <dir>/SUMMARY.json and a markdown table on stdout
"""
import argparse
import glob
import json
import os
import sys

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1"))
import p15_score as P  # noqa: E402
from ncrna_annot.eval import ruler  # noqa: E402

RUNS = ("g1b-42-ke", "g1b-123-ke", "g1b-456-ke")
COLS = [("auprc_lnc", "AUPRC_lnc"), ("n_lnc", "lnc 数"), ("lnc_L1_sn", "lnc L1 Sn"), ("lnc_L1_pr", "lnc L1 严格 Pr"),
        ("lnc_L2_sn", "lnc L2 Sn"), ("lnc_L3_sn", "lnc L3 Sn"), ("lnc_L3_pr", "lnc L3 Pr"), ("lnc_L4_sn", "lnc L4 Sn"),
        ("n_mrna", "mRNA 数"), ("mrna_L1_sn", "mRNA L1 Sn"), ("mrna_L2_sn", "mRNA L2 Sn"), ("mrna_L3_sn", "mRNA L3 Sn"),
        ("mrna_L3_pr", "mRNA L3 Pr"), ("g1", "护栏 1"), ("g2", "护栏 2")]


def row(core, fin):
    sl, sm = fin["lnc"], fin["mrna"]
    return {"auprc_lnc": core["auprc_lnc"], "n_lnc": fin["n_pred_lnc"], "lnc_L1_sn": sl["L1_locus"]["sn"],
            "lnc_L1_pr": sl["L1_locus"]["pr_strict"], "lnc_L2_sn": sl["L2_intron_chain"]["sn"],
            "lnc_L3_sn": sl["L3_junction"]["recall"], "lnc_L3_pr": sl["L3_junction"]["precision"],
            "lnc_L4_sn": sl["L4_exon"]["recall"], "n_mrna": fin["n_pred_mrna"], "mrna_L1_sn": sm["L1_locus"]["sn"],
            "mrna_L2_sn": sm["L2_intron_chain"]["sn"], "mrna_L3_sn": sm["L3_junction"]["recall"],
            "mrna_L3_pr": sm["L3_junction"]["precision"], "g1": fin["guards"]["lnc_genes_reported_as_mrna"],
            "g2": fin["guards"]["mrna_preds_hitting_lnc_only"]}


def fmt(v):
    return "–" if v is None else (f"{v:.3f}" if isinstance(v, float) else str(v))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="outputs/PHASE1/SPLICE_DIAG")
    ap.add_argument("--config", default="configs/p15_v19b.json")
    a = ap.parse_args()
    d = f"{ROOT}/{a.dir}"
    cfg = json.load(open(f"{ROOT}/{a.config}"))
    tg, tau = cfg["targets"]["chr18"], cfg["decode"]["tau"]
    truth = ruler.read_gtf(f"{ROOT}/{tg['truth']}", tg["truth_chrom"])
    rows = {}
    for r in RUNS:
        rd = f"{ROOT}/{cfg['out_dir']}/chr18/{r}"
        genes = P.load_genes(f"{rd}/genes.json")
        fin, _ = P.score_preds(truth, P.final_preds(genes, tau), boot=0)
        rows[(r, "0")] = row(P.core(genes, np.load(f"{rd}/post.npy", mmap_mode="r"), truth, tau), fin)
    for f in sorted(glob.glob(f"{d}/*/result.json")):
        res = json.load(open(f))
        key = "oracle" if res["oracle"] else f"{res['w']:g}"
        if key == "0":
            continue
        rows[(res["run"], key)] = row(res["core"], res["final"])
    variants = sorted({k for _, k in rows}, key=lambda k: (k == "oracle", float(k) if k != "oracle" else 0))
    out = {"rows": {f"{r}|{k}": v for (r, k), v in rows.items()}, "mean": {}}
    lines = ["| 变体 | run | " + " | ".join(c for _, c in COLS) + " |", "|" + "---|" * (len(COLS) + 2)]
    for k in variants:
        for r in RUNS:
            if (r, k) in rows:
                lines.append(f"| w={k} | {r} | " + " | ".join(fmt(rows[(r, k)][c]) for c, _ in COLS) + " |")
        have = [rows[(r, k)] for r in RUNS if (r, k) in rows]
        if len(have) == 3:
            m = {c: P.mean([h[c] for h in have]) for c, _ in COLS}
            out["mean"][k] = m
            lines.append(f"| w={k} | **3 seed 平均** | " + " | ".join(fmt(m[c]) for c, _ in COLS) + " |")
    json.dump(out, open(f"{d}/SUMMARY.json", "w"), indent=1)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
