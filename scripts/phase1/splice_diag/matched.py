#!/usr/bin/env python3
"""chr18 splice diagnostic, step 4 (POST HOC, written after the pre-registered readout of 2026-11-03): the splice
potential changes how many lnc genes clear tau* = 0.6, so compare at a matched operating point instead.

For every variant and run, the lnc genes (all MEAX genes >= 200 nt, genes.json) are ranked by confidence and the top
n are kept, n = the number the same run keeps at tau* without the potential (P1.5 batch-2 chr18 read); mRNA genes
are all kept, as in the final predictions.  Ruler lnc / mRNA L1-L4 without bootstrap.

  matched.py [--dir outputs/PHASE1/SPLICE_DIAG] -> <dir>/MATCHED.json and a markdown table on stdout
"""
import argparse
import glob
import json
import os
import sys

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1"))
import p15_score as P  # noqa: E402
from ncrna_annot.eval import ruler  # noqa: E402

RUNS = ("g1b-42-ke", "g1b-123-ke", "g1b-456-ke")
COLS = [("n_lnc", "lnc 数"), ("conf_cut", "置信度下限"), ("lnc_L1_sn", "lnc L1 Sn"), ("lnc_L1_pr", "lnc L1 严格 Pr"),
        ("lnc_L2_sn", "lnc L2 Sn"), ("lnc_L3_sn", "lnc L3 Sn"), ("lnc_L3_pr", "lnc L3 Pr"), ("lnc_L4_sn", "lnc L4 Sn"),
        ("mrna_L2_sn", "mRNA L2 Sn"), ("mrna_L3_sn", "mRNA L3 Sn"), ("mrna_L3_pr", "mRNA L3 Pr")]


def top_n(genes, n):
    lnc = sorted([g for g in genes if g["type"] == "l"], key=lambda g: -g["conf"])[:n]
    keep = [g for g in genes if g["type"] != "l"] + lnc
    preds = [{"id": str(i), "strand": g["strand"], "cls": "lnc" if g["type"] == "l" else "mrna",
              "exons": [tuple(x) for x in g["exons"]]} for i, g in enumerate(keep)]
    return preds, (lnc[-1]["conf"] if lnc else None)


def row(truth, preds, cut):
    fin, _ = P.score_preds(truth, preds, boot=0)
    sl, sm = fin["lnc"], fin["mrna"]
    return {"n_lnc": fin["n_pred_lnc"], "conf_cut": cut, "lnc_L1_sn": sl["L1_locus"]["sn"], "lnc_L1_pr": sl["L1_locus"]["pr_strict"],
            "lnc_L2_sn": sl["L2_intron_chain"]["sn"], "lnc_L3_sn": sl["L3_junction"]["recall"],
            "lnc_L3_pr": sl["L3_junction"]["precision"], "lnc_L4_sn": sl["L4_exon"]["recall"],
            "mrna_L2_sn": sm["L2_intron_chain"]["sn"], "mrna_L3_sn": sm["L3_junction"]["recall"],
            "mrna_L3_pr": sm["L3_junction"]["precision"]}


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
    base, nref, rows = {}, {}, {}
    for r in RUNS:
        g = P.load_genes(f"{ROOT}/{cfg['out_dir']}/chr18/{r}/genes.json")
        nref[r] = sum(1 for x in g if x["type"] == "l" and x["conf"] >= tau)
        rows[(r, "0")] = row(truth, *top_n(g, nref[r]))
    for f in sorted(glob.glob(f"{d}/*/genes.json")):
        res = json.load(open(os.path.dirname(f) + "/result.json")) if os.path.exists(os.path.dirname(f) + "/result.json") else None
        if res is None:
            continue
        key = "oracle" if res["oracle"] else f"{res['w']:g}"
        if key == "0":
            continue
        rows[(res["run"], key)] = row(truth, *top_n(P.load_genes(f), nref[res["run"]]))
    variants = sorted({k for _, k in rows}, key=lambda k: (k == "oracle", float(k) if k != "oracle" else 0))
    out = {"n_ref": nref, "rows": {f"{r}|{k}": v for (r, k), v in rows.items()}, "mean": {}}
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
    json.dump(out, open(f"{d}/MATCHED.json", "w"), indent=1)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
