#!/usr/bin/env python3
"""Score the official GENATATOR-leaderboard predictions on T2T chr20 (NC_060944.1) with ruler v1.3.

Truth: the leaderboard's own RefSeq chr20.gff (Space revision recorded in data/baselines/t2t_chr20).
Also counts lncRNA genes touched by ANY prediction (any type), i.e. what a tool "sees".
"""
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ruler  # noqa: E402

ROOT = os.path.expanduser("~/ncrna-annot")
B = f"{ROOT}/data/baselines/t2t_chr20"
TRUTH = f"{B}/leaderboard_data/ground_truth/chr20.gff"
OUT = f"{ROOT}/outputs/PHASE0/0.3_fp_audit/t2t_chr20"
CHROM = "NC_060944.1"
names = json.load(open(f"{B}/preds/model_name_mapping.json"))
UNSTRANDED = ("segment_nt", "NTv3")
os.makedirs(OUT, exist_ok=True)

truth = ruler.read_gff3_truth(TRUTH, CHROM)
lnc = {g: v for g, v in truth.items() if v["cls"] == "lnc"}
pc = {g: v for g, v in truth.items() if v["cls"] == "mrna"}
rows = {}
for fname, model in names.items():
    unstr = fname.startswith(UNSTRANDED)
    out = f"{OUT}/{fname.rsplit('.', 1)[0]}.json"
    cmd = [sys.executable, f"{ROOT}/scripts/phase0/ruler.py", "score", "--truth", TRUTH,
           "--truth-format", "gff3", "--chrom", CHROM, "--pred", f"{B}/preds/predictions/{fname}",
           "--boot", "1000", "--label", model, "--out", out] + (["--ignore-strand"] if unstr else [])
    subprocess.run(cmd, check=True, stderr=subprocess.DEVNULL)
    d = json.load(open(out))
    preds = ruler.read_pred_gff3(f"{B}/preds/predictions/{fname}", CHROM)
    t_l, t_p = ({k: dict(v, strand=".") for k, v in x.items()} if unstr else x for x in (lnc, pc))
    il, ip = ruler.ExonIndex(t_l), ruler.ExonIndex(t_p)
    seen, as_coding = set(), set()
    for p in preds:
        st = "." if unstr else p["strand"]
        h = il.hits(st, p["exons"])
        seen |= h
        if p["cls"] == "mrna" and h and not ip.hits(st, p["exons"]):
            as_coding |= h
    s = d["scores"]["new"]
    rows[model] = {
        "file": fname, "stranded": not unstr, "n_pred": d["n_pred_tx"],
        "lnc_L1_sn": s["lnc"]["L1_locus"]["sn"], "lnc_L1_sn_ci95": s["lnc"]["L1_locus"]["sn_ci95"],
        "lnc_L1_pr_strict": s["lnc"]["L1_locus"]["pr_strict"], "lnc_L2_sn": s["lnc"]["L2_intron_chain"]["sn"],
        "lnc_L3_recall": s["lnc"]["L3_junction"]["recall"], "lnc_L6_f1": s["lnc"]["L6_base"]["f1"],
        "mrna_L1_sn": s["mrna"]["L1_locus"]["sn"], "mrna_L1_pr_strict": s["mrna"]["L1_locus"]["pr_strict"],
        "mrna_L2_sn": s["mrna"]["L2_intron_chain"]["sn"], "mrna_L2_pr": s["mrna"]["L2_intron_chain"]["pr_strict"],
        "mrna_pred_hits_lnc_only": s["mrna"]["L1_locus"]["pred_hits_other_class_only"].get("lnc", 0),
        "lnc_genes_touched_by_any_pred": len(seen), "lnc_genes_reported_as_coding": len(as_coding)}
summary = {"truth": TRUTH, "chrom": CHROM, "n_lnc_genes": len(lnc), "n_mrna_genes": len(pc), "models": rows}
json.dump(summary, open(f"{OUT}/T2T_PANEL_SUMMARY.json", "w"), indent=1)
print("truth lnc genes", len(lnc), "mrna genes", len(pc))
f = lambda x: "-" if x is None else f"{x:.3f}"
print("model\tstr\tlncL1sn\tlncL1pr\tlncL2sn\tlnc_seen\tlnc_as_coding\tmrnaL1sn\tmrnaL1pr\tmrnaL2sn\tmrnaL2pr\tn_pred")
for m, r in rows.items():
    print(f"{m}\t{int(r['stranded'])}\t{f(r['lnc_L1_sn'])}\t{f(r['lnc_L1_pr_strict'])}\t{f(r['lnc_L2_sn'])}\t"
          f"{r['lnc_genes_touched_by_any_pred']}\t{r['lnc_genes_reported_as_coding']}\t{f(r['mrna_L1_sn'])}\t"
          f"{f(r['mrna_L1_pr_strict'])}\t{f(r['mrna_L2_sn'])}\t{f(r['mrna_L2_pr'])}\t{r['n_pred']}")
