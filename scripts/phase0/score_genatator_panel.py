#!/usr/bin/env python3
"""0.1a/0.1b: score our GENATATOR runs (Modal volume -> ont) on the panel with ruler v1.3.

Fetches outputs/PHASE0/0.1_genatator/<tag> from the Modal volume when missing locally.
Per species: lnc + mrna classes, lnc genes touched by any prediction, protocol 0.1b ratio
(lnc Sn / same-species mrna Sn) for L1 and L2.
"""
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ruler  # noqa: E402

ROOT = os.path.expanduser("~/ncrna-annot")
OUT = f"{ROOT}/outputs/PHASE0/0.1_genatator"
MODAL = os.path.expanduser("~/envs/modal/bin/modal")
PANEL = {"human_chr20": ("data/human/gencode.v50.annotation.gtf.gz", "chr20", "data/human/gencode.v41.annotation.gtf.gz"),
         "chicken_5": ("data/chicken/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.116.gtf.gz", "5", None),
         "zebrafish_22": ("data/zebrafish/Danio_rerio.GRCz11.116.gtf.gz", "22", None),
         "fly_2L": ("data/fly/Drosophila_melanogaster.BDGP6.54.63.gtf.gz", "2L", None),
         "worm_I": ("data/worm/Caenorhabditis_elegans.WBcel235.63.gtf.gz", "I", None)}
env = dict(os.environ, MODAL_PROFILE="ji-wang-unige")
summary = {}
for tag, (gtf, chrom, old) in PANEL.items():
    d = f"{OUT}/{tag}"
    if not os.path.exists(f"{d}/modal_run.json"):
        os.makedirs(OUT, exist_ok=True)
        subprocess.run([MODAL, "volume", "get", "ncrna-annot", f"/outputs/PHASE0/0.1_genatator/{tag}", OUT],
                       env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if not os.path.exists(f"{d}/pred.gff3"):
        summary[tag] = {"status": "not finished"}
        continue
    run = json.load(open(f"{d}/modal_run.json"))
    cmd = [sys.executable, f"{ROOT}/scripts/phase0/ruler.py", "score", "--truth", f"{ROOT}/{gtf}", "--chrom", chrom,
           "--pred", f"{d}/pred.gff3", "--boot", "1000", "--label", tag, "--out", f"{d}/score.json"]
    if old:
        cmd += ["--old-truth", f"{ROOT}/{old}"]
    subprocess.run(cmd, check=True, stderr=subprocess.DEVNULL)
    sc = json.load(open(f"{d}/score.json"))
    truth = ruler.read_gtf(f"{ROOT}/{gtf}", chrom)
    preds = ruler.read_pred_gff3(f"{d}/pred.gff3", chrom)
    il = ruler.ExonIndex({g: v for g, v in truth.items() if v["cls"] == "lnc"})
    seen = set()
    for p in preds:
        seen |= il.hits(p["strand"], p["exons"])
    s = sc["scores"]["new"]
    L = lambda c, m, k: s[c][m][k]
    r = {"status": "rc=%s" % run["rc"], "elapsed_h": round(run["elapsed_s"] / 3600, 2), "n_pred": sc["n_pred_tx"],
         "truth_genes": sc["truth_genes_by_class"],
         "lnc_L1_sn": L("lnc", "L1_locus", "sn"), "lnc_L1_sn_ci95": L("lnc", "L1_locus", "sn_ci95"),
         "lnc_L1_pr": L("lnc", "L1_locus", "pr_strict"), "lnc_L2_sn": L("lnc", "L2_intron_chain", "sn"),
         "mrna_L1_sn": L("mrna", "L1_locus", "sn"), "mrna_L1_pr": L("mrna", "L1_locus", "pr_strict"),
         "mrna_L2_sn": L("mrna", "L2_intron_chain", "sn"), "mrna_L2_pr": L("mrna", "L2_intron_chain", "pr_strict"),
         "lnc_genes_touched_by_any_pred": len(seen),
         "mrna_pred_hits_lnc_only": s["mrna"]["L1_locus"]["pred_hits_other_class_only"].get("lnc", 0)}
    for m in ("L1", "L2"):
        a, b = r[f"lnc_{m}_sn"], r[f"mrna_{m}_sn"]
        r[f"ratio_lnc_over_mrna_{m}"] = round(a / b, 4) if a is not None and b else None
    if old:
        r["delta_lnc_L1_sn"] = sc["delta"]["lnc"]["L1_locus"]["sn"]
        r["old_truth_lnc_L1_sn"] = sc["scores"]["old"]["lnc"]["L1_locus"]["sn"]
    summary[tag] = r
json.dump(summary, open(f"{OUT}/PANEL_SUMMARY.json", "w"), indent=1)
f = lambda x: "-" if x is None else (f"{x:.3f}" if isinstance(x, float) else str(x))
print("tag\tstatus\tlnc_genes\tlncL1sn\tlncL1pr\tlncL2sn\tlnc_seen\tmrnaL1sn\tmrnaL1pr\tmrnaL2sn\tratioL1\tn_pred")
for t, r in summary.items():
    if "lnc_L1_sn" not in r:
        print(t, r["status"])
        continue
    print("\t".join(f(x) for x in (t, r["status"], r["truth_genes"]["lnc"], r["lnc_L1_sn"], r["lnc_L1_pr"], r["lnc_L2_sn"],
                                   r["lnc_genes_touched_by_any_pred"], r["mrna_L1_sn"], r["mrna_L1_pr"], r["mrna_L2_sn"],
                                   r["ratio_lnc_over_mrna_L1"])) + f"\t{r['n_pred']}")
