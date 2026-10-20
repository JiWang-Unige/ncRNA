#!/usr/bin/env python3
"""0.6: score small-RNA predictions on the panel; recall split by truth biotype (descriptive)."""
import json
import os
import subprocess
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ruler  # noqa: E402

ROOT = os.path.expanduser("~/ncrna-annot")
PANEL = {"human_chr20": ("data/human/gencode.v50.annotation.gtf.gz", "chr20"),
         "chicken_5": ("data/chicken/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.116.gtf.gz", "5"),
         "zebrafish_22": ("data/zebrafish/Danio_rerio.GRCz11.116.gtf.gz", "22"),
         "fly_2L": ("data/fly/Drosophila_melanogaster.BDGP6.54.63.gtf.gz", "2L"),
         "worm_I": ("data/worm/Caenorhabditis_elegans.WBcel235.63.gtf.gz", "I")}
OUT = f"{ROOT}/outputs/PHASE0/0.6_smallrna"
types = dict(l.rstrip("\n").split("\t", 1) for l in open(f"{ROOT}/data/rfam/rfam_types.tsv") if "\t" in l)
summary = {}
for tag, (gtf, chrom) in PANEL.items():
    d = f"{OUT}/{tag}"
    subprocess.run([sys.executable, f"{ROOT}/scripts/phase0/ruler.py", "score", "--truth", f"{ROOT}/{gtf}",
                    "--chrom", chrom, "--pred", d, "--pred-format", "smallrna", "--rfam-types",
                    f"{ROOT}/data/rfam/rfam_types.tsv", "--classes", "small", "--boot", "1000",
                    "--label", tag, "--out", f"{d}/score.json"], check=True)
    sc = json.load(open(f"{d}/score.json"))
    truth = ruler.read_gtf(f"{ROOT}/{gtf}", chrom)
    preds = ruler.read_pred_smallrna(d, chrom, types)
    idx = ruler.ExonIndex({g: t for g, t in truth.items() if t["cls"] == "small"})
    hit = set()
    for p in preds:
        hit |= idx.hits(p["strand"], p["exons"])
    by = Counter(t["gene_type"] for t in truth.values() if t["cls"] == "small")
    byhit = Counter(truth[g]["gene_type"] for g in hit)
    run = json.load(open(f"{d}/run.json"))
    s = sc["scores"]["new"]["small"]
    summary[tag] = {"n_truth": s["n_truth_genes"], "n_pred": s["n_pred_tx"],
                    "L1_sn": s["L1_locus"]["sn"], "L1_sn_ci95": s["L1_locus"]["sn_ci95"],
                    "L1_pr_strict": s["L1_locus"]["pr_strict"], "L1_pr_lenient": s["L1_locus"]["pr_lenient"],
                    "n_novel_pred": s["L1_locus"]["n_novel_pred"],
                    "L5_recall": s["L5_single_exon"]["recall"], "L5_precision": s["L5_single_exon"]["precision"],
                    "recall_by_biotype": {k: [byhit[k], v] for k, v in by.most_common()},
                    "pred_by_source": Counter(p["src"].split("|")[0] if p["src"].startswith("tRNA") else
                                              (p["src"].split("|")[1].split(";")[1].strip() if ";" in p["src"] else "?")
                                              for p in preds).most_common(8),
                    "cmscan_cpu_h": run["cmscan_cpu_h"], "bp": run["bp_scanned"]}
json.dump(summary, open(f"{OUT}/PANEL_SUMMARY.json", "w"), indent=1)
for k, v in summary.items():
    print(k, json.dumps({x: v[x] for x in ("n_truth", "n_pred", "L1_sn", "L1_pr_strict", "n_novel_pred", "cmscan_cpu_h")}))
    print("   by biotype [hit, total]:", v["recall_by_biotype"])
