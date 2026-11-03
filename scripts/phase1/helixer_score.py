#!/usr/bin/env python3
"""Rescore the Phase 0 Helixer hg38 chr20 prediction (vertebrate_v0.3, not re-run; docs/phase1_results/PREREG_BASELINES.md
section 3) with the P1.5 ruler, so that it sits in the same table as Tiberius and SegmentNT.  Human is in Helixer's
training / validation data: optimistic.

Every Helixer mRNA is one mRNA prediction; its exons are the exon rows (UTR included).  No threshold.
  helixer_score.py [--gff data/baselines/human/helixer_vertebrate_v0.3_chr20.gff3]
                   -> outputs/PHASE1/BASELINES/helixer/SCORE.json
"""
import argparse
import json
import os
import sys
from collections import defaultdict

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1"))
import p15_score as P  # noqa: E402
from ncrna_annot.eval import ruler  # noqa: E402


def read(path):
    ex, strand = defaultdict(list), {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(f) < 9 or f[0] != "chr20":
            continue
        attr = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        if f[2] == "mRNA":
            strand[attr["ID"]] = f[6]
        elif f[2] == "exon":
            ex[attr["Parent"]].append((int(f[3]), int(f[4])))
    return [{"id": t, "cls": "mrna", "strand": strand[t], "exons": ruler._merge(sorted(ex[t]))} for t in sorted(strand)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", default="data/baselines/human/helixer_vertebrate_v0.3_chr20.gff3")
    a = ap.parse_args()
    preds = read(f"{ROOT}/{a.gff}")
    truth = ruler.read_gtf(f"{ROOT}/data/human/gencode.v50.annotation.gtf.gz", "chr20")
    fin, _ = P.score_preds(truth, preds)
    out = f"{ROOT}/outputs/PHASE1/BASELINES/helixer"
    os.makedirs(out, exist_ok=True)
    json.dump(P.norm({"n_tx": len(preds), "gff": a.gff, "final": fin}), open(f"{out}/SCORE.json", "w"), indent=1)
    sm = fin["mrna"]
    print(json.dumps({"n_tx": len(preds), "mrna_L1": [sm["L1_locus"]["sn"], sm["L1_locus"]["pr_strict"]],
                      "mrna_L2_sn": sm["L2_intron_chain"]["sn"], "guards": fin["guards"]}))


if __name__ == "__main__":
    main()
