#!/usr/bin/env python3
"""Score the Tiberius hg38 chr20 prediction (docs/phase1_results/PREREG_BASELINES.md) with the P1.5 ruler.

Every Tiberius transcript is one mRNA prediction; its exons are the exon rows of the GTF (the CDS rows when a
transcript has no exon rows).  No threshold.
  tiberius_score.py [--gtf outputs/PHASE1/BASELINES/tiberius/chr20.gtf] -> SCORE.json next to it
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
    ex, cds, strand = defaultdict(list), defaultdict(list), {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(f) < 9 or f[0] != "chr20" or f[2] not in ("exon", "CDS"):
            continue
        tid = f[8].split('transcript_id "', 1)[1].split('"', 1)[0]
        (ex if f[2] == "exon" else cds)[tid].append((int(f[3]), int(f[4])))
        strand[tid] = f[6]
    return [{"id": t, "cls": "mrna", "strand": strand[t], "exons": ruler._merge(sorted(ex.get(t) or cds[t]))}
            for t in sorted(strand)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", default="outputs/PHASE1/BASELINES/tiberius/chr20.gtf")
    a = ap.parse_args()
    preds = read(f"{ROOT}/{a.gtf}")
    truth = ruler.read_gtf(f"{ROOT}/data/human/gencode.v50.annotation.gtf.gz", "chr20")
    fin, _ = P.score_preds(truth, preds)
    json.dump(P.norm({"n_tx": len(preds), "final": fin}), open(os.path.join(os.path.dirname(f"{ROOT}/{a.gtf}"), "SCORE.json"), "w"), indent=1)
    sm = fin["mrna"]
    print(json.dumps({"n_tx": len(preds), "mrna_L1": [sm["L1_locus"]["sn"], sm["L1_locus"]["pr_strict"]],
                      "mrna_L2_sn": sm["L2_intron_chain"]["sn"], "guards": fin["guards"]}))


if __name__ == "__main__":
    main()
