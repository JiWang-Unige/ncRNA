#!/usr/bin/env python3
"""T2T chr20 (NC_060944.1): the GENATATOR leaderboard's own metric engine on one prediction GFF3 (protocol section 3,
Strong: lnc_RNA segmentation-level recall / precision at k = 250).

Engine: data/baselines/t2t_chr20/metric (HF Space AIRI-Institute/genatator-ab-initio-annotation-leaderboard,
revision d32d4dc), called as the Space does (exon branch, stranded) but one transcript type at a time, as in the
Phase 0 check (outputs/PHASE0/0.3_fp_audit/t2t_chr20/official_metric_lnc_vs_mrna.json).

  p15_t2t_official.py --check                 reproduce the Phase 0 GENATATOR numbers exactly
  p15_t2t_official.py --pred x.gff3 --out x.json
"""
import argparse
import json
import os
import sys

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
B = f"{ROOT}/data/baselines/t2t_chr20"
TRUTH = f"{B}/leaderboard_data/ground_truth/chr20.gff"
GENATATOR = f"{B}/preds/predictions/Hs_NC_060944.1_filtered_colored.gff"
PHASE0 = f"{ROOT}/outputs/PHASE0/0.3_fp_audit/t2t_chr20/official_metric_lnc_vs_mrna.json"
TYPES = ("lnc_RNA", "mRNA")
KS = (0, 250)


def official(pred, truth=TRUTH, types=TYPES, ks=KS):
    sys.path.insert(0, f"{B}/metric")
    from backend.gff_io import gff_text_to_dataframe
    from gene_level_final_final_fix import GeneLevelEvaluator
    ev = GeneLevelEvaluator()
    p_df, t_df = gff_text_to_dataframe(open(pred).read()), gff_text_to_dataframe(open(truth).read())
    out = {}
    for t in types:
        r = ev.evaluate_gff_exon(pred_gff=p_df, true_gff=t_df, k_values=list(ks), use_strand=True, transcript_types=[t])
        for k in ks:
            out[f"{t}|k{k}"] = {lv: {m: r[k][lv][m] for m in ("precision", "recall", "f1")}
                                for lv in ("interval-level", "segmentation-level")}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--pred")
    ap.add_argument("--out")
    a = ap.parse_args()
    if a.check:
        ref = json.load(open(PHASE0))
        got = official(GENATATOR)
        diff = {k: (v, ref[f"GENATATOR|{k}"]) for k, v in got.items() if v != ref[f"GENATATOR|{k}"]}
        print(json.dumps({"reproduced": not diff, "diff": diff, "lnc_RNA|k250": got["lnc_RNA|k250"]}))
        sys.exit(1 if diff else 0)
    res = official(a.pred)
    json.dump(res, open(a.out, "w"), indent=1)
    print(json.dumps({"lnc_RNA|k250": res["lnc_RNA|k250"]}))


if __name__ == "__main__":
    main()
