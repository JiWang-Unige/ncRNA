#!/usr/bin/env python3
"""P1.6 zero-shot panel scoring (CPU, descriptive, no gate) of the p16_read.py outputs, and the unified GFF3.

  p16_score.py --config configs/p16_v19.json

Per panel chromosome and main run (as P1.5 chr20, p15_score.py): core = screen2.summarize on all MEAX genes
(AUPRC_lnc), final = ruler v1.3 L1-L6 for lnc and mRNA of the final predictions (lnc confidence >= tau*, every mRNA)
with 2,000-resample CIs, guardrail counts, ECE; 3-seed means; GENATATOR on the same chromosome from Phase 0 0.1a
(outputs/PHASE0/0.1_genatator/<target>/score.json) for reference.  Writes <out_dir>/<target>/SCORE.json and
<out_dir>/PANEL_SCORE.json (with this script's sha256).
Unified GFF3 (protocol v1 section 2; Phase 0 0.6 integration rules): <out_dir>/<target>/<run>/unified.gff3 =
pred.gff3 + the 0.6 Rfam / tRNAscan hits as ncRNA genes (overlapping clan losers and Rfam lncRNA / Cis-reg / Intron
families dropped, as ruler.read_pred_smallrna; Metazoa_SRP and U6 hits kept only when full length, trunc "no"; no
RepeatMasker filter on the panel, no miRNA component); a small RNA that overlaps a predicted lnc exon on the same
strand keeps both and carries Note=overlaps_lnc:<gene id>.
"""
import argparse
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p15_score as P  # noqa: E402
from ncrna_annot.eval import ruler  # noqa: E402

SMALL = f"{ROOT}/outputs/PHASE0/0.6_smallrna"
GENATATOR = f"{ROOT}/outputs/PHASE0/0.1_genatator"
RFAM_TYPES = f"{ROOT}/data/rfam/rfam_types.tsv"
REPEAT_PRONE = ("Metazoa_SRP", "U6")


def small_hits(target, seqid):
    if not os.path.isdir(f"{SMALL}/{target}"):   # no Phase 0 small-RNA layer for this target (S5 mouse)
        return []
    types = dict(l.rstrip("\n").split("\t", 1) for l in open(RFAM_TYPES) if "\t" in l)
    trunc = {}
    for i, line in enumerate(open(f"{SMALL}/{target}/cmscan.tsv")):
        f = line.rstrip("\n").split("\t")
        if i:
            trunc[f"rfam{i}"] = f[10]
    out = []
    for p in ruler.read_pred_smallrna(f"{SMALL}/{target}", seqid, types):
        fam = p["src"].split("|")[0]
        if fam in REPEAT_PRONE and trunc.get(p["id"]) != "no":
            continue
        out.append(p)
    return out


def unified_gff3(run_dir, seqid, small, out):
    lines = [l for l in open(f"{run_dir}/pred.gff3") if not l.startswith("#")]
    gene_kind = {}
    for l in lines:
        f = l.rstrip("\n").split("\t")
        if f[2] in ("lnc_RNA", "mRNA"):
            gene_kind[f[8].split("Parent=")[1].split(";")[0]] = f[2]
    lnc_ex = []
    for l in lines:
        f = l.rstrip("\n").split("\t")
        if f[2] == "exon":
            gid = f[8].split("Parent=")[1].split(".t1")[0]
            if gene_kind.get(gid) == "lnc_RNA":
                lnc_ex.append((f[6], int(f[3]), int(f[4]), gid))
    with open(out, "w") as fh:
        fh.write("##gff-version 3\n")
        fh.writelines(lines)
        for i, p in enumerate(small):
            (a, b), st = p["exons"][0], p["strand"]
            gid, src = f"ncrna_annot_small_{seqid}_{i}", p["src"].replace(";", ",").replace("=", ":")
            ov = sorted({g for s_, x, y, g in lnc_ex if s_ == st and x <= b and a <= y})
            note = f";Note=overlaps_lnc:{','.join(ov)}" if ov else ""
            fh.write(f"{seqid}\tncrna_annot_small\tncRNA_gene\t{a}\t{b}\t.\t{st}\t.\tID={gid};source_hit={src}{note}\n")
            fh.write(f"{seqid}\tncrna_annot_small\tncRNA\t{a}\t{b}\t.\t{st}\t.\tID={gid}.t1;Parent={gid}\n")
            fh.write(f"{seqid}\tncrna_annot_small\texon\t{a}\t{b}\t.\t{st}\t.\tID={gid}.t1.e0;Parent={gid}.t1\n")
    return len(small)


def genatator(target):
    p = f"{GENATATOR}/{target}/score.json"
    if not os.path.exists(p):
        return None
    s = json.load(open(p))["scores"]["new"]
    pick = lambda c: {"n_pred": s[c]["n_pred_tx"], "L1_sn": s[c]["L1_locus"]["sn"], "L1_pr_strict": s[c]["L1_locus"]["pr_strict"],
                      "L2_sn": s[c]["L2_intron_chain"]["sn"], "L3_recall": s[c]["L3_junction"]["recall"]}
    return {"lnc": pick("lnc"), "mrna": pick("mrna")}


def brief(r):
    fl, fm = r["final"]["lnc"], r["final"]["mrna"]
    return {"auprc_lnc": r["core"]["auprc_lnc"], "n_truth_lnc": fl["n_truth_genes"], "n_pred_lnc": r["final"]["n_pred_lnc"],
            "lnc_L1_sn": fl["L1_locus"]["sn"], "lnc_L1_pr_strict": fl["L1_locus"]["pr_strict"],
            "lnc_L2_sn": fl["L2_intron_chain"]["sn"], "lnc_L3_recall": fl["L3_junction"]["recall"],
            "lnc_L3_precision": fl["L3_junction"]["precision"], "n_truth_mrna": fm["n_truth_genes"],
            "n_pred_mrna": r["final"]["n_pred_mrna"], "mrna_L1_sn": fm["L1_locus"]["sn"],
            "mrna_L1_pr_strict": fm["L1_locus"]["pr_strict"], "mrna_L2_sn": fm["L2_intron_chain"]["sn"],
            "mrna_L3_recall": fm["L3_junction"]["recall"], **r["final"]["guards"], "ece": (r["ece"] or {}).get("ece")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    a = ap.parse_args()
    cfg_path = a.config if os.path.isabs(a.config) else f"{ROOT}/{a.config}"
    cfg, cfg_sha = json.load(open(cfg_path)), P.sha256(cfg_path)
    out_dir, tau = f"{ROOT}/{cfg['out_dir']}", cfg["decode"]["tau"]
    panel = {}
    for target, tg in cfg["targets"].items():
        if target == "chr18":
            continue
        t0 = time.time()
        truth = ruler.read_gtf(f"{ROOT}/{tg['truth']}", tg["truth_chrom"])
        small = small_hits(target, tg["seqid"])
        res = {"runs": {}, "truth_genes_by_class": {c: sum(g["cls"] == c for g in truth.values()) for c in ("lnc", "mrna", "small")}}
        for name in tg["runs"]:
            rd = f"{out_dir}/{target}/{name}"
            m = json.load(open(f"{rd}/read_meta.json"))
            assert m["config_sha256"] == cfg_sha, f"{target} {name}: read with another config"
            r, _ = P.score_run(rd, truth, tau)
            r["n_small_in_unified_gff3"] = unified_gff3(rd, tg["seqid"], small, f"{rd}/unified.gff3")
            res["runs"][name] = r
        rows = {n: brief(r) for n, r in res["runs"].items()}
        res["brief"] = rows
        res["mean"] = {k: P.mean([v[k] for v in rows.values()]) for k in next(iter(rows.values()))}
        res["genatator"] = genatator(target)
        res.update(config=a.config, config_sha256=cfg_sha, score_s=round(time.time() - t0))
        json.dump(res, open(f"{out_dir}/{target}/SCORE.json", "w"), indent=1, default=lambda x: sorted(x) if isinstance(x, set) else str(x))
        panel[target] = {"mean": res["mean"], "by_run": rows, "genatator": res["genatator"],
                         "truth_genes_by_class": res["truth_genes_by_class"], "n_small": len(small)}
        print(json.dumps({target: {k: (round(v, 4) if isinstance(v, float) else v) for k, v in res["mean"].items()}}), flush=True)
    json.dump({"config": a.config, "config_sha256": cfg_sha, "score_script_sha256": P.sha256(os.path.abspath(__file__)), "panel": panel}, open(f"{out_dir}/PANEL_SCORE.json", "w"), indent=1)


if __name__ == "__main__":
    main()
