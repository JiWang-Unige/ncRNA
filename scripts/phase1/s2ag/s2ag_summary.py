#!/usr/bin/env python3
"""S2-AG summary (RESULT_S2_AG.md): AlphaGenome next to the splice models already read with s1_analyze.py on the same
chromosome, and the RNA-seq segment baseline next to our models.  Reads JSON outputs only; nothing is recomputed.

  s2ag_summary.py --out docs/phase1_results/S2AG/SUMMARY.json
"""
import argparse
import json
import os

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
O = f"{ROOT}/outputs/PHASE1"
S1_FILES = {  # chromosome -> s1_analyze outputs holding the other models
    "chr20": [f"{O}/S1/S1.json", f"{O}/S1/S1R_chr20_pangolin.json"],
    "chr18": [f"{O}/S1/S1R_chr18.json"],
    "chr9": [f"{O}/V20/rep/S1_chr9.json"],
    "chr7": [f"{O}/S1/S1R_chr7.json"],
}
LAYERS = ("1-9", "10-99", "100-999", ">=1000")


def layer_auroc(r1, cls, layer):  # as v20_replicate.py: positive-weighted mean over the layer's 5 sample strata
    st, n = r1[cls]["auroc_by_stratum"], r1[cls]["n_pos_by_stratum"]
    cells = [(st[layer * 5 + q], n[layer * 5 + q]) for q in range(5) if st[layer * 5 + q] is not None and n[layer * 5 + q] > 0]
    w = sum(c[1] for c in cells)
    return sum(a * k for a, k in cells) / w if w else None


def row(m):
    r1, r3 = m["R1"], m.get("R3", {})
    g = r1["gap"]
    return {"auroc_lnc": r1["lnc"]["auroc"], "auroc_mrna": r1["mrna"]["auroc"],
            "gap_unmatched": g["unmatched"], "gap_matched": g["matched"], "gap_matched_ci95": g["matched_ci95"],
            "layer_auroc_lnc": [layer_auroc(r1, "lnc", i) for i in range(4)],
            "layer_auroc_mrna": [layer_auroc(r1, "mrna", i) for i in range(4)],
            "topk": {k: r3[k]["topk_acc"] for k in ("mrna", "lnc_v41", "lnc_new") if k in r3},
            "auprc_r3": {k: r3[k]["auprc"] for k in ("mrna", "lnc_v41", "lnc_new") if k in r3}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    res = {"layers": LAYERS, "gap_sign": "mRNA AUROC - lnc AUROC (positive = lnc sites harder)", "splice": {}}
    for c, files in S1_FILES.items():
        models = {}
        for f in files:
            for name, m in json.load(open(f))["models"].items():
                models.setdefault(name, row(m))
        ag = json.load(open(f"{O}/S2AG/S1AG_{c}.json"))
        models["alphagenome"] = row(ag["models"]["alphagenome"])
        res["splice"][c] = {"models": models, "R4": ag["R4"]}
    res["offsets"] = json.load(open(f"{O}/S2AG/ag_chr20.npy.json"))["offsets"]
    res["site_ap"] = {c: {k: v for k, v in json.load(open(f"{O}/S2AG/ag_{c}.npy.json"))["ap"].items() if k in ("lnc_ap", "mrna_ap")}
                      for c in S1_FILES}
    seg = json.load(open(f"{O}/S2AG/SEGMENTS.json"))
    res["segments"] = {"grid_chr18": seg["grid_chr18"], "pick": seg["pick"], "pick_rule": seg["pick_rule"], "reads": {}}
    for c, r in seg["reads"].items():
        if isinstance(r, str):
            res["segments"]["reads"][c] = r
            continue
        l1 = r["lnc"]["L1_locus"]
        res["segments"]["reads"][c] = {"n_pred": r["n_pred_lnc"], "sn": l1["sn"], "sn_ci95": l1["sn_ci95"],
                                       "pr_strict": l1["pr_strict"], "pr_strict_ci95": l1["pr_strict_ci95"],
                                       "pr_lenient": l1["pr_lenient"], "n_novel_pred": l1["n_novel_pred"],
                                       "pred_hits_other_class_only": l1["pred_hits_other_class_only"]}
    os.makedirs(os.path.dirname(f"{ROOT}/{a.out}"), exist_ok=True)
    json.dump(res, open(f"{ROOT}/{a.out}", "w"), indent=1)
    f3 = lambda x: "—" if x is None else f"{x:.3f}"
    for c, d in res["splice"].items():
        print(f"== {c}")
        for name, r in d["models"].items():
            print(f"  {name:12s} AUROC lnc/mRNA {f3(r['auroc_lnc'])}/{f3(r['auroc_mrna'])}  gap matched {r['gap_matched']:+.3f} "
                  f"[{r['gap_matched_ci95'][0]:+.3f}, {r['gap_matched_ci95'][1]:+.3f}]  lnc layers "
                  + " / ".join(f3(x) for x in r["layer_auroc_lnc"]) + "  mRNA layers " + " / ".join(f3(x) for x in r["layer_auroc_mrna"])
                  + "  topk " + " / ".join(f3(r["topk"].get(k)) for k in ("mrna", "lnc_v41", "lnc_new")))
    print("segments pick", res["segments"]["pick"], res["segments"]["pick_rule"])
    for g in res["segments"]["grid_chr18"]:
        print("  grid", g)
    for c, r in res["segments"]["reads"].items():
        print("  read", c, r)


if __name__ == "__main__":
    main()
