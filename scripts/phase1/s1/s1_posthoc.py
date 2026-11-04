#!/usr/bin/env python3
"""S1 post hoc check (written after S1.json was read; descriptive, labelled post hoc in RESULT_S1.md):
the SpliceAI-paper style metric (R3) by usage.  Each multi-exon v50 transcript is binned by the gtexv2 coverage_sum
of its weakest annotated intron (exact junction match, 0 if absent): < 100, 100-999, >= 1,000 reads; R3 (top-k
accuracy and auPRC over all positions of the group's transcripts, other annotated sites removed) per class x bin.
Variant "clean": positions that are a site of any gtexv2 junction but not of the group's transcripts are also removed
from the candidates (real but unannotated sites are no longer counted as negatives).

  s1_posthoc.py --scores spliceai=outputs/PHASE1/S1/sai_chr20.npy own=outputs/PHASE1/S1/own_chr20.npy \
      vipsania=outputs/PHASE1/VIPSANIA/chr20/post6.npy:2 --out outputs/PHASE1/S1/S1_posthoc.json
"""
import argparse
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s1_analyze as A  # noqa: E402
from s1_analyze import ROOT, COL, ruler  # noqa: E402

BINS = ((0, 100), (100, 1000), (1000, 10 ** 12))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", default="chr20")
    ap.add_argument("--codes", default="data/p1/chr20_codes.npy")
    ap.add_argument("--snaptron", default="data/bioval/raw/snaptron/chr20_*.tsv")
    ap.add_argument("--v50", default="data/human/gencode.v50.annotation.gtf.gz")
    ap.add_argument("--scores", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    L = len(np.load(f"{ROOT}/{a.codes}", mmap_mode="r"))
    v50 = ruler.read_gtf(f"{ROOT}/{a.v50}", a.chrom)
    jn = A.read_snaptron(sorted(glob.glob(f"{ROOT}/{a.snaptron}")), a.chrom)
    jcov = {(int(s), int(e), st): int(c) for s, e, st, c in zip(jn["s"], jn["e"], jn["strand"], jn["cov"])}
    ann50 = A.annotated_keys(v50, L)
    used = []
    for strand in ("+", "-"):
        sel = (jn["strand"] == strand) | ~np.isin(jn["strand"], ("+", "-"))
        (lc, lp), (rc, rp) = A.junction_sites(jn["s"][sel], jn["e"][sel], np.full(sel.sum(), strand))
        used += [lc * L + lp, rc * L + rp]
    used = np.unique(np.concatenate(used))
    groups = {}
    for g in v50.values():
        if g["cls"] not in ("lnc", "mrna"):
            continue
        for ex in g["tx"].values():
            ex = ruler._merge(sorted(ex))
            if len(ex) < 2:
                continue
            u = min(jcov.get((s, e, g["strand"]), 0) for s, e in ruler.introns(ex))
            b = next(i for i, (lo, hi) in enumerate(BINS) if lo <= u < hi)
            groups.setdefault((g["cls"], b), []).append((g["strand"], ex))
    scores = [A.load_scores(s, L) for s in a.scores]
    res = {"post_hoc": True, "bins_reads": [list(b) for b in BINS], "models": {}}
    for name, arr, c0 in scores:
        out = {}
        for (cls, b), txs in sorted(groups.items()):
            span, posk = np.zeros((4, L), bool), []
            for strand, ex in txs:
                s0, e0 = max(ex[0][0] - 1, 0), min(ex[-1][1], L)
                for c in ((0, 1) if strand == "+" else (2, 3)):
                    span[c, s0:e0] = True
                for s, e in ruler.introns(ex):
                    posk += [COL[(strand, "L")] * L + s - 2, COL[(strand, "R")] * L + e]
            cand = np.concatenate([c * L + np.flatnonzero(span[c]) for c in range(4)])
            del span
            posk = np.unique(np.array(posk, np.int64))
            cand = cand[~np.isin(cand, np.setdiff1d(ann50, posk))]
            rec = {"n_tx": len(txs)}
            for variant, cc in (("", cand), ("clean_", cand[~np.isin(cand, np.setdiff1d(used, posk))])):
                lab = np.isin(cc, posk)
                sc = np.asarray(arr[cc % L, c0 + cc // L], np.float32)
                k = int(lab.sum())
                top = np.argpartition(-sc, k - 1)[:k]
                rec.update({f"{variant}n_sites": k, f"{variant}n_candidates": int(len(cc)),
                            f"{variant}topk_acc": float(lab[top].mean()), f"{variant}auprc": A.average_precision(sc, lab)})
            rec["n_used_unannotated_removed"] = rec["n_candidates"] - rec["clean_n_candidates"]
            out[f"{cls}_bin{b}"] = rec
        res["models"][name] = out
        print(json.dumps({name: {k: (v["n_tx"], v["n_sites"], v["n_used_unannotated_removed"], round(v["topk_acc"], 3), round(v["auprc"], 3),
                                     round(v["clean_topk_acc"], 3), round(v["clean_auprc"], 3)) for k, v in out.items()}}), flush=True)
    json.dump(res, open(f"{ROOT}/{a.out}", "w"), indent=1)


if __name__ == "__main__":
    main()
