#!/usr/bin/env python3
"""chr18 splice diagnostic, step 5 (POST HOC, descriptive): where lnc structure is lost, baseline vs SpliceAI w vs
truth-site oracle (seed 42 unless --run).  Final lnc predictions (confidence >= tau*) only.

  per predicted lnc gene: exon count; single-exon fraction
  per truth lnc gene detected at L1 (same-strand exon overlap): has >= 1 intron predicted exactly (any transcript),
    full intron chain of some transcript (L2), and the truth exon count of its longest transcript
  per predicted intron of predicted lnc genes that hit truth lnc: exact / donor-only / acceptor-only / neither
"""
import argparse
import json
import os
import sys
from collections import Counter

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1"))
import p15_score as P  # noqa: E402
from ncrna_annot.eval import ruler  # noqa: E402


def analyse(genes, truth, tau):
    lnc_t = {k: v for k, v in truth.items() if v["cls"] == "lnc"}
    idx = ruler.ExonIndex(lnc_t)
    preds = [g for g in genes if g["type"] == "l" and g["conf"] >= tau]
    t_introns = {gid: {i for ex in g["tx"].values() for i in ruler.introns(ruler._merge(ex))} for gid, g in lnc_t.items()}
    t_chains = {gid: {ruler.introns(ruler._merge(ex)) for ex in g["tx"].values() if len(ex) > 1} for gid, g in lnc_t.items()}
    t_left = {gid: {a for a, _ in s} for gid, s in t_introns.items()}
    t_right = {gid: {b for _, b in s} for gid, s in t_introns.items()}
    hit, any_exact, chain, kinds = set(), set(), set(), Counter()
    for p in preds:
        h = idx.hits(p["strand"], p["exons"])
        hit |= h
        pi = ruler.introns([tuple(x) for x in p["exons"]])
        for gid in h:
            if set(pi) & t_introns[gid]:
                any_exact.add(gid)
            if pi and pi in t_chains[gid]:
                chain.add(gid)
        if h:
            for a, b in pi:
                L = any(a in t_left[g] for g in h)
                R = any(b in t_right[g] for g in h)
                kinds["exact" if any((a, b) in t_introns[g] for g in h) else "both_ends_elsewhere" if L and R else
                      "left_only" if L else "right_only" if R else "neither"] += 1
    n_ex = Counter(min(len(p["exons"]), 5) for p in preds)
    t_ex = Counter(min(max(len(ruler._merge(ex)) for ex in lnc_t[g]["tx"].values()), 5) for g in hit)
    multi_hit = [g for g in hit if any(len(ex) > 1 for ex in lnc_t[g]["tx"].values())]
    return {"n_pred_lnc": len(preds), "pred_exon_count": {str(k): v for k, v in sorted(n_ex.items())},
            "pred_single_exon_frac": n_ex[1] / max(len(preds), 1),
            "truth_detected": len(hit), "truth_detected_multi_exon": len(multi_hit),
            "truth_detected_exon_count_longest_tx": {str(k): v for k, v in sorted(t_ex.items())},
            "detected_multi_with_any_exact_intron": len(any_exact & set(multi_hit)) / max(len(multi_hit), 1),
            "detected_multi_with_full_chain": len(chain & set(multi_hit)) / max(len(multi_hit), 1),
            "pred_introns_on_hit_genes": dict(kinds)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="outputs/PHASE1/SPLICE_DIAG")
    ap.add_argument("--config", default="configs/p15_v19b.json")
    ap.add_argument("--run", default="g1b-42-ke")
    a = ap.parse_args()
    cfg = json.load(open(f"{ROOT}/{a.config}"))
    tg, tau = cfg["targets"]["chr18"], cfg["decode"]["tau"]
    truth = ruler.read_gtf(f"{ROOT}/{tg['truth']}", tg["truth_chrom"])
    out = {"w=0": analyse(P.load_genes(f"{ROOT}/{cfg['out_dir']}/chr18/{a.run}/genes.json"), truth, tau)}
    for v in ("w0.5", "w1", "oracle"):
        f = f"{ROOT}/{a.dir}/{a.run}_{v}/genes.json"
        if os.path.exists(f):
            out[v] = analyse(P.load_genes(f), truth, tau)
    json.dump(out, open(f"{ROOT}/{a.dir}/WHY_{a.run}.json", "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
