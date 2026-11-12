#!/usr/bin/env python3
"""B1 (docs/phase1_results/PLAN_B1.md): strict one-to-one locus evaluation of existing predictions (CPU, no model).

A predicted gene (span = leftmost to rightmost exon boundary) and a truth gene (span over all transcripts) of the same
class and strand are a candidate pair when |p & g| / |g| >= theta and |p & g| / |p| >= theta; a maximum-cardinality
bipartite matching makes the pairs one-to-one.  S50: theta 0.5 (primary), S25: theta 0.25.  Sn = matched / truth
genes, Pr = matched / predictions.  The ruler's L1 (same-strand exon overlap >= 1 bp, many-to-many) is given next
to it, with descriptive numbers on merging, span ratio and end offsets.

  b1_strict.py [--out outputs/PHASE1/B1]   -> <out>/B1.json
"""
import argparse
import json
import os
import statistics
import sys
import time
from collections import defaultdict

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1/v20"))
import p15_score as P  # noqa: E402
import tiberius_score as T  # noqa: E402
from ncrna_annot.eval import ruler  # noqa: E402

V50 = f"{ROOT}/data/human/gencode.v50.annotation.gtf.gz"
THETAS = (0.5, 0.25)


def truth_spans(truth, cls):
    out = []
    for gid, g in truth.items():
        if g["cls"] == cls:
            s = min(e[0] for t in g["tx"].values() for e in t)
            e = max(e[1] for t in g["tx"].values() for e in t)
            out.append((gid, g["strand"], s, e))
    return out


def collapse(preds):
    """Transcripts of one predicted gene -> one prediction (union of exons); key = 'gene' field when present."""
    by = defaultdict(list)
    for i, p in enumerate(preds):
        by[(p.get("gene", f"_{i}"), p["cls"], p["strand"])].append(p)
    return [{"id": k[0], "cls": k[1], "strand": k[2],
             "exons": ruler._merge(sorted(x for p in ps for x in p["exons"]))} for k, ps in by.items()]


def matching(adj, n_right):
    """Maximum bipartite matching (Kuhn, iterative).  adj[i] = right nodes of left node i."""
    mr = [-1] * n_right
    ml = [-1] * len(adj)
    for i in sorted(range(len(adj)), key=lambda x: len(adj[x])):
        if not adj[i]:
            continue
        seen, stack, parent = set(), [i], {}
        found = None
        while stack and found is None:
            u = stack.pop()
            for v in adj[u]:
                if v in seen:
                    continue
                seen.add(v)
                parent[v] = u
                if mr[v] == -1:
                    found = v
                    break
                stack.append(mr[v])
        while found is not None:
            u = parent[found]
            prev = ml[u]
            mr[found], ml[u] = u, found
            found = prev if prev != -1 else None
    return ml


def strict(truth, preds, cls):
    ts = truth_spans(truth, cls)
    ps = [(p["strand"], p["exons"][0][0], p["exons"][-1][1]) for p in preds if p["cls"] == cls]
    order = sorted(range(len(ts)), key=lambda i: ts[i][2])
    starts = [ts[i][2] for i in order]
    maxlen = max((t[3] - t[2] for t in ts), default=0)
    res = {"n_truth": len(ts), "n_pred": len(ps)}
    cand = []  # (pred, truth, ov, Rg, Rp)
    for j, (st, s, e) in enumerate(ps):
        lo = np.searchsorted(starts, s - maxlen - 1)
        hi = np.searchsorted(starts, e, side="right")
        for k in range(lo, hi):
            i = order[k]
            _, gst, gs, ge = ts[i]
            if st != "." and gst != "." and st != gst:
                continue
            ov = min(e, ge) - max(s, gs) + 1
            if ov > 0:
                cand.append((j, i, ov, ov / (ge - gs + 1), ov / (e - s + 1)))
    for th in THETAS:
        adj = [[] for _ in ps]
        for j, i, ov, rg, rp in cand:
            if rg >= th and rp >= th:
                adj[j].append(i)
        ml = matching(adj, len(ts))
        m = [(j, i) for j, i in enumerate(ml) if i != -1]
        iou = []
        for j, i in m:
            s, e, gs, ge = ps[j][1], ps[j][2], ts[i][2], ts[i][3]
            ov = min(e, ge) - max(s, gs) + 1
            iou.append(ov / ((e - s + 1) + (ge - gs + 1) - ov))
        res[f"S{int(th * 100)}"] = {"matched": len(m), "sn": len(m) / len(ts) if ts else None,
                                    "pr": len(m) / len(ps) if ps else None,
                                    "median_iou": statistics.median(iou) if iou else None}
    # descriptive: span ratio / end offsets against the same-class truth gene with the largest overlap
    best = {}
    for j, i, ov, rg, rp in cand:
        if j not in best or ov > best[j][1]:
            best[j] = (i, ov)
    ratio, tss, tes = [], [], []
    for j, (i, _) in best.items():
        st, s, e = ps[j]
        _, gst, gs, ge = ts[i]
        ratio.append((e - s + 1) / (ge - gs + 1))
        minus = (gst if st == "." else st) == "-"
        tss.append(abs((e - ge) if minus else (s - gs)))
        tes.append(abs((s - gs) if minus else (e - ge)))
    res["span"] = {"n_pred_overlapping_truth_span": len(best),
                   "median_span_ratio": statistics.median(ratio) if ratio else None,
                   "median_tss_offset": statistics.median(tss) if tss else None,
                   "median_tes_offset": statistics.median(tes) if tes else None}
    return res


def l1(truth, preds, cls):
    idx = ruler.ExonIndex(truth)
    hit, nhit, multi, k = set(), 0, 0, []
    cp = [p for p in preds if p["cls"] == cls]
    for p in cp:
        strands = ("+", "-") if p["strand"] == "." else (p["strand"],)
        same = {x for s in strands for x in idx.hits(s, p["exons"]) if truth[x]["cls"] == cls}
        hit |= same
        if same:
            nhit += 1
            k.append(len(same))
            multi += len(same) >= 2
    nt = sum(g["cls"] == cls for g in truth.values())
    return {"sn": len(hit) / nt if nt else None, "pr_strict": nhit / len(cp) if cp else None,
            "frac_hit_preds_with_2plus_truth": multi / nhit if nhit else None,
            "mean_truth_per_hit_pred": sum(k) / len(k) if k else None, "max_truth_per_hit_pred": max(k, default=None)}


def score(truth, preds, classes=("lnc", "mrna")):
    return {c: dict(strict(truth, preds, c), L1=l1(truth, preds, c)) for c in classes}


def gff_preds(path, chrom):
    ps = ruler.read_pred_gff3(path, chrom)
    for p in ps:
        p["gene"] = p["attrs"].get("Parent", p["id"])
    return collapse(ps)


def mean_runs(runs, cls):
    out = {}
    for key in ("S50", "S25", "L1"):
        for f in ("sn", "pr", "pr_strict"):
            v = [r[cls][key].get(f) for r in runs if r[cls][key].get(f) is not None]
            if v:
                out[f"{key}_{f}"] = sum(v) / len(v)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="outputs/PHASE1/B1")
    a = ap.parse_args()
    t0 = time.time()
    os.makedirs(f"{ROOT}/{a.out}", exist_ok=True)
    res = {"plan": "docs/phase1_results/PLAN_B1.md", "thetas": THETAS, "sets": {}}

    c19 = json.load(open(f"{ROOT}/configs/p15_v19.json"))
    tau = c19["decode"]["tau"]
    truths = {}
    for tgt in ("chr18", "chr20"):
        tg = c19["targets"][tgt]
        truths[tgt] = truth = ruler.read_gtf(f"{ROOT}/{tg['truth']}", tg["truth_chrom"])
        runs = {}
        for run in c19["main"]["runs"].values():
            genes = P.load_genes(f"{ROOT}/{c19['out_dir']}/{tgt}/{run}/genes.json")
            runs[run] = score(truth, P.final_preds(genes, tau))
        res["sets"][f"v19_{tgt}"] = {"tau": tau, "runs": runs,
                                     "mean": {c: mean_runs(list(runs.values()), c) for c in ("lnc", "mrna")}}
        print(tgt, json.dumps(res["sets"][f"v19_{tgt}"]["mean"]), flush=True)

    # v2.0 chr9: baseline tau and calibrated tau (as v20_score.py)
    import v20_calib as C
    c20 = json.load(open(f"{ROOT}/configs/v20.json"))
    cal, d = c20["calibration"], c20["decode"]
    tab = json.load(open(f"{ROOT}/{cal['hexamer_table']}"))
    tg = c20["targets"]["chr9"]
    truth = ruler.read_gtf(f"{ROOT}/{tg['truth']}", tg["truth_chrom"])
    codes = np.load(f"{ROOT}/{tg['codes']}", mmap_mode="r")
    rb, rv = {}, {}
    for run in tg["runs"]:
        rd = f"{ROOT}/{c20['out_dir']}/chr9/{run}"
        genes = P.load_genes(f"{rd}/genes.json")
        post = np.load(f"{rd}/post.npy", mmap_mode="r")
        rows = C.gene_rows(genes, post, codes, tab)
        X = np.array([[r[f] for f in cal["model"]["features"]] for r in rows], float)
        p = C.apply(cal["model"], X) if cal["keep"] else X[:, cal["model"]["features"].index("conf")]
        sc = {r["k"]: float(s) for r, s in zip(rows, p)}
        gv = [dict(g, conf=sc.get(k, g["conf"])) for k, g in enumerate(genes)]
        rb[run] = score(truth, P.final_preds(genes, d["tau"]))
        rv[run] = score(truth, P.final_preds(gv, cal["tau_v20"]))
    for name, runs, t in (("v20base_chr9", rb, d["tau"]), ("v20_chr9", rv, cal["tau_v20"])):
        res["sets"][name] = {"tau": t, "runs": runs, "mean": {c: mean_runs(list(runs.values()), c) for c in ("lnc", "mrna")}}
        print(name, json.dumps(res["sets"][name]["mean"]), flush=True)

    # chr20 baselines and the v41 -> v50 yardstick
    t20 = truths["chr20"]
    base = {
        "genatator_chr20": gff_preds(P.G_CHR20, "chr20"),
        "segnt_t0.5_both_chr20": gff_preds(f"{ROOT}/outputs/PHASE1/BASELINES/segnt/pred_t0.5_both.gff3", "chr20"),
        "tiberius_chr20": collapse([dict(p, gene=p["id"].rsplit("_t", 1)[0])
                                    for p in T.read(f"{ROOT}/outputs/PHASE1/BASELINES/tiberius/chr20.gtf")]),
    }
    v41 = ruler.read_gtf(P.V41, "chr20")
    base["gencode_v41_as_pred_chr20"] = collapse([
        {"gene": gid, "cls": g["cls"], "strand": g["strand"], "exons": ruler._merge(sorted(t))}
        for gid, g in v41.items() if g["cls"] in ("lnc", "mrna") for t in g["tx"].values()])
    for name, preds in base.items():
        res["sets"][name] = score(t20, preds)
        print(name, json.dumps({c: {k: (v if not isinstance(v, dict) else {kk: vv for kk, vv in v.items() if kk in ("sn", "pr", "pr_strict")})
                                    for k, v in res["sets"][name][c].items() if k in ("n_pred", "S50", "S25", "L1")}
                                for c in ("lnc", "mrna")}), flush=True)
    res["seconds"] = round(time.time() - t0, 1)
    json.dump(P.norm(res), open(f"{ROOT}/{a.out}/B1.json", "w"), indent=1)
    print("written", f"{a.out}/B1.json", res["seconds"], "s")


if __name__ == "__main__":
    main()
