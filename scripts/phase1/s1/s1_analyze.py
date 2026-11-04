#!/usr/bin/env python3
"""S1 (PLAN_S1.md, card subjournal-portfolio 2026-11-04): splice-site predictability by usage, class and annotation
era on hg38 chr20, across splice models.  Descriptive; nothing here feeds a model choice.

Sites (independent of GENCODE): every Snaptron gtexv2 junction gives a left and a right site in the sai_score.py
convention (intron (s, e) 1-based -> 0-based s - 2 and e; columns acc+, don+, acc-, don-).
  positives  = sites of canonical junctions (+: GT..AG, -: CT..AC), usage = sum of coverage_sum over the site's
               junctions, samples = max samples_count;
  negatives  = positions with the same canonical dinucleotide that are in no gtexv2 junction (any motif) and are no
               v41 / v50 annotated site;
  class      = same-strand v50 gene body (first to last exon base): lnc only / mRNA only; both / neither are counted
               but left out of R1-R2;
  era        = positive site is a v41 site / only a v50 site / unannotated.
Readouts (PLAN_S1.md section 4): R1 AUROC (bootstrap by gene, 2,000, seed 0) and AP at a fixed positive:negative
ratio per class, gap mRNA - lnc unmatched and usage-matched (direct standardisation over usage layer x samples
quintile); R2 = R1 by era; R3 SpliceAI-paper style top-k accuracy / auPRC over all positions of v50 transcripts
(mRNA, lnc transcripts present in v41, lnc transcripts new since v41); R4 usage of v50 introns by class and era.

  s1_analyze.py --scores spliceai=outputs/PHASE1/S1/sai_chr20.npy own=outputs/PHASE1/S1/own_chr20.npy \
      vipsania=outputs/PHASE1/VIPSANIA/chr20/post6.npy:2 --out outputs/PHASE1/S1/S1.json
(name=path[:first column]; the four columns from there on are acc+, don+, acc-, don-.)
"""
import argparse
import glob
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
from ncrna_annot.eval import ruler  # noqa: E402

COLS = ("acc+", "don+", "acc-", "don-")
COL = {("+", "L"): 1, ("+", "R"): 0, ("-", "L"): 2, ("-", "R"): 3}
STRAND_OF_COL = ("+", "+", "-", "-")
A, C, G, T = 0, 1, 2, 3
USAGE_EDGES = (1, 10, 100, 1000)          # layers 1-9, 10-99, 100-999, >= 1000 reads
R4_EDGES = (0, 1, 10, 100, 1000)          # 0, 1-9, 10-99, 100-999, >= 1000
N_BOOT, NEG_SUB, AP_RATIO, AP_REPEATS = 2000, 300_000, 100, 20


def canonical_mask(codes):
    """[L, 4] bool: the canonical dinucleotide of each column's site is present (see module docstring)."""
    c = np.asarray(codes)
    L = len(c)
    m = np.zeros((L, 4), bool)
    nxt = lambda a, b: np.concatenate([(c[1:-1] == a) & (c[2:] == b), [False, False]])     # bases p+1, p+2
    prv = lambda a, b: np.concatenate([[False, False], (c[:-2] == a) & (c[1:-1] == b)])    # bases p-2, p-1
    m[:, 0], m[:, 1], m[:, 2], m[:, 3] = prv(A, G), nxt(G, T), nxt(C, T), prv(A, C)
    return m


def junction_sites(s, e, strand):
    """0-based (col, pos) of the left and right site of introns (s, e) on `strand` arrays."""
    plus = strand == "+"
    lc = np.where(plus, COL[("+", "L")], COL[("-", "L")])
    rc = np.where(plus, COL[("+", "R")], COL[("-", "R")])
    return (lc, s - 2), (rc, e)


def read_snaptron(paths, chrom):
    rows = {}
    for p in paths:
        with open(p) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            ix = {k: i for i, k in enumerate(head)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[ix["chromosome"]] != chrom:
                    continue
                rows[f[ix["snaptron_id"]]] = (int(f[ix["start"]]), int(f[ix["end"]]), f[ix["strand"]],
                                              f[ix["left_motif"]].upper(), f[ix["right_motif"]].upper(),
                                              int(f[ix["samples_count"]]), int(f[ix["coverage_sum"]]))
    j = list(rows.values())
    return {"s": np.array([r[0] for r in j], np.int64), "e": np.array([r[1] for r in j], np.int64),
            "strand": np.array([r[2] for r in j]), "lm": np.array([r[3] for r in j]), "rm": np.array([r[4] for r in j]),
            "samples": np.array([r[5] for r in j], np.int64), "cov": np.array([r[6] for r in j], np.int64)}


def annotated_keys(truth, L):
    keys = []
    for g in truth.values():
        for ex in g["tx"].values():
            for s, e in ruler.introns(ruler._merge(sorted(ex))):
                keys += [COL[(g["strand"], "L")] * L + s - 2, COL[(g["strand"], "R")] * L + e]
    return np.unique(np.array(keys, np.int64))


def body_index(truth, L, cls):
    """{strand: int32 [L] index of a same-strand `cls` gene whose body covers the base, -1 if none}, gene ids."""
    ids = sorted(g for g, v in truth.items() if v["cls"] == cls)
    out = {s: np.full(L, -1, np.int32) for s in ("+", "-")}
    for i, gid in enumerate(ids):
        g = truth[gid]
        s = min(a for ex in g["tx"].values() for a, _ in ex) - 1
        e = max(b for ex in g["tx"].values() for _, b in ex)
        seg = out[g["strand"]][max(s, 0):min(e, L)]
        seg[seg < 0] = i
    return out, ids


def wauroc_parts(pos_s, neg_sorted):
    """indices for weighted AUROC: for each positive the count of negatives strictly below / tied (searchsorted)."""
    return np.searchsorted(neg_sorted, pos_s, "left"), np.searchsorted(neg_sorted, pos_s, "right")


def wauroc(lo, hi, wpos, wneg_sorted, groups=None, n_groups=1):
    """weighted AUROC (ties 1/2); with `groups`, one AUROC per positive group against the same negatives."""
    cw = np.concatenate([[0.0], np.cumsum(wneg_sorted)])
    contrib = wpos * (cw[lo] + 0.5 * (cw[hi] - cw[lo]))
    W = cw[-1]
    if groups is None:
        return float(contrib.sum() / (wpos.sum() * W))
    num = np.bincount(groups, contrib, n_groups)
    den = np.bincount(groups, wpos, n_groups) * W
    with np.errstate(invalid="ignore", divide="ignore"):
        return num / den


def average_precision(scores, is_pos):
    order = np.argsort(-scores, kind="stable")
    hit = is_pos[order]
    n = hit.sum()
    if n == 0:
        return None
    prec = np.cumsum(hit) / np.arange(1, len(hit) + 1)
    return float(prec[hit].sum() / n)


def load_scores(spec, L):
    name, rest = spec.split("=", 1)
    path, _, c0 = rest.partition(":")
    a = np.load(f"{ROOT}/{path}", mmap_mode="r")
    c0 = int(c0 or 0)
    assert a.shape[0] >= L and a.shape[1] >= c0 + 4, (spec, a.shape)
    return name, a, c0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", default="chr20")
    ap.add_argument("--codes", default="data/p1/chr20_codes.npy")
    ap.add_argument("--snaptron", default="data/bioval/raw/snaptron/chr20_*.tsv")
    ap.add_argument("--v50", default="data/human/gencode.v50.annotation.gtf.gz")
    ap.add_argument("--v41", default="data/human/gencode.v41.annotation.gtf.gz")
    ap.add_argument("--scores", nargs="+", required=True)
    ap.add_argument("--judge", default="own", help="model used by the PLAN_S1 section 5 rule")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    t0 = time.time()
    rng = np.random.default_rng(0)
    codes = np.load(f"{ROOT}/{a.codes}", mmap_mode="r")
    L = len(codes)
    v50 = ruler.read_gtf(f"{ROOT}/{a.v50}", a.chrom)
    v41 = ruler.read_gtf(f"{ROOT}/{a.v41}", a.chrom)
    jn = read_snaptron(sorted(glob.glob(f"{ROOT}/{a.snaptron}")), a.chrom)
    canon = canonical_mask(codes)
    meta = {"n_junctions": int(len(jn["s"])), "L": L}

    # ---- every gtexv2 junction site (any motif; unknown strand -> both) is excluded from the negatives
    used = []
    for strand in ("+", "-"):
        sel = (jn["strand"] == strand) | ~np.isin(jn["strand"], ("+", "-"))
        (lc, lp), (rc, rp) = junction_sites(jn["s"][sel], jn["e"][sel], np.full(sel.sum(), strand))
        used += [lc * L + lp, rc * L + rp]
    used = np.unique(np.concatenate(used))
    ann50, ann41 = annotated_keys(v50, L), annotated_keys(v41, L)
    # ---- positives: canonical junctions
    plus = (jn["strand"] == "+") & (jn["lm"] == "GT") & (jn["rm"] == "AG")
    minus = (jn["strand"] == "-") & (jn["lm"] == "CT") & (jn["rm"] == "AC")
    ok = plus | minus
    meta["n_canonical_junctions"] = int(ok.sum())
    (lc, lp), (rc, rp) = junction_sites(jn["s"][ok], jn["e"][ok], jn["strand"][ok])
    col_, pos_ = np.concatenate([lc, rc]), np.concatenate([lp, rp])
    cov = np.concatenate([jn["cov"][ok]] * 2)
    smp = np.concatenate([jn["samples"][ok]] * 2)
    inside = (pos_ >= 0) & (pos_ < L)
    keys, cov, smp = col_[inside] * L + pos_[inside], cov[inside], smp[inside]
    pk, inv = np.unique(keys, return_inverse=True)
    pcov = np.bincount(inv, cov).astype(np.int64)
    psmp = np.zeros(len(pk), np.int64)
    np.maximum.at(psmp, inv, smp)
    pcol, ppos = pk // L, pk % L
    meta["convention_check_frac_canonical"] = float(canon[ppos, pcol].mean())
    assert meta["convention_check_frac_canonical"] > 0.99, meta  # coordinates / columns are right
    # ---- negatives
    nk = []
    for c in range(4):
        pos = np.flatnonzero(canon[:, c])
        nk.append(c * L + pos)
    nk = np.concatenate(nk)
    nk = nk[~np.isin(nk, used) & ~np.isin(nk, ann50) & ~np.isin(nk, ann41)]
    ncol, npos = nk // L, nk % L
    # ---- class and gene by body
    lb, lids = body_index(v50, L, "lnc")
    mb, mids = body_index(v50, L, "mrna")

    def classify(col, pos):
        lg = np.where(col < 2, lb["+"][pos], lb["-"][pos])
        mg = np.where(col < 2, mb["+"][pos], mb["-"][pos])
        cls = np.full(len(pos), "none", object)
        cls[(lg >= 0) & (mg < 0)] = "lnc"
        cls[(mg >= 0) & (lg < 0)] = "mrna"
        cls[(mg >= 0) & (lg >= 0)] = "both"
        gene = np.where(cls == "lnc", lg, np.where(cls == "mrna", mg, -1))
        return cls, gene
    pcls, pgene = classify(pcol, ppos)
    ncls, ngene = classify(ncol, npos)
    era = np.where(np.isin(pk, ann41), "v41", np.where(np.isin(pk, ann50), "v42-50", "unannotated"))
    layer = np.searchsorted(USAGE_EDGES, pcov, "right") - 1          # 0..3
    q = np.quantile(psmp[np.isin(pcls, ("lnc", "mrna"))], [0.2, 0.4, 0.6, 0.8])
    squint = np.searchsorted(q, psmp, "right")                        # 0..4
    stratum = layer * 5 + squint
    meta["samples_quintile_edges"] = q.tolist()
    meta["counts"] = {
        "positives": {k: int((pcls == k).sum()) for k in ("lnc", "mrna", "both", "none")},
        "negatives": {k: int((ncls == k).sum()) for k in ("lnc", "mrna", "both", "none")},
        "positives_by_era": {c: {e: int(((pcls == c) & (era == e)).sum()) for e in ("v41", "v42-50", "unannotated")} for c in ("lnc", "mrna")},
        "positives_by_layer": {c: np.bincount(layer[pcls == c], minlength=4).tolist() for c in ("lnc", "mrna")},
        "genes": {"lnc": len(lids), "mrna": len(mids)}}
    # fixed negative subsample per class for the bootstrap (seed 0)
    negsub = {}
    for c in ("lnc", "mrna"):
        idx = np.flatnonzero(ncls == c)
        negsub[c] = np.sort(rng.choice(idx, min(NEG_SUB, len(idx)), replace=False))
    pool_ratio = min(AP_RATIO, *[(ncls == c).sum() / max((pcls == c).sum(), 1) for c in ("lnc", "mrna")])
    meta["ap_ratio"] = float(pool_ratio)
    boot_mult = {c: [rng.multinomial(n, np.full(n, 1 / n)) for _ in range(N_BOOT)]
                 for c, n in (("lnc", len(lids)), ("mrna", len(mids)))}

    res = {"meta": meta, "models": {}}
    for spec in a.scores:
        name, arr, c0 = load_scores(spec, L)
        psc = np.asarray(arr[ppos, c0 + pcol], np.float32)
        nsc = np.asarray(arr[npos, c0 + ncol], np.float32)
        m = {}
        # ---------------- R1 / R2
        def block(sel_extra=None):
            out, parts = {}, {}
            for c in ("lnc", "mrna"):
                psel = pcls == c
                if sel_extra is not None:
                    psel &= sel_extra
                ni = negsub[c]
                order = np.argsort(nsc[ni], kind="stable")
                nsorted, ngs = nsc[ni][order], ngene[ni][order]
                pi = np.flatnonzero(psel)
                lo, hi = wauroc_parts(psc[pi], nsorted)
                parts[c] = (pi, lo, hi, nsorted, ngs)
                if len(pi) == 0:
                    out[c] = None
                    continue
                full = wauroc(lo, hi, np.ones(len(pi)), np.ones(len(ngs)))
                strat = wauroc(lo, hi, np.ones(len(pi)), np.ones(len(ngs)), stratum[pi], 20)
                # AP at the common fixed ratio, 20 negative draws
                allneg = np.flatnonzero(ncls == c)
                k = int(min(len(allneg), round(pool_ratio * len(pi))))
                aps = []
                for _ in range(AP_REPEATS):
                    nn = rng.choice(allneg, k, replace=False)
                    aps.append(average_precision(np.concatenate([psc[pi], nsc[nn]]),
                                                 np.concatenate([np.ones(len(pi), bool), np.zeros(k, bool)])))
                out[c] = {"n_pos": int(len(pi)), "auroc": full, "ap_fixed_ratio": float(np.mean(aps)),
                          "auroc_by_stratum": [None if np.isnan(x) else float(x) for x in strat],
                          "n_pos_by_stratum": np.bincount(stratum[pi], minlength=20).tolist()}
            if out["lnc"] is None or out["mrna"] is None:
                return out
            nl, nm = np.array(out["lnc"]["n_pos_by_stratum"]), np.array(out["mrna"]["n_pos_by_stratum"])
            good = (nl >= 5) & (nm >= 5)
            w = np.where(good, nl + nm, 0) / max((nl + nm)[good].sum(), 1)

            def deltas(mult=None):
                au = {}
                for c in ("lnc", "mrna"):
                    pi, lo, hi, nsorted, ngs = parts[c]
                    wp = np.ones(len(pi)) if mult is None else mult[c][pgene[pi]].astype(float)
                    wn = np.ones(len(ngs)) if mult is None else mult[c][ngs].astype(float)
                    full = wauroc(lo, hi, wp, wn)
                    st = wauroc(lo, hi, wp, wn, stratum[pi], 20)
                    au[c] = (full, float(np.nansum(np.where(good, st, 0) * w)))
                return au["mrna"][0] - au["lnc"][0], au["mrna"][1] - au["lnc"][1], au
            d_un, d_m, au = deltas()
            boots = np.array([deltas({c: boot_mult[c][b] for c in ("lnc", "mrna")})[:2] for b in range(N_BOOT)])
            out["gap"] = {"unmatched": d_un, "unmatched_ci95": np.nanpercentile(boots[:, 0], [2.5, 97.5]).tolist(),
                          "matched": d_m, "matched_ci95": np.nanpercentile(boots[:, 1], [2.5, 97.5]).tolist(),
                          "matched_auroc": {"lnc": au["lnc"][1], "mrna": au["mrna"][1]},
                          "strata_used": int(good.sum()), "pos_frac_in_used_strata": float(((nl + nm) * good).sum() / max((nl + nm).sum(), 1))}
            return out
        m["R1"] = block()
        m["R2"] = {e: block(era == e) for e in ("v41", "v42-50", "unannotated")}
        # ---------------- R3 SpliceAI-paper style
        r3 = {}
        tx41 = {t.split(".")[0] for g in v41.values() for t in g["tx"]}
        groups = {"mrna": [], "lnc_v41": [], "lnc_new": []}
        for g in v50.values():
            if g["cls"] not in ("lnc", "mrna"):
                continue
            for tid, ex in g["tx"].items():
                ex = ruler._merge(sorted(ex))
                if len(ex) < 2:
                    continue
                key = "mrna" if g["cls"] == "mrna" else ("lnc_v41" if tid.split(".")[0] in tx41 else "lnc_new")
                groups[key].append((g["strand"], ex))
        for key, txs in groups.items():
            span_mask, posk = np.zeros((4, L), bool), []
            for strand, ex in txs:
                s0, e0 = max(ex[0][0] - 1, 0), min(ex[-1][1], L)
                for c in ((0, 1) if strand == "+" else (2, 3)):
                    span_mask[c, s0:e0] = True
                for s, e in ruler.introns(ex):
                    posk += [COL[(strand, "L")] * L + s - 2, COL[(strand, "R")] * L + e]
            cand = np.concatenate([c * L + np.flatnonzero(span_mask[c]) for c in range(4)])
            del span_mask
            posk = np.unique(np.array(posk, np.int64))
            other = np.setdiff1d(ann50, posk)
            cand = cand[~np.isin(cand, other)]
            lab = np.isin(cand, posk)
            sc = np.asarray(arr[cand % L, c0 + cand // L], np.float32)
            k = int(lab.sum())
            top = np.argpartition(-sc, k - 1)[:k]
            r3[key] = {"n_tx": len(txs), "n_sites": k, "n_positions": int(len(cand)),
                       "topk_acc": float(lab[top].mean()), "auprc": average_precision(sc, lab)}
        m["R3"] = r3
        res["models"][name] = m
        print(json.dumps({"model": name, "R1_gap": m["R1"].get("gap"), "R3": r3}), flush=True)

    # ---------------- R4 usage of v50 introns (model-free)
    jkey = {(int(s), int(e), st): int(c) for s, e, st, c in zip(jn["s"], jn["e"], jn["strand"], jn["cov"])}
    in41 = {(s, e, g["strand"]) for g in v41.values() for ex in g["tx"].values() for s, e in ruler.introns(ruler._merge(sorted(ex)))}
    r4 = {}
    for c in ("lnc", "mrna"):
        intr = {(s, e, g["strand"]) for g in v50.values() if g["cls"] == c
                for ex in g["tx"].values() for s, e in ruler.introns(ruler._merge(sorted(ex)))}
        for era_name, sel in (("all", intr), ("v41", intr & in41), ("new", intr - in41)):
            u = np.array([jkey.get(x, 0) for x in sel], np.int64)
            lay = np.searchsorted(R4_EDGES, u, "right") - 1
            r4[f"{c}_{era_name}"] = {"n": int(len(u)), "frac_by_layer": (np.bincount(lay, minlength=5) / max(len(u), 1)).tolist(),
                                     "median_reads": float(np.median(u)) if len(u) else None}
    res["R4"] = {"layers": ["0", "1-9", "10-99", "100-999", ">=1000"], **r4}

    # ---------------- PLAN_S1 section 5
    j = res["models"].get(a.judge, {}).get("R1", {}).get("gap")
    verdict = None
    if j:
        lo_m = j["matched_ci95"][0]
        if (j["matched_ci95"][0] <= 0 <= j["matched_ci95"][1]) or abs(j["matched"]) < abs(j["unmatched"]) / 3:
            verdict = "gap mainly explained by usage"
        elif lo_m > 0 and j["matched"] >= j["unmatched"] / 2:
            verdict = "intrinsic lnc gap"
        else:
            verdict = "partly explained by usage"
    sai = res["models"].get("spliceai", {}).get("R3")
    lit = None
    if sai:
        lit = ("literature high score comes from old annotations"
               if sai["lnc_v41"]["topk_acc"] >= 0.7 and sai["lnc_new"]["topk_acc"] < 0.5 else "not explained by annotation era alone")
    res["verdict"] = {"judge": a.judge, "R1": verdict, "R3_spliceai": lit}
    res["seconds"] = round(time.time() - t0)
    os.makedirs(os.path.dirname(f"{ROOT}/{a.out}"), exist_ok=True)
    json.dump(res, open(f"{ROOT}/{a.out}", "w"), indent=1)
    print(json.dumps({"verdict": res["verdict"], "R4": res["R4"], "counts": meta["counts"], "seconds": res["seconds"]}), flush=True)


if __name__ == "__main__":
    main()
