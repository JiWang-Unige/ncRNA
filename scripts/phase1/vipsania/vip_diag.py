#!/usr/bin/env python3
"""chr18 diagnostic: does coding evidence that does not depend on our exon structure (Vipsania or Tiberius) help the
v2.0 calibration tell v1.9's lnc true positives from false positives?  Card vipsania-ref option 1 (2026-11-03);
design and decision rule in docs/phase1_results/REVIEW_VIPSANIA.md section 5 (GitHub 2edca99, written before any run).
chr18 only (development chromosome).  chr9 is not touched; chr20 is only scored as a descriptive baseline (part 4).

Implementation choices fixed here before the Vipsania / Tiberius outputs were read:

Rows      every lnc prediction (type "l") in outputs/PHASE1/P15/chr18/<run>/genes.json (MEAX gamma 3, before tau*),
          runs g1b-42-ke / g1b-123-ke / g1b-456-ke (configs/p15_v19.json).
Target    y = 1 when the prediction has a same-strand exon overlap with a v50 lnc gene (the "strict" L1 hit of
          p15_score.ece; truth outputs/PHASE1/DIAG_L2/chr18_v50.gtf as in the P1.5 config).
B         conf (v1.9 gene confidence), mrna_mass (mean over the predicted exon bases of the summed posterior of the
          four same-strand mRNA states mE/mL/mI/mR in post.npy), log mature length, exon count, log(1 + longest ORF)
          and longest ORF / transcript length, Fickett TESTCODE (CPAT tables, cpmodule/fickett.py) of the whole
          spliced transcript.  Longest ORF: ATG to the first in-frame TAA/TAG/TGA on the transcript strand, stop
          included, 0 when none; codons with N never match.  The hexamer score is left out (it needs a table from the
          training chromosomes, i.e. v2.0 development).
V         Vipsania (outputs/PHASE1/VIPSANIA/chr18/vipsania.gff3, CDS rows; post6.npy):
          v_cov_span  fraction of the prediction's span (first exon start .. last exon end) covered by same-strand CDS
          v_cov_exon  fraction of the predicted exon bases covered by same-strand CDS
          v_any_span  1 when any same-strand CDS overlaps the span
          v_post      mean same-strand coding posterior (post6 column cod+ / cod-) over the predicted exon bases
T         Tiberius (outputs/PHASE1/VIPSANIA/chr18/tiberius.gtf, CDS rows): t_cov_span, t_cov_exon, t_any_span.
Model     logistic regression, L2 with C = 1 (penalty 0.5 * |w|^2, intercept not penalised), features standardised with
          the training half's mean / sd; scipy L-BFGS.
Folds     chr18 halves by the prediction's first exon start (< 40,000,000 / >= 40,000,000, the v1.9 dev split); fit on
          one half, average precision (own_splice.average_precision) of y on the other, both directions, per run.
Rule      source X effective when, over the 3 runs, mean(delta AP of B+X vs B, averaged over the two halves) >= 0.02
          AND at least 5 of the 6 (run x half) deltas are > 0.  Both effective: take the larger mean delta; if they
          differ by < 0.01 take T (no new dependency).  B+V+T is used only if it beats the chosen single source by
          >= 0.01 mean delta.  Effective -> revise DRAFT_V20 section 3.3; none -> the draft stays.
Sanity    first: Vipsania chr18 mRNA L1 Sn >= 0.5 with the P1.5 ruler (every Vipsania transcript one mRNA prediction,
          exons = exon rows, else CDS rows).  If it fails, the verdict is "not readable" and usage is checked first.
Report only (no decision): lnc / mRNA site AP of the Vipsania acceptor / donor posteriors with own_splice.site_ap
          (SpliceAI 0.140, own S10k 0.161); the base offset of the Vipsania columns (-3..+3) is chosen per column by
          the mRNA AP, never by lnc; guard counts of Vipsania and Tiberius on chr18; chr20 Vipsania scored like the
          other chr20 baselines (RESULT_BASELINES.md) -> outputs/PHASE1/VIPSANIA/chr20/SCORE.json.

  vip_diag.py   (generanno env, CPU)  -> outputs/PHASE1/VIPSANIA/{DIAG.json, VERDICT.json, SITE_AP.json, chr18/SANITY.json}
"""
import json
import math
import os
import sys
from collections import defaultdict

import numpy as np
from scipy.optimize import minimize

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1/splice_diag"))
import own_splice  # noqa: E402
import p15_score as P  # noqa: E402
import sai_score  # noqa: E402
from ncrna_annot.eval import ruler  # noqa: E402
from ncrna_annot.grammar import g1  # noqa: E402

D = f"{ROOT}/outputs/PHASE1/VIPSANIA"
RUNS = ("g1b-42-ke", "g1b-123-ke", "g1b-456-ke")
SPLIT = 40_000_000
B_FEATS = ["conf", "mrna_mass", "log_len", "n_exon", "log_orf", "orf_frac", "fickett"]
V_FEATS = ["v_cov_span", "v_cov_exon", "v_any_span", "v_post"]
T_FEATS = ["t_cov_span", "t_cov_exon", "t_any_span"]
SETS = {"B": B_FEATS, "B+V": B_FEATS + V_FEATS, "B+T": B_FEATS + T_FEATS, "B+V+T": B_FEATS + V_FEATS + T_FEATS}

# ---------------------------------------------------------------- sequence features
FPOS = {"A": [.94, .68, .84, .93, .58, .68, .45, .34, .20, .22], "C": [.80, .70, .70, .81, .66, .48, .51, .33, .30, .23],
        "G": [.90, .88, .74, .64, .53, .48, .27, .16, .08, .08], "T": [.97, .97, .91, .68, .69, .44, .54, .20, .09, .09]}
FPOS_W = {"A": .26, "C": .18, "G": .31, "T": .33}
FPOS_P = [1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 0.0]
FCON = {"A": [.28, .49, .44, .55, .62, .49, .67, .65, .81, .21], "C": [.82, .64, .51, .64, .59, .59, .43, .44, .39, .31],
        "G": [.40, .54, .47, .64, .64, .73, .41, .41, .33, .29], "T": [.28, .24, .39, .40, .55, .75, .56, .69, .51, .58]}
FCON_W = {"A": .11, "C": .12, "G": .15, "T": .14}
FCON_P = [.33, .31, .29, .27, .25, .23, .21, .17, 0]


def _look(value, para, prob, w):
    for i, v in enumerate(para):
        if value >= v:
            return prob[i] * w


def fickett(s):
    """CPAT cpmodule/fickett.py fickett_value."""
    if len(s) < 2:
        return 0.0
    n, sc = len(s), 0.0
    ph = [s[i::3] for i in range(3)]
    for b in "ACGT":
        sc += _look(s.count(b) / n, FCON_P, FCON[b], FCON_W[b])
        c = [p.count(b) for p in ph]
        sc += _look(max(c) / (min(c) + 1.0), FPOS_P, FPOS[b], FPOS_W[b])
    return sc


STOPS = {"TAA", "TAG", "TGA"}


def longest_orf(s):
    best = 0
    for f in range(3):
        start = None
        for i in range(f, len(s) - 2, 3):
            c = s[i:i + 3]
            if start is None and c == "ATG":
                start = i
            elif start is not None and c in STOPS:
                best = max(best, i + 3 - start)
                start = None
    return best


COMP = str.maketrans("ACGTN", "TGCAN")


def transcript(codes, g):
    s = "".join("ACGTN"[min(int(c), 4)] for a, b in g["exons"] for c in codes[a - 1:b])
    return s.translate(COMP)[::-1] if g["strand"] == "-" else s


# ---------------------------------------------------------------- external predictions
def read_models(path, chrom):
    """{tid: (strand, exons, cds)} from a GFF3 or GTF (exons / CDS merged, 1-based inclusive)."""
    ex, cds, strand, parent = defaultdict(list), defaultdict(list), {}, {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(f) < 9 or f[0] != chrom or f[2] not in ("exon", "CDS"):
            continue
        if "transcript_id \"" in f[8]:
            tid = f[8].split('transcript_id "', 1)[1].split('"', 1)[0]
        else:
            a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            tid = a["Parent"].split(",")[0]
        (ex if f[2] == "exon" else cds)[tid].append((int(f[3]), int(f[4])))
        strand[tid] = f[6]
    return {t: (strand[t], ruler._merge(sorted(ex.get(t) or cds[t])), ruler._merge(sorted(cds.get(t, []))))
            for t in sorted(strand)}


def cds_mask(models, L):
    m = {"+": np.zeros(L + 2, bool), "-": np.zeros(L + 2, bool)}
    for s, _, cds in models.values():
        for a, b in cds:
            m[s][a:b + 1] = True  # index = 1-based position
    return m


def cov_feats(g, mask, pfx):
    s, ex = g["strand"], g["exons"]
    a, b = ex[0][0], ex[-1][1]
    span = mask[s][a:b + 1]
    exb = np.concatenate([mask[s][x:y + 1] for x, y in ex])
    return {f"{pfx}_cov_span": float(span.mean()), f"{pfx}_cov_exon": float(exb.mean()),
            f"{pfx}_any_span": float(span.any())}


# ---------------------------------------------------------------- logistic
def fit_predict(Xtr, ytr, Xte, C=1.0):
    mu, sd = Xtr.mean(0), Xtr.std(0)
    sd[sd == 0] = 1.0
    A, Z = (Xtr - mu) / sd, (Xte - mu) / sd
    A1 = np.hstack([A, np.ones((len(A), 1))])

    def f(w):
        z = A1 @ w
        loss = np.sum(np.logaddexp(0, z) - ytr * z) + 0.5 / C * np.dot(w[:-1], w[:-1])
        p = 1 / (1 + np.exp(-z))
        g = A1.T @ (p - ytr)
        g[:-1] += w[:-1] / C
        return loss, g

    w = minimize(f, np.zeros(A1.shape[1]), jac=True, method="L-BFGS-B", options={"maxiter": 2000}).x
    return 1 / (1 + np.exp(-(np.hstack([Z, np.ones((len(Z), 1))]) @ w)))


def auroc(s, y):
    from scipy.stats import rankdata
    r = rankdata(s)
    n1 = y.sum()
    n0 = len(y) - n1
    return float((r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)) if n1 and n0 else None


# ---------------------------------------------------------------- parts
def rows_for_run(run, truth, codes, vmask, tmask, vpost):
    rd = f"{ROOT}/outputs/PHASE1/P15/chr18/{run}"
    genes = P.load_genes(f"{rd}/genes.json")
    post = np.load(f"{rd}/post.npy", mmap_mode="r")
    lidx = ruler.ExonIndex({k: v for k, v in truth.items() if v["cls"] == "lnc"})
    rows = []
    for g in genes:
        if g["type"] != "l":
            continue
        s = g["strand"]
        ms = [g1.IDX[f"m{p}{s}"] for p in g1.PARTS]
        exb = np.concatenate([np.arange(a - 1, b) for a, b in g["exons"]])
        tx = transcript(codes, g)
        orf = longest_orf(tx)
        r = {"start": g["exons"][0][0], "y": 1.0 if lidx.hits(s, g["exons"]) else 0.0, "conf": g["conf"],
             "mrna_mass": float(np.asarray(post[exb][:, ms], np.float64).sum(1).mean()),
             "log_len": math.log(len(tx)), "n_exon": float(len(g["exons"])), "log_orf": math.log1p(orf),
             "orf_frac": orf / len(tx), "fickett": fickett(tx),
             "v_post": float(np.asarray(vpost[exb, 0 if s == "+" else 1], np.float64).mean())}
        r.update(cov_feats(g, vmask, "v"))
        r.update(cov_feats(g, tmask, "t"))
        rows.append(r)
    return rows


def evaluate(rows):
    y = np.array([r["y"] for r in rows])
    half = np.array([r["start"] >= SPLIT for r in rows])
    out = {"n": len(rows), "n_pos": int(y.sum()), "halves": {}}
    for h, name in ((False, "first"), (True, "second")):
        tr, te = half != h, half == h  # evaluate on half h, fit on the other one
        res = {"n_test": int(te.sum()), "pos_test": int(y[te].sum())}
        for sname, feats in SETS.items():
            X = np.array([[r[f] for f in feats] for r in rows], np.float64)
            p = fit_predict(X[tr], y[tr], X[te])
            res[sname] = {"ap": own_splice.average_precision(p, y[te] > 0), "auroc": auroc(p, y[te])}
        out["halves"][name] = res
    return out


def verdict(per_run):
    dl = {}
    for x in ("B+V", "B+T", "B+V+T"):
        d = [h[x]["ap"] - h["B"]["ap"] for r in per_run.values() for h in r["halves"].values()]
        per = [np.mean([h[x]["ap"] - h["B"]["ap"] for h in r["halves"].values()]) for r in per_run.values()]
        dl[x] = {"mean_delta": float(np.mean(per)), "per_run": [float(v) for v in per], "deltas": [float(v) for v in d],
                 "n_pos": int(sum(v > 0 for v in d)), "effective": bool(np.mean(per) >= 0.02 and sum(v > 0 for v in d) >= 5)}
    ev, et = dl["B+V"]["effective"], dl["B+T"]["effective"]
    if ev and et:
        mv, mt = dl["B+V"]["mean_delta"], dl["B+T"]["mean_delta"]
        choice = "V" if mv - mt >= 0.01 else "T"
    else:
        choice = "V" if ev else "T" if et else None
    if choice and dl["B+V+T"]["mean_delta"] - dl[f"B+{choice}"]["mean_delta"] >= 0.01:
        choice = "V+T"
    return {"sources": dl, "choice": choice, "effective": choice is not None}


def site_ap(truth, codes, vpost):
    """Report only. Vipsania columns -> own_splice.COLS order (acc+, don+, acc-, don-), offset by the mRNA AP."""
    import splice_diag
    ts = {c: sai_score.truth_sites(truth, c) for c in ("mrna", "lnc")}
    base = np.asarray(vpost[:, 2:6], np.float32)  # acc+, don+, acc-, don-
    ok = splice_diag.site_allowed(codes)
    best = {}
    for col in range(4):  # same candidate set as own_splice.site_ap for the mRNA class
        cand = np.flatnonzero(ok[:, col])
        pos = ts["mrna"][col]
        keep = cand[~np.isin(cand, np.setdiff1d(ts["lnc"][col], pos))]
        lab = np.isin(keep, pos)
        scores = {off: own_splice.average_precision(np.roll(base[:, col], off)[keep], lab) for off in range(-3, 4)}
        best[col] = max(scores, key=lambda k: scores[k] if scores[k] is not None else -1)
    sc = np.stack([np.roll(base[:, c], best[c]) for c in range(4)], 1)
    res = own_splice.site_ap(sc, codes, lambda cls, c: ts[cls][c])
    res["offsets"] = {own_splice.COLS[c]: best[c] for c in range(4)}
    res["compare"] = {"spliceai_lnc_ap": 0.140, "own_s10k_lnc_ap": 0.161}
    return res


def score_models(models, truth):
    preds = [{"id": t, "cls": "mrna", "strand": s, "exons": ex} for t, (s, ex, _) in models.items()]
    fin, _ = P.score_preds(truth, preds)
    sm = fin["mrna"]
    return fin, {"n_tx": len(preds), "mrna_L1": [sm["L1_locus"]["sn"], sm["L1_locus"]["pr_strict"]],
                 "mrna_L2_sn": sm["L2_intron_chain"]["sn"], "guards": fin["guards"]}


def main():
    truth = ruler.read_gtf(f"{ROOT}/outputs/PHASE1/DIAG_L2/chr18_v50.gtf", "chr18")
    codes = np.asarray(np.load(f"{ROOT}/data/p1/chr18_codes.npy", mmap_mode="r"))
    L = len(codes)
    vmod = read_models(f"{D}/chr18/vipsania.gff3", "chr18")
    tmod = read_models(f"{D}/chr18/tiberius.gtf", "chr18")
    # sanity first
    fin, short = score_models(vmod, truth)
    json.dump(P.norm({"vipsania": fin, "tiberius": score_models(tmod, truth)[0]}),
              open(f"{D}/chr18/SANITY.json", "w"), indent=1)
    t_short = score_models(tmod, truth)[1]
    print(json.dumps({"sanity_vipsania": short, "tiberius": t_short}), flush=True)
    readable = short["mrna_L1"][0] >= 0.5
    vpost = np.load(f"{D}/chr18/post6.npy", mmap_mode="r")
    assert len(vpost) == L, (len(vpost), L)
    vmask, tmask = cds_mask(vmod, L), cds_mask(tmod, L)
    per_run, desc = {}, {}
    for run in RUNS:
        rows = rows_for_run(run, truth, codes, vmask, tmask, vpost)
        per_run[run] = evaluate(rows)
        y = np.array([r["y"] for r in rows])
        desc[run] = {f: {"pos_mean": float(np.mean([r[f] for r in rows if r["y"]])),
                         "neg_mean": float(np.mean([r[f] for r in rows if not r["y"]]))}
                     for f in B_FEATS + V_FEATS + T_FEATS}
        desc[run]["n"], desc[run]["pos_rate"] = len(rows), float(y.mean())
        print(json.dumps({"run": run, **{h: {k: v["ap"] for k, v in r.items() if isinstance(v, dict)}
                                         for h, r in per_run[run]["halves"].items()}}), flush=True)
    vd = verdict(per_run)
    vd["readable"] = readable
    if not readable:
        vd["effective"], vd["choice"] = None, None
    json.dump({"per_run": per_run, "feature_means": desc, "sanity": {"vipsania": short, "tiberius": t_short}},
              open(f"{D}/DIAG.json", "w"), indent=1)
    json.dump(vd, open(f"{D}/VERDICT.json", "w"), indent=1)
    print(json.dumps({k: vd[k] for k in ("readable", "effective", "choice")}
                     | {x: round(v["mean_delta"], 4) for x, v in vd["sources"].items()}), flush=True)
    # report only
    ap = site_ap(truth, codes, vpost)
    json.dump(ap, open(f"{D}/SITE_AP.json", "w"), indent=1)
    print(json.dumps({"site_ap_lnc": ap["lnc_ap"], "site_ap_mrna": ap["mrna_ap"], "offsets": ap["offsets"]}), flush=True)
    t20 = ruler.read_gtf(f"{ROOT}/data/human/gencode.v50.annotation.gtf.gz", "chr20")
    v20 = read_models(f"{D}/chr20/vipsania.gff3", "chr20")
    fin20, s20 = score_models(v20, t20)
    json.dump(P.norm({"n_tx": len(v20), "gff": "outputs/PHASE1/VIPSANIA/chr20/vipsania.gff3", "final": fin20}),
              open(f"{D}/chr20/SCORE.json", "w"), indent=1)
    print(json.dumps({"chr20": s20}), flush=True)


if __name__ == "__main__":
    main()
