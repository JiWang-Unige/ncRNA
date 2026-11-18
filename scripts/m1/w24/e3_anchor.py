#!/usr/bin/env python3
"""M1 weeks 2-4, E3 (evidence-anchored units, chr18 development): evidence proposes transcription units, the sequence model
selects them.  For each ssCAGE TSS, the K nearest downstream same-strand QuantSeq sites within Dmax give candidate spans;
each span is scored by the v2.0 DNA-only posteriors of the run on that strand:
  lnc = mean over the span of P(lE) + P(lL) + P(lI) + P(lR),   mrna = the same for the m states,
  score = lnc - mrna.
Per TSS the best-scoring candidate is kept; units are then accepted greedily by score (no same-strand overlap) while
score >= theta.  Variant 'tib' also drops units overlapping a Tiberius coding gene (as the pairing baseline).
  e3_anchor.py [--frac 1.0 --seed 0] -> outputs/M1/W24/dev/E3_chr18_f<frac>_s<seed>.json
"""
import argparse
import json
import os

import numpy as np

import w24_lib as L
from dev_baselines import tib_hit

DMAX = (20_000, 50_000, 100_000)
KS = (1, 3)
THETAS = (-1.0, 0.0, 0.05, 0.1, 0.2, 0.3, 0.5)


def cums(post):
    from ncrna_annot.grammar import g1
    out = {}
    for k in ("l", "m"):
        for st in "+-":
            idx = [g1.IDX[f"{k}{p}{st}"] for p in ("E", "L", "I", "R")]
            v = np.zeros(len(post) + 1, np.float64)
            for a in range(0, len(post), 5_000_000):
                b = min(a + 5_000_000, len(post))
                v[a + 1:b + 1] = np.asarray(post[a:b][:, idx], np.float32).sum(1)
            out[(k, st)] = np.cumsum(v)
    return out


def span_mean(cs, k, st, s, e):
    c = cs[(k, st)]
    return (c[e] - c[s - 1]) / (e - s + 1)


def candidates(ev, dmax, K):
    out = []
    for st in "+-":
        T, P = ev["tss"][st], ev["pas"][st]
        if not len(T) or not len(P):
            continue
        for t, cnt in T:
            t = int(t)
            if st == "+":
                j = np.searchsorted(P, t, side="right")
                opts = [int(x) for x in P[j:j + K] if x - t <= dmax]
                spans = [(t, x) for x in opts]
            else:
                j = np.searchsorted(P, t, side="left")
                opts = [int(x) for x in P[max(0, j - K):j][::-1] if t - x <= dmax]
                spans = [(x, t) for x in opts]
            spans = [(s, e) for s, e in spans if e - s + 1 >= 200]
            if spans:
                out.append((st, int(cnt), spans))
    return out


def select(cands, cs, theta, tib=None):
    units = []
    for st, cnt, spans in cands:
        best = max(((span_mean(cs, "l", st, s, e) - span_mean(cs, "m", st, s, e), s, e) for s, e in spans))
        units.append((best[0], st, best[1], best[2], cnt))
    units.sort(key=lambda u: -u[0])
    return _greedy(units, theta, tib)


def select_evidence(cands, sig, rule, tib):
    """Evidence-only comparator on the same candidates: the end is the nearest site ('nearest'), the strongest QuantSeq
    site among the K ('strongest') or the farthest of the K ('farthest'); units accepted greedily by ssCAGE tags."""
    units = []
    for st, cnt, spans in cands:
        if rule == "nearest":
            s, e = spans[0]
        elif rule == "farthest":
            s, e = spans[-1]
        else:
            s, e = max(spans, key=lambda se: sig.get(((se[1] if st == "+" else se[0]), st), 0.0))
        units.append((float(cnt), st, s, e, cnt))
    units.sort(key=lambda u: -u[0])
    return _greedy(units, -1e9, tib)


def select_hybrid(cands, cs, sig, rule, theta, tib):
    """End chosen by the evidence rule (as select_evidence), units accepted greedily by ssCAGE tags, then kept only if the
    model's span score (lnc - mrna posterior) >= theta: the sequence model as a filter only."""
    units = []
    for st, cnt, spans in cands:
        if rule == "farthest":
            s, e = spans[-1]
        else:
            s, e = max(spans, key=lambda se: sig.get(((se[1] if st == "+" else se[0]), st), 0.0))
        sc = span_mean(cs, "l", st, s, e) - span_mean(cs, "m", st, s, e)
        units.append((float(cnt), st, s, e, cnt, sc))
    units.sort(key=lambda u: -u[0])
    acc, out = {"+": [], "-": []}, []
    for _, st, s, e, cnt, sc in units:
        if any(a <= e and s <= b for a, b in acc[st]) or (tib is not None and tib_hit(tib, st, s, e)):
            continue
        acc[st].append((s, e))
        if sc >= theta:
            out.append({"strand": st, "exons": [(s, e)], "score": sc, "tss_count": cnt})
    return out


def _greedy(units, theta, tib):
    acc, out = {"+": [], "-": []}, []
    for sc, st, s, e, cnt in units:
        if sc < theta:
            break
        if any(a <= e and s <= b for a, b in acc[st]):
            continue
        if tib is not None and tib_hit(tib, st, s, e):
            continue
        acc[st].append((s, e))
        out.append({"strand": st, "exons": [(s, e)], "score": sc, "tss_count": cnt})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frac", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    import simple_det as D
    C = "chr18"
    truth = L.load_truth(C)
    tib = D.read_tiberius(C)
    ev = L.evidence(C, a.frac, a.seed)
    sig = L.quant_signal(C)
    cfg = json.load(open(f"{L.ROOT}/configs/v20.json"))
    res = {"frac": a.frac, "seed": a.seed, "grid": {}}
    for run in cfg["targets"]["chr18"]["runs"]:
        post = np.load(f"{L.ROOT}/outputs/PHASE1/V20/read/{C}/{run}/post.npy", mmap_mode="r")
        cs = cums(post)
        for dm in DMAX:
            for K in KS:
                cands = candidates(ev, dm, K)
                for th in THETAS:
                    for var, tb in (("model", None), ("model+tib", tib)):
                        r = L.score_unknown(truth, select(cands, cs, th, tb), C)
                        r.pop("S50_ids")
                        res["grid"].setdefault(f"d{dm}_K{K}_th{th}_{var}", {})[run] = r
                if K > 1:
                    for rule in ("farthest", "strongest"):
                        for th in (-1.0, -0.05, 0.0, 0.05, 0.1):
                            r = L.score_unknown(truth, select_hybrid(cands, cs, sig, rule, th, tib), C)
                            r.pop("S50_ids")
                            res["grid"].setdefault(f"d{dm}_K{K}_hy-{rule}_th{th}+tib", {})[run] = r
                if run == cfg["targets"]["chr18"]["runs"][0]:   # evidence-only comparators do not depend on the run
                    for rule in ("nearest", "strongest", "farthest"):
                        r = L.score_unknown(truth, select_evidence(cands, sig, rule, tib), C)
                        r.pop("S50_ids")
                        res["grid"][f"d{dm}_K{K}_ev-{rule}+tib"] = {"evidence_only": r}
        del cs
    summ = {}
    for key, rr in res["grid"].items():
        summ[key] = {m: float(np.mean([r[m] for r in rr.values() if r[m] is not None])) for m in ("S50_sn", "S50_pr_known", "S50_pr", "n_pred")}
    res["mean"] = summ
    os.makedirs(f"{L.W}/dev", exist_ok=True)
    json.dump(res, open(f"{L.W}/dev/E3_{C}_f{a.frac}_s{a.seed}.json", "w"), indent=1)
    for key, v in sorted(summ.items(), key=lambda kv: -kv[1]["S50_sn"])[:25]:
        print(key, {m: round(x, 3) for m, x in v.items()})


if __name__ == "__main__":
    main()
