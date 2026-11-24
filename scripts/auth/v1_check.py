#!/usr/bin/env python3
"""AUTH v1 step 2 (docs/auth/v1/FEASIBILITY_V1.md section 5; author approved 2026-11-24): descriptive check of GM12878 and
H1 against WTC11.  No training, no scoring of any learned model, no threshold.

Cells and roles
- wtc11      : the AUTH inputs as frozen (LRGASP ssCAGE count.bed clusters; caches in outputs/AUTH, read only).
- wtc11_dpi  : WTC11 ssCAGE CTSS counted in the FANTOM5 DPI peak set -- isolates the effect of the peak definition.
- gm12878    : FANTOM5 hCAGE CTSS (3 reps) in DPI peaks; QuantSeq ENCSR275KDL (1 biological rep, 2 technical: frac 0.5 / 1.0;
               ENCFF448CNL lists the same sites and column 9 values as ENCFF259EZO, so only ENCFF259EZO is read);
               truth = PacBio Sequel II ENCSR838WFC rep1 (in the role of WTC11 rep2) x rep2 (role rep3, the deeper one).
- h1         : FANTOM5 hCAGE CTSS (CD34+ HSC time course day00, 3 reps) in DPI peaks; LRGASP QuantSeq (3 reps);
               truth = LRGASP PacBio rep1 (role rep2) x rep3 (role rep3, the deeper one); rep2 (0.25 M FLNC) is left out.
Rules copied unchanged: truth = auth_truth.build (W24 rule, no exclusion zone); 'observed' (uncertain) = transcripts of
the two truth replicates with count >= 1, not 'Genomic'; evidence filters, candidates, features, labels and R0 =
auth_lib (CAGE pooled tags >= 2, QuantSeq frac >= 0.67, 20 nearest sites within 50 kb, v41 coding mask).
DPI clusters: counts = CTSS tags inside the peak per sample; TSS = leftmost CTSS with the highest pooled count.

Chromosomes: truth COUNTS on all autosomes + X (train / test / secondary); evidence, coverage, candidates, features and R0
on the 14 training chromosomes only (auth_lib guard).  Test-chromosome loci are never written.
  v1_check.py -> outputs/AUTH_V1/CHECK.json, outputs/AUTH_V1/<cell>/{truth,evidence,pairs}/
"""
import gzip
import json
import math
import os
import re
import sys
import time
from collections import defaultdict

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, f"{ROOT}/scripts/auth")
import auth_lib as A                                  # noqa: E402
import auth_truth as AT                               # noqa: E402
import wtc11_truth as WT                              # noqa: E402

L = A.L
TRAIN, TEST, SECONDARY, ALL = AT.TRAIN, AT.TEST, AT.SECONDARY, AT.ALL
V1 = f"{ROOT}/data/v1"
OUTV = f"{ROOT}/outputs/AUTH_V1"
AUTH_OUT = f"{ROOT}/outputs/AUTH"
DPI = f"{V1}/fantom5/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz"
WE = f"{ROOT}/data/wtc11/encode"
CELLS = {
    "wtc11": {"truth": {"rep2": (f"{WE}/ENCFF527KWD.tsv", f"{WE}/ENCFF798DWV.gtf.gz"),
                        "rep3": (f"{WE}/ENCFF486EKE.tsv", f"{WE}/ENCFF422CGK.gtf.gz")},
              "ctss": [f"{L.GEO}/{s}.ctss.bed.gz" for s in L.CAGE_SAMPLES],
              "quant": [f"{WE}/{a}.bed.gz" for a in L.QUANT]},
    "gm12878": {"truth": {"rep2": (f"{V1}/gm12878/encode/ENCFF079JHG.tsv", f"{V1}/gm12878/encode/ENCFF871ODI.gtf.gz"),
                          "rep3": (f"{V1}/gm12878/encode/ENCFF386OOG.tsv", f"{V1}/gm12878/encode/ENCFF993PKH.gtf.gz")},
                "ctss": [f"{V1}/fantom5/gm12878_rep{i}.ctss.bed.gz" for i in (1, 2, 3)],
                "quant": [f"{V1}/gm12878/encode/ENCFF259EZO.bed.gz"]},
    "h1": {"truth": {"rep2": (f"{V1}/h1/encode/ENCFF580WQM.tsv", f"{V1}/h1/encode/ENCFF721HVH.gtf.gz"),
                     "rep3": (f"{V1}/h1/encode/ENCFF014CCY.tsv", f"{V1}/h1/encode/ENCFF236CGS.gtf.gz")},
           "alt_rep": (f"{V1}/h1/encode/ENCFF563JZI.tsv", f"{V1}/h1/encode/ENCFF814SOM.gtf.gz"),
           "ctss": [f"{V1}/fantom5/h1_rep{i}.ctss.bed.gz" for i in (1, 2, 3)],
           "quant": [f"{V1}/h1/encode/ENCFF970YYF.bed.gz", f"{V1}/h1/encode/ENCFF138ZWV.bed.gz", f"{V1}/h1/encode/ENCFF215NMJ.bed.gz"]},
}
CELLS["wtc11_dpi"] = dict(CELLS["wtc11"])


# ------------------------------------------------------------------ truth (auth_truth.build with other files)
def read_reps(files):
    WT.CHROMS = ALL
    reps = {}
    for r, (q, g) in files.items():
        qq = WT.read_quant(q)
        tx = WT.read_gtf(g, set(qq))
        reps[r] = {}
        for tid, t in tx.items():
            k = WT.key_of(tid, t, qq[tid]["tnov"])
            if k is not None:
                reps[r][k] = {**t, **qq[tid], "tid": tid}
    return reps


def truth_counts(reps, out_dir=None):
    n = {}
    for c in ALL:
        loci = AT.build(reps, c)
        n[c] = sum(x["k"] == "l" for x in loci)
        if out_dir and c in TRAIN:
            os.makedirs(f"{out_dir}/truth", exist_ok=True)
            json.dump({"chrom": c, "exclusion_zone": [], "loci": loci, "rule": "auth_truth.build (v1_check.py)"},
                      open(f"{out_dir}/truth/{c}.json", "w"))
    return {"all_autosomes_X": sum(n[c] for c in ALL if c != "chrY"), "train": sum(n[c] for c in TRAIN),
            "test": sum(n[c] for c in TEST), "secondary": sum(n[c] for c in SECONDARY),
            "test_per_chrom": {c: n[c] for c in TEST}}


def observed_spans(files, out_dir):
    """auth_lib.observed with the cell's two truth replicates -> <out>/evidence/<c>.observed_rep23.json (train only)."""
    spans = defaultdict(set)
    for q, g in files.values():
        want = set()
        with open(q) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                f = dict(zip(hdr, line.rstrip("\n").split("\t")))
                if f["transcript_novelty"] != "Genomic" and int(f[hdr[-1]]) >= 1:
                    want.add(f["annot_transcript_id"])
        with gzip.open(g, "rt") as fh:
            for line in fh:
                f = line.split("\t", 8)
                if f[0] in TRAIN and len(f) == 9 and f[2] == "transcript" and \
                        re.search(r'transcript_id "([^"]+)"', f[8]).group(1) in want:
                    spans[f[0]].add((f[6], int(f[3]), int(f[4])))
    os.makedirs(f"{out_dir}/evidence", exist_ok=True)
    for c in TRAIN:
        json.dump(sorted(spans[c]), open(f"{out_dir}/evidence/{c}.observed_rep23.json", "w"))


# ------------------------------------------------------------------ evidence
def load_ctss(paths):
    """{chrom: {strand: (positions, counts[n, samples])}} on the training chromosomes; genome-wide tags per sample."""
    raw = {c: {st: [defaultdict(int) for _ in paths] for st in "+-"} for c in TRAIN}
    lib = []
    for k, p in enumerate(paths):
        tot = 0
        with gzip.open(p, "rt") as fh:
            for line in fh:
                f = line.split("\t")
                n = int(f[4])
                if not n:
                    continue
                tot += n
                if f[0] in raw:
                    raw[f[0]][f[5].rstrip()][k][int(f[2])] += n
        lib.append(tot)
    out = {}
    for c in TRAIN:
        out[c] = {}
        for st in "+-":
            pos = np.array(sorted(set().union(*[d.keys() for d in raw[c][st]])), np.int64)
            M = np.zeros((len(pos), len(paths)), np.int64)
            for k, d in enumerate(raw[c][st]):
                if d:
                    ks = np.fromiter(d.keys(), np.int64, len(d))
                    M[np.searchsorted(pos, ks), k] = np.fromiter(d.values(), np.int64, len(d))
            out[c][st] = (pos, M)
    return out, lib


def load_dpi():
    pk = {c: {"+": [], "-": []} for c in TRAIN}
    with gzip.open(DPI, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if f[0] in pk:
                pk[f[0]][f[5]].append((int(f[1]) + 1, int(f[2])))
    return {c: {st: sorted(v) for st, v in d.items()} for c, d in pk.items()}


def dpi_clusters(ctss, peaks):
    out = []
    for st in "+-":
        pos, M = ctss[st]
        pooled = M.sum(1)
        for s, e in peaks[st]:
            i0, i1 = np.searchsorted(pos, s, "left"), np.searchsorted(pos, e, "right")
            if i1 <= i0:
                continue
            cnt = M[i0:i1].sum(0)
            if not cnt.sum():
                continue
            j = int(np.argmax(pooled[i0:i1]))
            out.append({"start": s, "end": e, "strand": st, "tss": int(pos[i0 + j]), "counts": cnt.tolist(),
                        "mode_n": int(pooled[i0 + j])})
    return sorted(out, key=lambda x: (x["start"], x["end"], x["strand"]))


def quant(paths, c):
    """w24_lib._quant_sites + quant_signal for other files (max over files of column 7 frac and column 9 signal)."""
    fr, sig = {}, {}
    for p in paths:
        with gzip.open(p, "rt") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[0] != c:
                    continue
                k = (int(f[3].split(":")[1]), f[5])
                fr[k] = max(fr.get(k, 0.0), float(f[6]))
                sig[k] = max(sig.get(k, 0.0), float(f[8]))
    return ([{"pos": p_, "strand": st, "frac": v} for (p_, st), v in sorted(fr.items())],
            [[p_, st, v] for (p_, st), v in sorted(sig.items())])


# ------------------------------------------------------------------ per-cell summary on the training chromosomes
def summarise(cell, out_dir, pc, cage_lib):
    A.OUT = out_dir
    P, ts, ocv, cov, res, ev_n = {}, {}, {}, [], [], {"tss": 0, "pas": 0}
    ceil3 = ceil20 = nt = 0
    for c in TRAIN:
        tr = A.load_truth(c)
        ev, _, _ = A.evidence(c)
        ev_n["tss"] += sum(len(ev["tss"][st]) for st in "+-")
        ev_n["pas"] += sum(len(ev["pas"][st]) for st in "+-")
        cv = L.coverage(tr, ev)
        cov.append(cv)
        Pc = A.build_pairs(c, pc[c])
        P[c], ts[c], ocv[c] = Pc, tr["ts"], A.obs_cover(A.observed([c])[c])
        hit = Pc["label"] == 2
        ceil20 += len(set(Pc["tidx"][hit].tolist()))
        ceil3 += len(set(Pc["tidx"][hit & (Pc["rank"] < A.K_R0)].tolist()))
        nt += len(tr["ts"])
        res.append(A.score_units(Pc, A.method_units(Pc, "R0"), ts[c], ocv[c]))
    n = sum(x["n"] for x in cov)
    X = np.concatenate([P[c]["X"] for c in TRAIN])
    lab = np.concatenate([P[c]["label"] for c in TRAIN])
    cnt = np.concatenate([P[c]["cnt"] for c in TRAIN]).astype(np.float64)
    tssk = np.concatenate([P[c]["strand"] * 10 ** 9 + P[c]["row"] + i * 10 ** 10 for i, c in enumerate(TRAIN)])
    first = np.r_[True, tssk[1:] != tssk[:-1]]
    r0 = A.metrics(res)
    return {"truth_lnc_train": nt, "evidence_train": ev_n,
            "coverage_train": {k: sum(x[k] * x["n"] for x in cov) / n for k in ("tss", "pas", "both")},
            "candidate_ceiling_train": {"K3": ceil3 / nt, "K20": ceil20 / nt},
            "pairs_train": {"n": int(len(lab)), "true": int((lab == 2).sum()), "uncertain": int((lab == 1).sum()),
                            "wrong": int((lab == 0).sum())},
            "R0_train": r0, "cage_lib_tags": cage_lib,
            "cage_cpm_median_per_tss": float(np.median(cnt[first] / sum(cage_lib) * 1e6)) if cage_lib else None,
            "_X": X, "_lab": lab, "_cpm": np.log1p(cnt / sum(cage_lib) * 1e6) if cage_lib else None}


def feature_table(S):
    ref = S["wtc11"]
    mu, sd = ref["_X"].mean(0), ref["_X"].std(0) + 1e-12
    out = {}
    for cell, s in S.items():
        X, y = s["_X"], s["_lab"]
        tw = y != 1
        f = {}
        for j, name in enumerate(A.FEATURES):
            v = X[:, j]
            f[name] = {"median": {"true": float(np.median(v[y == 2])) if (y == 2).any() else None,
                                  "uncertain": float(np.median(v[y == 1])) if (y == 1).any() else None,
                                  "wrong": float(np.median(v[y == 0])) if (y == 0).any() else None},
                       "auc_true_vs_wrong": float(A.auc(v[tw], (y[tw] == 2).astype(int))) if tw.any() else None,
                       "smd_vs_wtc11": float((v.mean() - mu[j]) / sd[j])}
        if s["_cpm"] is not None:
            v = s["_cpm"]
            f["log_cage_cpm"] = {"median": {k: float(np.median(v[y == i])) for k, i in (("true", 2), ("uncertain", 1), ("wrong", 0))},
                                 "auc_true_vs_wrong": float(A.auc(v[tw], (y[tw] == 2).astype(int)))}
        out[cell] = f
    return out


def main():
    t0 = time.time()
    res = {"doc": __doc__.split("  v1_check.py")[0], "cells": {}}
    pc = A.v41_pc(TRAIN)
    peaks = load_dpi()
    res["dpi_peaks_train"] = sum(len(v) for d in peaks.values() for v in d.values())
    S = {}
    for cell, cfg in CELLS.items():
        out_dir = AUTH_OUT if cell == "wtc11" else f"{OUTV}/{cell}"
        r = {}
        if cell in ("gm12878", "h1"):
            reps = read_reps(cfg["truth"])
            r["truth_observed_tx"] = {k: len(v) for k, v in reps.items()}
            r["truth_lnc_counts"] = truth_counts(reps, out_dir)
            if "alt_rep" in cfg:        # H1: the other two pairings (counts only)
                alt = read_reps({"rep2": cfg["alt_rep"]})["rep2"]
                r["truth_lnc_alt_pairs"] = {
                    "rep1_x_rep2": truth_counts({"rep2": alt, "rep3": reps["rep2"]})["all_autosomes_X"],
                    "rep2_x_rep3": truth_counts({"rep2": alt, "rep3": reps["rep3"]})["all_autosomes_X"]}
            observed_spans(cfg["truth"], out_dir)
        elif cell == "wtc11_dpi":
            for c in TRAIN:          # same truth / observed files as the frozen WTC11 run
                os.makedirs(f"{out_dir}/truth", exist_ok=True)
                os.makedirs(f"{out_dir}/evidence", exist_ok=True)
                for src, dst in ((f"{AUTH_OUT}/truth/{c}.json", f"{out_dir}/truth/{c}.json"),
                                 (f"{AUTH_OUT}/evidence/{c}.observed_rep23.json", f"{out_dir}/evidence/{c}.observed_rep23.json")):
                    if not os.path.exists(dst):
                        os.symlink(src, dst)
        else:
            r["truth_lnc_counts"] = json.load(open(f"{AUTH_OUT}/TRUTH_COUNTS.json"))["got"]
        lib = []
        if cell != "wtc11":
            if not os.path.exists(f"{out_dir}/seq"):
                os.symlink(f"{AUTH_OUT}/seq", f"{out_dir}/seq")
            ctss, lib = load_ctss(cfg["ctss"])
            for c in TRAIN:
                p = f"{out_dir}/evidence/{c}.json"
                if not os.path.exists(p):
                    q, sig = quant(cfg["quant"], c)
                    json.dump({"cage": dpi_clusters(ctss[c], peaks[c]), "quant": q, "signal": sig}, open(p, "w"))
            del ctss
        else:
            lib = load_ctss(cfg["ctss"])[1]
        for c in TRAIN:
            assert os.path.exists(f"{out_dir}/evidence/{c}.json") and os.path.exists(f"{out_dir}/truth/{c}.json"), (cell, c)
        s = summarise(cell, out_dir, pc, lib)
        S[cell] = s
        r.update({k: v for k, v in s.items() if not k.startswith("_")})
        res["cells"][cell] = r
        print(cell, json.dumps({k: r[k] for k in ("coverage_train", "candidate_ceiling_train", "pairs_train")}), flush=True)
    A.OUT = AUTH_OUT
    res["features_train"] = feature_table(S)
    res["seconds"] = round(time.time() - t0, 1)
    json.dump(res, open(f"{OUTV}/CHECK.json", "w"), indent=1)
    print(json.dumps({c: {k: v for k, v in r.items() if k in ("truth_lnc_counts", "R0_train")} for c, r in res["cells"].items()}, indent=1))


if __name__ == "__main__":
    main()
