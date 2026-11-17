#!/usr/bin/env python3
"""M1 experiment 3, stage B (docs/m1/PLAN_M1_AG_B.md): AlphaGenome CAGE tracks of cell types in which a start is
expressed (matched) vs tracks of cell types in which it has no tag (unmatched), same assay (hCAGE), same start.

Starts: rows of outputs/M1/AG/tss.json with a fold model and FANTOM5 CAGE evidence (same-strand robust peak whose
representative position is within 200 bp of the TSS, as end_metrics / DEDEND).  A start's expression in a FANTOM5
library = tags of those peaks / mapped tags * 1e6 (CPM).  Library -> AlphaGenome curie: the library's FF sample term
and its FF ancestors, relations derives_from (CL / UBERON) and is_a (CL / UBERON / EFO); libraries of the categories
'primary cells' and 'tissues' only (no time courses, perturbations; cell lines do not map to the AlphaGenome EFO
terms).  Group expression = mean CPM over the group's libraries.  matched: >= 1 CPM; unmatched: 0 tags in every
library of the group.  Per start: Loc1k over matched hCAGE tracks (fraction of tracks whose argmax is within 1 kb)
minus the same over unmatched tracks.
  ag_stageB.py prep  -> outputs/M1/AG/stageB/{expr.npz, PREP_B.json}   (FANTOM5 side only; reads no AlphaGenome output)
  ag_stageB.py read  -> outputs/M1/AG/stageB/RESULT_B.json             (needs stageB/FROZEN = sha256 of expr.npz)
"""
import csv
import gzip
import hashlib
import json
import os
import re
import sys
from collections import Counter, defaultdict

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "scripts/m1/ag"))
import f5_lib as F  # noqa: E402

O = f"{ROOT}/outputs/M1/AG"
B = f"{O}/stageB"
F5 = f"{ROOT}/data/fantom5_stageB"
BED = f"{ROOT}/data/bioval/raw/fantom5_cage_peaks.bed.gz"
COUNTS = f"{F5}/hg38_fair+new_CAGE_peaks_phase1and2_counts.osc.txt.gz"
RAW = f"{O}/colab_auto/raw"
CHROMS = ("chr8", "chr9", "chr18", "chr20")
GROUPS = {"main_chr8_chr9": ("chr8", "chr9"), "chr18": ("chr18",), "chr20": ("chr20",)}
EV_BP, TOL, DECOY, N_BOOT = 200, 1000, 10_000, 2000
CPM_MATCH = 1.0
CATS = ("primary cells", "tissues")


def sha_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def tracks(strand):
    return list(csv.DictReader(open(f"{RAW}/tracks_cage_{strand}.csv")))


def lib_curies(T, ff, seen=None, depth=0):
    seen = set() if seen is None else seen
    if ff in seen or ff not in T or depth > 8:
        return set()
    seen.add(ff)
    t = T[ff]
    out = {r[1] for r in t["rel"] if r[0] == "derives_from"} | {x for x in t["is_a"] if x.split(":")[0] in ("EFO", "CL", "UBERON")}
    for p in t["is_a"]:
        if p.startswith("FF:"):
            out |= lib_curies(T, p, seen, depth + 1)
    return out


def prep():
    os.makedirs(B, exist_ok=True)
    rows = [r for r in json.load(open(f"{O}/tss.json")) if r["chrom"] in CHROMS and r["model"] and r["TSS_ev"]]
    hc = [t for t in tracks("+") if t["Assay title"] == "hCAGE"]
    assert sorted(t["ontology_curie"] for t in hc) == sorted(t["ontology_curie"] for t in tracks("-") if t["Assay title"] == "hCAGE")
    curies = [t["ontology_curie"] for t in hc]
    assert len(set(curies)) == len(curies)
    # peaks within 200 bp of each start (same strand, representative position = thick midpoint, as S3.load_peaks)
    pk = defaultdict(lambda: {"+": ([], []), "-": ([], [])})
    with gzip.open(BED, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] in CHROMS and f[5] in "+-":
                pk[f[0]][f[5]][0].append((int(f[6]) + int(f[7])) // 2)
                pk[f[0]][f[5]][1].append(f[3])
    for c in pk:
        for s in "+-":
            pos, nm = pk[c][s]
            o = np.argsort(pos)
            pk[c][s] = (np.asarray(pos, np.int64)[o], [nm[i] for i in o])
    peaks_of = {}
    for r in rows:
        pos, nm = pk[r["chrom"]][r["s"]]
        i, j = np.searchsorted(pos, r["TSS"] - EV_BP), np.searchsorted(pos, r["TSS"] + EV_BP, side="right")
        peaks_of[r["id"]] = nm[i:j]
    need = {p for v in peaks_of.values() for p in v}
    # libraries
    T = F.read_obo(f"{F5}/ff-phase2-170801.obo.txt")
    sd = F.read_xlsx(f"{F5}/HumanSamples2.0.sdrf.xlsx")
    sdrf = {r[0]: {"ff": r[1], "cat": r[4], "desc": r[2]} for r in sd[1:]}
    counts, mapped, cols = {}, None, None
    with gzip.open(COUNTS, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            key = line[:line.find("\t")]
            if key not in need and key not in ("00Annotation", "01STAT:MAPPED"):
                continue
            f = line.rstrip("\n").split("\t")
            if key == "00Annotation":
                cols = f[1:]
            elif key == "01STAT:MAPPED":
                mapped = np.asarray(f[1:], np.float64)
            else:
                counts[key] = np.asarray(f[1:], np.int64)
    ext = [re.search(r"\.(\d+-\d+[A-Z]\d+)\.hg38", c) for c in cols]
    lib_ok, lib_cur, why = [], [], Counter()
    for k, m in enumerate(ext):
        s = sdrf.get(m.group(1)) if m else None
        if s is None:
            why["not_in_sdrf"] += 1
            lib_cur.append(set())
            continue
        if s["cat"] not in CATS:
            why[f"category:{s['cat']}"] += 1
            lib_cur.append(set())
            continue
        cs = lib_curies(T, s["ff"]) & set(curies)
        why["used" if cs else "no_ag_curie"] += 1
        lib_cur.append(cs)
        if cs:
            lib_ok.append(k)
    g_libs = [[k for k in lib_ok if c in lib_cur[k]] for c in curies]
    expr = np.zeros((len(rows), len(curies)), np.float32)
    zero = np.zeros((len(rows), len(curies)), bool)
    for i, r in enumerate(rows):
        tags = sum((counts[p] for p in peaks_of[r["id"]] if p in counts), np.zeros(len(cols), np.int64))
        cpm = tags / mapped * 1e6
        for g, ks in enumerate(g_libs):
            if ks:
                expr[i, g] = cpm[ks].mean()
                zero[i, g] = not tags[ks].any()
    has = np.array([bool(ks) for ks in g_libs])
    matched = (expr >= CPM_MATCH) & has
    unmatched = zero & has
    ids = np.array([r["id"] for r in rows])
    np.savez_compressed(f"{B}/expr.npz", ids=ids, curies=np.array(curies), expr=expr, zero=zero, has=has,
                        n_libs=np.array([len(k) for k in g_libs]))
    by = {r["id"]: r for r in rows}
    prep_out = {"plan": "docs/m1/PLAN_M1_AG_B.md", "inputs": {p: sha_file(p) for p in (COUNTS, BED, f"{F5}/HumanSamples2.0.sdrf.xlsx",
                                                                                      f"{F5}/ff-phase2-170801.obo.txt", f"{O}/tss.json")},
                "libraries": dict(why), "hcage_tracks": len(curies), "groups_with_libs": int(has.sum()),
                "libs_per_group": dict(sorted(Counter(int(x) for x in [len(k) for k in g_libs] if x).items())),
                "peaks_needed": len(need), "peaks_found": len(counts), "sets": {}}
    for g, cs in GROUPS.items():
        for kind in ("l", "m"):
            sel = [i for i, r in enumerate(rows) if r["chrom"] in cs and r["k"] == kind]
            ev = [i for i in sel if matched[i].any() and unmatched[i].any()]
            prep_out["sets"][f"{g}|{kind}"] = {
                "n_cage_starts": len(sel), "n_evaluable": len(ev),
                "n_no_peak_found": sum(1 for i in sel if not peaks_of[rows[i]["id"]]),
                "matched_tracks_median": float(np.median(matched[ev].sum(1))) if ev else None,
                "unmatched_tracks_median": float(np.median(unmatched[ev].sum(1))) if ev else None,
                "n_test_only": sum(1 for i in ev if by[ids[i]]["test_only"])}
    json.dump(prep_out, open(f"{B}/PREP_B.json", "w"), indent=1)
    print(json.dumps(prep_out, indent=1))
    print("expr.npz sha256", sha_file(f"{B}/expr.npz"))


def boot_paired(d, seed=0):
    d = np.asarray(d, float)
    if not len(d):
        return None
    rng = np.random.default_rng(seed)
    m = [d[rng.integers(0, len(d), len(d))].mean() for _ in range(N_BOOT)]
    return {"n": len(d), "diff": float(d.mean()), "ci95": np.percentile(m, [2.5, 97.5]).tolist()}


def read():
    want = open(f"{B}/FROZEN").read().strip()
    assert want == sha_file(f"{B}/expr.npz"), "expr.npz differs from the frozen one"
    z = np.load(f"{B}/expr.npz")
    ids, curies, expr, zero, has = list(z["ids"]), list(z["curies"]), z["expr"], z["zero"], z["has"]
    rows = {r["id"]: r for r in json.load(open(f"{O}/tss.json"))}
    lnc_all = defaultdict(list)
    for r in rows.values():
        if r["k"] == "l":
            for g, cs in GROUPS.items():
                if r["chrom"] in cs:
                    lnc_all[g].append(r["cpg"])
    q2 = {g: float(np.quantile(v, 0.5)) for g, v in lnc_all.items()}
    col = {s: {t["ontology_curie"]: j for j, t in enumerate(tracks(s)) if t["Assay title"] == "hCAGE"} for s in "+-"}
    arg = {}
    for m in ("fold-0", "fold-1", "fold-2", "fold-3"):
        a = np.load(f"{RAW}/{m}.npz")
        for j, rid in enumerate(a["rid"]):
            arg[str(rid)] = (a["cage_arg"][j].astype(int), int(a["agg_cage_mean_arg"][j]))
    per = []
    for i, rid in enumerate(ids):
        r = rows[rid]
        if rid not in arg:
            continue
        ca, mean_arg = arg[rid]
        mt = [col[r["s"]][curies[g]] for g in np.flatnonzero((expr[i] >= 1.0) & has)]
        ut = [col[r["s"]][curies[g]] for g in np.flatnonzero(zero[i] & has)]
        if not mt or not ut:
            continue
        hi = [col[r["s"]][curies[g]] for g in np.flatnonzero((expr[i] >= 5.0) & has)]
        lo = [col[r["s"]][curies[g]] for g in np.flatnonzero((expr[i] >= 1.0) & (expr[i] < 5.0) & has)]
        f = lambda t, d=0: float(np.mean(np.abs(ca[t] - d) <= TOL)) if t else None  # noqa: E731
        per.append({"id": rid, "chrom": r["chrom"], "k": r["k"], "test_only": r["test_only"], "cpg": r["cpg"],
                    "matched": f(mt), "unmatched": f(ut), "n_m": len(mt), "n_u": len(ut),
                    "matched_hi": f(hi), "matched_lo": f(lo),
                    "matched_decoy": (f(mt, DECOY) + f(mt, -DECOY)) / 2, "unmatched_decoy": (f(ut, DECOY) + f(ut, -DECOY)) / 2,
                    "mean_track": float(abs(mean_arg) <= TOL), "any_matched": float(np.abs(ca[mt]).min() <= TOL)})
    res = {"plan": "docs/m1/PLAN_M1_AG_B.md", "frozen_expr_sha256": want, "groups": {}}
    for g, cs in GROUPS.items():
        out = {"lnc_cpg_median": q2[g]}
        for kind, label in (("l", "lnc_cage"), ("m", "mrna_cage")):
            base = [p for p in per if p["chrom"] in cs and p["k"] == kind]
            subsets = {"all": base, "test_only": [p for p in base if p["test_only"]]}
            if kind == "l":
                subsets["low_cpg"] = [p for p in base if p["cpg"] <= q2[g]]
                subsets["high_cpg"] = [p for p in base if p["cpg"] > q2[g]]
            for sn, ss in subsets.items():
                blk = {"n": len(ss)}
                for key in ("matched", "unmatched", "matched_hi", "matched_lo", "matched_decoy", "unmatched_decoy",
                            "mean_track", "any_matched"):
                    v = [p[key] for p in ss if p[key] is not None]
                    blk[key] = float(np.mean(v)) if v else None
                blk["matched_minus_unmatched"] = boot_paired([p["matched"] - p["unmatched"] for p in ss])
                blk["matched_minus_mean_track"] = boot_paired([p["matched"] - p["mean_track"] for p in ss])
                blk["matched_minus_decoy"] = boot_paired([p["matched"] - p["matched_decoy"] for p in ss])
                pr = [p for p in ss if p["matched_hi"] is not None and p["matched_lo"] is not None]
                blk["hi_minus_lo_cpm"] = boot_paired([p["matched_hi"] - p["matched_lo"] for p in pr])
                blk["median_n_matched"] = float(np.median([p["n_m"] for p in ss])) if ss else None
                blk["median_n_unmatched"] = float(np.median([p["n_u"] for p in ss])) if ss else None
                out[f"{label}|{sn}"] = blk
        res["groups"][g] = out
    x = res["groups"]["main_chr8_chr9"]["lnc_cage|all"]["matched_minus_unmatched"]
    if x["diff"] >= 0.10 and x["ci95"][0] > 0:
        o = "C_celltype_dependent"
    elif x["diff"] <= 0.03 and x["ci95"][1] < 0.10:
        o = "not_dependent"
    else:
        o = "grey"
    res["verdict"] = {"primary": x, "outcome": o}
    json.dump(res, open(f"{B}/RESULT_B.json", "w"), indent=1)
    print(json.dumps(res["verdict"], indent=1))


if __name__ == "__main__":
    {"prep": prep, "read": read}[sys.argv[1]]()
