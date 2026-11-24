#!/usr/bin/env python3
"""AUTH v1 (docs/auth/v1/PLAN_V1.md, author approved 2026-11-24): shared code for the WTC11-DPI development and the
one-time H1 transfer test.

Evidence (plan 5.1): CAGE CTSS counted per sample inside the FANTOM5 hg38 DPI peaks (v1_check.dpi_clusters: pooled tags,
TSS = leftmost CTSS with the highest pooled count); QuantSeq sites as in AUTH (v1_check.quant).  Candidates, labels,
greedy, v41 mask and the 18 features are auth_lib unchanged; the features are cached raw and normalised on load:
  log_cage  -> log(1 + CPM), CPM = pooled tags / all CAGE tags of the cell (genome-wide, all samples) x 1e6
  n_samples -> fraction of the cell's CAGE samples with tags
  log_anti  -> log(1 + CPM) of the antisense tags within 500 bp
Cells: 'wtc11' = WTC11 ssCAGE (6 samples) in DPI peaks, AUTH truth (rep2 x rep3), outputs/AUTH_V1/wtc11_dpi;
       'h1'    = FANTOM5 H1 day00 (3 samples) in DPI peaks, LRGASP QuantSeq, truth rep1 x rep3, outputs/AUTH_V1/h1.
H1 test set (author 2026-11-24): chr2, 6, 10, 16, 22 + chr9, 12, 13, 20.
"""
import gzip
import json
import os
import re
import sys
from collections import defaultdict

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, f"{ROOT}/scripts/auth")
import auth_lib as A                                  # noqa: E402
import auth_truth as AT                               # noqa: E402
import v1_check as VC                                 # noqa: E402

TRAIN, TEST, SECONDARY = AT.TRAIN, AT.TEST, AT.SECONDARY
H1_TEST = TEST + SECONDARY
OUTV = f"{ROOT}/outputs/AUTH_V1"
CELLS = {
    "wtc11": {"out": f"{OUTV}/wtc11_dpi", "ctss": VC.CELLS["wtc11"]["ctss"], "quant": VC.CELLS["wtc11"]["quant"],
              "truth": VC.CELLS["wtc11"]["truth"], "check_key": "wtc11_dpi"},
    "h1": {"out": f"{OUTV}/h1", "ctss": VC.CELLS["h1"]["ctss"], "quant": VC.CELLS["h1"]["quant"],
           "truth": VC.CELLS["h1"]["truth"], "check_key": "h1"},
}
NORM_COLS = ("log_cage", "n_samples", "log_anti")


def use(cell):
    """Point auth_lib at the cell's output directory (evidence, truth, observed, pairs, seq caches)."""
    A.OUT = CELLS[cell]["out"]


def check_values():
    return json.load(open(f"{OUTV}/CHECK.json"))["cells"]


# ------------------------------------------------------------------ evidence on any chromosome set
def load_ctss(paths, chroms):
    """v1_check.load_ctss for a given chromosome set: {chrom: {strand: (positions, counts[n, samples])}}, genome-wide tags."""
    raw = {c: {st: [defaultdict(int) for _ in paths] for st in "+-"} for c in chroms}
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
    for c in chroms:
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


def load_dpi(chroms):
    pk = {c: {"+": [], "-": []} for c in chroms}
    with gzip.open(VC.DPI, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if f[0] in pk:
                pk[f[0]][f[5]].append((int(f[1]) + 1, int(f[2])))
    return {c: {st: sorted(v) for st, v in d.items()} for c, d in pk.items()}


def evidence_dicts(cell, chroms):
    """(evidence dicts {chrom: {'cage','quant','signal'}}, CAGE tags per sample genome-wide) -- nothing written."""
    cfg = CELLS[cell]
    for c in chroms:
        A.guard(c)
    ctss, lib = load_ctss(cfg["ctss"], chroms)
    peaks = load_dpi(chroms)
    ev = {}
    for c in chroms:
        q, sig = VC.quant(cfg["quant"], c)
        ev[c] = {"cage": VC.dpi_clusters(ctss[c], peaks[c]), "quant": q, "signal": sig}
    return ev, lib


def write_evidence(cell, chroms):
    """Write <out>/evidence/<c>.json for chromosomes that have none; returns the CAGE tags per sample."""
    out = CELLS[cell]["out"]
    os.makedirs(f"{out}/evidence", exist_ok=True)
    ev, lib = evidence_dicts(cell, chroms)
    for c in chroms:
        p = f"{out}/evidence/{c}.json"
        assert not os.path.exists(p), p
        json.dump(ev[c], open(p + ".tmp", "w"))
        os.replace(p + ".tmp", p)
    return lib


def write_observed(cell, chroms):
    """v1_check.observed_spans for a given chromosome set (the cell's two truth replicates)."""
    out = CELLS[cell]["out"]
    spans = defaultdict(set)
    for q, g in CELLS[cell]["truth"].values():
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
                if f[0] in chroms and len(f) == 9 and f[2] == "transcript" and \
                        re.search(r'transcript_id "([^"]+)"', f[8]).group(1) in want:
                    spans[f[0]].add((f[6], int(f[3]), int(f[4])))
    os.makedirs(f"{out}/evidence", exist_ok=True)
    for c in chroms:
        p = f"{out}/evidence/{c}.observed_rep23.json"
        assert not os.path.exists(p), p
        json.dump(sorted(spans[c]), open(p, "w"))


def write_truth(cell, chroms):
    """auth_truth.build with the cell's replicates (rep3 role = the deeper one) -> <out>/truth/<c>.json; lnc counts."""
    out = CELLS[cell]["out"]
    reps = VC.read_reps(CELLS[cell]["truth"])
    os.makedirs(f"{out}/truth", exist_ok=True)
    n = {}
    for c in chroms:
        loci = AT.build(reps, c)
        p = f"{out}/truth/{c}.json"
        assert not os.path.exists(p), p
        json.dump({"chrom": c, "exclusion_zone": [], "loci": loci, "rule": "auth_truth.build (v1_lib.write_truth)"}, open(p, "w"))
        n[c] = sum(x["k"] == "l" for x in loci)
    return n


# ------------------------------------------------------------------ normalised pairs
def normalise(P, lib_total, n_samples):
    X = P["X"].copy()
    cnt, anti = P["cnt"].astype(np.float64), P["anti"].astype(np.float64)
    assert np.allclose(X[:, A.FI["log_cage"]], np.log1p(cnt)) and np.allclose(X[:, A.FI["log_anti"]], np.log1p(anti))
    X[:, A.FI["log_cage"]] = np.log1p(cnt * 1e6 / lib_total)
    X[:, A.FI["n_samples"]] = X[:, A.FI["n_samples"]] / n_samples
    X[:, A.FI["log_anti"]] = np.log1p(anti * 1e6 / lib_total)
    return X


def load_pairs(cell, c, pc, lib):
    """auth_lib.build_pairs from the cell's caches (built if missing), with the normalised feature matrix."""
    use(cell)
    for k in (f"evidence/{c}.json", f"evidence/{c}.observed_rep23.json", f"truth/{c}.json"):
        assert os.path.exists(f"{A.OUT}/{k}"), (cell, k)        # never let auth_lib fall back to the WTC11 LRGASP files
    P = dict(A.build_pairs(c, pc[c]))
    P["X"] = normalise(P, float(sum(lib)), len(lib))
    return P


def ensure_seq_link(cell):
    out = CELLS[cell]["out"]
    if not os.path.exists(f"{out}/seq"):
        os.symlink(f"{ROOT}/outputs/AUTH/seq", f"{out}/seq")


def sha256(path):
    return A.sha256(path)
