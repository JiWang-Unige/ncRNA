#!/usr/bin/env python3
"""AUTH (docs/auth/PLAN_AUTH.md sections 5-6): shared code -- evidence, candidate pairs, labels, features, the selection
rules (model, R0-R6, P1), scoring (W24 S50 one-to-one matcher + the W24 'unknown' rule) and the paired 1 Mb block bootstrap.

Evidence   = w24_lib.evidence(chrom, 1.0, 0): ssCAGE clusters with pooled tags >= 2, TSS = the cluster's ctss mode;
             QuantSeq sites supported by >= 2 of 3 replicates.  Rebuilt here in one pass with two extra cluster fields
             (mode_n, the pooled ctss tags at the mode) and cached in outputs/AUTH/evidence/ (W24 caches untouched).
Candidates = each TSS x the K_WIDE = 20 nearest downstream same-strand QuantSeq sites within 50 kb, spans < 200 bp
             dropped (= e3_anchor.candidates(ev, 50_000, 20)); 'rank' = position among the 20 nearest (0 = nearest).
Labels     = 2 true: reciprocal span overlap >= 0.5 with a truth lnc locus (same strand);
             1 uncertain: not true, same-strand overlap with a transcript seen in rep2 or rep3 (w24_lib.observed_any rule);
             0 wrong: neither.
Features   = FEATURES below (section 5.2; no long reads, no lncRNA annotation, no sequence model).
Selection  = per TSS one pair chosen by a rule, units accepted greedily by a key (no same-strand overlap; a unit that
             overlaps a GENCODE v41 protein_coding gene on its strand is dropped and does not block: e3_anchor._greedy),
             then the top fraction q of the accepted units (per chromosome) for the methods with a q.
Scoring    = S50 (one-to-one maximum matching, reciprocal span overlap >= 0.5) against truth lnc loci; a prediction that
             is not matched but overlaps a rep2 / rep3 transcript on its strand is 'unknown' (w24_lib.score_unknown).
             recall = matched loci / loci; precision (known) = TP / (n - unknown); strict precision = TP / n.

Chromosome guard: evidence, candidates, features, labels and truth files are open for the training chromosomes only;
auth_score.py --mode test calls open_test() after it has verified outputs/AUTH/FROZEN.
"""
import bisect
import gzip
import json
import math
import os
import re
import subprocess
import sys
from collections import defaultdict

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
for _d in ("scripts/m1/w24", "scripts/phase1/s3", "scripts/auth"):
    sys.path.insert(0, f"{ROOT}/{_d}")
import w24_lib as L                                  # noqa: E402  (also puts src / phase1 / m1 on the path)
import b1_strict as B                                # noqa: E402
import s3_features as S3                             # noqa: E402
import simple_det as D                               # noqa: E402
from e3_anchor import candidates, select_evidence    # noqa: E402
from auth_truth import TRAIN, TEST, SECONDARY        # noqa: E402

OUT = f"{ROOT}/outputs/AUTH"
DMAX, K_WIDE, K_R0 = 50_000, 20, 3
MIN_LEN = 200
A_MAX = 0.6                       # R3: downstream 20 nt A content above this = internal priming
PAS_HEX = ("AATAAA", "ATTAAA")     # R3 / feature: hexamer fully inside 10-35 nt upstream of the site
PAS_RC = ("TTTATT", "TTTAAT")
WIN_ANTI, BIDIR = 500, 0.5        # R4: antisense ssCAGE tags within 500 bp >= half of the TSS's tags
WIN_GC = 500                      # +-500 bp around the TSS (CpG island: GC >= 0.5 and o/e >= 0.6, M1 experiment 2)
BLOCK, N_BOOT = 1_000_000, 2000
LOGIT_C = 1.0                     # simple_det.fit_logit default
FEATURES = ("log_cage", "n_samples", "log_width", "mode_share", "log_anti", "cgi", "gc500",
            "quant_nrep", "log_quant_signal", "a_down20", "pas_hex", "end_rank",
            "log_len", "log_cage_inside", "log_quant_inside", "pc_same", "pc_anti", "repeat_frac")
FI = {f: i for i, f in enumerate(FEATURES)}
STR = "+-"

_OPEN = set(TRAIN)


def open_test():
    _OPEN.update(TEST)
    _OPEN.update(SECONDARY)


def guard(c):
    assert c in _OPEN, f"{c} is locked until the freeze (training chromosomes only)"


# ------------------------------------------------------------------ interval helper
def cover(iv):
    """Merged 1-based closed intervals (a, b, C) for vectorised bp-coverage / overlap queries."""
    iv = np.asarray(iv, np.int64).reshape(-1, 2)
    if not len(iv):
        z = np.zeros(0, np.int64)
        return z, z, z
    iv = iv[np.argsort(iv[:, 0], kind="stable")]
    a, b = [int(iv[0, 0])], [int(iv[0, 1])]
    for s, e in iv[1:]:
        if s <= b[-1] + 1:
            b[-1] = max(b[-1], int(e))
        else:
            a.append(int(s))
            b.append(int(e))
    a, b = np.array(a, np.int64), np.array(b, np.int64)
    return a, b, np.concatenate([[0], np.cumsum(b - a + 1)])[:-1]


def cov_bp(cv, s, e):
    """bp of [s, e] (1-based closed, arrays) covered by the merged intervals."""
    a, b, C = cv
    s, e = np.asarray(s, np.int64), np.asarray(e, np.int64)
    if not len(a):
        return np.zeros(np.broadcast(s, e).shape, np.int64)

    def F(x):
        j = np.searchsorted(a, x, side="right") - 1
        jj = np.maximum(j, 0)
        return np.where(j >= 0, C[jj] + np.minimum(b[jj], x) - a[jj] + 1, 0)
    return F(e) - F(s - 1)


# ------------------------------------------------------------------ truth, sequence, annotation
def load_truth(c):
    guard(c)
    d = json.load(open(f"{OUT}/truth/{c}.json"))
    lnc = [x for x in d["loci"] if x["k"] == "l"]
    ts = sorted([(x["id"], x["strand"], x["start"], x["end"]) for x in lnc], key=lambda t: t[2])
    return {"ts": ts, "lnc": {x["id"]: x for x in lnc}, "zone": []}


def extract_seq(chroms):
    """Soft-masked chromosomes from data/human/hg38.fa.gz -> outputs/AUTH/seq/<c>.fa (lines copied as they are)."""
    want = [c for c in chroms if not os.path.exists(f"{OUT}/seq/{c}.fa")]
    for c in want:
        guard(c)
    if not want:
        return
    os.makedirs(f"{OUT}/seq", exist_ok=True)
    p = subprocess.Popen(["zcat", f"{ROOT}/data/human/hg38.fa.gz"], stdout=subprocess.PIPE, bufsize=1 << 24)
    fh, cur = None, None
    for line in p.stdout:
        if line.startswith(b">"):
            if fh:
                fh.close()
                os.replace(f"{OUT}/seq/{cur}.fa.tmp", f"{OUT}/seq/{cur}.fa")
                fh = None
            cur = line[1:].split()[0].decode()
            if cur in want:
                fh = open(f"{OUT}/seq/{cur}.fa.tmp", "wb")
                fh.write(line)
        elif fh:
            fh.write(line)
    if fh:
        fh.close()
        os.replace(f"{OUT}/seq/{cur}.fa.tmp", f"{OUT}/seq/{cur}.fa")
    assert p.wait() == 0
    missing = [c for c in want if not os.path.exists(f"{OUT}/seq/{c}.fa")]
    assert not missing, missing


def load_seq(c):
    guard(c)
    with open(f"{OUT}/seq/{c}.fa", "rb") as fh:
        assert fh.readline().startswith(b">" + c.encode())
        return b"".join(x.rstrip() for x in fh)


def chrom_len(c):
    p = f"{OUT}/seq/{c}.len"
    if not os.path.exists(p):
        open(p, "w").write(str(len(load_seq(c))))
    return int(open(p).read())


def v41_pc(chroms):
    """GENCODE v41 protein_coding gene spans {chrom: {strand: Nx2}} (as auth_c.v41_pc)."""
    out = {c: {"+": [], "-": []} for c in chroms}
    with gzip.open(f"{ROOT}/data/human/gencode.v41.annotation.gtf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if f[0] in out and f[2] == "gene" and 'gene_type "protein_coding"' in f[8]:
                out[f[0]][f[6]].append((int(f[3]), int(f[4])))
    return {c: {st: np.array(sorted(v), np.int64).reshape(-1, 2) for st, v in d.items()} for c, d in out.items()}


def observed(chroms):
    """w24_lib.observed_any for several chromosomes in one pass (spans of rep2 / rep3 transcripts with count >= 1, not
    'Genomic'), cached in outputs/AUTH/evidence/<c>.observed_rep23.json."""
    for c in chroms:
        guard(c)
    todo = {c for c in chroms if not os.path.exists(f"{OUT}/evidence/{c}.observed_rep23.json")}
    if todo:
        os.makedirs(f"{OUT}/evidence", exist_ok=True)
        spans = defaultdict(set)
        for q, g in (("ENCFF527KWD.tsv", "ENCFF798DWV.gtf.gz"), ("ENCFF486EKE.tsv", "ENCFF422CGK.gtf.gz")):
            want = set()
            with open(f"{L.ENC}/{q}") as fh:
                hdr = fh.readline().rstrip("\n").split("\t")
                for line in fh:
                    f = dict(zip(hdr, line.rstrip("\n").split("\t")))
                    if f["transcript_novelty"] != "Genomic" and int(f[hdr[-1]]) >= 1:
                        want.add(f["annot_transcript_id"])
            with gzip.open(f"{L.ENC}/{g}", "rt") as fh:
                for line in fh:
                    f = line.split("\t", 8)
                    if f[0] in todo and len(f) == 9 and f[2] == "transcript" and \
                            re.search(r'transcript_id "([^"]+)"', f[8]).group(1) in want:
                        spans[f[0]].add((f[6], int(f[3]), int(f[4])))
        for c in todo:
            json.dump(sorted(spans[c]), open(f"{OUT}/evidence/{c}.observed_rep23.json", "w"))
    return {c: [tuple(x) for x in json.load(open(f"{OUT}/evidence/{c}.observed_rep23.json"))] for c in chroms}


def obs_cover(spans):
    return {st: cover([(s, e) for st2, s, e in spans if st2 == st]) for st in STR}


# ------------------------------------------------------------------ evidence
def _cage_ext(chrom):
    """w24_lib._cage_clusters (same clusters and TSS = ctss mode, leftmost on ties) plus mode_n."""
    cl = {}
    for k, smp in enumerate(L.CAGE_SAMPLES):
        with gzip.open(f"{L.GEO}/{smp}.count.bed.gz", "rt") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[0] != chrom:
                    continue
                key = (int(f[1]), int(f[2]), f[5])
                cl.setdefault(key, np.zeros(len(L.CAGE_SAMPLES), np.int64))[k] = int(f[6])
    ctss = defaultdict(int)
    for smp in L.CAGE_SAMPLES:
        with gzip.open(f"{L.GEO}/{smp}.ctss.bed.gz", "rt") as fh:
            for line in fh:
                if not line.startswith(chrom + "\t"):
                    continue
                f = line.rstrip("\n").split("\t")
                if int(f[4]):
                    ctss[(int(f[2]), f[5])] += int(f[4])
    out = []
    for (s, e, st), cnt in sorted(cl.items()):
        best, mode = -1, (s + 1 + e) // 2
        for x in range(s + 1, e + 1):
            c = ctss.get((x, st), 0)
            if c > best:
                best, mode = c, x
        out.append({"start": s + 1, "end": e, "strand": st, "tss": mode, "counts": cnt.tolist(), "mode_n": best})
    return out


def evidence_raw(c):
    guard(c)
    p = f"{OUT}/evidence/{c}.json"
    if not os.path.exists(p):
        os.makedirs(f"{OUT}/evidence", exist_ok=True)
        sig = L.quant_signal(c)
        raw = {"cage": _cage_ext(c), "quant": L._quant_sites(c), "signal": [[p_, st, v] for (p_, st), v in sorted(sig.items())]}
        json.dump(raw, open(p + ".tmp", "w"))
        os.replace(p + ".tmp", p)
    return json.load(open(p))


def evidence(c):
    """(ev, rows, raw): ev exactly as w24_lib.evidence(c, 1.0, 0); rows[strand] = sorted (tss, tags, cluster index)."""
    raw = evidence_raw(c)
    rows = {"+": [], "-": []}
    for i, cl in enumerate(raw["cage"]):
        k = int(sum(cl["counts"]))
        if k >= L.CAGE_MIN:
            rows[cl["strand"]].append((cl["tss"], k, i))
    pas = {"+": [], "-": []}
    for q in raw["quant"]:
        if q["frac"] >= L.QUANT_MIN_FRAC:
            pas[q["strand"]].append(q["pos"])
    rows = {st: sorted(v) for st, v in rows.items()}
    ev = {"tss": {st: np.array([r[:2] for r in v], np.int64).reshape(-1, 2) for st, v in rows.items()},
          "pas": {st: np.array(sorted(v), np.int64) for st, v in pas.items()}}
    return ev, rows, raw


# ------------------------------------------------------------------ candidate pairs, labels, features
def label_pairs(st_arr, s, e, ts, ocv):
    """(label, truth index): 2 true / 1 uncertain / 0 wrong (see module doc)."""
    lab = np.zeros(len(s), np.int8)
    tid = np.full(len(s), -1, np.int64)
    for j, (_, gst, gs, ge) in enumerate(ts):
        m = st_arr == STR.index(gst)
        ov = np.minimum(e, ge) - np.maximum(s, gs) + 1
        ok = m & (ov > 0) & (ov / (ge - gs + 1) >= 0.5) & (ov / (e - s + 1) >= 0.5) & (tid == -1)
        tid[ok] = j
    lab[tid >= 0] = 2
    for k, st in enumerate(STR):
        m = (st_arr == k) & (lab == 0)
        lab[m] = np.where(cov_bp(ocv[st], s[m], e[m]) > 0, 1, 0)
    return lab, tid


def build_pairs(c, pc, truth=None, obs=None):
    """All candidate pairs of one chromosome with labels and features -> dict of arrays (cached in outputs/AUTH/pairs)."""
    guard(c)
    path = f"{OUT}/pairs/{c}.npz"
    if os.path.exists(path):
        z = np.load(path)
        return {k: z[k] for k in z.files}
    ev, rows, raw = evidence(c)
    seq = load_seq(c)
    up = seq.upper()
    low = np.frombuffer(seq, np.uint8) >= 97
    rcum = np.concatenate([[0], np.cumsum(low, dtype=np.int64)])
    del low
    qfrac = {(q["pos"], q["strand"]): q["frac"] for q in raw["quant"]}
    qsig = {(p_, st): v for p_, st, v in raw["signal"]}
    qall = {st: np.array(sorted(q["pos"] for q in raw["quant"] if q["strand"] == st), np.int64) for st in STR}
    craw = {st: [(cl["tss"], sum(cl["counts"])) for cl in raw["cage"] if cl["strand"] == st] for st in STR}
    craw = {st: (np.array([t for t, _ in sorted(v)], np.int64), np.concatenate([[0], np.cumsum([n for _, n in sorted(v)])]))
            for st, v in craw.items()}
    cols = defaultdict(list)
    for k, st in enumerate(STR):
        T, P = ev["tss"][st], ev["pas"][st]
        if not len(T) or not len(P):
            continue
        Tp = T[:, 0]
        opp = craw["-" if st == "+" else "+"]
        for i, (t, cnt) in enumerate(T):
            t, cnt = int(t), int(cnt)
            if st == "+":
                j = np.searchsorted(P, t, side="right")
                opts = [int(x) for x in P[j:j + K_WIDE] if x - t <= DMAX]
            else:
                j = np.searchsorted(P, t, side="left")
                opts = [int(x) for x in P[max(0, j - K_WIDE):j][::-1] if t - x <= DMAX]
            cl = raw["cage"][rows[st][i][2]]
            # TSS-level
            a0, a1 = np.searchsorted(opp[0], t - WIN_ANTI, "left"), np.searchsorted(opp[0], t + WIN_ANTI, "right")
            anti = int(opp[1][a1] - opp[1][a0])
            w = up[max(t - 1 - WIN_GC, 0):t - 1 + WIN_GC + 1].decode()
            nn = len(w) - w.count("N")
            gc = (w.count("G") + w.count("C")) / nn if nn else 0.0
            cgi = float(gc >= 0.5 and S3.cpg_oe(w) >= 0.6)
            tssf = (math.log1p(cnt), float(sum(x > 0 for x in cl["counts"])), math.log1p(cl["end"] - cl["start"] + 1),
                    min(cl["mode_n"] / cnt, 1.0), math.log1p(anti), cgi, gc)
            for r, x in enumerate(opts):
                s, e = (t, x) if st == "+" else (x, t)
                if e - s + 1 < MIN_LEN:
                    continue
                if st == "+":
                    wd = up[x:x + 20].decode()
                    a_dn = wd.count("A") / max(len(wd), 1)
                    wu = up[max(x - 36, 0):x - 10].decode()
                    pas = float(any(h in wu for h in PAS_HEX))
                    n_cin = np.searchsorted(Tp, x, "right") - np.searchsorted(Tp, t, "right")
                    n_qin = np.searchsorted(qall[st], x, "left") - np.searchsorted(qall[st], t, "right")
                else:
                    wd = up[max(x - 21, 0):x - 1].decode()
                    a_dn = wd.count("T") / max(len(wd), 1)
                    wu = up[x + 9:x + 35].decode()
                    pas = float(any(h in wu for h in PAS_RC))
                    n_cin = np.searchsorted(Tp, t, "left") - np.searchsorted(Tp, x, "left")
                    n_qin = np.searchsorted(qall[st], t, "left") - np.searchsorted(qall[st], x, "right")
                cols["strand"].append(k)
                cols["row"].append(i)
                cols["tss"].append(t)
                cols["cnt"].append(cnt)
                cols["anti"].append(anti)
                cols["s"].append(s)
                cols["e"].append(e)
                cols["rank"].append(r)
                cols["endpos"].append(x)
                cols["X"].append(tssf + (float(round(qfrac[(x, st)] * 3)), math.log1p(qsig.get((x, st), 0.0)), a_dn, pas,
                                         float(r), math.log(e - s + 1), math.log1p(int(n_cin)), math.log1p(int(n_qin))))
    P = {k: np.array(v, np.float64 if k == "X" else np.int64) for k, v in cols.items()}
    P["X"] = P["X"].reshape(-1, len(FEATURES) - 3)
    s, e = P["s"], P["e"]
    ln = (e - s + 1).astype(np.float64)
    pcs = np.zeros(len(s))
    pca = np.zeros(len(s))
    for k, st in enumerate(STR):
        m = P["strand"] == k
        pcs[m] = cov_bp(cover(pc[st]), s[m], e[m]) / ln[m]
        pca[m] = cov_bp(cover(pc["-" if st == "+" else "+"]), s[m], e[m]) / ln[m]
    rep = (rcum[e] - rcum[s - 1]) / ln
    P["X"] = np.column_stack([P["X"], pcs, pca, rep])
    assert P["X"].shape[1] == len(FEATURES)
    truth = truth or load_truth(c)
    obs = obs or observed([c])[c]
    P["label"], P["tidx"] = label_pairs(P["strand"], s, e, truth["ts"], obs_cover(obs))
    os.makedirs(f"{OUT}/pairs", exist_ok=True)
    np.savez(path + ".tmp.npz", **P)
    os.replace(path + ".tmp.npz", path)
    return P


# ------------------------------------------------------------------ selection
def groups(P):
    """Start / end index of each TSS's pairs (pairs are stored TSS by TSS, nearest end first)."""
    key = P["strand"] * 10 ** 9 + P["row"]
    cut = np.nonzero(np.diff(key))[0] + 1
    return np.concatenate([[0], cut]), np.concatenate([cut, [len(key)]])


def choose(P, rule, score=None, ok_pair=None, ok_tss=None):
    """One pair index per TSS (or none).  rule: 'farthest3' (farthest of the 3 nearest = ev_farthest), 'nearest1' (the
    nearest site, W24 'pair'), 'score' (highest score among all 20; first on ties)."""
    g0, g1 = groups(P)
    out = []
    for a, b in zip(g0, g1):
        if ok_tss is not None and not ok_tss[a]:
            continue
        idx = np.arange(a, b)
        if ok_pair is not None:
            idx = idx[ok_pair[a:b]]
        if rule == "farthest3":
            idx = idx[P["rank"][idx] < K_R0]
            if len(idx):
                out.append(int(idx[np.argmax(P["rank"][idx])]))
        elif rule == "nearest1":
            idx = idx[P["rank"][idx] == 0]
            if len(idx):
                out.append(int(idx[0]))
        else:
            if len(idx):
                out.append(int(idx[np.argmax(score[idx])]))
    return np.array(out, np.int64)


def greedy(P, idx, key, masked=None):
    """e3_anchor._greedy: by key (descending, stable), accept units with no same-strand overlap with accepted ones; a
    masked unit (coding overlap) is dropped and does not block.  Returns accepted pair indices in acceptance order."""
    order = sorted(range(len(idx)), key=lambda u: -key[u])
    acc = {0: ([], []), 1: ([], [])}
    out = []
    for u in order:
        i = int(idx[u])
        A, Bn = acc[int(P["strand"][i])]
        s, e = int(P["s"][i]), int(P["e"][i])
        j = bisect.bisect_right(A, e) - 1
        if j >= 0 and Bn[j] >= s:
            continue
        if masked is not None and masked[i]:
            continue
        A.insert(j + 1, s)
        Bn.insert(j + 1, e)
        out.append(i)
    return np.array(out, np.int64)


def top_q(acc, q):
    return acc[:int(math.ceil(q * len(acc) - 1e-9))] if q < 1 else acc


def masks(P, tib=None):
    """Same-strand coding overlap per pair: v41 (feature pc_same > 0) or a Tiberius gene-span array."""
    if tib is None:
        return P["X"][:, FI["pc_same"]] > 0
    m = np.zeros(len(P["s"]), bool)
    for k, st in enumerate(STR):
        sel = P["strand"] == k
        m[sel] = cov_bp(cover(tib[st]), P["s"][sel], P["e"][sel]) > 0
    return m


def end_ok(P):
    """R3 end quality: downstream 20 nt A content <= 0.6 and a PAS hexamer 10-35 nt upstream."""
    X = P["X"]
    return (X[:, FI["a_down20"]] <= A_MAX) & (X[:, FI["pas_hex"]] > 0)


def tss_ok_bidir(P):
    """R4: keep TSSs whose antisense ssCAGE tags within 500 bp are < half of their own tags."""
    return P["anti"] < BIDIR * P["cnt"]


def method_units(P, name, score=None, q=1.0, masked=None):
    """Accepted unit indices (acceptance order) of one method on one chromosome.
    R0 ev_farthest | R1 nearest | R2 R0 + top q by tags | R3 end quality, then farthest | R4 no bidirectional TSS |
    R5 CpG-island TSS | R6 R2 + R3 + R4 | model (top q by score) | P1 (end by score, greedy by tags, no q)."""
    masked = masks(P) if masked is None else masked
    cnt = P["cnt"].astype(np.float64)
    if name in ("R0", "R2"):
        idx = choose(P, "farthest3")
    elif name == "R1":
        idx = choose(P, "nearest1")
    elif name == "R3":
        idx = choose(P, "farthest3", ok_pair=end_ok(P))
    elif name == "R4":
        idx = choose(P, "farthest3", ok_tss=tss_ok_bidir(P))
    elif name == "R5":
        idx = choose(P, "farthest3", ok_tss=P["X"][:, FI["cgi"]] > 0)
    elif name == "R6":
        idx = choose(P, "farthest3", ok_pair=end_ok(P), ok_tss=tss_ok_bidir(P))
    elif name in ("model", "P1"):
        idx = choose(P, "score", score=score)
    else:
        raise ValueError(name)
    key = score[idx] if name == "model" else cnt[idx]
    acc = greedy(P, idx, key, masked)
    return top_q(acc, q) if name in ("R2", "R6", "model") else acc


# ------------------------------------------------------------------ scoring
def s50_adj(ts, ss, ee, sts):
    """Adjacency of simple_det.s50 (reciprocal span overlap >= 0.5, same strand) for predictions given as arrays."""
    starts = [t[2] for t in ts]
    maxlen = max(t[3] - t[2] for t in ts) if ts else 0
    adj = []
    for s, e, st in zip(ss, ee, sts):
        lo, hi = np.searchsorted(starts, s - maxlen - 1), np.searchsorted(starts, e, side="right")
        a = []
        for i in range(lo, hi):
            _, gst, gs, ge = ts[i]
            if gst != st:
                continue
            ov = min(e, ge) - max(s, gs) + 1
            if ov > 0 and ov / (ge - gs + 1) >= 0.5 and ov / (e - s + 1) >= 0.5:
                a.append(i)
        adj.append(a)
    return adj


def score_spans(ts, ocv, ss, ee, sts, adj=None):
    """Per-prediction labels (2 TP / 1 unknown / 0 wrong), matched truth flags, counts."""
    ss, ee = np.asarray(ss, np.int64), np.asarray(ee, np.int64)
    n = len(ss)
    adj = s50_adj(ts, ss, ee, sts) if adj is None else adj
    ml = B.matching(adj, len(ts)) if n and ts else [-1] * n
    ov = np.zeros(n, bool)
    sa = np.array([STR.index(x) for x in sts], np.int64) if n else np.zeros(0, np.int64)
    for k, st in enumerate(STR):
        m = sa == k
        if m.any():
            ov[m] = cov_bp(ocv[st], ss[m], ee[m]) > 0
    tp = np.array([x != -1 for x in ml], bool)
    lab = np.where(tp, 2, np.where(ov, 1, 0)).astype(np.int8)
    matched = np.zeros(len(ts), bool)
    for x in ml:
        if x != -1:
            matched[x] = True
    return {"lab": lab, "matched": matched, "n_truth": len(ts), "n_pred": n, "tp": int(tp.sum()),
            "unk": int((lab == 1).sum()), "wrong": int((lab == 0).sum()), "tp_not_observed": int((tp & ~ov).sum())}


def score_units(P, acc, ts, ocv):
    sts = [STR[int(x)] for x in P["strand"][acc]]
    return score_spans(ts, ocv, P["s"][acc], P["e"][acc], sts)


def metrics(rs):
    """Pooled metrics over a list of score_spans results."""
    nt, n, tp, unk = (sum(r[k] for r in rs) for k in ("n_truth", "n_pred", "tp", "unk"))
    return {"recall": tp / nt if nt else None, "prec_known": tp / (n - unk) if n - unk > 0 else None,
            "prec_strict": tp / n if n else None, "n_pred": n, "tp": tp, "unknown": unk, "wrong": n - unk - tp, "n_truth": nt}


def q_curve(P, acc, ts, ocv, qs):
    """Counts for every q in qs on one chromosome (the accepted list is ranked; prefixes are re-matched one-to-one)."""
    sts = [STR[int(x)] for x in P["strand"][acc]]
    ss, ee = P["s"][acc], P["e"][acc]
    adj = s50_adj(ts, ss, ee, sts)
    out = []
    for q in qs:
        k = int(math.ceil(q * len(acc) - 1e-9)) if q < 1 else len(acc)
        out.append(score_spans(ts, ocv, ss[:k], ee[:k], sts[:k], adj[:k]))
    return out


# ------------------------------------------------------------------ model
def fit(X, y, C=LOGIT_C):
    """Standardise on the fitting rows, then simple_det.fit_logit (L2, L-BFGS)."""
    mu, sd = X.mean(0), X.std(0)
    sd = np.where(sd > 0, sd, 1.0)
    w = D.fit_logit((X - mu) / sd, y.astype(np.float64), C)
    return {"mu": mu.tolist(), "sd": sd.tolist(), "w": w.tolist()}


def logit(model, X):
    w = np.asarray(model["w"])
    return ((X - np.asarray(model["mu"])) / np.asarray(model["sd"])) @ w[1:] + w[0]


def sigmoid(z):
    return 1 / (1 + np.exp(-np.clip(z, -50, 50)))


def auc(score, y):
    """Mann-Whitney AUC (ties averaged)."""
    from scipy.stats import rankdata
    r = rankdata(score)
    n1 = int(y.sum())
    n0 = len(y) - n1
    return float((r[y.astype(bool)].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)) if n1 and n0 else None


def calibration(prob, y, bins=10):
    out = []
    for b in range(bins):
        lo, hi = b / bins, (b + 1) / bins
        m = (prob >= lo) & ((prob < hi) if b < bins - 1 else (prob <= hi))
        out.append({"bin": [lo, hi], "n": int(m.sum()), "mean_prob": float(prob[m].mean()) if m.any() else None,
                    "frac_true": float(y[m].mean()) if m.any() else None})
    return out


# ------------------------------------------------------------------ block bootstrap
def block_ids(chroms, lengths):
    blocks = [(c, b) for c in chroms for b in range((lengths[c] - 1) // BLOCK + 1)]
    return {blk: i for i, blk in enumerate(blocks)}


def block_table(res_by_chrom, ts_by_chrom, spans_by_chrom, bid):
    """Per-block counts (n_truth, matched, n_pred, tp, unknown) of one method; blocks by 1-based start // 1 Mb."""
    T = np.zeros((len(bid), 5))
    for c, r in res_by_chrom.items():
        for j, t in enumerate(ts_by_chrom[c]):
            b = bid[(c, (t[2] - 1) // BLOCK)]
            T[b, 0] += 1
            T[b, 1] += r["matched"][j]
        for s, lab in zip(spans_by_chrom[c], r["lab"]):
            b = bid[(c, (int(s) - 1) // BLOCK)]
            T[b, 2] += 1
            T[b, 3] += lab == 2
            T[b, 4] += lab == 1
    return T


def boot_weights(nb, seed=0, n=N_BOOT):
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, nb, (n, nb))
    W = np.zeros((n, nb))
    for i in range(n):
        W[i] = np.bincount(idx[i], minlength=nb)
    return W


def boot_metrics(W, T):
    s = W @ T
    with np.errstate(invalid="ignore", divide="ignore"):
        return {"recall": s[:, 1] / s[:, 0], "prec_known": s[:, 3] / (s[:, 2] - s[:, 4]), "prec_strict": s[:, 3] / s[:, 2]}


def point(T):
    s = T.sum(0)
    return {"recall": s[1] / s[0] if s[0] else None, "prec_known": s[3] / (s[2] - s[4]) if s[2] - s[4] else None,
            "prec_strict": s[3] / s[2] if s[2] else None, "n_pred": float(s[2])}


def paired(W, T1, T2):
    """Point difference (method 1 - method 2) and 95 % percentile CI from the same block resamples."""
    b1, b2 = boot_metrics(W, T1), boot_metrics(W, T2)
    p1, p2 = point(T1), point(T2)
    out = {}
    for k in ("recall", "prec_known", "prec_strict"):
        d = b1[k] - b2[k]
        d = d[np.isfinite(d)]
        out[k] = {"diff": (p1[k] or 0) - (p2[k] or 0), "ci95": np.percentile(d, [2.5, 97.5]).tolist() if len(d) else None}
    return out


def sha256(path):
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()
