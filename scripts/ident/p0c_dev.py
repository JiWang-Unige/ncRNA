#!/usr/bin/env python3
"""IDENT P0c development checks (docs/ident/P0c/DEV_CHECKS_PLAN.md, criteria written before running).

Dev chromosomes chr16, chr22.  rep1 (ENCFF373TKM) split by read-name hash into D 20% / A 80%; V = ENCFF388HXU.
Intron-chain unit: graph nodes are admitted introns (>= 2 QC'd D reads, canonical motif), ordered 5'->3' on the
transcript strand; each D read's intron sequence is a hard subpath observation; open boundaries.

  p0c_dev.py extract     -> per-chromosome read tables (QC'd) for rep1 and V
  p0c_dev.py build       -> loci, G-strict / G-expanded graphs, states (check 2), quant states (check 3)
  p0c_dev.py qo          -> q_o models and A-sampling outcomes (check 1)
  p0c_dev.py report      -> DEV_CHECKS.json
"""
import collections
import hashlib
import json
import os
import pickle
import random
import re
import subprocess
import sys
import time
from multiprocessing import get_context

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ident_core as C        # noqa: E402

ROOT = os.path.expanduser("~/ncrna-annot")
ENC = f"{ROOT}/data/wtc11/encode"
FASTA = f"{ROOT}/data/x_audit/genomes/hg38.fa"
SAMTOOLS = os.path.expanduser("~/envs/lrtools/bin/samtools")
OUT = f"{ROOT}/outputs/IDENT/P0c_dev"
CHROMS = ("chr16", "chr22")
BAMS = {"rep1": "ENCFF373TKM.bam", "V": "ENCFF388HXU.bam"}
D_FRAC = 0.20
MIN_INTRON, MAX_CLIP, MIN_MAPQ = 50, 30, 30
MIN_D_READS, LOCUS_BUDGET = 20, 300.0
CANON = {("GT", "AG"), ("GC", "AG"), ("AT", "AC")}
RC = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(RC)[::-1]


# ------------------------------------------------------------------ extract
def load_chrom(ch):
    """Chromosome extracted from data/human/hg38.fa.gz into OUT/genome/<ch>.fa (one-off, by zcat | awk)."""
    with open(f"{OUT}/genome/{ch}.fa") as fh:
        return "".join(l.strip() for l in fh if not l.startswith(">")).upper()


def parse_reads(bam, ch, genome):
    """-> list of (qname, strand, introns tuple in transcript order, blocks, annotated_flags) for QC'd reads."""
    p = subprocess.Popen([SAMTOOLS, "view", "-q", str(MIN_MAPQ), "-F", "0x904", bam, ch], stdout=subprocess.PIPE, text=True)
    out, stats = [], collections.Counter()
    for line in p.stdout:
        f = line.rstrip("\n").split("\t")
        stats["seen"] += 1
        flag, pos, cig = int(f[1]), int(f[3]), f[5]
        tags = dict((t[:2], t[5:]) for t in f[11:])
        ops = re.findall(r"(\d+)([MIDNSHX=])", cig)
        if (ops[0][1] == "S" and int(ops[0][0]) > MAX_CLIP) or (ops[-1][1] == "S" and int(ops[-1][0]) > MAX_CLIP):
            stats["clip"] += 1
            continue
        g, introns, blocks, bstart = pos - 1, [], [], pos - 1
        short = False
        for n, op in ops:
            n = int(n)
            if op in "M=XD":
                g += n
            elif op == "N":
                if n < MIN_INTRON:
                    short = True
                blocks.append((bstart, g)); introns.append((g, g + n)); g += n; bstart = g
        blocks.append((bstart, g))
        if short:
            stats["short_intron"] += 1
            continue
        rstrand = "-" if flag & 16 else "+"
        ts = tags.get("ts")
        strand = rstrand if ts in (None, "+") else ("-" if rstrand == "+" else "+")
        bad = False
        for a, b in introns:
            don, acc = genome[a:a + 2], genome[b - 2:b]
            mot = (don, acc) if strand == "+" else (rc(acc), rc(don))
            if mot not in CANON:
                bad = True
                break
        if bad:
            stats["noncanonical"] += 1
            continue
        jm = tags.get("jM", "")
        ann = [int(v) >= 20 for v in jm.split(",")[1:]] if jm.startswith("c,") or jm.startswith("B:c") else []
        if strand == "-":
            introns = introns[::-1]
        out.append((f[0], strand, tuple(introns), tuple(blocks), ann))
        stats["kept"] += 1
    p.wait()
    return out, stats


def is_D(qname):
    return int(hashlib.md5(qname.encode()).hexdigest()[:8], 16) / 0xFFFFFFFF < D_FRAC


def cmd_extract():
    os.makedirs(OUT, exist_ok=True)
    meta = {}
    for ch in CHROMS:
        genome = load_chrom(ch)
        for lab, bam in BAMS.items():
            reads, stats = parse_reads(f"{ENC}/{bam}", ch, genome)
            ann = [a for r in reads for a in r[4]]
            meta[f"{ch}_{lab}"] = dict(stats, annotated_junction_frac=(sum(ann) / len(ann)) if ann else None)
            pickle.dump(reads, open(f"{OUT}/{ch}_{lab}.pkl", "wb"))
            print(ch, lab, dict(stats), flush=True)
    names = {lab: set() for lab in BAMS}
    for ch in CHROMS:
        for lab in BAMS:
            names[lab] |= {r[0].rsplit("/", 1)[0] for r in pickle.load(open(f"{OUT}/{ch}_{lab}.pkl", "rb"))}
    meta["movie_zmw_overlap_rep1_V"] = len(names["rep1"] & names["V"])
    json.dump(meta, open(f"{OUT}/extract_meta.json", "w"), indent=1)
    print(json.dumps(meta, indent=1))


# ------------------------------------------------------------------ loci and graphs
class UF:
    def __init__(self):
        self.p = {}

    def find(self, x):
        self.p.setdefault(x, x)
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        self.p[self.find(a)] = self.find(b)


def build_loci(ch):
    rep1 = pickle.load(open(f"{OUT}/{ch}_rep1.pkl", "rb"))
    D = [r for r in rep1 if is_D(r[0])]
    A = [r for r in rep1 if not is_D(r[0])]
    sup = collections.Counter((r[1], i) for r in D for i in r[2])
    admitted = {k for k, v in sup.items() if v >= 2}
    Dq = [r for r in D if r[2] and all((r[1], i) in admitted for i in r[2])]
    dropped = sum(1 for r in D if r[2] and not all((r[1], i) in admitted for i in r[2]))
    uf = UF()
    for strand, iv in admitted:                     # union overlapping introns (sweep)
        pass
    by_strand = collections.defaultdict(list)
    for strand, iv in admitted:
        by_strand[strand].append(iv)
    for strand, ivs in by_strand.items():
        ivs.sort()
        cur_end, cur = -1, None
        for a, b in ivs:
            if a < cur_end:
                uf.union((strand, (a, b)), cur)
            else:
                uf.find((strand, (a, b)))
            if b > cur_end:
                cur_end, cur = b, (strand, (a, b))
    for r in Dq:                                    # reads link their introns
        for i in r[2][1:]:
            uf.union((r[1], r[2][0]), (r[1], i))
    loci = collections.defaultdict(lambda: {"D": [], "A": [], "introns": set()})
    for k in admitted:
        loci[uf.find(k)]["introns"].add(k[1])
    for r in Dq:
        loci[uf.find((r[1], r[2][0]))]["D"].append(r)
    intron_locus = {k[1]: (k[0], uf.find(k)) for k in admitted}
    for r in A:
        if r[2]:
            hit = {intron_locus.get(i) for i in r[2]} - {None}
            hit = {h for h in hit if h[0] == r[1]}
            if len(hit) == 1:
                loci[next(iter(hit))[1]]["A"].append(r)
    out = []
    for key, L in loci.items():
        L["strand"] = key[0]
        L["key"] = f"{ch}:{key[0]}:{min(a for a, _ in L['introns'])}"
        out.append(L)
    return out, {"D": len(D), "D_qc": len(Dq), "D_dropped_unadmitted": dropped, "A": len(A), "admitted_introns": len(admitted)}


def node_order(introns, strand):
    return sorted(introns, key=lambda iv: iv[0], reverse=(strand == "-"))


def make_graph(L, expanded=False):
    """-> (Locus, idx map intron->node, hard observations)."""
    strand = L["strand"]
    order = node_order(L["introns"], strand)
    idx = {iv: k for k, iv in enumerate(order)}
    edges, starts, ends, hard = set(), collections.Counter(), collections.Counter(), set()
    for r in L["D"]:
        ch = tuple(idx[i] for i in r[2])
        hard.add(("sub", ch))
        starts[ch[0]] += 1
        ends[ch[-1]] += 1
        for a, b in zip(ch, ch[1:]):
            edges.add((a, b))
    if expanded:
        blocks = [b for r in L["D"] for b in r[3]]
        for i in order:
            for j in order:
                if idx[i] >= idx[j]:
                    continue
                lo, hi = (i[1], j[0]) if strand == "+" else (j[1], i[0])
                if lo >= hi:
                    continue
                if any(bs <= lo and be >= hi for bs, be in blocks):
                    edges.add((idx[i], idx[j]))
    sources = {f"s{v}": (v, False) for v in starts}
    sinks = {f"t{v}": (v, False) for v in ends}
    loc = C.Locus(len(order), edges, sources, sinks)
    return loc, idx, hard


def chain_of(path):
    return path[1]


def nested_flags(paths):
    chains = [chain_of(p) for p in paths]
    out = []
    for k, c in enumerate(chains):
        L = len(c)
        out.append(any(len(d) > L and any(d[i:i + L] == c for i in range(len(d) - L + 1)) for d in chains))
    return out


def states_inf_fast(ev):
    """delta = infinity, structure only: NECESSARY iff some hard observation is produced only by this path; else
    POSSIBLE (proposition in ROUND8).  Paths producing nothing are still POSSIBLE (open boundaries)."""
    P = ev.locus.paths
    st = ["POSSIBLE"] * len(P)
    for o in ev.coverage_obs():
        prod = [k for k, p in enumerate(P) if C.produces(p, o)]
        if not prod:
            return None
        if len(prod) == 1:
            st[prod[0]] = "NECESSARY"
    return st


def analyse_locus(args):
    L, ch = args
    t0 = time.time()
    res = {"key": L["key"], "nD": len(L["D"]), "nA": len(L["A"]), "graphs": {}}
    for gname, expanded in (("strict", False), ("expanded", True)):
        loc, idx, hard = make_graph(L, expanded)
        g = {"n_paths": len(loc.paths), "too_many": loc.too_many}
        if loc.too_many or len(loc.paths) < 2:
            res["graphs"][gname] = g
            continue
        ev = C.Evidence(loc, hard=hard)
        g["chains"] = [chain_of(p) for p in loc.paths]
        g["order"] = node_order(L["introns"], L["strand"])
        g["nested"] = nested_flags(loc.paths)
        g["inf"] = states_inf_fast(ev)
        deadline = t0 + LOCUS_BUDGET
        t1 = C.three_state(ev, 1, deadline=deadline, time_limit=20.0)
        g["d1_model"], g["kstar"] = t1["model"], t1.get("kstar")
        g["d1"] = [t1["paths"][k]["state"] for k in range(len(loc.paths))] if t1["model"] == "OK" else None
        if gname == "strict":
            # check 3: same-library junction counts as quantitative evidence, delta = infinity
            y = collections.Counter()
            for r in L["D"]:
                for i in r[2]:
                    y[("sub", (idx[i],))] += 1
            evq = C.Evidence(loc, hard=hard, quant=[("short2", {o: float(v) for o, v in y.items()}, 0.30, 2.0)])
            if time.time() < deadline:
                tq = C.three_state(evq, None, deadline=deadline, time_limit=20.0)
                g["q_model"] = tq["model"]
                g["qinf"] = [tq["paths"][k]["state"] for k in range(len(loc.paths))] if tq["model"] == "OK" else None
            else:
                g["q_model"] = "UNRESOLVED"
            g["ev"] = ev
        res["graphs"][gname] = g
    res["seconds"] = round(time.time() - t0, 1)
    res["timed_out"] = res["seconds"] > LOCUS_BUDGET
    return res


def cmd_build(procs=10):
    allres, meta = [], {}
    for ch in CHROMS:
        loci, m = build_loci(ch)
        meta[ch] = dict(m, loci=len(loci))
        elig = [L for L in loci if len(L["D"]) >= MIN_D_READS]
        meta[ch]["loci_ge20D"] = len(elig)
        with get_context("spawn").Pool(procs) as pool:
            for r in pool.imap_unordered(analyse_locus, [(L, ch) for L in elig], chunksize=1):
                allres.append(r)
        pickle.dump(elig, open(f"{OUT}/{ch}_loci.pkl", "wb"))
        print(ch, meta[ch], flush=True)
    pickle.dump(allres, open(f"{OUT}/build.pkl", "wb"))
    json.dump(meta, open(f"{OUT}/build_meta.json", "w"), indent=1)


# ------------------------------------------------------------------ check 1: q_o
def trunc_dist(loci):
    """P(k introns covered | candidate has n introns): empirical read intron counts, truncated at n."""
    ks = collections.Counter(len(r[2]) for L in loci for r in L["D"])
    tot = sum(ks.values())
    return {k: v / tot for k, v in ks.items()}


def p_obs_given_path(o, chain, pk):
    n, L = len(chain), len(o)
    occ = [i for i in range(n - L + 1) if chain[i:i + L] == o]
    if not occ:
        return 0.0
    prob = 0.0
    tail = sum(v for k, v in pk.items() if k >= n)
    for k in range(1, n + 1):
        w = pk.get(k, 0.0) if k < n else tail
        if w == 0 or k < L:
            continue
        nwin = n - k + 1
        good = sum(1 for s in range(nwin) if any(s <= i and i + L <= s + k for i in occ))
        prob += w * good / nwin
    return prob


def mixture(chains, reads_chains, pk, model):
    m = len(chains)
    if model == "Q0":
        return np.full(m, 1.0 / m)
    lik = np.array([[p_obs_given_path_exact(rc_, c, pk) for c in chains] for rc_ in reads_chains])
    if model == "Q1":
        w = np.zeros(m)
        for row in lik:
            nz = row > 0
            if nz.any():
                w[nz] += 1.0 / nz.sum()
        return (w + 1e-3) / (w + 1e-3).sum()
    pi = np.full(m, 1.0 / m)                       # Q2: EM under the truncation model
    for _ in range(100):
        r = lik * pi
        s = r.sum(1, keepdims=True)
        s[s == 0] = 1
        pi = (r / s).sum(0) + 1e-3
        pi /= pi.sum()
    return pi


def p_obs_given_path_exact(read_chain, chain, pk):
    """Probability that a read from `chain` has exactly intron chain `read_chain` (window equality)."""
    n, L = len(chain), len(read_chain)
    if L > n:
        return 0.0
    tail = sum(v for k, v in pk.items() if k >= n)
    w = pk.get(L, 0.0) if L < n else tail
    nwin = n - L + 1
    return w * sum(1 for s in range(nwin) if chain[s:s + L] == read_chain) / nwin


def qo_locus(args):
    res, L, pk, reps, seed = args
    g = res["graphs"].get("strict", {})
    if not g.get("d1") or g.get("d1_model") != "OK":
        return []
    ev, order = g["ev"], g["order"]
    idx = {iv: k for k, iv in enumerate(order)}
    chains = g["chains"]
    reads_chains = [tuple(idx[i] for i in r[2]) for r in L["D"]]
    pis = {mname: mixture(chains, reads_chains, pk, mname) for mname in ("Q0", "Q1", "Q2")}
    A_chains = []
    for r in L["A"]:
        if all(i in idx for i in r[2]):
            A_chains.append(tuple(idx[i] for i in r[2]))
        else:
            A_chains.append(None)                  # out-of-graph read: consumes budget, produces nothing in-graph
    if not A_chains:
        return []
    b = max(1, round(0.10 * len(L["D"])))
    rng = random.Random(seed)
    out = []
    for k, st in enumerate(g["d1"]):
        if st != "POSSIBLE" or g["nested"][k]:
            continue
        dg = C.diagnose(ev, k, g["kstar"] + 1, protocols=("long",), time_limit=20.0)
        if dg["status"] != "SUFFICIENT":
            out.append({"key": res["key"], "path": k, "status": dg["status"]})
            continue
        obs = [o[1] for o in dg["obs"]]
        qo = {mn: [float(sum(pi[j] * p_obs_given_path(o, chains[j], pk) for j in range(len(chains)))) for o in obs] for mn, pi in pis.items()}
        for rep in range(reps):
            samp = [A_chains[i] for i in rng.sample(range(len(A_chains)), min(b, len(A_chains)))]
            hit = [any(c is not None and any(c[i:i + len(o)] == o for i in range(len(c) - len(o) + 1)) for c in samp) for o in obs]
            out.append({"key": res["key"], "path": k, "status": "SUFFICIENT", "rep": rep, "b": len(samp), "n_obs": len(obs),
                        "hit": hit, "qo": qo, "all_hit": all(hit)})
    return out


def cmd_qo(procs=12, reps=5):
    allres = {r["key"]: r for r in pickle.load(open(f"{OUT}/build.pkl", "rb"))}
    jobs = []
    for ch in CHROMS:
        loci = pickle.load(open(f"{OUT}/{ch}_loci.pkl", "rb"))
        pk = trunc_dist(loci)
        for k, L in enumerate(loci):
            if L["key"] in allres:
                jobs.append((allres[L["key"]], L, pk, reps, k))
    rows = []
    with get_context("spawn").Pool(procs) as pool:
        for r in pool.imap_unordered(qo_locus, jobs, chunksize=1):
            rows.extend(r)
    pickle.dump(rows, open(f"{OUT}/qo.pkl", "wb"))
    print("qo rows", len(rows))


# ------------------------------------------------------------------ report
def brier_auc(rows, model):
    y, p, keys = [], [], []
    for r in rows:
        if r.get("status") != "SUFFICIENT":
            continue
        for h, q in zip(r["hit"], r["qo"][model]):
            y.append(1.0 if h else 0.0); p.append(1 - (1 - q) ** r["b"]); keys.append(r["key"])
    return np.array(y), np.array(p), np.array(keys)


def auc(y, p):
    pos, neg = p[y == 1], p[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return None
    order = np.argsort(np.r_[pos, neg])
    ranks = np.empty(len(order)); ranks[order] = np.arange(1, len(order) + 1)
    return float((ranks[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def cmd_report():
    build = pickle.load(open(f"{OUT}/build.pkl", "rb"))
    meta = json.load(open(f"{OUT}/build_meta.json"))
    ext = json.load(open(f"{OUT}/extract_meta.json"))
    rep = {"extract": ext, "build": meta}
    # V exact chains (genomic introns), per strand
    Vchains = collections.Counter()
    for ch in CHROMS:
        for r in pickle.load(open(f"{OUT}/{ch}_V.pkl", "rb")):
            Vchains[(r[1], r[2])] += 1
    # check 2
    c2 = {}
    for gname in ("strict", "expanded"):
        elig = [r for r in build if r["graphs"].get(gname, {}).get("n_paths", 0) >= 2 or r["graphs"].get(gname, {}).get("too_many")]
        done = [r for r in elig if r["graphs"][gname].get("d1") is not None]
        nonnested_possible_loci = [r for r in done if any(s == "POSSIBLE" and not n for s, n in zip(r["graphs"][gname]["d1"], r["graphs"][gname]["nested"]))]
        tot_paths = sum(len(r["graphs"][gname]["chains"]) for r in done)
        nested = sum(sum(r["graphs"][gname]["nested"]) for r in done)
        poss_inf = sum(sum(1 for s in r["graphs"][gname]["inf"] if s == "POSSIBLE") for r in done)
        poss_d1 = sum(sum(1 for s in r["graphs"][gname]["d1"] if s == "POSSIBLE") for r in done)
        c2[gname] = {"eligible_loci": len(elig), "computed": len(done), "computed_frac": len(done) / max(1, len(elig)),
                     "too_many_paths": sum(1 for r in elig if r["graphs"][gname].get("too_many")),
                     "d1_model_counts": collections.Counter(r["graphs"][gname].get("d1_model", "NA") for r in elig),
                     "paths": tot_paths, "nested_frac": nested / max(1, tot_paths),
                     "possible_inf": poss_inf, "possible_d1": poss_d1,
                     "loci_with_nonnested_possible_d1": len(nonnested_possible_loci)}
    # strict NECESSARY at d1 that become POSSIBLE in expanded (match by genomic chain)
    flip, base = 0, 0
    for r in build:
        s, e = r["graphs"].get("strict", {}), r["graphs"].get("expanded", {})
        if not s.get("d1") or not e.get("d1"):
            continue
        es = {tuple(e["order"][i] for i in c): st for c, st in zip(e["chains"], e["d1"])}
        for c, st in zip(s["chains"], s["d1"]):
            if st == "NECESSARY":
                base += 1
                if es.get(tuple(s["order"][i] for i in c)) == "POSSIBLE":
                    flip += 1
    c2["strict_d1_necessary"] = base
    c2["strict_d1_necessary_possible_in_expanded"] = flip
    c2["PASS"] = c2["strict"]["loci_with_nonnested_possible_d1"] >= 30 and c2["strict"]["computed_frac"] >= 0.80
    rep["check2"] = c2
    # check 3
    groups = {"evidence_forced": [], "parsimony_forced": [], "quant_forced": []}
    q_inf, q_elig = 0, 0
    for r in build:
        g = r["graphs"].get("strict", {})
        if not g.get("inf"):
            continue
        q_elig += 1
        if g.get("q_model") == "INFEASIBLE_MODEL":
            q_inf += 1
        strand = r["key"].split(":")[1]
        for k, c in enumerate(g["chains"]):
            gen = tuple(g["order"][i] for i in c)
            vs = Vchains.get((strand, gen), 0) > 0
            if g["inf"][k] == "NECESSARY":
                groups["evidence_forced"].append(vs)
            elif g.get("d1") and g["d1"][k] == "NECESSARY":
                groups["parsimony_forced"].append(vs)
            if g.get("qinf") and g["qinf"][k] == "NECESSARY" and g["inf"][k] != "NECESSARY":
                groups["quant_forced"].append(vs)
    c3 = {k: {"n": len(v), "V_support_rate": (sum(v) / len(v)) if v else None} for k, v in groups.items()}
    c3["quant_infeasible_loci"] = q_inf
    c3["quant_infeasible_frac"] = q_inf / max(1, q_elig)
    c3["PASS"] = (c3["quant_forced"]["n"] >= 20 and c3["quant_forced"]["V_support_rate"] is not None
                  and c3["parsimony_forced"]["V_support_rate"] is not None
                  and c3["quant_forced"]["V_support_rate"] >= c3["parsimony_forced"]["V_support_rate"]
                  and c3["quant_infeasible_frac"] <= 0.20)
    rep["check3"] = c3
    # check 1
    if os.path.exists(f"{OUT}/qo.pkl"):
        rows = pickle.load(open(f"{OUT}/qo.pkl", "rb"))
        c1 = {"diagnosis_status": collections.Counter(r["status"] for r in rows if r.get("rep", 0) == 0)}
        y0, p0, keys = brier_auc(rows, "Q0")
        c1["pairs"] = int(len(y0)); c1["base_rate"] = float(y0.mean()) if len(y0) else None
        b0 = (y0 - p0) ** 2
        rng = np.random.default_rng(0)
        uk = np.unique(keys)
        for mn in ("Q0", "Q1", "Q2"):
            y, p, _ = brier_auc(rows, mn)
            bm = (y - p) ** 2
            skill = 1 - bm.mean() / b0.mean() if b0.mean() > 0 else None
            boots = []
            for _ in range(2000):
                pick = rng.choice(uk, len(uk))
                msk = np.concatenate([np.nonzero(keys == k)[0] for k in pick])
                if b0[msk].mean() > 0:
                    boots.append(1 - bm[msk].mean() / b0[msk].mean())
            c1[mn] = {"brier": float(bm.mean()), "skill_vs_Q0": skill, "skill_ci": [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))] if boots else None,
                      "auc": auc(y, p), "mean_pred": float(p.mean())}
        bins = np.linspace(0, 1, 6)
        for mn in ("Q1", "Q2"):
            y, p, _ = brier_auc(rows, mn)
            c1[mn]["calibration"] = [(float(bins[i]), float(bins[i + 1]), int(((p >= bins[i]) & (p < bins[i + 1] + (i == 4))).sum()),
                                      float(y[(p >= bins[i]) & (p < bins[i + 1] + (i == 4))].mean()) if ((p >= bins[i]) & (p < bins[i + 1] + (i == 4))).any() else None) for i in range(5)]
        c1["PASS"] = c1["pairs"] >= 50 and any(c1[mn]["skill_ci"] and c1[mn]["skill_ci"][0] > 0 and (c1[mn]["auc"] or 0) >= 0.70 for mn in ("Q1", "Q2"))
        rep["check1"] = c1
    json.dump(rep, open(f"{OUT}/DEV_CHECKS.json", "w"), indent=1, default=str)
    print(json.dumps(rep, indent=1, default=str))


if __name__ == "__main__":
    {"extract": cmd_extract, "build": cmd_build, "qo": cmd_qo, "report": cmd_report}[sys.argv[1]]()
