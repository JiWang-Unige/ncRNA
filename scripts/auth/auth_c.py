#!/usr/bin/env python3
"""AUTH plan (docs/auth/PLAN_AUTH.md section 2): descriptive checks C-1 and C-2 on the W24 test chromosomes (chr9 / chr12 /
chr20 / chr13), which M1 weeks 2-4 already read once.  No gate; nothing here changes the AUTH gates or its test split.

The W24 predictions are rebuilt with the frozen W24 code (w24_score / e3_anchor / dev_baselines; 100 % evidence, seed 0)
and the pooled counts are checked against outputs/M1/W24/RESULT.json before anything is classified.

C-1  What the 'uncertain' predictions are.  Each ev_farthest prediction (after the exclusion zone / half cuts) is
       TP    S50-matched (one-to-one, reciprocal span overlap >= 0.5) to an evaluable truth lnc locus;
       UNC   not TP, same-strand span overlap with a transcript observed in rep2 or rep3 (count >= 1, not 'Genomic');
       WRONG neither.
     UNC is split by the first rule that applies:
       U-a  overlaps an evaluable truth lnc locus (same strand): 'taken' (reciprocal >= 0.5 with a locus, but the
            one-to-one match went elsewhere), else by ends against the most-overlapped locus: 5' end > 1 kb from every
            truth transcript 5' end ('5p'), 3' end > 1 kb from every 3' end ('3p'), 'both', or 'span' (ends within 1 kb
            but S50 still fails);
       U-b  overlaps an observed transcript of a V29 protein_coding gene (a coding gene the Tiberius mask missed);
       U-c  overlaps an observed lnc-type transcript (V29 lnc gene type, or novel gene Antisense / Intergenic) outside
            the evaluable truth: reason 'both' (in rep2 and rep3, locus excluded / cut / < 200 bp), else 'single'
            (one replicate), 'ism' (ISM only) or 'mono' (novel single-exon);
       U-d  other V29 gene types only (pseudogenes, small-RNA hosts, ...).
     For U-a with a 5' end within 1 kb we also ask whether a better end existed among the candidate QuantSeq sites (the
     K = 3 used by ev_farthest, or the 20 nearest within 50 kb): the room for the 'choose the end' half of AUTH.
     The same labels are reported for the predictions the sequence filter removed (ev_farthest minus ours, per v2.0
     run, averaged), and for the v41 protein_coding mask that AUTH will use instead of Tiberius.
C-2  Truth loci found by ours (mean over runs >= 0.5) but not by rep1 IsoQuant (W24: 31): rep1 reads in the span (same
     strand, primary), rep1 / rep2 / rep3 TALON counts of the locus transcripts, IsoQuant at the site ('absent', 'masked'
     = removed by the Tiberius mask or the 200 bp floor, 'mismatch' = gene span fails reciprocal 0.5, 'taken'), whether
     one IsoQuant transcript (not the gene span W24 scored) would pass, StringTie, length, exons, novelty, CpG island at
     the 5' end (+-500 bp GC >= 0.5 and o/e >= 0.6), compared with the truth loci IsoQuant did find.
  auth_c.py -> outputs/AUTH/C/{RESULT_C.json, c1_predictions.tsv, c2_loci.tsv}
"""
import gzip
import json
import os
import re
import subprocess
import sys
import time
from collections import Counter, defaultdict

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, f"{ROOT}/scripts/m1/w24")
sys.path.insert(0, f"{ROOT}/scripts/phase1/s3")
import w24_lib as L                     # noqa: E402  (also puts src / phase1 / m1 on the path)
import w24_score as WS                  # noqa: E402
import b1_strict as B                   # noqa: E402
from e3_anchor import candidates, select_evidence, select_hybrid   # noqa: E402
from dev_baselines import tib_hit       # noqa: E402
from wtc11_truth import LNC_TYPES       # noqa: E402
import s3_features as S3                # noqa: E402

OUT = f"{ROOT}/outputs/AUTH/C"
ENC = L.ENC
SAM = os.path.expanduser("~/envs/lrtools/bin/samtools")
REPS = {"rep1": ("ENCFF803XYM.tsv", "ENCFF370SDY.gtf.gz"), "rep2": ("ENCFF527KWD.tsv", "ENCFF798DWV.gtf.gz"),
        "rep3": ("ENCFF486EKE.tsv", "ENCFF422CGK.gtf.gz")}
OFF = 1_000
DMAX, K, K_WIDE = WS.DMAX, WS.K, 20


# ------------------------------------------------------------------ inputs
def observed(chroms):
    """{rep: list of observed transcripts (count >= 1, not 'Genomic')} and {rep: total count} for CPM."""
    out, tot = {}, {}
    for r, (q, g) in REPS.items():
        want, n_all = {}, 0
        with open(f"{ENC}/{q}") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                f = dict(zip(hdr, line.rstrip("\n").split("\t")))
                n = int(f[hdr[-1]])
                n_all += n
                if n >= 1 and f["transcript_novelty"] != "Genomic":
                    want[f["annot_transcript_id"]] = (n, f["transcript_novelty"], f["gene_novelty"])
        tot[r] = n_all
        tx = {}
        with gzip.open(f"{ENC}/{g}", "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split("\t", 8)
                if f[0] not in chroms or f[2] not in ("transcript", "exon"):
                    continue
                tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
                if tid not in want:
                    continue
                t = tx.setdefault(tid, {"exons": []})
                if f[2] == "transcript":
                    gt = re.search(r'gene_type "([^"]+)"', f[8])
                    t.update(chrom=f[0], strand=f[6], s=int(f[3]), e=int(f[4]), gene_type=gt.group(1) if gt else None,
                             gene_id=re.search(r'gene_id "([^"]+)"', f[8]).group(1))
                else:
                    t["exons"].append((int(f[3]), int(f[4])))
        lst = []
        for tid, t in tx.items():
            if "strand" not in t:
                continue
            t["exons"].sort()
            n, tnov, gnov = want[tid]
            lst.append(dict(t, tid=tid, rep=r, count=n, tnov=tnov, gnov=gnov, key=key_of(tid, t, tnov), kl=klass(t, gnov)))
        out[r] = lst
    return out, tot


def key_of(tid, t, tnov):
    """Cross-replicate identity as wtc11_truth.key_of (ISM kept here, flagged separately)."""
    if tnov == "Known":
        return ("K", tid)
    if len(t["exons"]) < 2:
        return None
    return ("N", t["chrom"], t["strand"], tuple((a[1], b[0]) for a, b in zip(t["exons"], t["exons"][1:])))


def klass(t, gnov):
    if gnov in ("Antisense", "Intergenic"):
        return "l"
    if t["gene_type"] == "protein_coding":
        return "m"
    if t["gene_type"] in LNC_TYPES:
        return "l"
    return "o"


def v41_pc(chroms):
    out = defaultdict(lambda: {"+": [], "-": []})
    with gzip.open(f"{ROOT}/data/human/gencode.v41.annotation.gtf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if f[0] in chroms and f[2] == "gene" and 'gene_type "protein_coding"' in f[8]:
                out[f[0]][f[6]].append((int(f[3]), int(f[4])))
    return {c: {st: np.array(sorted(v), np.int64).reshape(-1, 2) for st, v in d.items()} for c, d in out.items()}


def s50_pairs(ts, preds, th=0.5):
    """simple_det.s50 that also returns, per prediction, the matched truth index (-1) and the reciprocal-overlap set."""
    starts = [t[2] for t in ts]
    maxlen = max(t[3] - t[2] for t in ts)
    adj = []
    for p in preds:
        st, s, e = p["strand"], p["exons"][0][0], p["exons"][-1][1]
        lo, hi = np.searchsorted(starts, s - maxlen - 1), np.searchsorted(starts, e, side="right")
        a = []
        for i in range(lo, hi):
            _, gst, gs, ge = ts[i]
            if gst != st:
                continue
            ov = min(e, ge) - max(s, gs) + 1
            if ov > 0 and ov / (ge - gs + 1) >= th and ov / (e - s + 1) >= th:
                a.append(i)
        adj.append(a)
    return B.matching(adj, len(ts)), adj


def recip(s, e, gs, ge):
    ov = min(e, ge) - max(s, gs) + 1
    return ov > 0 and ov / (ge - gs + 1) >= 0.5 and ov / (e - s + 1) >= 0.5


def ov_any(arr, s, e):
    return bool(len(arr)) and bool(((arr[:, 0] <= e) & (s <= arr[:, 1])).any())


def ends(st, s, e):
    """(5' end, 3' end) of a span on a strand."""
    return (s, e) if st == "+" else (e, s)


def cpg_island(seq, pos1):
    w = seq[max(pos1 - 1 - 500, 0):pos1 - 1 + 501].upper()
    n = len(w) - w.count("N")
    gc = (w.count("G") + w.count("C")) / n if n else 0.0
    oe = S3.cpg_oe(w)
    return {"gc": round(gc, 3), "oe": round(oe, 3), "cgi": bool(gc >= 0.5 and oe >= 0.6)}


def reads(bam, c, st, s, e):
    flt = ["-F", "2324"] if st == "+" else ["-f", "16", "-F", "2308"]
    same = int(subprocess.run([SAM, "view", "-c", *flt, bam, f"{c}:{s}-{e}"], capture_output=True, text=True, check=True).stdout)
    tot = int(subprocess.run([SAM, "view", "-c", "-F", "2308", bam, f"{c}:{s}-{e}"], capture_output=True, text=True, check=True).stdout)
    return same, tot


# ------------------------------------------------------------------ main
def main():
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)
    chroms = WS.TEST
    cfg = json.load(open(f"{ROOT}/configs/v20.json"))
    tab = json.load(open(f"{ROOT}/{cfg['calibration']['hexamer_table']}"))
    runs = cfg["targets"]["chr18"]["runs"]
    obs, tot = observed(set(chroms))
    print("observed", {r: len(v) for r, v in obs.items()}, "totals", tot, round(time.time() - t0), "s", flush=True)
    keys23 = {r: {t["key"] for t in obs[r] if t["key"] is not None and t["tnov"] != "ISM"} for r in ("rep2", "rep3")}
    pc41 = v41_pc(set(chroms))
    W24 = json.load(open(f"{L.W}/RESULT.json"))

    rows, c2rows = [], []
    removed = defaultdict(lambda: defaultdict(Counter))      # run -> chrom -> label counts of removed ev_farthest preds
    kept = defaultdict(lambda: defaultdict(Counter))
    ceiling = Counter()
    check = Counter()
    for c in chroms:
        tc = time.time()
        truth = L.load_truth(c)
        cuts = WS.cuts_of(c)
        truth["ts"] = [t for t in truth["ts"] if not any(t[2] - 1 < x < t[3] for x in cuts)]
        truth["lnc"] = {i: v for i, v in truth["lnc"].items() if i in {t[0] for t in truth["ts"]}}
        truth["zone"] = list(truth["zone"]) + [(x - 1, x + 1) for x in cuts]
        ts = truth["ts"]
        tsid = {t[0]: j for j, t in enumerate(ts)}
        tib = WS.read_tib(f"{L.W}/tiberius/{c}.gtf" if c in ("chr12", "chr13") else f"{ROOT}/outputs/PHASE1/BASELINES/tiberius/{c}.gtf", c)
        sig = L.quant_signal(c)
        ev = L.evidence(c, 1.0, 0)
        cands = candidates(ev, DMAX, K)
        wide = {(st, sp[0][0] if st == "+" else sp[0][1]): sp for st, _, sp in candidates(ev, DMAX, K_WIDE)}
        keep = lambda P: [p for p in P if not L.in_zone(truth["zone"], p["exons"][0][0], p["exons"][-1][1])]   # noqa: E731
        evf = keep(select_evidence(cands, sig, "farthest", tib))
        # candidate-set ceiling (no one-to-one, no mask): truth loci with some candidate span at reciprocal >= 0.5
        for kk, sets in (("k3", cands), ("k20", [(st, 0, sp) for (st, _), sp in wide.items()])):
            arr = {st: np.array([se for st2, _, sp in sets if st2 == st for se in sp], np.int64).reshape(-1, 2) for st in "+-"}
            for _, st, gs, ge in ts:
                a = arr[st]
                ov = np.minimum(a[:, 1], ge) - np.maximum(a[:, 0], gs) + 1
                ok = (ov > 0) & (ov / (ge - gs + 1) >= 0.5) & (ov / (a[:, 1] - a[:, 0] + 1) >= 0.5)
                ceiling[f"{kk}_hit"] += bool(ok.any())
        ceiling["n_truth"] += len(ts)

        # observed rep2 / rep3 transcripts on this chromosome, by strand
        o23 = defaultdict(list)
        for r in ("rep2", "rep3"):
            for t in obs[r]:
                if t["chrom"] == c:
                    o23[t["strand"]].append(t)
        o23_arr = {st: np.array([(t["s"], t["e"]) for t in v], np.int64).reshape(-1, 2) for st, v in o23.items()}

        ml, adj = s50_pairs(ts, evf)
        n_tp = sum(1 for m in ml if m != -1)
        check["evf_n"] += len(evf)
        check["evf_tp"] += n_tp
        lab = {}
        for j, p in enumerate(evf):
            st, s, e = p["strand"], p["exons"][0][0], p["exons"][-1][1]
            p5, p3 = ends(st, s, e)
            r = {"chrom": c, "strand": st, "start": s, "end": e, "len": e - s + 1, "cage": p["tss_count"],
                 "v41pc_same": ov_any(pc41[c][st], s, e), "v41pc_anti": ov_any(pc41[c]["-" if st == "+" else "+"], s, e)}
            if ml[j] != -1:
                r.update(label="TP", sub="", truth_id=ts[ml[j]][0])
            else:
                a = o23_arr.get(st, np.zeros((0, 2), np.int64))
                hit = np.nonzero((a[:, 0] <= e) & (s <= a[:, 1]))[0] if len(a) else []
                ovt = [o23[st][i] for i in hit]
                tl = [t for t in ts if t[1] == st and t[2] <= e and s <= t[3]]
                if not ovt:            # W24's 'uncertain' rule looks at observed rep2 / rep3 spans only
                    anti = o23_arr.get("-" if st == "+" else "+", np.zeros((0, 2), np.int64))
                    r.update(label="WRONG", sub="antisense_obs" if ov_any(anti, s, e) else "none")
                elif tl:
                    r["label"] = "UNC"
                    if adj[j]:
                        r.update(sub="U-a:taken", truth_id=ts[adj[j][0]][0])
                    else:
                        best = max(tl, key=lambda t: min(e, t[3]) - max(s, t[2]) + 1)
                        tx = truth["lnc"][best[0]]["tx"]
                        t5 = [ends(st, x["exons"][0][0], x["exons"][-1][1])[0] for x in tx]
                        t3 = [ends(st, x["exons"][0][0], x["exons"][-1][1])[1] for x in tx]
                        d5 = min(abs(p5 - x) for x in t5)
                        j3 = min(range(len(t3)), key=lambda i: abs(p3 - t3[i]))
                        d3s = (p3 - t3[j3]) * (1 if st == "+" else -1)      # > 0: our 3' end lies downstream
                        k5, k3 = d5 > OFF, abs(d3s) > OFF
                        sub = "both" if k5 and k3 else "5p" if k5 else "3p" if k3 else "span"
                        r.update(sub=f"U-a:{sub}", truth_id=best[0], d5=d5, d3=d3s, len_ratio=round((e - s + 1) / (best[3] - best[2] + 1), 3))
                        if not k5:     # 5' end right: would another candidate end have matched this locus?
                            for kk, sp in (("alt_k3", next((sp for st2, _, sp in cands if st2 == st and (sp[0][0] if st == "+" else sp[0][1]) == p5), [])),
                                           ("alt_k20", wide.get((st, p5), []))):
                                r[kk] = any(recip(a_, b_, best[2], best[3]) for a_, b_ in sp)
                else:
                    kl = {t["kl"] for t in ovt}
                    if "m" in kl:
                        r.update(label="UNC", sub="U-b")
                    elif "l" in kl:
                        lt = [t for t in ovt if t["kl"] == "l"]
                        if any(t["key"] is not None and t["tnov"] != "ISM" and t["key"] in keys23["rep2"] and t["key"] in keys23["rep3"] for t in lt):
                            why = "both"
                        elif any(t["key"] is not None and t["tnov"] != "ISM" for t in lt):
                            why = "single"
                        elif any(t["tnov"] == "ISM" for t in lt):
                            why = "ism"
                        else:
                            why = "mono"
                        r.update(label="UNC", sub=f"U-c:{why}")
                    else:
                        r.update(label="UNC", sub="U-d", types=";".join(sorted({str(t["gene_type"]) for t in ovt})))
            lab[(st, s, e)] = r
            rows.append(r)

        # ours per run (W24 'ours' = ev_farthest units kept by the v2.0 span score >= 0)
        per_ours = defaultdict(list)
        for run in runs:
            cs, _ = WS.run_inputs(c, run, cfg, tab)
            ours = keep(select_hybrid(cands, cs, sig, "farthest", WS.THETA, tib))
            del cs
            ko = {(p["strand"], p["exons"][0][0], p["exons"][-1][1]) for p in ours}
            assert ko <= set(lab), "ours must be a subset of ev_farthest"
            check[f"ours_n|{run}"] += len(ours)
            for k_, r in lab.items():
                tag = r["label"] if r["label"] != "UNC" else r["sub"].split(":")[0]
                (kept if k_ in ko else removed)[run][c][tag] += 1
                if k_ not in ko and ":" in r["sub"]:     # sub-labels (U-a:*, U-c:*) on top of the top-level tag
                    removed[run][c][r["sub"]] += 1
            mo, _ = s50_pairs(ts, ours)
            m_ids = {ts[i][0] for i in mo if i != -1}
            for t in ts:
                per_ours[t[0]].append(t[0] in m_ids)
            print(c, run, "ours", len(ours), "tp", len(m_ids), round(time.time() - tc), "s", flush=True)

        # C-2: long-read baselines at the truth loci
        lr = {}
        for tool in ("isoquant", "stringtie"):
            raw = WS.longread(c, tool)
            filt = keep([p for p in raw if not tib_hit(tib, p["strand"], p["exons"][0][0], p["exons"][-1][1])])
            mm, aa = s50_pairs(ts, filt)
            lr[tool] = {"filt": filt, "ids": {ts[i][0] for i in mm if i != -1}, "adj": aa}
        gtf = f"{L.W}/longread/{c}/isoquant/rep1/rep1.transcript_models.gtf"
        iso_tx = defaultdict(list)
        for line in open(gtf):
            f = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(f) < 9 or f[0] != c or f[2] != "transcript" or f[6] not in "+-":
                continue
            iso_tx[f[6]].append((int(f[3]), int(f[4])))
        iso_tx = {k_: np.array(sorted(v), np.int64).reshape(-1, 2) for k_, v in iso_tx.items()}
        raw_iso = defaultdict(list)
        for p in WS.longread(c, "isoquant"):
            raw_iso[p["strand"]].append((p["exons"][0][0], p["exons"][-1][1]))
        seq = S3.read_fasta(f"{ROOT}/data/human/{c}.fa")
        bam = f"{L.W}/longread/{c}/rep1.bam"
        k1 = defaultdict(int)
        for t in obs["rep1"]:
            if t["chrom"] == c and t["key"] is not None:
                k1[t["key"]] += t["count"]
        for t in ts:
            gid, st, s, e = t
            loc = truth["lnc"][gid]
            ours_m = float(np.mean(per_ours[gid]))
            iso_m = gid in lr["isoquant"]["ids"]
            p5 = ends(st, s, e)[0]
            same, total = reads(bam, c, st, s, e)
            keys = set()
            for x in loc["tx"]:
                kk = ("K", x["id3"]) if x["nov"] == "Known" else (
                    ("N", c, st, tuple((a[1], b[0]) for a, b in zip(x["exons"], x["exons"][1:]))) if len(x["exons"]) > 1 else None)
                if kk is not None:
                    keys.add(kk)
            a_tx = iso_tx.get(st, np.zeros((0, 2), np.int64))
            ov_tx = a_tx[(a_tx[:, 0] <= e) & (s <= a_tx[:, 1])] if len(a_tx) else a_tx
            gs = [g for g in raw_iso[st] if g[0] <= e and s <= g[1]]
            fs = [(p["exons"][0][0], p["exons"][-1][1]) for p in lr["isoquant"]["filt"] if p["strand"] == st
                  and p["exons"][0][0] <= e and s <= p["exons"][-1][1]]
            if iso_m:
                site = "found"
            elif not len(ov_tx):
                site = "absent"
            elif not fs:
                site = "masked"
            elif any(recip(a_, b_, s, e) for a_, b_ in fs):
                site = "taken"
            else:
                site = "mismatch"
            c2rows.append({"chrom": c, "id": gid, "strand": st, "start": s, "end": e, "len": e - s + 1,
                           "max_exons": max(len(x["exons"]) for x in loc["tx"]), "n_tx": loc["n_tx"], "novel": loc["novel"],
                           "ours": ours_m, "isoquant": iso_m, "stringtie": gid in lr["stringtie"]["ids"], "iso_site": site,
                           "iso_tx_recip": bool(any(recip(int(a_), int(b_), s, e) for a_, b_ in ov_tx)),
                           "iso_gene_len_ratio": round(max((g[1] - g[0] + 1) for g in gs) / (e - s + 1), 3) if gs else None,
                           "rep1_reads_same": same, "rep1_reads_all": total, "talon_rep1": sum(k1.get(k_, 0) for k_ in keys),
                           "talon_rep2": loc["count2"], "talon_rep3": loc["count3"], **cpg_island(seq, p5)})
        del seq
        print(c, "done", round(time.time() - tc), "s", flush=True)

    # ---------------------------------------------------------------- checks against W24 RESULT.json
    pw = W24["pooled"]
    assert int(check["evf_n"]) == int(round(pw["ev_farthest|f1.0_s0"]["n_pred"])), (check["evf_n"], pw["ev_farthest|f1.0_s0"])
    assert abs(check["evf_tp"] - pw["ev_farthest|f1.0_s0"]["S50_sn"] * W24["n_truth_lnc"]) < 0.5
    n_ours = np.mean([check[f"ours_n|{r}"] for r in runs])
    assert abs(n_ours - pw["ours|f1.0_s0"]["n_pred"]) < 0.5, (n_ours, pw["ours|f1.0_s0"])
    only = [r for r in c2rows if r["ours"] >= 0.5 and not r["isoquant"]]
    assert len(only) == W24["ours_not_isoquant"]["n"], len(only)

    # ---------------------------------------------------------------- summaries
    lab_c = Counter(r["label"] for r in rows)
    sub_c = Counter(r["sub"] for r in rows if r["label"] == "UNC")
    top = Counter(r["sub"].split(":")[0] for r in rows if r["label"] == "UNC")
    ua = [r for r in rows if r["sub"].startswith("U-a:") and r["sub"] != "U-a:taken"]
    ua_ok5 = [r for r in ua if r["d5"] <= OFF]
    rem = {}
    for run in runs:
        tot_r = Counter()
        for c in chroms:
            tot_r.update(removed[run][c])
        rem[run] = dict(tot_r)
    keys_all = sorted({k_ for v in rem.values() for k_ in v})
    rem_mean = {k_: float(np.mean([rem[r].get(k_, 0) for r in runs])) for k_ in keys_all}
    kept_mean = {}
    for k_ in ("TP", "U-a", "U-b", "U-c", "U-d", "WRONG"):
        kept_mean[k_] = float(np.mean([sum(kept[r][c][k_] for c in chroms) for r in runs]))

    def grp(sel):
        if not sel:
            return {"n": 0}
        f = lambda k_: float(np.median([x[k_] for x in sel]))   # noqa: E731
        return {"n": len(sel), "median_rep1_reads_same": f("rep1_reads_same"), "median_talon_rep1": f("talon_rep1"),
                "median_talon_rep23": float(np.median([x["talon_rep2"] + x["talon_rep3"] for x in sel])),
                "median_len": f("len"), "frac_cgi": float(np.mean([x["cgi"] for x in sel])),
                "frac_novel": float(np.mean([x["novel"] for x in sel])), "frac_spliced": float(np.mean([x["max_exons"] > 1 for x in sel])),
                "frac_rep1_reads_lt5": float(np.mean([x["rep1_reads_same"] < 5 for x in sel])),
                "frac_rep1_reads_0": float(np.mean([x["rep1_reads_same"] == 0 for x in sel])),
                "rep1_reads_same_q25_q75": np.percentile([x["rep1_reads_same"] for x in sel], [25, 75]).tolist(),
                "frac_talon_rep1_0": float(np.mean([x["talon_rep1"] == 0 for x in sel]))}

    res = {"plan": "docs/auth/PLAN_AUTH.md section 2", "chroms": list(chroms), "runs": runs,
           "check_vs_W24": {"ev_farthest_n": int(check["evf_n"]), "ev_farthest_tp": int(check["evf_tp"]), "ours_n_mean": float(n_ours),
                            "ours_not_isoquant": len(only), "ok": True},
           "C1": {"ev_farthest_labels": dict(lab_c), "unc_top": dict(top), "unc_sub": dict(sub_c),
                  "unc_w24_estimate": int(round(pw["ev_farthest|f1.0_s0"]["n_pred"] - check["evf_tp"] / pw["ev_farthest|f1.0_s0"]["S50_pr_known"])),
                  "U_a_ends": {"n": len(ua), "d5_le_1kb": len(ua_ok5),
                               "d3_downstream": sum(r["d3"] > OFF for r in ua), "d3_upstream": sum(r["d3"] < -OFF for r in ua),
                               "median_len_ratio": float(np.median([r["len_ratio"] for r in ua])) if ua else None,
                               "alt_k3_would_match": sum(bool(r.get("alt_k3")) for r in ua_ok5),
                               "alt_k20_would_match": sum(bool(r.get("alt_k20")) for r in ua_ok5)},
                  "v41pc_same_strand": {lb: sum(r["v41pc_same"] for r in rows if (r["label"] if r["label"] != "UNC" else r["sub"].split(":")[0]) == lb)
                                        for lb in ("TP", "U-a", "U-b", "U-c", "U-d", "WRONG")},
                  "wrong_sub": dict(Counter(r["sub"] for r in rows if r["label"] == "WRONG")),
                  "after_v41pc_mask_approx": dict(Counter((r["label"] if r["label"] != "UNC" else r["sub"].split(":")[0])
                                                          for r in rows if not r["v41pc_same"])),
                  "removed_by_sequence_filter_mean": rem_mean, "removed_by_run": rem, "kept_by_sequence_filter_mean": kept_mean,
                  "candidate_ceiling": {"n_truth": ceiling["n_truth"], "k3": ceiling["k3_hit"], "k20": ceiling["k20_hit"]}},
           "C2": {"ours_not_isoquant": grp(only),
                  "iso_found": grp([r for r in c2rows if r["isoquant"]]),
                  "missed_by_both": grp([r for r in c2rows if r["ours"] < 0.5 and not r["isoquant"]]),
                  "ours_not_isoquant_site": dict(Counter(r["iso_site"] for r in only)),
                  "ours_not_isoquant_iso_tx_recip": sum(r["iso_tx_recip"] for r in only),
                  "ours_not_isoquant_stringtie": sum(r["stringtie"] for r in only),
                  "n_truth": len(c2rows)},
           "seconds": round(time.time() - t0, 1)}
    json.dump(res, open(f"{OUT}/RESULT_C.json", "w"), indent=1)
    cols = ["chrom", "strand", "start", "end", "len", "cage", "label", "sub", "truth_id", "d5", "d3", "len_ratio", "alt_k3",
            "alt_k20", "types", "v41pc_same", "v41pc_anti"]
    with open(f"{OUT}/c1_predictions.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(k_, "")) for k_ in cols) + "\n")
    cols2 = list(c2rows[0])
    with open(f"{OUT}/c2_loci.tsv", "w") as fh:
        fh.write("\t".join(cols2) + "\n")
        for r in c2rows:
            fh.write("\t".join(str(r[k_]) for k_ in cols2) + "\n")
    print(json.dumps({k_: v for k_, v in res.items() if k_ in ("check_vs_W24", "C1", "C2")}, indent=1, default=str))


if __name__ == "__main__":
    main()
