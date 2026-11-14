#!/usr/bin/env python3
"""Selective unified GFF3 (docs/phase1_results/PLAN_UGFF.md).

  ugff.py dev                 chr18 development: guardrail effect, intron-precision curve -> t_int, CpG cut-off,
                              start check  -> outputs/PHASE1/UGFF/RULES.json
  ugff.py build --chrom C     C in chr18 / chr20 / chr9 -> outputs/PHASE1/UGFF/C/{unified.gff3, SUMMARY.json}

No model is re-read: lnc layers come from the genes.json / post.npy of reads already made (chr20 v1.9, chr18 / chr9
v2.0).  Coding layer = Tiberius GTF (outputs/PHASE1/BASELINES/tiberius/<chrom>.gtf).  CPU only.
"""
import argparse
import json
import math
import os
import sys
import time
from collections import defaultdict

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
for p in ("src", "scripts/phase1", "scripts/phase1/v20", "scripts/phase1/vipsania", "scripts/phase1/s3"):
    sys.path.insert(0, os.path.join(ROOT, p))
import p15_score as P  # noqa: E402
import v20_calib as C  # noqa: E402
from ncrna_annot.eval import ruler  # noqa: E402
from ncrna_annot.grammar import g1  # noqa: E402

OUT = f"{ROOT}/outputs/PHASE1/UGFF"
TIB = f"{ROOT}/outputs/PHASE1/BASELINES/tiberius"
GRID = [round(0.50 + 0.05 * i, 2) for i in range(10)]
V20 = ["v20-g1b-42-ke", "v20-g1b-123-ke", "v20-g1b-456-ke"]
CH = {
    "chr18": dict(cfg="configs/v20.json", rd="outputs/PHASE1/V20/read/chr18", runs=V20, model="v2.0", small=None),
    "chr9": dict(cfg="configs/v20.json", rd="outputs/PHASE1/V20/read/chr9", runs=V20, model="v2.0", small=None),
    "chr20": dict(cfg="configs/p15_v19.json", rd="outputs/PHASE1/P15/chr20", runs=["g1b-42-ke", "g1b-123-ke", "g1b-456-ke"],
                  model="v1.9", small="human_chr20"),
}


# ------------------------------------------------------------------ inputs
def load_chrom(chrom):
    c = CH[chrom]
    cfg = json.load(open(f"{ROOT}/{c['cfg']}"))
    tg = cfg["targets"][chrom]
    truth = ruler.read_gtf(f"{ROOT}/{tg['truth']}", tg["truth_chrom"])
    codes = np.load(f"{ROOT}/{tg['codes']}", mmap_mode="r")
    cal = cfg.get("calibration") if c["model"] == "v2.0" else None
    if cal:
        assert cal["keep"] and P.sha256(f"{ROOT}/{cal['hexamer_table']}") == cal["hexamer_table_sha256"]
        tab = json.load(open(f"{ROOT}/{cal['hexamer_table']}"))
        tau = cal["tau_v20"]
    else:
        tau = cfg["decode"]["tau"]
    runs = {}
    for run in c["runs"]:
        rd = f"{ROOT}/{c['rd']}/{run}"
        genes = P.load_genes(f"{rd}/genes.json")
        post = np.load(f"{rd}/post.npy", mmap_mode="r")
        if cal:
            rows = C.gene_rows(genes, post, codes, tab)
            X = np.array([[r[f] for f in cal["model"]["features"]] for r in rows], float)
            sc = dict(zip((r["k"] for r in rows), C.apply(cal["model"], X)))
            for k, g in enumerate(genes):
                if k in sc:
                    g["conf"] = float(sc[k])
        lnc = [dict(g, k=k) for k, g in enumerate(genes) if g["type"] == "l" and g["conf"] >= tau]
        runs[run] = {"lnc": lnc, "post": post}
    return {"truth": truth, "codes": codes, "tau": tau, "runs": runs, "cfg": c}


def read_tiberius(chrom):
    tx = {}
    for line in open(f"{TIB}/{chrom}.gtf"):
        f = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(f) < 9 or f[0] != chrom or f[2] not in ("CDS", "exon"):
            continue
        tid = f[8].split('transcript_id "', 1)[1].split('"', 1)[0]
        gid = f[8].split('gene_id "', 1)[1].split('"', 1)[0]
        t = tx.setdefault(tid, {"id": tid, "gene": gid, "strand": f[6], "exons": [], "cds": []})
        if f[2] == "exon":
            t["exons"].append((int(f[3]), int(f[4])))
        else:
            t["cds"].append((int(f[3]), int(f[4]), f[7]))
    for t in tx.values():
        t["exons"] = ruler._merge(sorted(t["exons"] or [(a, b) for a, b, _ in t["cds"]]))
        t["cds"].sort()
    return tx


class IvIndex:
    """Same-strand interval overlap on plain (strand, a, b) intervals."""
    BIN = 100_000

    def __init__(self, ivs):
        self.bins = defaultdict(list)
        for s, a, b, key in ivs:
            for x in range(a // self.BIN, b // self.BIN + 1):
                self.bins[(s, x)].append((a, b, key))

    def hits(self, s, exons):
        out = set()
        for a, b in exons:
            for x in range(a // self.BIN, b // self.BIN + 1):
                out |= {k for ga, gb, k in self.bins.get((s, x), ()) if ga <= b and a <= gb}
        return out


# ------------------------------------------------------------------ per-gene quantities
INTRON = {s: [g1.IDX[f"l{p}{s}"] for p in ("L", "I", "R")] for s in g1.STRANDS}


def intron_conf(g, post):
    ex = g["exons"]
    return [float(np.asarray(post[b1:a2 - 1][:, INTRON[g["strand"]]], np.float64).sum(1).mean())
            for (_, b1), (a2, _) in zip(ex, ex[1:])]


def five_prime(strand, exons):
    return exons[0][0] if strand == "+" else exons[-1][1]


def cpg_at(codes, p):
    s = "".join("ACGTN"[min(int(c), 4)] for c in codes[max(0, p - 501):p + 500])
    c, g, cg = s.count("C"), s.count("G"), s.count("CG")
    return cg * len(s) / (c * g) if c and g else 0.0


def truth_sets(truth):
    lnc_int, any_int, lnc_chain = set(), set(), set()
    for g in truth.values():
        for ex in g["tx"].values():
            itr = ruler.introns(ruler._merge(sorted(ex)))
            any_int |= {(g["strand"],) + i for i in itr}
            if g["cls"] == "lnc":
                lnc_int |= {(g["strand"],) + i for i in itr}
                if itr:
                    lnc_chain.add((g["strand"], itr))
    return lnc_int, any_int, lnc_chain


def guard(lnc, cds_idx, tidx, truth):
    keep, drop = [], defaultdict(int)
    for g in lnc:
        if cds_idx.hits(g["strand"], g["exons"]):
            h = tidx.hits(g["strand"], g["exons"])
            cl = {truth[x]["cls"] for x in h}
            drop["lnc" if "lnc" in cl else "mrna_only" if cl == {"mrna"} else "other" if cl else "none"] += 1
        else:
            keep.append(g)
    return keep, dict(drop, total=sum(drop.values()))


def structure_stats(genes, ts, t):
    """Written introns / chains under threshold t (None = forced: every decoded transcript written)."""
    lnc_int, any_int, lnc_chain = ts
    n_tx = n_int = n_lnc = n_any = n_chain = 0
    for g in genes:
        if len(g["exons"]) < 2 or (t is not None and min(g["iconf"]) < t):
            continue
        itr = ruler.introns(g["exons"])
        n_tx += 1
        n_int += len(itr)
        n_lnc += sum((g["strand"],) + i in lnc_int for i in itr)
        n_any += sum((g["strand"],) + i in any_int for i in itr)
        n_chain += (g["strand"], itr) in lnc_chain
    return {"n_tx_multi": n_tx, "n_introns": n_int, "intron_pr_lnc": n_lnc / n_int if n_int else None,
            "intron_pr_any": n_any / n_int if n_int else None, "exact_lnc_chains": n_chain}


def prepare(chrom):
    d = load_chrom(chrom)
    truth, codes = d["truth"], d["codes"]
    tib = read_tiberius(chrom)
    cds_idx = IvIndex([(t["strand"], a, b, t["id"]) for t in tib.values() for a, b, _ in t["cds"]])
    tidx = ruler.ExonIndex(truth)
    for run, r in d["runs"].items():
        r["kept"], r["guard"] = guard(r["lnc"], cds_idx, tidx, truth)
        for g in r["kept"]:
            g["iconf"] = intron_conf(g, r["post"])
    d["tib"], d["tidx"] = tib, tidx
    return d


def lnc_preds(genes):
    return [{"id": str(g["k"]), "strand": g["strand"], "cls": "lnc", "exons": [tuple(x) for x in g["exons"]]} for g in genes]


def tib_preds(tib):
    return [{"id": t["id"], "strand": t["strand"], "cls": "mrna", "exons": t["exons"]} for t in tib.values()]


def l1(truth, preds, cls="lnc"):
    s = ruler.score_class(cls, truth, preds, ruler.ExonIndex(truth), 0)
    return {"sn": s["L1_locus"]["sn"], "pr_strict": s["L1_locus"]["pr_strict"], "n_pred": s["n_pred_tx"],
            **({"L2": {k: v for k, v in s["L2_intron_chain"].items() if not isinstance(v, (list, dict))}}
               if cls == "mrna" else {})}


def mean(xs):
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else None


# ------------------------------------------------------------------ dev (chr18)
def dev():
    t0 = time.time()
    d = prepare("chr18")
    truth, codes = d["truth"], d["codes"]
    ts = truth_sets(truth)
    res = {"chrom": "chr18", "tau": d["tau"], "runs": {}, "curve": {}}
    for run, r in d["runs"].items():
        res["runs"][run] = {"n_kept": len(r["lnc"]), "guard_dropped": r["guard"],
                            "lnc_L1_before": l1(truth, lnc_preds(r["lnc"])), "lnc_L1_after": l1(truth, lnc_preds(r["kept"])),
                            "forced": structure_stats(r["kept"], ts, None)}
    for t in GRID:
        per = {run: structure_stats(r["kept"], ts, t) for run, r in d["runs"].items()}
        res["curve"][str(t)] = {"runs": per, "mean_intron_pr_lnc": mean([p["intron_pr_lnc"] for p in per.values()]),
                                "min_introns": min(p["n_introns"] for p in per.values())}
    ok = [t for t in GRID if (res["curve"][str(t)]["mean_intron_pr_lnc"] or 0) >= 0.50 and res["curve"][str(t)]["min_introns"] >= 20]
    res["t_int"] = ok[0] if ok else None
    # CpG cut-off: upper quartile of promoter CpG o/e at annotated v50 lnc gene 5' ends (chr18)
    starts = []
    for g in truth.values():
        if g["cls"] == "lnc":
            ex = sorted(e for t in g["tx"].values() for e in t)
            starts.append(ex[0][0] if g["strand"] == "+" else max(b for _, b in ex))
    res["cpg_cut"] = float(np.quantile([cpg_at(codes, p) for p in starts], 0.75))
    res["n_truth_lnc_starts"] = len(starts)
    # start check: predictions hitting a v50 lnc gene, 5' end within 1 kb of a 5' end of a hit gene's transcripts
    lidx = ruler.ExonIndex({k: v for k, v in truth.items() if v["cls"] == "lnc"})
    chk = {"CpG_rich": [0, 0], "other": [0, 0]}
    for r in d["runs"].values():
        for g in r["kept"]:
            h = lidx.hits(g["strand"], g["exons"])
            if not h:
                continue
            p5 = five_prime(g["strand"], g["exons"])
            ends = [five_prime(truth[x]["strand"], ruler._merge(sorted(t))) for x in h for t in truth[x]["tx"].values()]
            cl = "CpG_rich" if cpg_at(codes, p5) >= res["cpg_cut"] else "other"
            chk[cl][0] += 1
            chk[cl][1] += any(abs(p5 - e) <= 1000 for e in ends)
    res["start_check"] = {k: {"n_hit_preds": n, "within_1kb": w, "share": w / n if n else None} for k, (n, w) in chk.items()}
    res["seconds"] = round(time.time() - t0)
    os.makedirs(OUT, exist_ok=True)
    json.dump(P.norm(res), open(f"{OUT}/RULES.json", "w"), indent=1)
    print(json.dumps({"t_int": res["t_int"], "cpg_cut": res["cpg_cut"], "start_check": res["start_check"],
                      "curve": {t: v["mean_intron_pr_lnc"] for t, v in res["curve"].items()}}))


# ------------------------------------------------------------------ build
def build(chrom):
    t0 = time.time()
    rules = json.load(open(f"{OUT}/RULES.json"))
    t_int, cut = rules["t_int"], rules["cpg_cut"]
    d = prepare(chrom)
    truth, codes, c = d["truth"], d["codes"], d["cfg"]
    ship = c["runs"][0]
    genes = d["runs"][ship]["kept"]
    others = IvIndex([(g["strand"], a, b, (run, i)) for run in c["runs"] for i, g in enumerate(d["runs"][run]["kept"])
                      for a, b in g["exons"]])
    recs = []  # (start, end, order, lines)
    src = "ncrna_annot"
    assert len({t["gene"] for t in d["tib"].values()}) == len(d["tib"]), "Tiberius gene with several transcripts"
    for t in sorted(d["tib"].values(), key=lambda t: (t["exons"][0][0], t["id"])):
        gid, tid = f"tiberius_{chrom}_{t['gene']}", f"tiberius_{chrom}_{t['id']}"
        a, b, s = t["exons"][0][0], t["exons"][-1][1], t["strand"]
        ln = [f"{chrom}\tTiberius\tgene\t{a}\t{b}\t.\t{s}\t.\tID={gid};biotype=protein_coding",
              f"{chrom}\tTiberius\tmRNA\t{a}\t{b}\t.\t{s}\t.\tID={tid};Parent={gid}"]
        ln += [f"{chrom}\tTiberius\texon\t{x}\t{y}\t.\t{s}\t.\tID={tid}.e{j};Parent={tid}" for j, (x, y) in enumerate(t["exons"])]
        ln += [f"{chrom}\tTiberius\tCDS\t{x}\t{y}\t.\t{s}\t{ph}\tID={tid}.cds{j};Parent={tid}" for j, (x, y, ph) in enumerate(t["cds"])]
        recs.append((a, b, 0, ln))
    n_res = 0
    lnc_ivs = []
    for i, g in enumerate(sorted(genes, key=lambda g: g["exons"][0][0])):
        ex, s = g["exons"], g["strand"]
        a, b = ex[0][0], ex[-1][1]
        gid = f"{src}_{chrom}_L{i}"
        sup = len({run for run, _ in others.hits(s, ex)})
        p5 = five_prime(s, ex)
        oe = cpg_at(codes, p5)
        resolved = t_int is not None and len(ex) > 1 and min(g["iconf"]) >= t_int
        n_res += resolved
        att = [f"ID={gid}", "biotype=lncRNA", f"confidence={g['conf']:.4f}", f"support={sup}/3",
               f"structure={'resolved' if resolved else 'unresolved'}", f"n_decoded_exons={len(ex)}"]
        if g["iconf"]:
            att.append(f"min_intron_conf={min(g['iconf']):.4f}")
        att += [f"promoter_cpg_oe={oe:.3f}", f"start_class={'CpG_rich' if oe >= cut else 'other'}", "end_precision=low",
                f"model=ncrna_annot_{c['model']}"]
        if not resolved:
            att.append("decoded_exons=" + ",".join(f"{x}-{y}" for x, y in ex))
        ln = [f"{chrom}\t{src}\tlncRNA_gene\t{a}\t{b}\t{g['conf']:.4f}\t{s}\t.\t{';'.join(att)}"]
        if resolved:
            ln.append(f"{chrom}\t{src}\tlnc_RNA\t{a}\t{b}\t{g['conf']:.4f}\t{s}\t.\tID={gid}.t1;Parent={gid}")
            ln += [f"{chrom}\t{src}\texon\t{x}\t{y}\t.\t{s}\t.\tID={gid}.t1.e{j};Parent={gid}.t1" for j, (x, y) in enumerate(ex)]
        recs.append((a, b, 1, ln))
        lnc_ivs += [(s, x, y, gid) for x, y in ex]
    n_small = 0
    if c["small"]:
        import p16_score as P16
        lidx = IvIndex(lnc_ivs)
        for i, p in enumerate(P16.small_hits(c["small"], chrom)):
            (a, b), s = (p["exons"][0][0], p["exons"][-1][1]), p["strand"]
            gid, hit = f"{src}_small_{chrom}_{i}", p["src"].replace(";", ",").replace("=", ":")
            ov = sorted(lidx.hits(s, p["exons"]))
            note = f";Note=overlaps_lnc:{','.join(ov)}" if ov else ""
            recs.append((a, b, 2, [f"{chrom}\t{src}\tncRNA_gene\t{a}\t{b}\t.\t{s}\t.\tID={gid};biotype=small_ncRNA;source_hit={hit}{note}",
                                   f"{chrom}\t{src}\tncRNA\t{a}\t{b}\t.\t{s}\t.\tID={gid}.t1;Parent={gid}",
                                   f"{chrom}\t{src}\texon\t{a}\t{b}\t.\t{s}\t.\tID={gid}.t1.e0;Parent={gid}.t1"]))
            n_small += 1
    recs.sort(key=lambda r: (r[0], r[2], r[1]))
    os.makedirs(f"{OUT}/{chrom}", exist_ok=True)
    path = f"{OUT}/{chrom}/unified.gff3"
    with open(path, "w") as fh:
        fh.write(f"##gff-version 3\n##sequence-region {chrom} 1 {len(codes)}\n")
        fh.write(f"#! coding: Tiberius mammalia_nosofttmasking_v2 (CDS only, no UTR)\n")
        fh.write(f"#! lncRNA: ncrna_annot {c['model']} run {ship}, threshold {d['tau']}; Tiberius-CDS guardrail; "
                 f"transcripts only when every intron confidence >= {t_int} (PLAN_UGFF.md)\n")
        fh.write(f"#! small ncRNA: {'Rfam/tRNAscan (Phase 0 0.6, P1.6 rules)' if c['small'] else 'not run for this chromosome'}\n")
        for r in recs:
            fh.write("\n".join(r[3]) + "\n")
    val = validate(path, chrom)
    rt = ruler.read_pred_gff3(path, chrom)
    ts = truth_sets(truth)
    summ = {"chrom": chrom, "model": c["model"], "shipped_run": ship, "tau": d["tau"], "t_int": t_int, "cpg_cut": cut,
            "n_coding_tx": len(d["tib"]), "n_lnc_genes": len(genes), "n_lnc_resolved_tx": n_res, "n_small": n_small,
            "validation": val, "roundtrip": {"n_lnc_RNA_read": sum(p["cls"] == "lnc" for p in rt),
                                              "n_mRNA_read": sum(p["cls"] == "mrna" for p in rt)},
            "coding_tiberius": l1(truth, tib_preds(d["tib"]), "mrna"),
            "guard": {run: r["guard"] for run, r in d["runs"].items()},
            "lnc_L1_before_guard": {run: l1(truth, lnc_preds(r["lnc"])) for run, r in d["runs"].items()},
            "lnc_L1_after_guard": {run: l1(truth, lnc_preds(r["kept"])) for run, r in d["runs"].items()},
            "structure": {"forced": structure_stats(genes, ts, None),
                          "selective": structure_stats(genes, ts, t_int) if t_int is not None else {"n_tx_multi": 0}},
            "start_class_counts": None, "seconds": None}
    sc = defaultdict(int)
    for line in open(path):
        if "\tlncRNA_gene\t" in line:
            sc[line.split("start_class=")[1].split(";")[0]] += 1
    summ["start_class_counts"] = dict(sc)
    summ["sha256"] = P.sha256(path)
    summ["seconds"] = round(time.time() - t0)
    json.dump(P.norm(summ), open(f"{OUT}/{chrom}/SUMMARY.json", "w"), indent=1)
    print(json.dumps({k: summ[k] for k in ("chrom", "n_coding_tx", "n_lnc_genes", "n_lnc_resolved_tx", "n_small", "validation")}))


def validate(path, chrom):
    ids, parents, span, bad, prev, n = {}, [], {}, [], 0, 0
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        n += 1
        if len(f) != 9 or f[0] != chrom or int(f[3]) > int(f[4]):
            bad.append(("format", line[:80]))
            continue
        a = dict(x.split("=", 1) for x in f[8].split(";"))
        if a["ID"] in ids:
            bad.append(("dup_id", a["ID"]))
        ids[a["ID"]] = (int(f[3]), int(f[4]), f[6])
        if "Parent" in a:
            parents.append((a["ID"], a["Parent"]))
        else:
            if int(f[3]) < prev:
                bad.append(("unsorted", a["ID"]))
            prev = int(f[3])
    for i, p in parents:
        if p not in ids:
            bad.append(("no_parent", i))
        elif not (ids[p][0] <= ids[i][0] and ids[i][1] <= ids[p][1] and ids[p][2] == ids[i][2]):
            bad.append(("outside_parent", i))
    return {"n_rows": n, "n_ids": len(ids), "problems": len(bad), "first": bad[:5]}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["dev", "build"])
    ap.add_argument("--chrom", default="chr18", choices=list(CH))
    a = ap.parse_args()
    dev() if a.cmd == "dev" else build(a.chrom)
