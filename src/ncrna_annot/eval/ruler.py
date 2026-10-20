#!/usr/bin/env python3
"""Phase 0 ruler v1 (PROTOCOL.md section 2): score transcript predictions against GENCODE.

Standard library only. Coordinates are 1-based closed (GFF convention).

  ruler.py score --truth gencode.v50.gtf.gz [--old-truth gencode.v41.gtf.gz] --chrom chr20 \
                 --pred pred.gff3 --pred-format gff3 --out metrics.json
  ruler.py score ... --pred gencode.v41.gtf.gz --pred-format gtf     # old-vs-new self-test

Prediction classes: GFF3 transcript type lnc_RNA -> lnc, mRNA -> mrna. A GTF used as a
prediction maps gene_type lncRNA -> lnc and protein_coding -> mrna.
"""
import argparse
import gzip
import json
import random
import sys
from collections import defaultdict

SMALL = {"snRNA", "snoRNA", "miRNA", "pre_miRNA", "misc_RNA", "scaRNA", "vault_RNA", "Mt_tRNA",
         "Mt_rRNA", "rRNA", "ribozyme", "sRNA", "scRNA", "Y_RNA", "tRNA", "piRNA", "SRP_RNA",
         "RNase_P_RNA", "RNase_MRP_RNA", "telomerase_RNA"}
LNC = {"lncRNA", "lincRNA", "antisense", "antisense_RNA", "sense_overlapping", "sense_intronic",
       "3prime_overlapping_ncRNA", "bidirectional_promoter_lncRNA", "macro_lncRNA"}
LEN_SPLIT = {"ncRNA", "processed_transcript"}  # FlyBase/WormBase generic ncRNA, legacy zebrafish
BIN = 100_000
RULER_VERSION = "ruler-v1.3-20261020"  # v1.3: GFF3 truth, Tiberius attrs, --ignore-strand (metrics unchanged)


def gene_class(gene_type, max_tx_len=None):
    """PROTOCOL v1.2 class map. Generic ncRNA types split at 200 nt (longest spliced tx)."""
    if gene_type in LNC:
        return "lnc"
    if gene_type == "protein_coding":
        return "mrna"
    if gene_type in SMALL:
        return "small"
    if gene_type in LEN_SPLIT:
        return "lnc" if (max_tx_len or 0) >= 200 else "small"
    if "pseudogene" in gene_type:
        return "pseudo"
    return "other"


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def _attrs_gtf(s):
    out = {}
    for part in s.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        k, _, v = part.partition(" ")
        v = v.strip('"')
        if k == "tag":
            out.setdefault("tags", set()).add(v)
        else:
            out[k] = v
    return out


def _attrs_gff3(s):
    return dict(kv.split("=", 1) for kv in s.strip().split(";") if "=" in kv)


def read_gtf(path, chrom, basic_only=False):
    """GENCODE GTF -> {gene_id: {strand, cls, gene_type, tx: {tx_id: [(s, e), ...]}}}."""
    genes = {}
    with _open(path) as fh:
        for line in fh:
            if not line.startswith(chrom + "\t"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] != "exon":
                continue
            a = _attrs_gtf(f[8])
            if basic_only and "basic" not in a.get("tags", ()):
                continue
            gtype = a.get("gene_type") or a.get("gene_biotype", "")  # GENCODE / Ensembl
            g = genes.setdefault(a["gene_id"], {"strand": f[6], "gene_type": gtype,
                                                "cls": gene_class(gtype), "tx": {}})
            g["tx"].setdefault(a["transcript_id"], []).append((int(f[3]), int(f[4])))
    for g in genes.values():
        for t in g["tx"].values():
            t.sort()
        if g["gene_type"] in LEN_SPLIT:
            longest = max(sum(e - s + 1 for s, e in _merge(t)) for t in g["tx"].values())
            g["cls"] = gene_class(g["gene_type"], longest)
    return genes


def _tx_id(col9):
    """ID= (GFF3) or a bare id (Tiberius gene/mRNA rows)."""
    a = _attrs_gff3(col9)
    return a.get("ID") or col9.strip().split(";")[0].strip()


def _exon_parents(col9):
    """Parent= (GFF3) or transcript_id "x" (GTF-style, incl. doubled quotes from CSV export)."""
    a = _attrs_gff3(col9)
    if "Parent" in a:
        return a["Parent"].split(",")
    for part in col9.split(";"):
        part = part.strip().replace('""', '"').lstrip('"')
        if part.startswith("transcript_id"):
            return [part.split(None, 1)[1].strip().strip('"')]
    return []


def read_pred_gff3(path, chrom):
    """GFF3 with gene / mRNA|lnc_RNA / exon features -> list of predicted transcripts."""
    tx, exons = {}, defaultdict(list)
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.startswith(chrom + "\t"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if f[2] in ("mRNA", "lnc_RNA", "transcript"):
                tid = _tx_id(f[8])
                tx[tid] = {"id": tid, "strand": f[6],
                           "cls": "lnc" if f[2] == "lnc_RNA" else "mrna", "attrs": _attrs_gff3(f[8])}
            elif f[2] == "exon":
                for p in _exon_parents(f[8]):
                    exons[p].append((int(f[3]), int(f[4])))
    out = []
    for tid, t in tx.items():
        if exons.get(tid):
            t["exons"] = _merge(sorted(exons[tid]))
            out.append(t)
    return out


def read_gff3_truth(path, chrom):
    """RefSeq-style GFF3 truth (gene_biotype on gene rows) -> same structure as read_gtf."""
    gtype, tx_gene, tx_strand, ex = {}, {}, {}, defaultdict(list)
    with _open(path) as fh:
        for line in fh:
            if not line.startswith(chrom + "\t"):
                continue
            f = line.rstrip("\n").split("\t")
            a = _attrs_gff3(f[8])
            if f[2] == "gene":
                gtype[a["ID"]] = a.get("gene_biotype", "")
            elif f[2] == "exon":
                for p in a.get("Parent", "").split(","):
                    ex[p].append((int(f[3]), int(f[4])))
            elif "Parent" in a and "ID" in a:
                tx_gene[a["ID"]], tx_strand[a["ID"]] = a["Parent"].split(",")[0], f[6]
    genes = {}
    for tid, exons in ex.items():
        gid = tx_gene.get(tid)
        if gid is None or gid not in gtype:
            continue
        g = genes.setdefault(gid, {"strand": tx_strand[tid], "gene_type": gtype[gid],
                                   "cls": gene_class(gtype[gid]), "tx": {}})
        g["tx"][tid] = sorted(exons)
    for g in genes.values():
        if g["gene_type"] in LEN_SPLIT:
            g["cls"] = gene_class(g["gene_type"], max(sum(e - s + 1 for s, e in _merge(t)) for t in g["tx"].values()))
    return genes


def _unstrand_truth(genes):
    for g in genes.values():
        g["strand"] = "."
    return genes


def read_pred_smallrna(path, chrom, rfam_types=None):
    """Directory with cmscan.tsv (smallrna_modal.py) and trnascan.txt -> 'small' predictions.

    Rfam hits of lncRNA families and of non-gene types (Cis-reg, Intron) are dropped when an
    Rfam family-type table (rfam_acc -> type string) is given; overlapping clan losers ('=')
    are always dropped.
    """
    out = []
    for i, line in enumerate(open(f"{path}/cmscan.tsv")):
        f = line.rstrip("\n").split("\t")
        if i == 0 or f[0] != chrom or f[9] == "=":
            continue
        typ = (rfam_types or {}).get(f[5], "")
        if rfam_types is not None and ("lncRNA" in typ or typ.startswith(("Cis-reg", "Intron"))):
            continue
        out.append({"id": f"rfam{i}", "strand": f[3], "cls": "small", "src": f"{f[4]}|{typ}",
                    "exons": [(int(f[1]), int(f[2]))]})
    for i, line in enumerate(open(f"{path}/trnascan.txt")):
        f = line.split()
        if len(f) < 9 or not f[1].isdigit() or f[0] != chrom:
            continue
        b, e = int(f[2]), int(f[3])
        out.append({"id": f"trna{i}", "strand": "+" if b <= e else "-", "cls": "small",
                    "src": f"tRNA-{f[4]}", "exons": [(min(b, e), max(b, e))]})
    return out


def preds_from_gtf(path, chrom):
    out = []
    for gid, g in read_gtf(path, chrom).items():
        if g["cls"] not in ("lnc", "mrna", "small"):
            continue
        for tid, ex in g["tx"].items():
            out.append({"id": tid, "strand": g["strand"], "cls": g["cls"], "exons": _merge(ex)})
    return out


def _merge(iv):
    out = []
    for s, e in iv:
        if out and s <= out[-1][1] + 1:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def introns(ex):
    return tuple((ex[i][1] + 1, ex[i + 1][0] - 1) for i in range(len(ex) - 1))


class ExonIndex:
    """Same-strand exon-overlap lookup: which genes have an exon overlapping a query exon."""

    def __init__(self, genes):
        self.bins = defaultdict(list)
        for gid, g in genes.items():
            for ex in g["tx"].values():
                for s, e in ex:
                    for b in range(s // BIN, e // BIN + 1):
                        self.bins[(g["strand"], b)].append((s, e, gid))

    def hits(self, strand, exons):
        found = set()
        for s, e in exons:
            for b in range(s // BIN, e // BIN + 1):
                for gs, ge, gid in self.bins.get((strand, b), ()):
                    if gs <= e and s <= ge:
                        found.add(gid)
        return found


def _prf(tp_pred, n_pred, tp_truth, n_truth):
    p = tp_pred / n_pred if n_pred else None
    r = tp_truth / n_truth if n_truth else None
    f = 2 * p * r / (p + r) if p and r else (0.0 if p is not None and r is not None else None)
    return {"precision": p, "recall": r, "f1": f, "n_pred": n_pred, "n_truth": n_truth}


def _boot(vec, n=1000, seed=0):
    if not vec or n <= 0:
        return None
    rng, k = random.Random(seed), len(vec)
    ms = sorted(sum(vec[rng.randrange(k)] for _ in range(k)) / k for _ in range(n))
    return [ms[int(0.025 * n)], ms[int(0.975 * n) - 1]]


def _bases(ivs_by_strand):
    return {st: _merge(sorted(iv)) for st, iv in ivs_by_strand.items()}


def _overlap_len(a, b):
    i = j = tot = 0
    while i < len(a) and j < len(b):
        s, e = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if s <= e:
            tot += e - s + 1
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return tot


def score_class(cls, truth, preds, index, boot_n):
    """All ruler metrics for one class ('lnc' or 'mrna'). truth = all genes of the version."""
    tgenes = {gid: g for gid, g in truth.items() if g["cls"] == cls}
    cpreds = [p for p in preds if p["cls"] == cls]
    res = {"n_truth_genes": len(tgenes), "n_pred_tx": len(cpreds)}

    # L1 locus + novel + type confusion
    hit_genes, pred_hit, pred_novel, wrong_cls = set(), [], [], defaultdict(int)
    for p in cpreds:
        g = index.hits(p["strand"], p["exons"])
        same = {x for x in g if truth[x]["cls"] == cls}
        hit_genes |= same
        pred_hit.append(1 if same else 0)
        pred_novel.append(0 if g else 1)
        if g and not same:
            wrong_cls[sorted(truth[x]["cls"] for x in g)[0]] += 1
    sn_vec = [1 if gid in hit_genes else 0 for gid in tgenes]
    strict = [h for h in pred_hit]
    lenient = [h for h, nv in zip(pred_hit, pred_novel) if not nv]
    res["L1_locus"] = {
        "sn": sum(sn_vec) / len(sn_vec) if sn_vec else None, "sn_ci95": _boot(sn_vec, boot_n),
        "pr_strict": sum(strict) / len(strict) if strict else None, "pr_strict_ci95": _boot(strict, boot_n),
        "pr_lenient": sum(lenient) / len(lenient) if lenient else None,
        "n_novel_pred": sum(pred_novel), "pred_hits_other_class_only": dict(wrong_cls)}

    # L2 intron chain (multi-exon)
    chain2gene = {}
    for gid, g in tgenes.items():
        for ex in g["tx"].values():
            if len(ex) > 1:
                chain2gene[(g["strand"], introns(_merge(ex)))] = gid
    multi_genes = set(chain2gene.values())
    got, pm, pm_len = set(), [], []
    for p, nv in zip(cpreds, pred_novel):
        if len(p["exons"]) < 2:
            continue
        gid = chain2gene.get((p["strand"], introns(p["exons"])))
        if gid:
            got.add(gid)
        pm.append(1 if gid else 0)
        if not nv:
            pm_len.append(1 if gid else 0)
    sn2 = [1 if gid in got else 0 for gid in multi_genes]
    res["L2_intron_chain"] = {
        "sn": sum(sn2) / len(sn2) if sn2 else None, "sn_ci95": _boot(sn2, boot_n),
        "pr_strict": sum(pm) / len(pm) if pm else None, "pr_strict_ci95": _boot(pm, boot_n),
        "pr_lenient": sum(pm_len) / len(pm_len) if pm_len else None,
        "n_truth_multi_genes": len(sn2), "n_pred_multi_tx": len(pm)}

    # L3 junctions, L4 exons (first/last: inner boundary only), L5 single-exon, L6 bases
    tj, pj, te, pe = set(), set(), set(), set()
    tsingle, psingle = [], []
    tb, pb = defaultdict(list), defaultdict(list)

    def exon_keys(st, ex):
        if len(ex) < 2:
            return set()
        keys = {(st, "L", ex[0][1]), (st, "R", ex[-1][0])}
        keys |= {(st, "I", s, e) for s, e in ex[1:-1]}
        return keys

    for g in tgenes.values():
        for ex in g["tx"].values():
            ex = _merge(ex)
            tj |= {(g["strand"],) + i for i in introns(ex)}
            te |= exon_keys(g["strand"], ex)
            if len(ex) == 1:
                tsingle.append((g["strand"],) + ex[0])
            tb[g["strand"]] += ex
    for p in cpreds:
        pj |= {(p["strand"],) + i for i in introns(p["exons"])}
        pe |= exon_keys(p["strand"], p["exons"])
        if len(p["exons"]) == 1:
            psingle.append((p["strand"],) + p["exons"][0])
        pb[p["strand"]] += p["exons"]
    res["L3_junction"] = _prf(len(tj & pj), len(pj), len(tj & pj), len(tj))
    res["L4_exon"] = _prf(len(te & pe), len(pe), len(te & pe), len(te))

    def recip(a, b):
        ov = min(a[2], b[2]) - max(a[1], b[1]) + 1
        return a[0] == b[0] and ov > 0 and ov >= 0.5 * (a[2] - a[1] + 1) and ov >= 0.5 * (b[2] - b[1] + 1)

    tsingle, psingle = sorted(set(tsingle)), sorted(set(psingle))
    t_hit = sum(any(recip(t, p) for p in psingle if p[1] <= t[2] and t[1] <= p[2]) for t in tsingle)
    p_hit = sum(any(recip(t, p) for t in tsingle if p[1] <= t[2] and t[1] <= p[2]) for p in psingle)
    res["L5_single_exon"] = _prf(p_hit, len(psingle), t_hit, len(tsingle))

    tbm, pbm = _bases(tb), _bases(pb)
    inter = sum(_overlap_len(tbm.get(s, []), pbm.get(s, [])) for s in set(tbm) | set(pbm))
    tlen = sum(e - s + 1 for v in tbm.values() for s, e in v)
    plen = sum(e - s + 1 for v in pbm.values() for s, e in v)
    res["L6_base"] = _prf(inter, plen, inter, tlen)
    return res


def delta_genes(new, old):
    """lncRNA genes in `new` with no same-strand exon overlap with any gene in `old`."""
    idx = ExonIndex(old)
    return {gid for gid, g in new.items() if g["cls"] == "lnc"
            and not any(idx.hits(g["strand"], _merge(ex)) for ex in g["tx"].values())}


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("score")
    s.add_argument("--truth", required=True, help="GENCODE GTF (T_new)")
    s.add_argument("--truth-format", choices=["gtf", "gff3"], default="gtf")
    s.add_argument("--ignore-strand", action="store_true",
                   help="strand-agnostic matching (for unstranded predictors, e.g. SegmentNT/NTv3)")
    s.add_argument("--old-truth", help="GENCODE GTF (T_old); adds old-truth scores and the delta set")
    s.add_argument("--chrom", required=True)
    s.add_argument("--pred", required=True)
    s.add_argument("--pred-format", choices=["gff3", "gtf", "smallrna"], default="gff3")
    s.add_argument("--rfam-types", help="TSV rfam_acc<TAB>type (Rfam family.txt cols 1 and 19)")
    s.add_argument("--classes", default="lnc,mrna", help="comma list of lnc,mrna,small")
    s.add_argument("--basic-only", action="store_true")
    s.add_argument("--boot", type=int, default=1000)
    s.add_argument("--label", default="")
    s.add_argument("--out", required=True)
    a = ap.parse_args()

    classes = a.classes.split(",")
    if a.pred_format == "gff3":
        preds = read_pred_gff3(a.pred, a.chrom)
    elif a.pred_format == "gtf":
        preds = preds_from_gtf(a.pred, a.chrom)
    else:
        types = None
        if a.rfam_types:
            types = dict(l.rstrip("\n").split("\t", 1) for l in open(a.rfam_types) if "\t" in l)
        preds = read_pred_smallrna(a.pred, a.chrom, types)
    truths = {"new": read_gtf(a.truth, a.chrom, a.basic_only) if a.truth_format == "gtf"
              else read_gff3_truth(a.truth, a.chrom)}
    if a.old_truth:
        truths["old"] = read_gtf(a.old_truth, a.chrom, a.basic_only)
    if a.ignore_strand:
        for t in truths.values():
            _unstrand_truth(t)
        for p in preds:
            p["strand"] = "."
    out = {"ruler": RULER_VERSION, "label": a.label, "chrom": a.chrom, "pred": a.pred,
           "truth": {"new": a.truth, "old": a.old_truth}, "basic_only": a.basic_only,
           "ignore_strand": a.ignore_strand,
           "n_pred_tx": {c: sum(p["cls"] == c for p in preds) for c in classes},
           "truth_genes_by_class": {c: sum(g["cls"] == c for g in truths["new"].values())
                                    for c in ("lnc", "mrna", "small", "pseudo", "other")},
           "scores": {}}
    for name, truth in truths.items():
        idx = ExonIndex(truth)
        out["scores"][name] = {c: score_class(c, truth, preds, idx, a.boot) for c in classes}
    if "old" in truths:
        new, old = truths["new"], truths["old"]
        dset = delta_genes(new, old)
        sub_truth = {gid: g for gid, g in new.items() if gid in dset or g["cls"] != "lnc"}
        idx = ExonIndex(sub_truth)
        out["delta"] = {"n_genes": len(dset),
                        "lnc": score_class("lnc", sub_truth, preds, idx, a.boot)}
        # incompleteness: lnc preds that are FP under old truth but TP under new truth
        io, inw = ExonIndex(old), ExonIndex(new)
        lp = [p for p in preds if p["cls"] == "lnc"]
        fp_old = [p for p in lp if not any(old[g]["cls"] == "lnc" for g in io.hits(p["strand"], p["exons"]))]
        rescued = [p for p in fp_old if any(new[g]["cls"] == "lnc" for g in inw.hits(p["strand"], p["exons"]))]
        out["fp_old_to_tp_new"] = {"n_fp_old": len(fp_old), "n_rescued": len(rescued),
                                   "frac": len(rescued) / len(fp_old) if fp_old else None}
    json.dump(out, open(a.out, "w"), indent=1)
    for c in classes:
        l1 = out["scores"]["new"][c]["L1_locus"]
        l2 = out["scores"]["new"][c]["L2_intron_chain"]
        print(f"[{a.label}] {c} L1 sn={l1['sn']} pr_strict={l1['pr_strict']} | L2 sn={l2['sn']} pr={l2['pr_strict']}",
              file=sys.stderr)


if __name__ == "__main__":
    main()
