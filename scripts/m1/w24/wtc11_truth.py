#!/usr/bin/env python3
"""M1 weeks 2-4 (docs/m1/PLAN_M1_W24.md section 5): WTC11 validation truth from the LRGASP PacBio replicates.

Rules (written before any prediction on these chromosomes is read):
- Replicate roles: rep1 (ENCFF373TKM BAM) is the input of the full-transcriptome baselines; truth = transcripts observed in
  BOTH rep2 (ENCFF527KWD / ENCFF798DWV) and rep3 (ENCFF486EKE / ENCFF422CGK).
- Observed in a replicate: listed in its TALON quantification with count >= 1 and transcript_novelty in
  {Known, NIC, NNC, Antisense, Intergenic} (ISM fragments and 'Genomic' are dropped).
- Identity across replicates (TALON ran per replicate, novel IDs differ): Known -> GENCODE transcript id; other ->
  (chrom, strand, intron chain); novel single-exon transcripts are dropped (no reliable cross-replicate identity).
- Coordinates: outer bounds of the two replicates' models (exons from rep3's model, the deeper replicate).
- Class: lnc if the GENCODE V29 gene_type is a long non-coding type, or the gene is novel (gene_novelty Antisense /
  Intergenic); m if gene_type protein_coding; everything else is ignored.
- Loci: truth transcripts of one class clustered by same-strand exon overlap; locus span >= 200 bp.
- Exclusion zone (low-cost option, author 2026-11-17): every v2.0 training window (data/p1v20/train_windows.tsv) +- 24,576 bp,
  merged; a locus is evaluable only if its whole span lies outside.
  wtc11_truth.py -> outputs/M1/W24/truth/<chrom>.json, outputs/M1/W24/TRUTH_STATS.json
"""
import gzip
import json
import os
import re
from collections import Counter, defaultdict

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
E = f"{ROOT}/data/wtc11/encode"
OUT = f"{ROOT}/outputs/M1/W24"
REPS = {"rep2": ("ENCFF527KWD.tsv", "ENCFF798DWV.gtf.gz"), "rep3": ("ENCFF486EKE.tsv", "ENCFF422CGK.gtf.gz")}
CHROMS = ("chr18", "chr9", "chr12", "chr20", "chr13")   # dev chr18; test chr9 + chr12 + chr20 + chr13 (author, option A)
TEST = ("chr9", "chr12", "chr20", "chr13")
KEEP = {"Known", "NIC", "NNC", "Antisense", "Intergenic"}
LNC_TYPES = {"lincRNA", "antisense", "processed_transcript", "sense_intronic", "sense_overlapping", "3prime_overlapping_ncRNA",
             "bidirectional_promoter_lncRNA", "macro_lncRNA", "non_coding", "lncRNA"}
WIN, PAD = 24_576, 24_576


def read_quant(path):
    out = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = dict(zip(hdr, line.rstrip("\n").split("\t")))
            if f["transcript_novelty"] in KEEP and int(f[hdr[-1]]) >= 1:
                out[f["annot_transcript_id"]] = {"count": int(f[hdr[-1]]), "tnov": f["transcript_novelty"], "gnov": f["gene_novelty"]}
    return out


def read_gtf(path, want):
    tx = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if f[0] not in CHROMS or f[2] not in ("transcript", "exon"):
                continue
            tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            if tid not in want:
                continue
            if f[2] == "transcript":
                gt = re.search(r'gene_type "([^"]+)"', f[8])
                tx.setdefault(tid, {"exons": []}).update(chrom=f[0], strand=f[6], gene_type=gt.group(1) if gt else None,
                                                          gene_id=re.search(r'gene_id "([^"]+)"', f[8]).group(1))
            else:
                tx.setdefault(tid, {"exons": []})["exons"].append((int(f[3]), int(f[4])))
    for t in tx.values():
        t["exons"].sort()
    return tx


def key_of(tid, t, tnov):
    if tnov == "Known":
        return ("K", tid)
    if len(t["exons"]) < 2:
        return None
    introns = tuple((a[1], b[0]) for a, b in zip(t["exons"], t["exons"][1:]))
    return ("N", t["chrom"], t["strand"], introns)


def klass(t, gnov):
    if gnov in ("Antisense", "Intergenic"):
        return "l"
    if t["gene_type"] == "protein_coding":
        return "m"
    if t["gene_type"] in LNC_TYPES:
        return "l"
    return None


def exclusion():
    iv = defaultdict(list)
    for line in open(f"{ROOT}/data/p1v20/train_windows.tsv"):
        c, s, _ = line.rstrip("\n").split("\t")
        s = int(s)
        iv[c].append((max(0, s - PAD), s + WIN + PAD))
    out = {}
    for c, v in iv.items():
        m = []
        for a, b in sorted(v):
            if m and a <= m[-1][1]:
                m[-1][1] = max(m[-1][1], b)
            else:
                m.append([a, b])
        out[c] = m
    return out


def main():
    os.makedirs(f"{OUT}/truth", exist_ok=True)
    reps = {}
    for r, (q, g) in REPS.items():
        qq = read_quant(f"{E}/{q}")
        tx = read_gtf(f"{E}/{g}", set(qq))
        reps[r] = {}
        for tid, t in tx.items():
            k = key_of(tid, t, qq[tid]["tnov"])
            if k is None:
                continue
            reps[r][k] = {**t, **qq[tid], "tid": tid}
    both = set(reps["rep2"]) & set(reps["rep3"])
    ex = exclusion()
    stats = {"rules": __doc__.split("Rules")[1].split("  wtc11_truth.py")[0].strip(),
             "observed": {r: len(v) for r, v in reps.items()}, "in_both": len(both), "chroms": {}}
    for c in CHROMS:
        trs = []
        for k in both:
            a, b = reps["rep2"][k], reps["rep3"][k]
            if b["chrom"] != c:
                continue
            cl = klass(b, b["gnov"])
            if cl is None:
                continue
            ex_ = b["exons"]
            s, e = min(a["exons"][0][0], ex_[0][0]), max(a["exons"][-1][1], ex_[-1][1])
            ex_ = [(s, ex_[0][1])] + ex_[1:] if len(ex_) > 1 else [(s, e)]
            if len(ex_) > 1:
                ex_[-1] = (ex_[-1][0], e)
            trs.append({"k": cl, "strand": b["strand"], "exons": ex_, "start": s, "end": e, "nov": b["tnov"],
                        "gnov": b["gnov"], "count2": a["count"], "count3": b["count"], "id3": b["tid"], "gene_type": b["gene_type"]})
        loci = []
        for cl in ("l", "m"):
            for st in "+-":
                ts = sorted([t for t in trs if t["k"] == cl and t["strand"] == st], key=lambda t: t["start"])
                cur = None
                for t in ts:
                    if cur and t["start"] <= cur["end"] and any(x0 <= y1 and y0 <= x1 for x0, x1 in t["exons"] for y0, y1 in cur["exons_all"]):
                        cur["tx"].append(t)
                        cur["end"] = max(cur["end"], t["end"])
                        cur["exons_all"] += t["exons"]
                    else:
                        cur = {"k": cl, "strand": st, "start": t["start"], "end": t["end"], "tx": [t], "exons_all": list(t["exons"])}
                        loci.append(cur)
        zone = ex.get(c, [])
        out = []
        for i, L in enumerate(sorted(loci, key=lambda x: (x["start"], x["strand"]))):
            if L["end"] - L["start"] + 1 < 200:
                continue
            inside = any(a < L["end"] and L["start"] - 1 < b for a, b in zone)   # 1-based closed span vs 0-based half-open zone
            out.append({"id": f"{c}:w{i}", "k": L["k"], "strand": L["strand"], "start": L["start"], "end": L["end"],
                        "n_tx": len(L["tx"]), "count2": sum(t["count2"] for t in L["tx"]), "count3": sum(t["count3"] for t in L["tx"]),
                        "novel": all(t["gnov"] != "Known" for t in L["tx"]), "spliced": any(len(t["exons"]) > 1 for t in L["tx"]),
                        "excluded": inside,
                        "tx": [{"exons": t["exons"], "nov": t["nov"], "id3": t["id3"], "count2": t["count2"], "count3": t["count3"]} for t in L["tx"]]})
        json.dump({"chrom": c, "exclusion_zone": zone, "loci": out}, open(f"{OUT}/truth/{c}.json", "w"))
        zl = sum(b - a for a, b in zone)
        cs = {}
        for cl in ("l", "m"):
            sel = [x for x in out if x["k"] == cl]
            ok = [x for x in sel if not x["excluded"]]
            cs[cl] = {"loci": len(sel), "evaluable": len(ok), "excluded": len(sel) - len(ok),
                      "evaluable_novel": sum(x["novel"] for x in ok), "evaluable_spliced": sum(x["spliced"] for x in ok),
                      "median_count_eval": sorted(x["count2"] + x["count3"] for x in ok)[len(ok) // 2] if ok else None,
                      "median_count_excl": sorted(x["count2"] + x["count3"] for x in sel if x["excluded"])[
                          len([x for x in sel if x["excluded"]]) // 2] if len(sel) > len(ok) else None,
                      "median_len_eval": sorted(x["end"] - x["start"] + 1 for x in ok)[len(ok) // 2] if ok else None}
        stats["chroms"][c] = {"zone_bp": zl, "zone_intervals": len(zone), "classes": cs,
                              "tx_by_novelty_lnc": dict(Counter(t["nov"] for t in trs if t["k"] == "l"))}
    n_test = sum(stats["chroms"][c]["classes"]["l"]["evaluable"] for c in TEST)
    stats["power_rule"] = "test chromosomes (chr9 + chr12 + chr20 + chr13) evaluable lnc loci >= 150"
    stats["test_evaluable_lnc"] = n_test
    stats["power_ok"] = n_test >= 150
    json.dump(stats, open(f"{OUT}/TRUTH_STATS.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in stats.items() if k != "rules"}, indent=1))


if __name__ == "__main__":
    main()
