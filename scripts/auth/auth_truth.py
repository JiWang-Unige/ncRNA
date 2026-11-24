#!/usr/bin/env python3
"""AUTH (docs/auth/PLAN_AUTH.md section 4): genome-wide WTC11 truth under the W24 rule, without the exclusion zone.

Rule = scripts/m1/w24/wtc11_truth.py (its functions are imported, the locus loop is copied unchanged): transcripts seen
in both PacBio rep2 and rep3 (TALON count >= 1, novelty Known / NIC / NNC / Antisense / Intergenic; ISM and 'Genomic'
dropped; novel single-exon dropped), class lnc / mRNA, outer bounds of the two replicates, clustered per class and strand
by exon overlap, locus span >= 200 bp.  AUTH uses no DNA model, so there is no exclusion zone ('excluded' is always False).

Before any candidate is generated the lnc locus counts are checked against section 4 of the plan (all 1,337; train 842;
test 312; secondary 183).  Only counts are written for chromosomes outside the training set; locus files are written for
the training chromosomes only (the test / secondary loci are built by auth_score.py --mode test after the freeze).
  auth_truth.py            -> outputs/AUTH/TRUTH_COUNTS.json, outputs/AUTH/truth/<train chrom>.json
"""
import json
import os
import sys
import time

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, f"{ROOT}/scripts/m1/w24")
import wtc11_truth as WT  # noqa: E402

OUT = f"{ROOT}/outputs/AUTH"
TRAIN = ("chr1", "chr3", "chr4", "chr5", "chr7", "chr8", "chr11", "chr14", "chr15", "chr17", "chr19", "chr21", "chrX", "chr18")
TEST = ("chr2", "chr6", "chr10", "chr16", "chr22")
SECONDARY = ("chr9", "chr12", "chr13", "chr20")
ALL = tuple(f"chr{i}" for i in range(1, 23)) + ("chrX", "chrY")
EXPECTED = {"all_autosomes_X": 1337, "train": 842, "test": 312, "secondary": 183}


def read_reps(chroms):
    WT.CHROMS = tuple(chroms)          # read_gtf keeps these chromosomes only
    reps = {}
    for r, (q, g) in WT.REPS.items():
        qq = WT.read_quant(f"{WT.E}/{q}")
        tx = WT.read_gtf(f"{WT.E}/{g}", set(qq))
        reps[r] = {}
        for tid, t in tx.items():
            k = WT.key_of(tid, t, qq[tid]["tnov"])
            if k is None:
                continue
            reps[r][k] = {**t, **qq[tid], "tid": tid}
    return reps


def build(reps, c):
    """Loci of one chromosome: the body of wtc11_truth.main() for one chromosome, with an empty exclusion zone."""
    both = sorted(set(reps["rep2"]) & set(reps["rep3"]))     # sorted: deterministic order for start ties
    trs = []
    for k in both:
        a, b = reps["rep2"][k], reps["rep3"][k]
        if b["chrom"] != c:
            continue
        cl = WT.klass(b, b["gnov"])
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
    out = []
    for i, L in enumerate(sorted(loci, key=lambda x: (x["start"], x["strand"]))):
        if L["end"] - L["start"] + 1 < 200:
            continue
        out.append({"id": f"{c}:w{i}", "k": L["k"], "strand": L["strand"], "start": L["start"], "end": L["end"],
                    "n_tx": len(L["tx"]), "count2": sum(t["count2"] for t in L["tx"]), "count3": sum(t["count3"] for t in L["tx"]),
                    "novel": all(t["gnov"] != "Known" for t in L["tx"]), "spliced": any(len(t["exons"]) > 1 for t in L["tx"]),
                    "excluded": False,
                    "tx": [{"exons": t["exons"], "nov": t["nov"], "id3": t["id3"], "count2": t["count2"], "count3": t["count3"]} for t in L["tx"]]})
    return out


def write(c, loci):
    os.makedirs(f"{OUT}/truth", exist_ok=True)
    json.dump({"chrom": c, "exclusion_zone": [], "loci": loci, "rule": "wtc11_truth.py, no exclusion zone (auth_truth.py)"},
              open(f"{OUT}/truth/{c}.json", "w"))


def main():
    t0 = time.time()
    reps = read_reps(ALL)
    n = {}
    for c in ALL:
        loci = build(reps, c)
        n[c] = {"lnc": sum(x["k"] == "l" for x in loci), "mrna": sum(x["k"] == "m" for x in loci)}
        if c in TRAIN:
            write(c, loci)
        del loci
    got = {"all_autosomes_X": sum(n[c]["lnc"] for c in ALL if c != "chrY"), "train": sum(n[c]["lnc"] for c in TRAIN),
           "test": sum(n[c]["lnc"] for c in TEST), "secondary": sum(n[c]["lnc"] for c in SECONDARY)}
    res = {"rule": "wtc11_truth.py without exclusion zone", "expected_plan_section4": EXPECTED, "got": got,
           "match": got == EXPECTED, "lnc_per_chrom": {c: n[c]["lnc"] for c in ALL}, "chrY_lnc_not_in_any_set": n["chrY"]["lnc"],
           "mrna_total": sum(v["mrna"] for v in n.values()), "seconds": round(time.time() - t0, 1)}
    json.dump(res, open(f"{OUT}/TRUTH_COUNTS.json", "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
