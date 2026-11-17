"""Truth-only census: lnc loci per chromosome under the frozen-draft truth rule (both rep2 & rep3) and two looser rules,
inside / outside the v2.0 training-window exclusion zone.  No prediction is read."""
import json, os, sys
sys.path.insert(0, os.path.expanduser("~/ncrna-annot/scripts/m1/w24"))
import wtc11_truth as W
CH = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
W.CHROMS = tuple(CH)
reps = {}
for r, (q, g) in {"rep1": ("ENCFF803XYM.tsv", "ENCFF370SDY.gtf.gz"), **W.REPS}.items():
    qq = W.read_quant(f"{W.E}/{q}")
    tx = W.read_gtf(f"{W.E}/{g}", set(qq))
    reps[r] = {}
    for tid, t in tx.items():
        k = W.key_of(tid, t, qq[tid]["tnov"])
        if k is not None:
            reps[r][k] = {**t, **qq[tid]}
ex = W.exclusion()
rules = {"both_rep2_rep3": lambda k: k in reps["rep2"] and k in reps["rep3"],
         "either_rep2_rep3_ge2reads": lambda k: (reps["rep2"].get(k, {}).get("count", 0) + reps["rep3"].get(k, {}).get("count", 0)) >= 2,
         "any_rep_ge2reads": lambda k: sum(reps[r].get(k, {}).get("count", 0) for r in reps) >= 2}
allk = set().union(*[set(v) for v in reps.values()])
out = {}
for name, f in rules.items():
    per = {}
    for c in CH:
        loci = []
        ts = []
        for k in allk:
            if not f(k):
                continue
            t = next(reps[r][k] for r in ("rep3", "rep2", "rep1") if k in reps[r])
            if t["chrom"] != c or W.klass(t, t["gnov"]) != "l":
                continue
            ts.append(t)
        for st in "+-":
            cur = None
            for t in sorted([t for t in ts if t["strand"] == st], key=lambda t: t["exons"][0][0]):
                s, e = t["exons"][0][0], t["exons"][-1][1]
                if cur and s <= cur[1] and any(a <= y1 and y0 <= b for a, b in t["exons"] for y0, y1 in cur[2]):
                    cur[1] = max(cur[1], e); cur[2] += t["exons"]
                else:
                    cur = [s, e, list(t["exons"])]; loci.append(cur)
        loci = [L for L in loci if L[1] - L[0] + 1 >= 200]
        zone = ex.get(c, [])
        ok = [L for L in loci if not any(a < L[1] and L[0] - 1 < b for a, b in zone)]
        per[c] = [len(loci), len(ok)]
    out[name] = {"per_chrom_[all,outside_zone]": per, "total": [sum(v[0] for v in per.values()), sum(v[1] for v in per.values())]}
json.dump(out, open(os.path.expanduser("~/ncrna-annot/outputs/M1/W24/TRUTH_CENSUS.json"), "w"), indent=1)
for name, v in out.items():
    print(name, "total", v["total"])
    print("  ", " ".join(f"{c}:{a}/{b}" for c, (a, b) in v["per_chrom_[all,outside_zone]"].items()))
