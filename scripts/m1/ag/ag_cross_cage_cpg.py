"""Post hoc, descriptive (not in either plan): first-stage Loc1k on chr8+chr9 lnc starts split by FANTOM5 CAGE evidence
x low CpG (Q1+Q2), for AlphaGenome (mean over CAGE tracks), Puffin and the CpG control."""
import json, os, sys
import numpy as np
ROOT = os.path.expanduser("~/ncrna-annot")
for d in ("src", "scripts/phase1", "scripts/phase1/s3", "scripts/phase1/ends", "scripts/phase1/dedend", "scripts/m1/ag"):
    sys.path.insert(0, os.path.join(ROOT, d))
import dedend_metrics as DM
import s3_features as S3
import ag_fold_metrics as AM
O = f"{ROOT}/outputs/M1/AG"
rows = {r["id"]: r for r in json.load(open(f"{O}/tss.json"))}
for m in ("fold-0", "fold-1", "fold-2", "fold-3"):
    z = np.load(f"{O}/colab_auto/raw/{m}.npz")
    for j, rid in enumerate(z["rid"]):
        rows[str(rid)]["ag"] = float(abs(int(z["agg_cage_mean_arg"][j])) <= 1000)
lnc = [r for r in rows.values() if r["chrom"] in ("chr8", "chr9") and r["k"] == "l"]
q2 = float(np.quantile([r["cpg"] for r in lnc], 0.5))
rd = [r for r in lnc if "ag" in r]
for c in ("chr8", "chr9"):
    seq = S3.read_fasta(f"{ROOT}/data/human/{c}.fa")
    tr = np.load(f"{ROOT}/outputs/PHASE1/DEDEND/puffin/{c}.npy", mmap_mode="r")
    for r in rd:
        if r["chrom"] == c:
            r["cpg_ctrl"] = float(AM.cpg_dist(seq, r) <= 1000)
            r["puffin"] = float(DM.puffin_dist(tr, r) <= 1000)
out = {}
for ev in (True, False):
    for low in (True, False):
        ss = [r for r in rd if bool(r["TSS_ev"]) == ev and (r["cpg"] <= q2) == low]
        out[f"cage={ev}|lowcpg={low}"] = {"n": len(ss), **{k: round(float(np.mean([r[k] for r in ss])), 3) for k in ("ag", "puffin", "cpg_ctrl")}}
print(json.dumps(out, indent=1))
json.dump(out, open(f"{O}/stageB/CROSS_cage_x_cpg.json", "w"), indent=1)
