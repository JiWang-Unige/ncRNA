"""v1.9 development check (chr18 only): does checkpoint ranking by lnc AUPRC (exon-level MEA, gamma 1) agree between
the two halves of chr18?  half_split.py root name post out.json [gamma]  (genes / truth genes assigned to a half by start)."""
import json
import sys

import numpy as np

root, name, path, out = sys.argv[1:5]
GAMMA = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0
sys.path.insert(0, root + "/src")
from ncrna_annot.decode import chrom
from ncrna_annot.eval import ruler
from ncrna_annot.grammar import g1
from ncrna_annot.train import screen2 as s2

SPLIT = 40_000_000
post = np.load(path, mmap_mode="r")
codes = np.load(f"{root}/data/p1/chr18_codes.npy", mmap_mode="r")[:len(post)]
T = g1.G1Transitions()().detach().double().numpy()
truth = ruler.read_gtf(f"{root}/outputs/PHASE1/DIAG_L2/chr18_v50.gtf", "chr18")
genes = chrom.decode_blocks(post, codes, T, 20, 50, mode="mea", gamma=GAMMA, level="exon", log_gain=False)
tstart = {k: min(e[0] for t in v["tx"].values() for e in t) for k, v in truth.items()}
res = {"name": name, "gamma": GAMMA}
for h, keep in (("A", lambda s: s < SPLIT), ("B", lambda s: s >= SPLIT)):
    tr = {k: v for k, v in truth.items() if keep(tstart[k])}
    gs = [dict(g) for g in genes if keep(g["exons"][0][0])]
    r = s2.summarize(gs, post, tr, s2.CONF_GRID)
    res[h] = {"auprc": r["auprc_lnc"], "n": r["n_pred_lnc"], "n_truth": r["n_truth_lnc"],
              "lnc_L1": r["lnc_L1"], "mrna_L1": r["mrna_L1"], "conf_grid": r["conf_grid"]}
print(json.dumps({name: {h: {k: res[h][k] for k in ("auprc", "n", "n_truth", "lnc_L1")} for h in "AB"}}), flush=True)
json.dump(res, open(out, "w"), indent=1)
