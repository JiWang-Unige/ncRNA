"""chr18 diagnostic: is the posterior ensemble a faithful mean, and how do the members disagree?"""
import json, sys
import numpy as np
sys.path.insert(0, sys.argv[1] + "/src")
from ncrna_annot.eval import ruler
from ncrna_annot.grammar import g1
root, ens_p, mem = sys.argv[1], sys.argv[2], sys.argv[3:]
truth = ruler.read_gtf(f"{root}/outputs/PHASE1/DIAG_L2/chr18_v50.gtf", "chr18")
posts = {p.split("/")[-2]: np.load(p, mmap_mode="r") for p in mem}
posts["ENS"] = np.load(ens_p, mmap_mode="r")
L = len(posts["ENS"])
cls = np.zeros(L, np.int8)
for g in truth.values():
    ex = [e for t in g["tx"].values() for e in t]
    a, b = min(e[0] for e in ex), max(e[1] for e in ex)
    cls[max(a - 5001, 0):min(b + 5000, L)] = 3
for g in truth.values():
    if g["cls"] in ("lnc", "mrna"):
        for t in g["tx"].values():
            for a, b in t:
                cls[a - 1:b] = 1 if g["cls"] == "lnc" else 2
rng = np.random.default_rng(0)
ix = np.sort(rng.choice(L, 2_000_000, replace=False))
c = cls[ix]
P = {n: np.asarray(p[ix], np.float64) for n, p in posts.items()}
names = [n for n in P if n != "ENS" and "avg" not in n]
mean = sum(P[n] for n in names) / len(names)
grp = {"IG": [0], "lnc": [g1.IDX[f"l{q}{s}"] for s in g1.STRANDS for q in g1.PARTS],
       "mrna": [g1.IDX[f"m{q}{s}"] for s in g1.STRANDS for q in g1.PARTS]}
lE = [g1.IDX[f"lE{s}"] for s in g1.STRANDS]
def ap(score, y):
    o = np.argsort(-score); y = y[o]; tp = np.cumsum(y); return float((tp / np.arange(1, len(y) + 1))[y == 1].mean())
out = {"ens_minus_mean_maxabs": float(np.abs(P["ENS"] - mean).max()),
       "ens_minus_mean_p99": float(np.quantile(np.abs(P["ENS"] - mean), 0.99))}
am = {n: P[n].argmax(1) for n in P}
for n, p in P.items():
    s = p.sum(1)
    r = {"sum_q01_q99": [round(float(x), 4) for x in np.quantile(s, [0.01, 0.99])],
         "ap_lncE": round(ap(p[:, lE].sum(1), (c == 1).astype(int)), 4),
         "maxpost_median": round(float(np.median(p.max(1))), 4)}
    for cn, cv in (("farIG", 0), ("lncE", 1), ("mrnaE", 2)):
        m = c == cv
        r[cn] = {k: round(float(p[m][:, v].sum(1).mean()), 3) for k, v in grp.items()}
        r[cn]["argmax"] = {k: round(float(np.isin(am[n][m], v).mean()), 3) for k, v in grp.items()}
    out[n] = r
for i, a in enumerate(names):
    for b in names[i + 1:]:
        out[f"agree_{a}_{b}"] = round(float((am[a] == am[b]).mean()), 4)
# argmax segment counts on a contiguous 5 Mb window
a0 = 30_000_000
for n, p in posts.items():
    x = np.asarray(p[a0:a0 + 5_000_000], np.float32).argmax(1)
    out[n]["segments_5Mb"] = int((np.diff(x) != 0).sum())
print(json.dumps(out, indent=1))
