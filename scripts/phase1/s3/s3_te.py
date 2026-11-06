#!/usr/bin/env python3
"""S3-TE (PLAN_S3_TE.md): detection of TE-derived lncRNA loci, stratified, from existing reads only (descriptive).

TE = UCSC hg38 rmsk repClass in {LINE, SINE, LTR, DNA, Retroposon, RC} (no '?' classes), strand ignored.
Per v50 lnc gene: f = fraction of its merged exon bases (all transcripts) inside a TE; TSS-in-TE = the TSS base of
its Ensembl_canonical transcript (longest transcript if none is tagged) inside a TE.
Strata: f < 0.1 / 0.1 <= f < 0.5 / f >= 0.5; TSS in TE yes / no; "TE-derived" = f >= 0.5 or TSS in TE.
Reads: chr20 = v1.9 main runs (configs/p15_v19.json, tau 0.6); chr9 = v2.0 runs (configs/v20.json, calibrated
probability >= tau_v20, as v20_score.py).  Sn per stratum = three-run mean of the fraction of the stratum's genes hit
(strict L1, same-strand exon overlap).  CI: 2,000 resamples of the truth genes (numpy default_rng(0)), strata kept.

  s3_te.py --rmsk data/te/rmsk.txt.gz --out outputs/PHASE1/S3TE/S3TE.json
"""
import argparse
import gzip
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1/v20"))
import p15_score as P  # noqa: E402
from ncrna_annot.eval import ruler  # noqa: E402

TE_CLASSES = {"LINE", "SINE", "LTR", "DNA", "Retroposon", "RC"}
V50 = "data/human/gencode.v50.annotation.gtf.gz"
N_BOOT = 2000


def te_mask(path, chrom, L):
    m = np.zeros(L, bool)
    n = 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if f[5] != chrom or f[11] not in TE_CLASSES:
                continue
            m[int(f[6]):int(f[7])] = True   # rmsk genoStart is 0-based, genoEnd exclusive
            n += 1
    return m, n


def canonical(chrom):
    out = {}
    with gzip.open(f"{ROOT}/{V50}", "rt") as fh:
        for line in fh:
            if not line.startswith(chrom + "\t"):
                continue
            f = line.split("\t")
            if f[2] == "transcript" and 'tag "Ensembl_canonical"' in f[8]:
                out[f[8].split('gene_id "')[1].split('"')[0]] = f[8].split('transcript_id "')[1].split('"')[0]
    return out


def gene_table(truth, mask, canon):
    rows = {}
    for gid, g in truth.items():
        if g["cls"] != "lnc":
            continue
        merged = ruler._merge(sorted(e for ex in g["tx"].values() for e in ex))
        bases = np.concatenate([np.arange(a - 1, b) for a, b in merged])
        tid = canon.get(gid) if canon.get(gid) in g["tx"] else max(g["tx"], key=lambda t: sum(b - a + 1 for a, b in g["tx"][t]))
        ex = sorted(g["tx"][tid])
        tss = ex[0][0] - 1 if g["strand"] == "+" else ex[-1][1] - 1
        f = float(mask[bases].mean())
        rows[gid] = {"f": f, "tss_te": bool(mask[tss]), "len": int(len(bases)), "n_exon": len(ex),
                     "stratum": "lt0.1" if f < 0.1 else "0.1-0.5" if f < 0.5 else "ge0.5"}
        rows[gid]["te_derived"] = rows[gid]["f"] >= 0.5 or rows[gid]["tss_te"]
    return rows


def hits_chr20():
    cfg = json.load(open(f"{ROOT}/configs/p15_v19.json"))
    tau = cfg["decode"]["tau"]
    runs = sorted(cfg["main"]["runs"].values())
    out = {}
    for run in runs:
        genes = P.load_genes(f"{ROOT}/{cfg['out_dir']}/chr20/{run}/genes.json")
        out[run] = [p for p in P.final_preds(genes, tau) if p["cls"] == "lnc"]
    return out, {"config": "configs/p15_v19.json", "tau": tau}


def hits_chr9():
    import v20_calib as C
    cfg = json.load(open(f"{ROOT}/configs/v20.json"))
    cal = cfg["calibration"]
    tab = json.load(open(f"{ROOT}/{cal['hexamer_table']}"))
    codes = np.load(f"{ROOT}/{cfg['targets']['chr9']['codes']}", mmap_mode="r")
    out = {}
    for run in cfg["targets"]["chr9"]["runs"]:
        rd = f"{ROOT}/{cfg['out_dir']}/chr9/{run}"
        genes = P.load_genes(f"{rd}/genes.json")
        rows = C.gene_rows(genes, np.load(f"{rd}/post.npy", mmap_mode="r"), codes, tab)
        X = np.array([[r[f] for f in cal["model"]["features"]] for r in rows], float)
        p = C.apply(cal["model"], X) if cal["keep"] else X[:, cal["model"]["features"].index("conf")]
        sc = {r["k"]: float(s) for r, s in zip(rows, p)}
        gv = [dict(g, conf=sc.get(k, g["conf"])) for k, g in enumerate(genes)]
        out[run] = [q for q in P.final_preds(gv, cal["tau_v20"]) if q["cls"] == "lnc"]
    return out, {"config": "configs/v20.json", "tau_v20": cal["tau_v20"], "calibration_kept": cal["keep"]}


def analyse(chrom, preds_by_run, rmsk):
    truth = ruler.read_gtf(f"{ROOT}/{V50}", chrom)
    L = len(np.load(f"{ROOT}/data/p1/{chrom}_codes.npy", mmap_mode="r"))
    mask, n_te = te_mask(rmsk, chrom, L)
    tab = gene_table(truth, mask, canonical(chrom))
    gids = sorted(tab)
    lidx = ruler.ExonIndex({k: truth[k] for k in gids})
    H = []
    for run, preds in sorted(preds_by_run.items()):
        got = set()
        for p in preds:
            got |= set(lidx.hits(p["strand"], p["exons"]))
        H.append(np.array([g in got for g in gids], float))
    H = np.mean(H, 0)                                   # three-run mean hit per truth gene
    strata = {"lt0.1": np.array([tab[g]["stratum"] == "lt0.1" for g in gids]),
              "0.1-0.5": np.array([tab[g]["stratum"] == "0.1-0.5" for g in gids]),
              "ge0.5": np.array([tab[g]["stratum"] == "ge0.5" for g in gids]),
              "tss_te": np.array([tab[g]["tss_te"] for g in gids]),
              "tss_not_te": np.array([not tab[g]["tss_te"] for g in gids])}
    ted = np.array([tab[g]["te_derived"] for g in gids])
    rng = np.random.default_rng(0)
    boot = {k: [] for k in list(strata) + ["d_f", "d_tss"]}
    for _ in range(N_BOOT):
        i = rng.integers(0, len(gids), len(gids))
        sn = {k: H[i][m[i]].mean() if m[i].any() else np.nan for k, m in strata.items()}
        for k in strata:
            boot[k].append(sn[k])
        boot["d_f"].append(sn["ge0.5"] - sn["lt0.1"])
        boot["d_tss"].append(sn["tss_te"] - sn["tss_not_te"])
    ci = lambda x: np.nanpercentile(x, [2.5, 97.5]).tolist()
    res = {"n_te_intervals": n_te, "te_bp_frac": float(mask.mean()), "n_truth_lnc": len(gids), "strata": {}}
    for k, m in strata.items():
        res["strata"][k] = {"n": int(m.sum()), "sn": float(H[m].mean()) if m.any() else None, "sn_ci95": ci(boot[k]),
                            "miss": float(1 - H[m].mean()) if m.any() else None,
                            "median_len": float(np.median([tab[g]["len"] for g, x in zip(gids, m) if x])) if m.any() else None,
                            "median_n_exon": float(np.median([tab[g]["n_exon"] for g, x in zip(gids, m) if x])) if m.any() else None}
    res["delta"] = {"ge0.5_minus_lt0.1": {"sn": res["strata"]["ge0.5"]["sn"] - res["strata"]["lt0.1"]["sn"], "ci95": ci(boot["d_f"])},
                    "tss_te_minus_not": {"sn": res["strata"]["tss_te"]["sn"] - res["strata"]["tss_not_te"]["sn"], "ci95": ci(boot["d_tss"])}}
    miss_w = 1 - H                                       # expected missed fraction per gene (three-run mean)
    share_all = float(ted.mean())
    share_missed = float((miss_w * ted).sum() / miss_w.sum())
    res["te_derived"] = {"n": int(ted.sum()), "share_all": share_all, "share_missed": share_missed,
                         "enrichment_missed": share_missed / share_all if share_all else None}
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rmsk", default="data/te/rmsk.txt.gz")
    ap.add_argument("--out", default="outputs/PHASE1/S3TE/S3TE.json")
    a = ap.parse_args()
    t0 = time.time()
    rmsk = f"{ROOT}/{a.rmsk}"
    out = {"plan": "docs/phase1_results/PLAN_S3_TE.md", "rmsk_sha256": P.sha256(rmsk), "te_classes": sorted(TE_CLASSES), "chroms": {}}
    for chrom, fn in (("chr20", hits_chr20), ("chr9", hits_chr9)):
        preds, meta = fn()
        out["chroms"][chrom] = {"reads": meta, "runs": sorted(preds), **analyse(chrom, preds, rmsk)}
        print(chrom, json.dumps({k: v for k, v in out["chroms"][chrom].items() if k in ("strata", "delta", "te_derived")})[:1500], flush=True)
    d = {c: out["chroms"][c]["delta"] for c in out["chroms"]}
    both = lambda key: all(d[c][key]["sn"] < 0 and d[c][key]["ci95"][1] < 0 for c in d)
    one = lambda key: any(d[c][key]["sn"] < 0 and d[c][key]["ci95"][1] < 0 for c in d)
    out["verdict"] = {k: ("TE-derived lnc harder" if both(k) else "inconsistent" if one(k) else "no evidence")
                      for k in ("ge0.5_minus_lt0.1", "tss_te_minus_not")}
    out["seconds"] = round(time.time() - t0)
    os.makedirs(os.path.dirname(f"{ROOT}/{a.out}"), exist_ok=True)
    json.dump(out, open(f"{ROOT}/{a.out}", "w"), indent=1)
    print(json.dumps(out["verdict"]), flush=True)


if __name__ == "__main__":
    main()
