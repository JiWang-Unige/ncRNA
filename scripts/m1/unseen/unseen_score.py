#!/usr/bin/env python3
"""M1 audit option A (docs/m1/PLAN_M1_UNSEEN.md): score giraffe chr14 (parts a / b / c, pooled).

Ours: v2.0 runs read by p16_read.py under configs/m1_unseen.json, confidences replaced by the configs/v20.json
calibration (tau_v20), exactly as simple_det.load_ours does for chr9.  Detector: the frozen experiment-2 RULES.json
(sha256 checked against outputs/M1/SIMPLE/FROZEN) through simple_det.py's own functions; only the sequence loader and
the Tiberius directory are pointed at the giraffe parts.  Truth: RefSeq lncRNA genes lying entirely inside one part.
  unseen_score.py   -> outputs/M1/UNSEEN/SCORE.json
"""
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "scripts/m1"))
import simple_det as D  # noqa: E402  (adds src / phase1 paths)
import p15_score as P  # noqa: E402
import v20_calib as Cc  # noqa: E402
from ncrna_annot.eval import ruler  # noqa: E402

O = f"{ROOT}/outputs/M1/UNSEEN"
CFG = f"{ROOT}/configs/m1_unseen.json"
PARTS = ("giraffe_14a", "giraffe_14b", "giraffe_14c")
REF_LOWCPG_S50 = (0.024458420684835773 + 0.027153558052434454) / 2   # chr20, chr9 (RESULT_M1_SIMPLE.md)

_SEQ = {}


def giraffe_seq():
    if "s" not in _SEQ:
        with open(f"{ROOT}/data/giraffe/NC_137288.1.fa", "rb") as fh:
            fh.readline()
            _SEQ["s"] = b"".join(line.strip() for line in fh)
    return _SEQ["s"]


def load_part(name):
    """simple_det.load_seq replacement: the part's slice of the giraffe chromosome."""
    s0, e0 = PJ["parts"][name]
    raw = giraffe_seq()[s0:e0]
    arr = np.frombuffer(raw, np.uint8)
    lower = (arr >= 97) & (arr <= 122)
    lut = np.full(256, 4, np.uint8)
    for i, c in enumerate(b"ACGT"):
        lut[c] = i
        lut[c + 32] = i
    return raw.upper(), lut[arr], lower


def split_tiberius():
    """Tiberius GTF of the whole chromosome -> one GTF per part (rows starting in the part, shifted, clipped)."""
    rows = {p: [] for p in PARTS}
    for line in open(f"{O}/tiberius/NC_137288.1.gtf"):
        f = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(f) < 9:
            continue
        s, e = int(f[3]), int(f[4])
        for p in PARTS:
            a, b = PJ["parts"][p]
            if a < s <= b:
                f2 = list(f)
                f2[0], f2[3], f2[4] = p, str(s - a), str(min(e, b) - a)
                rows[p].append("\t".join(f2))
    for p in PARTS:
        open(f"{O}/tiberius/{p}.gtf", "w").write("\n".join(rows[p]) + "\n")
    return {p: len(v) for p, v in rows.items()}


def part_truth(full, name):
    a, b = PJ["parts"][name]
    out, crossing = {}, 0
    for gid, g in full.items():
        s = min(x[0] for t in g["tx"].values() for x in t)
        e = max(x[1] for t in g["tx"].values() for x in t)
        if a < s and e <= b:
            out[gid] = dict(g, tx={t: [(x - a, y - a) for x, y in ex] for t, ex in g["tx"].items()})
        elif a < s <= b < e:   # crosses this part's right cut (each crossing gene counted once)
            crossing += g["cls"] == "lnc"
    return out, crossing


def cpg_oe(seq, gene):
    s0 = min(x[0] for t in gene["tx"].values() for x in t)
    e0 = max(x[1] for t in gene["tx"].values() for x in t)
    t = s0 - 1 if gene["strand"] == "+" else e0 - 1
    w = seq[max(t - 500, 0):t + 501].upper()
    c, g, cg = w.count(b"C"), w.count(b"G"), w.count(b"CG")
    return cg * len(w) / (c * g) if c and g else 0.0


def main():
    t0 = time.time()
    frozen = open(f"{ROOT}/outputs/M1/SIMPLE/FROZEN").read().split()[0]
    rules_path = f"{ROOT}/outputs/M1/SIMPLE/RULES.json"
    assert D.sha256(rules_path) == frozen, "RULES.json differs from the frozen experiment-2 rules"
    rules = json.load(open(rules_path))
    cfg = json.load(open(CFG))
    cal, tau = cfg["calibration"], cfg["calibration"]["tau_v20"]
    assert cal["keep"] and P.sha256(f"{ROOT}/{cal['hexamer_table']}") == cal["hexamer_table_sha256"]
    tab = json.load(open(f"{ROOT}/{cal['hexamer_table']}"))
    D.load_seq = load_part
    D.TIB = f"{O}/tiberius"
    tib_rows = split_tiberius()
    full = ruler.read_gff3_truth(f"{ROOT}/data/giraffe/chr14.gff3", "NC_137288.1")
    runs = cfg["main"]["runs"].values()
    per_part, gids, vec_o, vec_d, cpg = {}, [], {"S50": [], "L1s": [], "L1e": []}, {"S50": [], "L1s": []}, []
    pooled_counts = {"ours": {r: {"n_pred": 0, "S50": 0, "L1s": 0, "L1e": 0, "L1s_hitpred": 0, "L1e_hitpred": 0} for r in runs},
                     "det": {m: {"n_pred": 0, "S50": 0, "L1s": 0, "L1s_hitpred": 0} for m in ("S50", "L1s")}}
    n_cross = 0
    for name in PARTS:
        truth, cross = part_truth(full, name)
        n_cross += cross
        ts, idx = D.truth_spans(truth), ruler.ExonIndex(truth)
        pg = [t[0] for t in ts]
        seq = giraffe_seq()[PJ["parts"][name][0]:PJ["parts"][name][1]]
        cpg += [cpg_oe(seq, truth[g]) for g in pg]
        codes = np.load(f"{ROOT}/data/p1/{name}_codes.npy", mmap_mode="r")
        hits_o = {"S50": [], "L1s": [], "L1e": []}
        res_o = {}
        for run in runs:
            rd = f"{ROOT}/{cfg['out_dir']}/{name}/{run}"
            genes = P.load_genes(f"{rd}/genes.json")
            post = np.load(f"{rd}/post.npy", mmap_mode="r")
            rows = Cc.gene_rows(genes, post, codes, tab)
            X = np.array([[r[f] for f in cal["model"]["features"]] for r in rows], float)
            sc = dict(zip((r["k"] for r in rows), Cc.apply(cal["model"], X)))
            fin = [{"id": f"o{k}", "cls": "lnc", "strand": g["strand"], "exons": [tuple(x) for x in g["exons"]]}
                   for k, g in enumerate(genes) if g["type"] == "l" and sc.get(k, -1) >= tau]
            mt, m50, h1 = D.metrics(truth, ts, idx, [D.as_span(g) for g in fin])
            he, ne = D.exon_l1(idx, truth, fin)
            res_o[run] = dict(mt, L1_exon_sn=len(he) / len(ts) if ts else None)
            hits_o["S50"].append(np.array([g in m50 for g in pg], float))
            hits_o["L1s"].append(np.array([g in h1 for g in pg], float))
            hits_o["L1e"].append(np.array([g in he for g in pg], float))
            pc = pooled_counts["ours"][run]
            pc["n_pred"] += len(fin)
            pc["S50"] += len(m50)
            pc["L1s"] += len(h1)
            pc["L1e"] += len(he)
            pc["L1s_hitpred"] += round((mt["L1s_pr"] or 0) * len(fin))
            pc["L1e_hitpred"] += ne
        n_match = int(round(np.mean([r["n_pred"] for r in res_o.values()])))
        F = D.features(name)
        tib = D.read_tiberius(name)
        res_d = {}
        for metric in ("S50", "L1s"):
            c = rules["selected"][metric]["config"]
            s = D.score_lr(F, rules["lr_model"]) if c["variant"] == "lr" else D.score_cgi(F)
            det = D.predict(s, F["K"], F["L"], c, n=n_match, tib=tib)
            dm, d50, d1 = D.metrics(truth, ts, idx, det)
            res_d[metric] = dm
            vec_d_m = {"S50": np.array([g in d50 for g in pg], float), "L1s": np.array([g in d1 for g in pg], float)}
            vec_d.setdefault(f"{metric}_cfg", []).append(vec_d_m)
            pd_ = pooled_counts["det"][metric]
            pd_["n_pred"] += len(det)
            pd_["S50"] += len(d50)
            pd_["L1s"] += len(d1)
            pd_["L1s_hitpred"] += round((dm["L1s_pr"] or 0) * len(det))
        for k in hits_o:
            vec_o[k].append(np.array(hits_o[k]))
        gids += pg
        per_part[name] = {"n_truth_lnc": len(ts), "n_lnc_crossing_excluded": cross, "n_matched": n_match,
                          "ours": res_o, "det": res_d, "tiberius_rows": tib_rows[name]}
        print(name, json.dumps(per_part[name]), flush=True)
    # pooled per-gene vectors: ours (runs x genes), detector per metric config
    o = {k: np.concatenate(v, axis=1) for k, v in vec_o.items()}
    d = {m: {k: np.concatenate([p[k] for p in vec_d[f"{m}_cfg"]]) for k in ("S50", "L1s")} for m in ("S50", "L1s")}
    cpg = np.array(cpg)
    qs = np.quantile(cpg, [0.25, 0.5, 0.75])
    q = np.searchsorted(qs, cpg, side="right")
    low = q <= 1
    nt = len(gids)
    pooled = {"n_truth_lnc": nt, "n_lnc_crossing_excluded": n_cross,
              "ours_mean": {k: float(o[k].mean()) for k in o},
              "ours_runs": {r: {k: v for k, v in pc.items()} for r, pc in pooled_counts["ours"].items()},
              "det": {m: {"S50_sn": float(d[m]["S50"].mean()), "L1s_sn": float(d[m]["L1s"].mean()),
                          "counts": pooled_counts["det"][m]} for m in ("S50", "L1s")}}
    cmp_ = {"S50": D.boot(d["S50"]["S50"], list(o["S50"])), "L1s": D.boot(d["L1s"]["L1s"], list(o["L1s"])),
            "L1s_det_vs_ours_exonL1": D.boot(d["L1s"]["L1s"], list(o["L1e"]))}
    lowc = {"S50": D.boot(d["S50"]["S50"][low], [v[low] for v in o["S50"]]),
            "L1s": D.boot(d["L1s"]["L1s"][low], [v[low] for v in o["L1s"]])}
    strata = {f"Q{i + 1}": {"n": int((q == i).sum()),
                            "S50_det": float(d["S50"]["S50"][q == i].mean()), "S50_ours": float(o["S50"].mean(0)[q == i].mean()),
                            "L1s_det": float(d["L1s"]["L1s"][q == i].mean()), "L1s_ours": float(o["L1s"].mean(0)[q == i].mean())}
              for i in range(4)}
    x = lowc["S50"]
    if x["diff_ci95"][0] > 0:
        verdict = "REPRODUCED"
    elif x["ours_minus_det"] <= REF_LOWCPG_S50 / 2:
        verdict = "NOT_REPRODUCED"
    else:
        verdict = "GREY"
    r = {m: cmp_[m]["ratio"] for m in ("S50", "L1s")}
    res = {"plan": "docs/m1/PLAN_M1_UNSEEN.md", "config_sha256": D.sha256(CFG), "rules_sha256": frozen,
           "parts": PJ["parts"], "per_part": per_part, "pooled": pooled, "compare": cmp_, "lowcpg_q1q2": lowc,
           "cpg_quartile_cuts": qs.tolist(), "cpg_strata": strata,
           "verdict_primary": verdict, "ref_lowcpg_s50_chr20_chr9_mean": REF_LOWCPG_S50,
           "exp2_main_rule_descriptive": {"R_S50": r["S50"], "R_L1": r["L1s"],
                                          "retract_condition_met": r["S50"] >= 0.8 and r["L1s"] >= 0.8,
                                          "keep_condition_met": r["L1s"] < 0.5 and cmp_["S50"]["diff_ci95"][0] > 0},
           "seconds": round(time.time() - t0, 1)}
    json.dump(P.norm(res), open(f"{O}/SCORE.json", "w"), indent=1)
    print(json.dumps({"verdict": verdict, "lowcpg": lowc, "ratios": r, "ours_mean": pooled["ours_mean"],
                      "det": {m: {k: v for k, v in pooled["det"][m].items() if k != "counts"} for m in pooled["det"]}}), flush=True)


PJ = json.load(open(f"{ROOT}/data/giraffe/PARTS.json"))

if __name__ == "__main__":
    main()
