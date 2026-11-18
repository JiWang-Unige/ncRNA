#!/usr/bin/env python3
"""M1 weeks 2-4: markdown tables from a w24_score.py output (reporting only; no scoring).
  result_tables.py <RESULT.json | SCORE_DEV_chr18.json>  -> markdown on stdout
"""
import json
import sys

from statistics import mean

METHODS = ("ours", "ev_farthest", "pair", "snap", "dna", "stringtie", "isoquant")
GRAD = ("ours", "ev_farthest", "pair", "snap")


def f(x, n=3):
    return "—" if x is None else f"{x:.{n}f}"


def main():
    d = json.load(open(sys.argv[1]))
    P, CI, D = d["pooled"], d["ci95_S50_sn"], d["paired"]
    chroms = list(d["info"]["chroms"])
    print(f"真值 lnc 位点数（合并）：{d['n_truth_lnc']}；染色体：{', '.join(chroms)}\n")
    print("**100% 证据（合并）**\n")
    print("| 方法 | S50 召回 [95% CI] | 精度（剔除不确定） | 严格精度 | 预测数 |")
    print("|---|---|---|---|---|")
    for m in METHODS:
        k = f"{m}|f1.0_s0"
        if k in P:
            p, c = P[k], CI[k]
            print(f"| {m} | {f(p['S50_sn'])} [{f(c[0])}, {f(c[1])}] | {f(p['S50_pr_known'])} | {f(p['S50_pr'])} | {p['n_pred']:.0f} |")
    print("\n**配对差（按位点自助，2,000 次）**\n")
    print("| 证据量 | ours − ev_farthest | ours − snap | ours − pair |")
    print("|---|---|---|---|")
    tags = sorted({k.split("|")[1] for k in D}, key=lambda t: (-float(t[1:].split("_")[0]), t))
    for t in tags:
        cells = []
        for o in ("ev_farthest", "snap", "pair"):
            x = D[f"ours-{o}|{t}"]
            cells.append(f"{x['diff']:+.3f} [{x['ci95'][0]:+.3f}, {x['ci95'][1]:+.3f}]")
        print(f"| {t} | " + " | ".join(cells) + " |")
    print("\n**证据量梯度（合并；30% / 10% 为 3 个种子的平均，括号内为种子间范围）**\n")
    print("| 方法 | 100% 召回 / 精度 | 30% 召回 / 精度 | 10% 召回 / 精度 |")
    print("|---|---|---|---|")
    for m in GRAD:
        cells = []
        for fr in ("1.0", "0.3", "0.1"):
            ks = [k for k in P if k.startswith(f"{m}|f{fr}_")]
            sn = [P[k]["S50_sn"] for k in ks]
            pk = [P[k]["S50_pr_known"] for k in ks if P[k]["S50_pr_known"] is not None]
            rng = f" ({min(sn):.3f}–{max(sn):.3f})" if len(sn) > 1 else ""
            cells.append(f"{mean(sn):.3f}{rng} / {f(float(mean(pk)) if pk else None)}")
        print(f"| {m} | " + " | ".join(cells) + " |")
    print("\n**每条染色体（100% 证据；S50 召回 / 精度（剔除不确定） / 预测数）**\n")
    print("| 染色体 | 真值 | 证据覆盖（起点 / 终点 / 两者） | " + " | ".join(METHODS) + " |")
    print("|---|---|---|" + "---|" * len(METHODS))
    for c in chroms:
        ci = d["info"]["chroms"][c]
        cov = ci["coverage_f1.0_s0"]
        cells = []
        for m in METHODS:
            r = ci["by_method"].get(f"{m}|f1.0_s0")
            cells.append("—" if r is None else f"{r['S50_sn']:.3f} / {f(r['S50_pr_known'], 2)} / {r['n_pred']:.0f}")
        print(f"| {c} | {ci['n_truth_lnc']} | {cov['tss']:.2f} / {cov['pas']:.2f} / {cov['both']:.2f} | " + " | ".join(cells) + " |")
    n = sum(d["info"]["chroms"][c]["coverage_f1.0_s0"]["n"] for c in chroms)
    both = sum(d["info"]["chroms"][c]["coverage_f1.0_s0"]["both"] * d["info"]["chroms"][c]["coverage_f1.0_s0"]["n"] for c in chroms)
    tss = sum(d["info"]["chroms"][c]["coverage_f1.0_s0"]["tss"] * d["info"]["chroms"][c]["coverage_f1.0_s0"]["n"] for c in chroms)
    pas = sum(d["info"]["chroms"][c]["coverage_f1.0_s0"]["pas"] * d["info"]["chroms"][c]["coverage_f1.0_s0"]["n"] for c in chroms)
    print(f"\n证据覆盖（合并，±100 bp）：起点 {tss / n:.3f}，终点 {pas / n:.3f}，两者 {both / n:.3f}（n = {n}）")
    for fr in ("0.3", "0.1"):
        vals = [d["info"]["chroms"][c][k]["both"] * d["info"]["chroms"][c][k]["n"] for c in chroms
                for k in d["info"]["chroms"][c] if k.startswith(f"coverage_f{fr}_")]
        if vals:
            print(f"- {int(float(fr) * 100)}% 证据：两者都有 {sum(vals) / (n * 3):.3f}（3 个种子平均）")
    if "ours_not_isoquant" in d:
        o = d["ours_not_isoquant"]
        print(f"\nours 定位到（≥ 2/3 模型）而 IsoQuant（rep1）没定位到：{o['n']} 个，占真值 {o['frac_of_truth']:.3f}")
    v = d["verdict"]
    g = v["gate"]
    print("\n**判读**\n")
    print(f"- T1：召回 {f(g['S50_sn'])}（≥ 0.30），精度（剔除不确定）{f(g['S50_pr_known'])}（≥ 0.40），"
          f"ours − snap {g['vs_snap']['diff']:+.3f} CI [{g['vs_snap']['ci95'][0]:+.3f}, {g['vs_snap']['ci95'][1]:+.3f}]（下限 > 0）→ "
          f"{'达到' if g['met'] else '未达到'}")
    e = D["ours-ev_farthest|f1.0_s0"]
    print(f"- T2：ours − ev_farthest 召回 {e['diff']:+.3f} CI [{e['ci95'][0]:+.3f}, {e['ci95'][1]:+.3f}]，"
          f"精度差 {v['precision_known_minus_ev']:+.3f} → {'有贡献' if v['sequence_prior_contributes'] else '无贡献'}")
    print(f"- 映射：{v['decision_hint']}")


if __name__ == "__main__":
    main()
