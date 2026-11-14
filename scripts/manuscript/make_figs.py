#!/usr/bin/env python3
"""Manuscript Figures 6 (locus detection) and 7 (end localisation), drawn only from result JSONs in
docs/phase1_results/.  Every plotted value is also written to docs/manuscript/figures/fig_data.json for audit.

  make_figs.py [--root <repo>]   -> docs/manuscript/figures/fig{6,7}.{pdf,png}, fig_data.json
"""
import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# reference palette (dataviz skill, light mode): categorical slots in fixed order + neutral inks
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
LNC, MRNA = BLUE, ORANGE
plt.rcParams.update({"font.size": 7, "axes.labelsize": 7, "axes.titlesize": 8, "xtick.labelsize": 6.5,
                     "ytick.labelsize": 6.5, "legend.fontsize": 6.5, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.titleweight": "bold", "axes.titlelocation": "left", "pdf.fonttype": 42,
                     "font.family": "DejaVu Sans", "lines.linewidth": 1.5})
DATA = {}


def J(root, p):
    return json.load(open(os.path.join(root, "docs/phase1_results", p)))


def ygrid(ax):
    ax.yaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def chance(ax, y=0.05, label="chance (±1 kb of ±20 kb)"):
    ax.axhline(y, color=MUTED, linewidth=0.8, linestyle=(0, (3, 2)))
    ax.text(ax.get_xlim()[1], y, " " + label, va="center", ha="left", color=MUTED, fontsize=6, clip_on=False)


# ---------------------------------------------------------------- Figure 7
def fig7(root, out):
    b1 = J(root, "B1/B1.json")["sets"]
    sets = [("Ours\nchr20", "v19_chr20", True), ("Ours\nchr9", "v20_chr9", True),
            ("v41 as\nprediction", "gencode_v41_as_pred_chr20", False), ("GENA-\nTATOR", "genatator_chr20", False),
            ("Segment-\nNT", "segnt_t0.5_both_chr20", False), ("Tiberius", "tiberius_chr20", False)]
    rows = []
    for lab, k, mean in sets:
        v = b1[k]["mean"] if mean else {c: {"L1_sn": b1[k][c]["L1"]["sn"], "S50_sn": b1[k][c]["S50"]["sn"]} for c in ("lnc", "mrna")}
        rows.append((lab, v["lnc"]["L1_sn"], v["lnc"]["S50_sn"], v["mrna"]["L1_sn"], v["mrna"]["S50_sn"]))
    DATA["fig7a"] = [dict(zip(("set", "lnc_L1", "lnc_S50", "mrna_L1", "mrna_S50"), r)) for r in rows]

    s11 = J(root, "S11/GATE.json")["summary"]
    cond = [("−1 kb", "I1"), ("transcript\ninterval", "A"), ("gene\nspan", "B"), ("+1 kb", "O1"), ("+5 kb", "O5"),
            ("+10 kb", "O10"), ("+20 kb", "O20")]
    DATA["fig7b"] = {c: s11[k]["chain_correct"] for c, k in cond}

    e18, e20, dd = J(root, "ENDS/ENDS.json")["sources"], J(root, "REP20/ENDS_chr20.json")["sources"], J(root, "DEDEND/DEDEND.json")
    q = lambda src: [[s["lnc_TSS"][f"cpg_q{i}"]["loc1k"] for i in range(1, 5)] for s in src.values()]  # noqa: E731
    c18, c20 = np.array(q(e18)), np.array(q(e20))
    pf = {g: [dd["puffin"][g]["lnc_TSS"][f"cpg_q{i}"]["loc1k"] for i in range(1, 5)] for g in dd["puffin"]}
    DATA["fig7c"] = {"ours_chr18_7models": c18.tolist(), "ours_chr20_7models": c20.tolist(), "puffin": pf}

    tes = lambda src, k: [s[f"{k}_TES"]["with_evidence"]["loc1k"] for s in src.values()]  # noqa: E731
    pa = {c: {k: dd["polya"][c]["primary"][f"{k}_TES"]["with_evidence"]["loc1k"] for k in ("lnc", "mrna")} for c in ("chr18", "chr20")}
    DATA["fig7d"] = {"ours_chr18": {k: tes(e18, k) for k in ("lnc", "mrna")}, "ours_chr20": {k: tes(e20, k) for k in ("lnc", "mrna")},
                     "polyaid": pa}

    fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.6), gridspec_kw={"hspace": 0.62, "wspace": 0.38})
    # a
    ax = axs[0, 0]
    x = np.arange(len(rows))
    w = 0.2
    for i, (key, col, hatch, lab) in enumerate([(1, LNC, None, "lncRNA, exon overlap (L1)"), (2, LNC, "////", "lncRNA, one-to-one span (S50)"),
                                                 (3, MRNA, None, "mRNA, L1"), (4, MRNA, "////", "mRNA, S50")]):
        vals = [r[key] for r in rows]
        ax.bar(x + (i - 1.5) * w, vals, w * 0.9, color=col if hatch is None else "white", edgecolor=col, hatch=hatch,
               linewidth=0.8, label=lab)
    ax.set_xticks(x, [r[0] for r in rows], fontsize=5.8)
    ax.set_ylabel("Recall of reference genes")
    ax.set_ylim(0, 1)
    ygrid(ax)
    ax.legend(frameon=False, ncol=2, loc="lower left", bbox_to_anchor=(0.0, 1.01), fontsize=5.8, handlelength=1.2)
    ax.set_title("a  Detected but not delimited", pad=30)
    # b
    ax = axs[0, 1]
    xs = np.arange(len(cond))
    vals = [DATA["fig7b"][c] for c, _ in cond]
    ax.plot(xs, vals, color=LNC, marker="o", markersize=4)
    for xi, v in zip(xs, vals):
        ax.text(xi + 0.12, v + 0.025, f"{v:.2f}", ha="left", fontsize=5.8, color=INK2)
    ax.set_xticks(xs, [c for c, _ in cond], fontsize=5.8)
    ax.set_ylabel("Exact lncRNA intron chain")
    ax.set_ylim(0, 0.45)
    ygrid(ax)
    ax.set_title("b  Structure needs sub-kilobase ends", pad=30)
    ax.text(0, -0.2, "← trimmed", transform=ax.get_xaxis_transform(), fontsize=5.6, color=MUTED, ha="center")
    ax.text(4.5, -0.2, "padded →", transform=ax.get_xaxis_transform(), fontsize=5.6, color=MUTED, ha="center")
    ax.text(0.98, 0.95, "GENATATOR segmentation,\n200 lncRNA genes, chr18", transform=ax.transAxes, ha="right", va="top", fontsize=5.8, color=MUTED)
    # c
    ax = axs[1, 0]
    xq = np.arange(1, 5)
    for arr, col, lab in [(c18, LNC, "Our features, chr18 (7 models)"), (c20, AQUA, "Our features, chr20 held-out (7)")]:
        ax.fill_between(xq, arr.min(0), arr.max(0), color=col, alpha=0.15, linewidth=0)
        ax.plot(xq, arr.mean(0), color=col, marker="o", markersize=3.5, label=lab)
    ax.plot(xq, pf["heldout_chr8_chr9"], color=ORANGE, marker="s", markersize=3.5, label="Puffin, chr8+9 (its held-out)")
    ax.set_xticks(xq, ["Q1\n(low CpG)", "Q2", "Q3", "Q4\n(CpG island)"])
    ax.set_ylabel("lncRNA starts located within 1 kb")
    ax.set_ylim(0, 0.7)
    ax.set_xlim(0.8, 4.2)
    ygrid(ax)
    chance(ax, label="chance")
    ax.legend(frameon=False, loc="upper left", fontsize=5.8)
    ax.set_title("c  Starts: only at CpG-island promoters")
    # d
    ax = axs[1, 1]
    groups = [("Our features\nchr18", DATA["fig7d"]["ours_chr18"]), ("Our features\nchr20", DATA["fig7d"]["ours_chr20"]),
              ("PolyaID\nchr18", {k: [pa["chr18"][k]] for k in ("lnc", "mrna")}), ("PolyaID\nchr20", {k: [pa["chr20"][k]] for k in ("lnc", "mrna")})]
    for gi, (lab, d) in enumerate(groups):
        for j, (k, col) in enumerate((("lnc", LNC), ("mrna", MRNA))):
            v = np.array(d[k])
            xpos = gi + (j - 0.5) * 0.32
            ax.bar(xpos, v.mean(), 0.28, color=col, label=("lncRNA ends" if k == "lnc" else "mRNA ends") if gi == 0 else None)
            if len(v) > 1:
                ax.scatter(np.full(len(v), xpos), v, s=6, color=INK2, zorder=3, linewidths=0)
    ax.set_xticks(range(len(groups)), [g[0] for g in groups], fontsize=5.8)
    ax.set_ylabel("polyA-supported ends within 1 kb")
    ax.set_ylim(0, 0.7)
    ax.set_xlim(-0.6, 3.6)
    ygrid(ax)
    chance(ax, label="chance")
    ax.legend(frameon=False, loc="upper left", fontsize=5.8)
    ax.set_title("d  Ends: not located, even by a polyA model")
    for f in ("pdf", "png"):
        fig.savefig(f"{out}/fig7.{f}", dpi=300, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- Figure 6
def fig6(root, out):
    b1, bm = J(root, "B1/B1.json")["sets"], J(root, "B1/B1_MORE.json")["sets"]
    sets = [("chr20", b1["v19_chr20"]["mean"], b1["genatator_chr20"]),
            ("T2T chr20", bm["v19_t2t_chr20"]["mean"], bm["genatator_t2t_chr20"]),
            ("chr9", b1["v20_chr9"]["mean"], None),
            ("mouse chr19", bm["v19_mouse_chr19"]["mean"], None),
            ("chicken 5", bm["v19_chicken_5"]["mean"], bm["genatator_chicken_5"]),
            ("zebrafish 22", bm["v19_zebrafish_22"]["mean"], bm["genatator_zebrafish_22"]),
            ("fly 2L", bm["v19_fly_2L"]["mean"], bm["genatator_fly_2L"]),
            ("worm I", bm["v19_worm_I"]["mean"], bm["genatator_worm_I"])]
    rows = []
    for lab, m, g in sets:
        rows.append({"set": lab, "ours_L1": m["lnc"]["L1_sn"], "ours_S50": m["lnc"]["S50_sn"], "ours_pr": m["lnc"]["L1_pr_strict"],
                     "gen_L1": g["lnc"]["L1"]["sn"] if g else None, "gen_S50": g["lnc"]["S50"]["sn"] if g else None})
    DATA["fig6a"] = rows

    bv = J(root, "BIOVAL/RESULT.json")
    strata = [("All", "all"), ("Inter-\ngenic", "S_inter"), ("Other", "S_other")]
    bvals = [(lab, bv["mean"][k]["primary"], bv["mean_ctrl"][k]["primary"]) for lab, k in strata]
    bvals.append(("≥2 of 3\nmodels", bv["consensus"]["primary"], bv["consensus"]["ctrl_mean"]))
    DATA["fig6b"] = bvals

    s3 = {c: J(root, f"S3/S3_{c}.json")["R2"] for c in ("chr18", "chr20")}
    DATA["fig6c"] = {c: {f: [s3[c][f][f"Q{i}"]["detected"] for i in range(1, 5)] for f in ("cpg_oe", "log_span")} for c in s3}

    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.6), gridspec_kw={"wspace": 0.42, "width_ratios": [1.5, 1.1, 1.0]})
    ax = axs[0]
    y = np.arange(len(rows))[::-1]
    for yi, r in zip(y, rows):
        ax.plot([r["ours_S50"], r["ours_L1"]], [yi, yi], color=GRID, linewidth=2, zorder=1)
        ax.scatter(r["ours_L1"], yi, s=22, color=LNC, zorder=3, label="Ours, exon overlap (L1)" if yi == y[0] else None)
        ax.scatter(r["ours_S50"], yi, s=22, facecolor="white", edgecolor=LNC, linewidth=1.2, zorder=3,
                   label="Ours, one-to-one span (S50)" if yi == y[0] else None)
        if r["gen_L1"] is not None:
            ax.scatter(r["gen_L1"], yi, s=14, marker="D", color=MUTED, zorder=2, label="GENATATOR, L1" if yi == y[0] else None)
    ax.set_yticks(y, [r["set"] for r in rows])
    ax.set_xlabel("lncRNA genes recovered")
    ax.set_xlim(0, 0.45)
    ax.xaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.axhline(4.5, color=MUTED, linewidth=0.5)
    ax.text(0.445, 4.62, "human (held-out)", ha="right", va="bottom", fontsize=5.6, color=MUTED)
    ax.text(0.445, 4.38, "other species (zero-shot)", ha="right", va="top", fontsize=5.6, color=MUTED)
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.45, -0.2), ncol=2, fontsize=5.6)
    ax.set_title("a  Held-out detection")
    ax = axs[1]
    x = np.arange(len(bvals))
    ax.bar(x - 0.18, [b[1] for b in bvals], 0.34, color=LNC, label="Predicted regions")
    ax.bar(x + 0.18, [b[2] for b in bvals], 0.34, color=GRID, edgecolor=MUTED, linewidth=0.6, label="Matched random")
    ax.set_xticks(x, [b[0] for b in bvals], fontsize=5.6)
    ax.set_ylabel("With transcription evidence")
    ax.set_ylim(0, 1)
    ygrid(ax)
    ax.legend(frameon=False, loc="upper left", fontsize=5.8)
    ax.set_title("b  Outside the reference")
    ax.set_xlabel("chr20 predictions outside GENCODE v50", fontsize=6)
    ax = axs[2]
    xq = np.arange(1, 5)
    for f, col, lab in (("cpg_oe", LNC, "promoter CpG o/e"), ("log_span", AQUA, "gene span")):
        for c, ls in (("chr18", (0, (2, 1.5))), ("chr20", "-")):
            ax.plot(xq, DATA["fig6c"][c][f], color=col, linestyle=ls, marker="o", markersize=3, label=f"{lab}, {c}")
    ax.set_xticks(xq, ["Q1", "Q2", "Q3", "Q4"])
    ax.set_ylabel("lncRNA genes detected (L1)")
    ax.set_ylim(0, 0.65)
    ygrid(ax)
    ax.legend(frameon=False, loc="upper left", fontsize=5.4)
    ax.set_title("c  What is detected")
    for f in ("pdf", "png"):
        fig.savefig(f"{out}/fig6.{f}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
    a = ap.parse_args()
    out = os.path.join(a.root, "docs/manuscript/figures")
    os.makedirs(out, exist_ok=True)
    fig6(a.root, out)
    fig7(a.root, out)
    json.dump(DATA, open(f"{out}/fig_data.json", "w"), indent=1, default=float)
    print("ok", out)


if __name__ == "__main__":
    main()
