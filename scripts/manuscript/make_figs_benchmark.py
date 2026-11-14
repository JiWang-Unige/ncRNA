#!/usr/bin/env python3
"""Manuscript Figures 1-5 (benchmark, sections 1-6), drawn only from the Phase 0 result JSONs copied into
docs/phase0_results/data/ (from 05_phase0/results).  Chromosome lengths for Fig. 1b are the only constants; they are
the values of RESULT_truth_density.md.  Every plotted value goes to docs/manuscript/figures/fig_data_benchmark.json.

  make_figs_benchmark.py [--root <repo>]   -> docs/manuscript/figures/fig{1..5}.{pdf,png}
"""
import argparse
import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

BLUE, ORANGE, AQUA, YELLOW, MAGENTA = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
LNC, MRNA = BLUE, ORANGE
plt.rcParams.update({"font.size": 7, "axes.labelsize": 7, "axes.titlesize": 8, "xtick.labelsize": 6.5,
                     "ytick.labelsize": 6.5, "legend.fontsize": 6.2, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.titleweight": "bold", "axes.titlelocation": "left", "pdf.fonttype": 42,
                     "font.family": "DejaVu Sans", "lines.linewidth": 1.5})
DATA = {}
LEN_MB = {"human_chr20": 64.4, "chicken_5": 59.5, "zebrafish_22": 39.1, "fly_2L": 23.5, "worm_I": 15.1}  # RESULT_truth_density.md
SPECIES = [("human_chr20", "Human chr20"), ("chicken_5", "Chicken 5"), ("zebrafish_22", "Zebrafish 22"),
           ("fly_2L", "Fly 2L"), ("worm_I", "Worm I")]


def J(root, p):
    return json.load(open(os.path.join(root, "docs/phase0_results/data", p)))


def ygrid(ax, axis="y"):
    (ax.yaxis if axis == "y" else ax.xaxis).grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def save(fig, out, name):
    for f in ("pdf", "png"):
        fig.savefig(f"{out}/{name}.{f}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def fig1(root, out):
    rows = []
    for f in glob.glob(os.path.join(root, "docs/phase0_results/data/0.2_ruler/genomewide/v41asPred_chr*.json")):
        d = json.load(open(f))["scores"]["new"]
        rows.append({"chrom": json.load(open(f))["chrom"], "lnc_sn": d["lnc"]["L1_locus"]["sn"], "lnc_n": d["lnc"]["n_truth_genes"],
                     "mrna_sn": d["mrna"]["L1_locus"]["sn"], "mrna_n": d["mrna"]["n_truth_genes"]})
    key = lambda r: (100 if r["chrom"] == "chrX" else int(r["chrom"][3:]))  # noqa: E731
    rows.sort(key=key)
    gw = {c: sum(r[f"{c}_sn"] * r[f"{c}_n"] for r in rows) / sum(r[f"{c}_n"] for r in rows) for c in ("lnc", "mrna")}
    ps = J(root, "0.1_genatator/PANEL_SUMMARY.json")
    dens = [{"species": lab, "lnc_per_mb": ps[k]["truth_genes"]["lnc"] / LEN_MB[k],
             "lnc_over_mrna": ps[k]["truth_genes"]["lnc"] / ps[k]["truth_genes"]["mrna"]} for k, lab in SPECIES]
    DATA["fig1"] = {"per_chrom": rows, "genome_wide_weighted": gw, "density": dens}

    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.5), gridspec_kw={"wspace": 0.3, "width_ratios": [1.7, 1]})
    ax = axs[0]
    x = np.arange(len(rows))
    ax.scatter(x, [r["mrna_sn"] for r in rows], s=12, color=MRNA, label="protein-coding loci", zorder=3)
    ax.scatter(x, [r["lnc_sn"] for r in rows], s=12, color=LNC, label="lncRNA loci", zorder=3)
    for c, col in (("lnc", LNC), ("mrna", MRNA)):
        ax.axhline(gw[c], color=col, linewidth=0.8, linestyle=(0, (3, 2)))
        ax.text(-0.4, gw[c] + (0.025 if c == "mrna" else -0.16), f"genome-wide {gw[c]:.3f}", ha="left", fontsize=5.8, color=INK2)
    ax.set_xticks(x, [r["chrom"][3:] for r in rows], fontsize=5.6)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("v50 loci recovered by v41")
    ax.set_ylim(0, 1.08)
    ygrid(ax)
    ax.legend(frameon=False, loc="center right", bbox_to_anchor=(1.0, 0.68))
    ax.set_title("a  lncRNA truth is half new; coding truth is stable")
    ax = axs[1]
    y = np.arange(len(dens))[::-1]
    ax.barh(y, [d["lnc_per_mb"] for d in dens], 0.6, color=LNC)
    for yi, d in zip(y, dens):
        ax.text(d["lnc_per_mb"] + 0.5, yi, f"lnc:mRNA {d['lnc_over_mrna']:.2f}", va="center", fontsize=5.8, color=INK2)
    ax.set_yticks(y, [d["species"] for d in dens])
    ax.set_xlabel("Annotated lncRNA genes per Mb")
    ax.set_xlim(0, 38)
    ygrid(ax, "x")
    ax.set_title("b  Truth density across lineages")
    save(fig, out, "fig1")


def fig2(root, out):
    t = J(root, "0.3_fp_audit/t2t_chr20/T2T_PANEL_SUMMARY.json")["models"]
    om = J(root, "0.3_fp_audit/t2t_chr20/official_metric_lnc_vs_mrna.json")
    rows = [{"model": k, "stranded": v["stranded"], "lnc": v["lnc_L1_sn"], "mrna": v["mrna_L1_sn"],
             "lnc_reported_as_coding": v["lnc_genes_reported_as_coding"], "n_lnc_pred": v["n_pred"]["lnc"]} for k, v in t.items()]
    off = {c: om[f"GENATATOR|{c}|k250"]["segmentation-level"] for c in ("lnc_RNA", "mRNA")}
    DATA["fig2"] = {"panel": rows, "official_k250_segmentation": off}
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.9), gridspec_kw={"wspace": 0.55, "width_ratios": [1.6, 0.7, 0.9]})
    ax = axs[0]
    order = [r for r in rows if r["stranded"]] + [r for r in rows if not r["stranded"]]
    y = np.arange(len(order))[::-1]
    for yi, r in zip(y, order):
        ax.plot([r["lnc"], r["mrna"]], [yi, yi], color=GRID, linewidth=2, zorder=1)
        ax.scatter(r["mrna"], yi, s=16, color=MRNA, zorder=3, label="mRNA loci" if yi == y[0] else None)
        ax.scatter(r["lnc"], yi, s=16, color=LNC, zorder=3, label="lncRNA loci" if yi == y[0] else None)
    ns = sum(r["stranded"] for r in rows)
    ax.axhline(y[ns - 1] - 0.5, color=MUTED, linewidth=0.5)
    ax.text(0.5, y[0] + 0.55, "gene annotators (stranded)", ha="center", fontsize=5.6, color=MUTED)
    ax.text(0.03, y[ns] + 0.0, "segmenters\n(unstranded)", ha="left", va="center", fontsize=5.6, color=MUTED)
    ax.set_yticks(y, [r["model"] for r in order], fontsize=5.8)
    ax.set_xlabel("Loci recovered (exon overlap, L1)")
    ax.set_xlim(0, 1.02)
    ygrid(ax, "x")
    ax.legend(frameon=False, loc="lower left", bbox_to_anchor=(0.0, -0.32), ncol=2)
    ax.set_title("a  T2T-CHM13 chr20 leaderboard")
    ax = axs[1]
    vals = [off["lnc_RNA"]["recall"], off["mRNA"]["recall"]]
    ax.bar([0, 1], vals, 0.6, color=[LNC, MRNA])
    for xi, v in enumerate(vals):
        ax.text(xi, v + 0.015, f"{v:.3f}", ha="center", fontsize=5.8, color=INK2)
    ax.set_xticks([0, 1], ["lncRNA", "mRNA"])
    ax.set_ylabel("Segmentation recall, k = 250 bp")
    ax.set_ylim(0, 0.6)
    ygrid(ax)
    ax.set_title("b  Official metric")
    ax.text(0.5, -0.22, "GENATATOR", transform=ax.transAxes, ha="center", fontsize=5.8, color=MUTED)
    ax = axs[2]
    st = [r for r in rows if r["stranded"]]
    y2 = np.arange(len(st))[::-1]
    ax.barh(y2, [r["lnc_reported_as_coding"] for r in st], 0.6, color=MAGENTA)
    for yi, r in zip(y2, st):
        ax.text(r["lnc_reported_as_coding"] + 1, yi, str(r["lnc_reported_as_coding"]), va="center", fontsize=5.8, color=INK2)
    ax.set_yticks(y2, [r["model"] for r in st], fontsize=5.8)
    ax.set_xlabel("lncRNA loci called coding")
    ax.set_xlim(0, 62)
    ygrid(ax, "x")
    ax.set_title("c  Misreported")
    save(fig, out, "fig2")


def fig3(root, out):
    ps = J(root, "0.1_genatator/PANEL_SUMMARY.json")
    rows = [{"species": lab, "lnc": ps[k]["lnc_L1_sn"], "mrna": ps[k]["mrna_L1_sn"], "n_lnc_pred": ps[k]["n_pred"]["lnc"]} for k, lab in SPECIES]
    DATA["fig3"] = rows
    fig, ax = plt.subplots(figsize=(3.4, 2.4))
    x = np.arange(len(rows))
    ax.plot(x, [r["mrna"] for r in rows], color=MRNA, marker="o", markersize=4, label="mRNA loci")
    ax.plot(x, [r["lnc"] for r in rows], color=LNC, marker="o", markersize=4, label="lncRNA loci")
    for xi, r in zip(x, rows):
        ax.text(xi, r["lnc"] + 0.04, f"{r['n_lnc_pred']} lnc\ncalls", ha="center", fontsize=5.2, color=MUTED)
    ax.set_xticks(x, [r["species"].replace(" ", "\n", 1) for r in rows], fontsize=5.8)
    ax.set_ylabel("Loci recovered (L1)")
    ax.set_ylim(0, 1.05)
    ygrid(ax)
    ax.legend(frameon=False, loc="upper right")
    ax.set_title("GENATATOR beyond its training clade")
    save(fig, out, "fig3")


def fig4(root, out):
    names = [("probe_generanno_cds", "GENERanno-1.2b (CDS)", BLUE), ("probe_generanno_base", "GENERanno-0.5b", AQUA),
             ("probe_ntv2", "NT-v2-500M", YELLOW), ("probe_caduceus_ps", "Caduceus-PS", MAGENTA)]
    res = {}
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.4), gridspec_kw={"wspace": 0.38})
    for f, lab, col in names:
        d = J(root, f"0.4_probe/{f}.json")
        layers = sorted(map(int, d["auprc"]))
        depth = np.array(layers) / max(layers)
        res[lab] = {}
        for ax, task in zip(axs, ("lnc_exon", "mrna_exon", "lnc_tss")):
            fold = np.array([d["auprc"][str(L)][task] for L in layers]) / d["prevalence"][task]
            res[lab][task] = {"prior": d["prevalence"][task], "auprc_by_layer": [d["auprc"][str(L)][task] for L in layers]}
            ax.plot(depth, fold, color=col, label=lab, linewidth=1.2)
    pri = J(root, "0.4_probe/probe_generanno_cds.json")["prevalence"]
    for ax, title in zip(axs, (f"a  lncRNA exon (prior {pri['lnc_exon']:.3f})", f"b  mRNA exon (prior {pri['mrna_exon']:.3f})", f"c  lncRNA TSS (prior {pri['lnc_tss']:.3f})")):
        ax.axhline(1, color=MUTED, linewidth=0.8, linestyle=(0, (3, 2)))
        ax.set_xlabel("Relative layer depth")
        ax.set_ylabel("AUPRC / prior")
        ax.set_ylim(bottom=0)
        ygrid(ax)
        ax.set_title(title)
    axs[0].legend(frameon=False, loc="upper left", fontsize=5.6)
    DATA["fig4"] = res
    save(fig, out, "fig4")


def fig5(root, out):
    s = J(root, "0.5_unknown/SUMMARY_0.5.json")
    e = J(root, "0.5_unknown/EXPLORATORY_auprc.json")
    arms = [("F", "Background\n(F)"), ("Uev", "Unknown,\nevidence (U-ev)"), ("Upu", "Unknown,\nall (U-pu)")]
    seeds = ("s42", "s123", "s456")
    rec = {a: [s["runs"][f"{a}_{sd}"]["delta_L1_sn"] for sd in seeds] for a, _ in arms}
    ap = {a: [e["runs"][f"{a}_{sd}"]["A_delta_vs_intergenic_ap"] for sd in seeds] for a, _ in arms}
    pr = {a: [s["runs"][f"{a}_{sd}"]["strict_pr_v50"] for sd in seeds] for a, _ in arms}
    DATA["fig5"] = {"delta_L1_sn": rec, "A_ap": ap, "strict_pr_v50": pr, "Uev_minus_F": s["Uev_minus_F"], "Upu_minus_F": s["Upu_minus_F"]}
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.3), gridspec_kw={"wspace": 0.45})
    for ax, d, lab, t in ((axs[0], rec, "Post-2022 loci recovered (L1)", "a  Recovery of new loci"),
                          (axs[1], ap, "AUPRC, new-locus exons", "b  Ranking (exploratory)"),
                          (axs[2], pr, "Strict precision (v50)", "c  Precision")):
        for i, (a, al) in enumerate(arms):
            v = np.array(d[a])
            col = [MUTED, BLUE, AQUA][i]   # truncated axis: means as ticks, not bars (bars must start at zero)
            ax.hlines(v.mean(), i - 0.25, i + 0.25, color=col, linewidth=2)
            ax.scatter(np.full(3, i), v, s=12, color=col, zorder=3, linewidths=0)
        ax.set_xticks(range(3), [al for _, al in arms], fontsize=5.6)
        lo = min(min(x) for x in d.values())
        hi = max(max(x) for x in d.values())
        ax.set_ylim(lo - (hi - lo) * 1.5, hi + (hi - lo) * 0.6)
        ax.set_ylabel(lab)
        ygrid(ax)
        ax.set_title(t)
    axs[0].text(0.02, 0.97, "U-ev − F: +{:.3f} (95% CI {:.3f}–{:.3f})".format(s["Uev_minus_F"]["delta_L1_sn_diff"], *s["Uev_minus_F"]["ci95"]),
                transform=axs[0].transAxes, fontsize=5.4, color=INK2, va="top")
    save(fig, out, "fig5")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
    a = ap.parse_args()
    out = os.path.join(a.root, "docs/manuscript/figures")
    os.makedirs(out, exist_ok=True)
    for f in (fig1, fig2, fig3, fig4, fig5):
        f(a.root, out)
    json.dump(DATA, open(f"{out}/fig_data_benchmark.json", "w"), indent=1, default=float)
    print("ok")


if __name__ == "__main__":
    main()
