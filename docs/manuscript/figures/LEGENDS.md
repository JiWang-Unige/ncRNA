# Figure legends (Figures 6 and 7, first versions, 2026-11-14)

> Drawn by `scripts/manuscript/make_figs.py` from the result JSONs only; every plotted value is in `fig_data.json`.
> Figures 1–5 (benchmark) are not drawn yet. The ablation panel planned for Figure 6 moved to a supplementary table
> (the per-seed deltas for the third seed live in the batch-2 tables of `RESULT_P15.md`, not in one JSON).

**Figure 6 | Sequence-only detection of lncRNA-associated transcribed regions.**
**a**, Recall of annotated lncRNA genes on held-out human chromosomes (hg38 chr20 and chr9, T2T-CHM13 chr20) and,
without any training on them, on chromosomes of other species. Filled circles, exon-overlap criterion (L1); open
circles, one-to-one span criterion (S50; each prediction and gene cover ≥ 50% of the other's span); diamonds,
GENATATOR under L1. Means of three models (chr9: calibrated v2.0 models). **b**, Share of chr20 predictions outside
GENCODE v50 that carry same-strand transcription evidence (spliced RNA-seq junctions, long-read or assembled
transcripts, or CAGE plus polyA at both ends) against position-matched random intervals; intergenic: no annotated
gene within 5 kb on either strand; last pair: regions called by at least two of three models. **c**, Detection
rate of lncRNA genes by quartile of promoter CpG observed/expected ratio and of gene span, chr18 (development) and
chr20 (held-out).

**Figure 7 | The ends of transcription units limit ab initio lncRNA annotation.**
**a**, Recall of lncRNA and mRNA genes under exon overlap (filled) and one-to-one span matching (hatched) for our
models (chr20, v1.9; chr9, v2.0), GENCODE v41 scored as a prediction of v50, GENATATOR, SegmentNT and Tiberius
(chr20). **b**, Exact lncRNA intron chains recovered by GENATATOR's interval segmentation model when given the
annotated transcript interval, the gene span, or the transcript interval trimmed or padded at both ends (200
multi-exon lncRNA genes, chr18; chr18 is in GENATATOR's training data, so values are upper bounds). **c**, Share of
lncRNA starts whose highest-scoring position within ±20 kb lies within 1 kb of the annotated start, by promoter
CpG quartile: our boundary logits and end probes (seven models; line, mean; band, range) on chr18 and on held-out
chr20, and Puffin on its own held-out chromosomes chr8 and chr9. **d**, The same read-out for polyA-supported
transcript ends: our seven models (bars, mean; dots, individual models) and PolyaID / PolyaStrength (strongest
site among positions classified as polyadenylation sites). Dashed line, chance (2 kb of a 40-kb window).
