# Figure legends (Figures 5 and 6, first versions, 2026-11-14; renumbered 2026-11-15)

> Drawn by `scripts/manuscript/make_figs.py` from the result JSONs only; every plotted value is in `fig_data.json`.
> The ablation panel planned for this figure moved to Supplementary Table 2
> (the per-seed deltas for the third seed live in the batch-2 tables of `RESULT_P15.md`, not in one JSON).

**Figure 5 | Sequence-only detection of lncRNA-associated transcribed regions.**
**a**, Recall of annotated lncRNA genes on held-out human chromosomes (hg38 chr20 and chr9, T2T-CHM13 chr20) and,
without any training on them, on chromosomes of other species. Filled circles, exon-overlap criterion (L1); open
circles, one-to-one span criterion (S50; each prediction and gene cover ≥ 50% of the other's span); diamonds,
GENATATOR under L1. Means of three models (chr9: calibrated v2.0 models). **b**, Share of chr20 predictions outside
GENCODE v50 that carry same-strand transcription evidence (spliced RNA-seq junctions, long-read or assembled
transcripts, or CAGE plus polyA at both ends) against position-matched random intervals; intergenic: no annotated
gene within 5 kb on either strand; last pair: regions called by at least two of three models. **c**, Detection
rate of lncRNA genes by quartile of promoter CpG observed/expected ratio and of gene span, chr18 (development) and
chr20 (held-out).

**Figure 6 | The ends of transcription units limit ab initio lncRNA annotation.**
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

# Figure legends (Figures 1–4 and Supplementary Fig. 1, first versions, 2026-11-14; renumbered 2026-11-15)

> Drawn by `scripts/manuscript/make_figs_benchmark.py` from the Phase 0 JSONs copied into `docs/phase0_results/data/`
> (from `05_phase0/results`); every plotted value is in `fig_data_benchmark.json`. Chromosome lengths in Fig. 1b are
> the values of `RESULT_truth_density.md`.

**Figure 1 | Reference lncRNA annotation is incomplete and uneven.** **a**, Share of GENCODE v50 loci on each
chromosome (chr1–22, X) that GENCODE v41 already contained (same-strand exon overlap), for lncRNA (blue) and
protein-coding (orange) genes; dashed lines, genome-wide means weighted by gene number. **b**, Annotated lncRNA
genes per Mb and lncRNA : mRNA gene ratio on the five chromosomes of the evaluation panel (human GENCODE v50;
others Ensembl 116 / Metazoa 63).

**Figure 2 | Existing ab initio annotators do not recover lncRNA genes.** **a**, Recovery of lncRNA and mRNA loci
(exon overlap) by the 13 public predictions of the T2T-CHM13 chr20 leaderboard (RefSeq truth: 434 lncRNA, 546 mRNA
genes); unstranded segmenters are scored strand-agnostically. **b**, The leaderboard's own segmentation recall at
its headline tolerance (k = 250 bp) for GENATATOR, lnc_RNA versus mRNA. **c**, lncRNA loci covered by a predicted
coding gene that hits no reference mRNA, per gene annotator.

**Figure 3 | Transfer beyond the training clade.** Locus recovery by GENATATOR (trained on 39 mammals) on human
chr20 and four non-mammalian chromosomes, mRNA and lncRNA; labels give the number of lncRNA transcripts predicted.

**Figure 4 | Weak but present sequence signal for lncRNA exons.** Linear-probe AUPRC divided by the class prior
across the layers of four frozen DNA language models (hg38 chr21 training, chr20 evaluation; 300,000 positions):
**a**, lncRNA exon; **b**, mRNA exon; **c**, lncRNA transcription start. Dashed line, prior.

**Supplementary Figure 1 | Treating evidence-supported unannotated sequence as unknown.** Pre-registered test
on frozen features with a minimal decoder, three seeds per arm (dots; short lines, mean): unannotated sequence treated
as background (F), as unknown where pre-2022 RNA-seq junctions support it (U-ev), or everywhere (U-pu). **a**,
Recovery of lncRNA loci first annotated after GENCODE v41 (n = 473, chr20); **b**, AUPRC for exons of those loci
(exploratory, post hoc); **c**, strict precision against v50.
