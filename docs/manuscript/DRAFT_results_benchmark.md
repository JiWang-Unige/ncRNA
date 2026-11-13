# Draft — benchmark sections (plan B) of the joint coding × lncRNA annotation paper

> Status: first draft, 2026-10-20; 2026-11-13: S50 numbers and the L1 caveat added (see `DRAFT_results_ends.md`). Numbers are from the Phase 0 reports (`05_phase0/results/`); every number
> has a pointer to its JSON. Not yet checked by a second reader. Figures are planned, not drawn.
> Language: English (manuscript); discussion notes for the authors are in Chinese at the end.

## Working title

*Non-coding genes are invisible to ab initio genome annotation: a benchmark under incomplete truth, and a
grammar-constrained annotator for coding and lncRNA genes*

## Results

### 1. lncRNA truth is incomplete and expanding fast, whereas coding truth is stable

We first asked how stable the reference used to train and evaluate annotators is. Treating GENCODE v41
(July 2022) as a "prediction" of GENCODE v50 (2026), the v41 lncRNA set recovers only **49.0%** of v50 lncRNA
loci genome-wide (17,960 → 34,677 genes on chr1–22, X; 15,604 loci entirely new), versus **99.8%** for
protein-coding loci (19,941 → 20,033) [Fig. 1a; `0.2_ruler/genomewide/*.json`]. Per chromosome, lncRNA locus
recovery ranges from 0.39 (chrX) to 0.58 (chr16); intron-chain recovery on chr20 is 0.41. Precision of v41
against v50 is ≥ 0.98, i.e. the change is almost purely additive [`RESULT_0.2_ruler.md`].

*Implication.* A model trained on any fixed lncRNA annotation is taught that about half of today's lncRNA loci
are background, and an evaluation against any fixed annotation scores correct novel predictions as false
positives. Both training and evaluation must therefore treat unannotated sequence as partially unknown.

Across lineages, lncRNA truth density differs by more than an order of magnitude (lncRNA : mRNA gene ratio
0.05 in *C. elegans* chr I to 1.74 in human chr20; 3.5 lncRNA genes / Mb on zebrafish chr22 under the legacy
GRCz11 build versus 11.7 on chicken chr5 and 14.8 on human chr20) [Fig. 1b; `RESULT_truth_density.md`].
Part of this reflects biology, part annotation effort; precision on sparsely annotated genomes mostly measures
the annotation gap.

### 2. Existing ab initio annotators do not recover lncRNA gene structures

We scored the 13 public predictions of the GENATATOR leaderboard (T2T-CHM13 chr20; RefSeq truth with 434
lncRNA and 546 mRNA genes) with a stranded ruler (Methods) [Fig. 2a; `0.3_fp_audit/t2t_chr20/T2T_PANEL_SUMMARY.json`]:

* The protein-coding annotators — Tiberius, ANNEVO, Helixer and AUGUSTUS — output no lncRNA. Each reports
  some lncRNA loci as protein-coding genes (3, 8, 4 and 52 loci respectively).
* GENATATOR, the only pipeline designed to output lncRNA, recovers **3.5%** of lncRNA loci (locus-level,
  any exon overlap) and **0.7%** of lncRNA intron chains, while recovering 93% of mRNA loci and 65% of mRNA
  intron chains. With the leaderboard's own metric at its headline tolerance (k = 250 bp), lnc_RNA
  segmentation recall is **0.7%** (precision 3.7%), versus 49.6% for mRNA [`official_metric_lnc_vs_mrna.json`].
* Unstranded segmenters (NTv3, SegmentNT) mark 7–51% of lncRNA loci as transcribed but produce no correct
  lncRNA intron chain (≤ 0.5%), with locus precision 0.09–0.23.

Running the released GENATATOR pipeline ourselves on hg38 chr20 against GENCODE v50 gave the same picture
(lncRNA locus Sn 0.028 [95% CI 0.019–0.040], intron chain 0.012; mRNA locus Sn 0.967, intron chain 0.708;
46 lncRNA loci reported only as mRNA) [`0.1_genatator/PANEL_SUMMARY.json`]. Under the stricter one-to-one span
criterion of section 12 (S50), its lncRNA recall is 0.006 and its mRNA recall 0.794 [`B1/B1.json`]. The bottleneck is interval discovery:
on chr21, only 33 of 590 lncRNA genes were touched by any predicted transcript [`RESULT_0.1a_0.3`].

For a representative coding annotator (Helixer, vertebrate model, hg38 chr20), 97% of lncRNA loci are not
touched at all, and **56% of its locus-level false positives fall on annotated lncRNA loci** (18 of 32)
[`0.3_fp_audit/helixer_chr20*.json`]. ANNEVO independently reports that 18.3% of its false positives overlap
annotated non-coding transcripts (different metric) [ANNEVO, Nat Methods 2026].

### 3. Transfer beyond the training clade collapses

GENATATOR is trained on 39 mammals. On four non-mammalian chromosomes (Ensembl 116), it emits almost no lncRNA
(0–13 transcripts; locus Sn ≤ 0.022), and its mRNA recovery falls with phylogenetic distance: locus Sn 0.84
(chicken chr5), 0.49 (zebrafish chr22), 0.21 (*Drosophila* 2L) and 0.006 (*C. elegans* chr I, 15 genes
predicted on the whole chromosome) [Fig. 3; `RESULT_0.1b_cross_lineage.md`].

### 4. Structured small RNAs are well served by covariance models; the gaps are elsewhere

Infernal/Rfam 15.1 and tRNAscan-SE recover nearly all tRNA, snRNA, rRNA and vertebrate snoRNA genes on the
five test chromosomes, but miss most miRNA in chicken (4/40), all piRNA and most WormBase short ncRNA in
*C. elegans*. Most "false positives" are annotation gaps (tRNA genes absent from Ensembl GTFs) or
repeat-derived hits (357 of 397 novel human chr20 hits are Metazoa_SRP matches to Alu) [`RESULT_0.6_smallrna.md`].
Because the reference small-RNA annotations are themselves produced by these tools, we use them as an
integration component, not as a method claim.

### 5. Sequence signal for lncRNA exons is weak but present

Linear probes on frozen DNA language models (hg38 chr21 → chr20) show that lncRNA exon positions are
separable from background, but far less than mRNA exons: best AUPRC 0.38 versus a 0.16 prior
(GENERanno-1.2b, CDS-annotation checkpoint), compared with ~9× the prior for mRNA exons; TSS/TES signals are
weak for every backbone (AUPRC < 0.09) [`RESULT_0.4_probe.md`]. A CDS-supervised backbone was *better*, not
worse, for lncRNA exons and splice sites than the self-supervised ones (NT-v2-500M, GENERanno-0.5b,
Caduceus-PS) — exon/splice features learned from coding supervision transfer to non-coding exons.

### 6. Treating evidence-supported unannotated sequence as unknown gives a small, consistent gain

In a pre-registered test (frozen features, minimal decoder), marginalising the labels of unannotated regions
that carry pre-2022 RNA-seq junction evidence improved recovery of loci first annotated after 2022 by 1.7
percentage points (95% CI 0.4–3.2; pre-declared threshold 5 points, not met) and their ranking (AUPRC +9–11%
relative, all seeds positive) [`RESULT_0.5_unknown.md`]. We therefore use unknown-handling as a component and
report it as an ablation rather than as the main contribution.

## Methods (evaluation framework)

**Ruler.** Predictions and references are reduced to transcripts (chromosome, strand, exons). Metrics per class
(lncRNA, mRNA, small RNA):
* L1 locus: exon overlap on the same strand; Sn over reference genes, Pr over predicted transcripts. L1 does not
  require a prediction to delimit a gene (one prediction may hit several genes); the one-to-one span criterion
  S50 / S25 (`DRAFT_results_ends.md`, Methods) is reported alongside it for every locus-level claim.
* L2 intron chain: exact match with any reference isoform.
* L3 junction and L4 exon: exact; terminal exons are matched on their inner boundary only.
* L5 single-exon: reciprocal overlap ≥ 0.5.
* L6 base-level F1.
* Type confusion between classes.

Strict precision counts predictions in unannotated sequence as errors; lenient precision excludes them
("novel"). Truth incompleteness is quantified as the fraction of old-reference false positives that become
true positives under a newer release. Confidence intervals come from 1,000 bootstrap resamples of genes.
Supported inputs: GENCODE/Ensembl GTF, RefSeq GFF3, and unstranded predictors (scored strand-agnostically,
as in the GENATATOR leaderboard). Unit tests check synthetic cases and the identity truth = prediction for
all classes. [`ruler.py` v1.3]

**Class mapping across annotations.** The lnc class is lncRNA, lincRNA, antisense and sense_overlapping, plus
generic ncRNA / processed_transcript genes whose longest spliced transcript is ≥ 200 nt. Shorter generic
ncRNA genes are mapped to small RNA. Transposable-element genes and TEC are excluded (protocol v1.2).

**Time-split truth.** GENCODE v41 is the "old" reference and v50 the "new" one; the Δ set is the v50 lncRNA
genes with no same-strand exon overlap with any v41 gene (473 on chr20, 15,604 genome-wide).

**Evidence mask (U-ev).** Snaptron srav3h (recount3, 2021) junctions kept if samples ≥ 100, coverage ≥ 500,
canonical motif, length ≤ 500 kb, and no overlap with annotated genes. A ±500 bp flank is added, then
annotated gene spans are subtracted.

## Figure plan

1. a) v41 → v50 growth per chromosome, lncRNA vs mRNA; b) truth density across lineages.
2. a) T2T chr20 leaderboard panel, lncRNA vs mRNA recovery for 13 predictors; b) the official-metric check;
   c) lncRNA loci misreported as mRNA.
3. Cross-lineage transfer of GENATATOR (mRNA and lncRNA).
4. Probe AUPRC per layer and backbone.
5. Unknown-handling: Δ recovery and AUPRC (supplementary).
6. Locus detection: see `DRAFT_results_locus.md` (sections 7–11).
7. End localisation: see `DRAFT_results_ends.md` (sections 12–15).

---

## 作者备注（中文）

- 所有数字都来自 Phase 0 的 JSON，路径见方括号。正式成稿前需要逐条复核；尤其是 GENATATOR 论文正文对 lncRNA 的表述，以及 ANNEVO 的 18.3% 的原始定义，要读原文确认。
- 第 2 节的措辞刻意避免"从未有人做过"，改为可量化的"端到端几乎无法恢复"。审稿人最可能的反驳是"GENATATOR 的 lncRNA 结果是在已知区间内做的"：我们的立场正是"端到端 gene finding"与"区间内分割"是两个不同的任务，需要在 Discussion 里讲清楚。
- 第 5、6 节是方法动机，可能移到方法部分或补充材料，取决于 Phase 1 的结果。
