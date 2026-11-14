# Sequence-only annotation finds lncRNA-associated transcribed regions but cannot delimit them: the ends of transcription units limit ab initio lncRNA annotation

> Assembled by `scripts/manuscript/assemble.py` from the section drafts; edit the drafts, not this file. Bracketed file names are provenance pointers for the authors and are removed at submission.

## Abstract

Ab initio gene finders now annotate protein-coding genes from sequence alone, but long non-coding RNA (lncRNA)
genes — as numerous in the human genome as coding genes — remain absent from their output. Here we ask whether a
sequence-only model can annotate lncRNA genes and, where it cannot, which step fails. Using an evaluation framework
built for an incomplete and rapidly growing reference, we developed a joint coding–lncRNA annotator that decodes a
gene grammar over frozen DNA language-model features. On chromosomes read once after all rules were frozen, it
detected 30% of lncRNA genes at a strict precision of 0.42–0.50, more than ten times the recall of existing tools,
and the detection transferred to mouse and non-mammalian vertebrates. Yet only about 7% of lncRNA genes were
delimited by a single prediction, against 55% of coding genes. Decomposing this gap shows that the ends of
transcription units, not exon structure, limit annotation: given ends within a kilobase, structure becomes partly
recoverable, whereas lncRNA starts can be located from sequence only at CpG-island promoters and lncRNA ends not at
all — also by dedicated models of transcription initiation and polyadenylation.

## Introduction

Most newly sequenced eukaryotic genomes are annotated without matched transcriptomic data for every tissue and
developmental stage, so gene models come in large part from ab initio predictors that read the genome sequence
alone. For protein-coding genes these predictors have improved markedly: deep-learning models combined with
hidden-Markov decoding now recover most human coding genes at gene level and transfer across distant lineages
(Tiberius, Helixer, ANNEVO) [verify exact figures before citing]. Long non-coding RNA (lncRNA) genes are now
more numerous than protein-coding genes in the human reference (34,677 against 20,033 in GENCODE v50, section 1),
yet they are absent from the output of these tools.
Their annotation still rests on transcript evidence — short- and long-read RNA-seq, 5′-end and 3′-end
sequencing — which is tissue-limited and, for lncRNAs, low in expression and often specific to a few cell types
(Derrien et al., 2012). As a result, lncRNA annotation outside a handful of model organisms is sparse and uneven,
and even in human it is far from complete.

Two properties of lncRNA genes make sequence-only annotation hard to build and hard to evaluate. First, the
reference itself is moving. Large long-read capture projects have roughly doubled the human lncRNA catalogue in
recent years (Lagarde et al., 2017; Kaur et al., 2024); between GENCODE v41 (2022) and v50 the set of lncRNA loci
grew such that v41 recovers only 49% of v50 loci, against 99.8% for coding loci (section 1). A model trained on
any fixed annotation is told that unannotated lncRNA genes are background, and an evaluation against it counts
correct novel predictions as errors. Second, the sequence signal is weak: lncRNA exons lack the codon structure
that drives coding gene finders, and their promoters are more often enhancer-like and less CpG-rich than those of
coding genes (Hon et al., 2017; Mattioli et al., 2019). Consistent with this, the only public pipeline designed to
output lncRNA gene models end to end recovers about 3% of lncRNA loci on a held-out human chromosome (section 2).

DNA language models offer a different starting point: features learned from genomes at scale that separate
lncRNA exons from background better than chance (section 5), and that have been fine-tuned to segment transcripts
into exons, introns and untranslated regions (Shmelev et al., 2026) or to predict RNA-seq coverage across tissues
(AlphaGenome) [verify]. These approaches answer a narrower question than annotation does. Segmentation models are
evaluated inside given transcript intervals, and coverage predictors output tracks, not genes. Whether a
sequence-only model can find and delimit lncRNA genes along a whole chromosome — and, if it cannot, which part of
the task fails — has not been measured.

Here we address this in three steps. We first build an evaluation framework that treats the reference as
incomplete: stranded locus-, structure- and base-level measures, a time-split truth (v41 against v50), evidence
masks that mark likely unannotated transcription as unknown, and a one-to-one span criterion that asks whether a
prediction delimits a gene rather than merely touches it. We then build a sequence-only annotator that decodes a
joint grammar of coding and lncRNA genes over frozen DNA language-model features and that marginalises unknown
sequence during training. On chromosomes read once after all rules were frozen, it detects about 30% of lncRNA
genes at a strict precision of 0.42–0.50, more than ten times the recall of existing tools at the same precision,
and the detection transfers to mouse and to non-mammalian vertebrates. Under the one-to-one span criterion,
however, only about 7% of lncRNA genes are delimited, against 55% of coding genes. Finally, we decompose this gap.
Exon structure becomes partly recoverable once the ends of the transcription unit are known to within a kilobase;
the model's posteriors carry no lncRNA end signal, for reasons that lie neither in tiling, decoding nor unreliable
reference ends; and with these features — and equally with dedicated models of transcription initiation and
polyadenylation — the start of an lncRNA gene can be placed only where its promoter is a
CpG island, while lncRNA transcript ends cannot be placed at all — observations that we replicated on a held-out
chromosome. The ends of transcription units, not exon
structure, are the rate-limiting step of ab initio lncRNA annotation.

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
  some lncRNA loci as protein-coding genes (3, 8, 4 and 52 loci respectively; GENATATOR 22) [Fig. 2c].
* GENATATOR, the only pipeline designed to output lncRNA, recovers **3.5%** of lncRNA loci (locus-level,
  any exon overlap) and **0.7%** of lncRNA intron chains, while recovering 93% of mRNA loci and 65% of mRNA
  intron chains. With the leaderboard's own metric at its headline tolerance (k = 250 bp), lnc_RNA
  segmentation recall is **0.7%** (precision 3.7%), versus 49.6% for mRNA [Fig. 2b; `official_metric_lnc_vs_mrna.json`].
* Unstranded segmenters (NTv3, SegmentNT) mark 7–51% of lncRNA loci as transcribed but produce no correct
  lncRNA intron chain (≤ 0.5%), with locus precision 0.09–0.23.

Running the released GENATATOR pipeline ourselves on hg38 chr20 against GENCODE v50 gave the same picture
(lncRNA locus Sn 0.028 [95% CI 0.019–0.040], intron chain 0.012; mRNA locus Sn 0.967, intron chain 0.708;
34 lncRNA genes covered by a predicted mRNA that hits no reference mRNA, and 46 predicted mRNAs that fall only on
lncRNA loci) [`0.1_genatator/PANEL_SUMMARY.json`]. Under the stricter one-to-one span
criterion of section 12 (S50), its lncRNA recall is 0.006 and its mRNA recall 0.794 [`B1/B1.json`]. The bottleneck is interval discovery:
on chr21, only 33 of 590 lncRNA genes were touched by any predicted transcript [`RESULT_0.1a_0.3`]. A later
release of the same group's models (GENATATORs; Shmelev et al., 2026) fine-tunes DNA language models for
segmentation into exon, intron and UTR classes and reports that lncRNA segmentation is harder than mRNA; it scores
segmentation, not the end-to-end discovery and delimitation of loci measured here, and we did not run it.

For a representative coding annotator (Helixer, vertebrate model, hg38 chr20), 97% of lncRNA loci are not
touched at all, and **56% of its locus-level false positives fall on annotated lncRNA loci** (18 of 32)
[`0.3_fp_audit/helixer_chr20*.json`]. ANNEVO independently reports that 18.3% of its false positives overlap
annotated non-coding transcripts (different metric) [ANNEVO, Nat Methods 2026].

### 3. Transfer beyond the training clade collapses

GENATATOR is trained on 39 mammals. On four non-mammalian chromosomes (Ensembl 116), it emits almost no lncRNA
(0–13 transcripts; locus Sn ≤ 0.022), and its mRNA recovery falls with phylogenetic distance: locus Sn 0.84
(chicken chr5), 0.49 (zebrafish chr22), 0.21 (*Drosophila* 2L) and 0.006 (*C. elegans* chr I, 15 genes
predicted on the whole chromosome), against 0.97 on human chr20 [Fig. 3; `RESULT_0.1b_cross_lineage.md`].

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
weak for every backbone (AUPRC < 0.09) [Fig. 4; `RESULT_0.4_probe.md`]. A CDS-supervised backbone was *better*, not
worse, for lncRNA exons and splice sites than the self-supervised ones (NT-v2-500M, GENERanno-0.5b,
Caduceus-PS) — exon/splice features learned from coding supervision transfer to non-coding exons.

### 6. Treating evidence-supported unannotated sequence as unknown gives a small, consistent gain

In a pre-registered test (frozen features, minimal decoder), marginalising the labels of unannotated regions
that carry pre-2022 RNA-seq junction evidence improved recovery of loci first annotated after 2022 by 1.7
percentage points (95% CI 0.4–3.2; pre-declared threshold 5 points, not met) and, in a post hoc analysis, their
ranking (AUPRC +9–11% relative, all seeds positive) [Supplementary Fig. 1; `RESULT_0.5_unknown.md`]. We therefore use unknown-handling as a component and
report it as an ablation rather than as the main contribution.

### 7. A grammar-constrained annotator detects lncRNA-associated transcribed regions on held-out chromosomes

We built an annotator that labels every base of a chromosome with one of 17 states of a joint grammar for
protein-coding and lncRNA genes (intergenic, plus exon and three intron sub-states for each class on each strand),
using frozen features of a DNA language model (GENERanno-1.2b, CDS-annotation checkpoint), a light per-base head and
a hidden-Markov conditional random field with gene-entry and -exit bonuses (Methods). Unannotated sequence carrying
RNA-seq junction evidence is treated as unknown during training rather than as background (section 6). The model
was trained on 3,000 windows of 24.6 kb (2.4% of the genome) from human chromosomes other than chr18 and chr20; chr18 was used for
development only, and chr20 (model version v1.9) and chr9 (v2.0) were each read once after the decoding rules
had been frozen (Methods; `PHASE1_PROTOCOL.md` §16–17).

On chr20, the three v1.9 models recover 30.6% of the 955 GENCODE v50 lncRNA genes (exon overlap, L1) at a strict
precision of 0.422, where every prediction outside the reference counts as an error; GENATATOR recovers 2.8% at the
same precision (gain +0.278, 95% CI 0.256–0.300) [Fig. 5a; `RESULT_P15.md`]. Under the one-to-one span measure the
recall is 6.4% (precision 9.4%) against 0.6% for GENATATOR [`RESULT_B1.md`]. On the T2T-CHM13 chr20 benchmark
(RefSeq truth) the picture is the same: L1 recall 36.9% against 3.5%, S50 recall 8.4% against 1.4%, with lower
strict precision than on hg38 (L1 0.236; GENATATOR 0.259) [`P15/t2t_chr20_SCORE.json`; `B1/B1_MORE.json`].

On chr9, a second held-out chromosome that no model, score or RNA-seq evidence had touched before the read, the
calibrated v2.0 models reach L1 recall 30.4% at strict precision 0.498, and S50 recall 6.9% at precision 10.8%
[`RESULT_V20_CHR9.md`, `RESULT_B1.md`]. Calibration — a logistic model fitted on chr18 to the decoded genes'
exon count, coding-state posterior mass, open-reading-frame and hexamer features — raised strict precision by
0.041 (95% CI 0.028–0.055) at a cost of 0.014 in recall, and reduced the expected calibration error of the reported
gene confidence from 0.295 to 0.041. The gains on chr9 matched those predicted on chr18 (+0.041 and −0.014).

Other sequence-based models scored with the same measures recover fewer lncRNA regions or do so at much lower
precision [Supplementary Table 1; `RESULT_BASELINES.md`, `RESULT_S2_AG.md`]:

* SegmentNT, the only deep-learning annotator with an lncRNA class that was not trained on chr20: L1 recall
  10.8% at strict precision 0.113 (S50 0.8%); 97% of its lncRNA calls carry no strand.
* AlphaGenome, whose predicted RNA-seq coverage we segmented with a threshold chosen on chr18: recall 52.1% at
  strict precision 0.097, with 10,003 calls (16 times ours). No threshold reached a strict precision of 0.43 on
  chr18 (maximum 0.090). These chromosomes are in AlphaGenome's training data, so the comparison favours it.
* Tiberius, Helixer and Vipsania output no lncRNA.

Our models report fewer lncRNA genes as protein-coding than any coding annotator: 10.0 lncRNA genes covered by a
predicted mRNA that hits no reference mRNA, against 15 for Tiberius, 14 for Helixer and 34 for GENATATOR. On the
coding side, however, our gene structures are far weaker than those of dedicated tools: mRNA L1 recall 0.861
(Tiberius 0.964), mRNA intron-chain recall 0.054 (Tiberius 0.557) [`RESULT_BASELINES.md`]. We therefore do not
claim competitive coding-gene annotation; in the unified output (Methods) coding structures are taken from
Tiberius.

Recall differs between the three models trained with different seeds (chr20 L1 0.214–0.402; S50 0.041–0.102),
mostly because they report different numbers of lncRNA genes (500–891); strict precision is stable (0.41–0.43).
Detection is higher for longer lncRNA genes and genes with more exons (single-exon genes 7.8%, five or more exons
64.0%), and nearly as high for genes first annotated after GENCODE v41 (26.1%) and for genes without detectable
homology to any training window (29.2%) as for all genes [`RESULT_P15.md` §5].

### 8. Each component of the grammar contributes; a richer grammar does not help

Pre-registered ablations, each trained with three seeds and read once on chr20, all lowered the area under the
precision–recall curve for lncRNA genes relative to the same-seed full model [Supplementary Table 2; `RESULT_P15.md` §7]:
removing the gene-entry and -exit bonuses (−0.056, three of three seeds lower), replacing unknown handling by
treating all unannotated sequence as background (−0.039, 3/3), and removing the coding-gene pathway so that the
model learns lncRNA genes alone (−0.061, 3/3). The last result means that learning what a coding gene is helps
the model delimit non-coding genes. Decoding the same posteriors without the grammar (per-base argmax) also lowered
it (−0.022, 3/3). Training on the 2022 annotation (GENCODE v41) instead of v50 lowered it in two seeds and raised it
in one (+0.040), so we cannot claim that the newest annotation is required. An 81-state grammar with explicit
duration distributions was much worse (0.010 against 0.145 for the same seed). The signs of all twelve ablation
runs agreed between the development chromosome and chr20.

### 9. Detection transfers to other species, but adding mouse training does not improve it

Without any training or tuning on other species, the human-trained models detect lncRNA regions on mouse chr19
(GENCODE vM39; L1 recall 28.0% at strict precision 0.372, S50 6.7%), close to the held-out human chromosome
[`RESULT_S5.md`; `B1/B1_MORE.json`]. On four non-mammalian chromosomes with sparser annotation, L1 recall is 31.8%
(chicken chr5), 34.8% (zebrafish chr22), 21.4% (*Drosophila* 2L) and 7.3% (*C. elegans* chr I), where GENATATOR
recovers 0–2.2% [Fig. 5a; `RESULT_P16.md`]. Strict precision is low on these chromosomes (0.246, 0.065, 0.089 and
0.008), partly because few lncRNA genes are annotated there (135 on zebrafish chr22, 150 on worm chr I); on
*C. elegans* the failure is real, as half to four-fifths of the lncRNA calls fall on coding genes. The coding
pathway transfers too (mRNA L1 recall 0.86, 0.84, 0.64 and 0.68, against 0.84, 0.49, 0.21 and 0.006 for GENATATOR).

In a pre-registered paired test, adding 3,000 mouse training windows to the human windows did not improve
zero-shot detection on six held-out chromosomes of chicken, zebrafish, fly, worm, rat and cow (mean recall
difference against a human-only control with the same number of windows +0.040, with the confidence interval
excluding zero on one of six sets and 1.49 times as many calls; against the original human model −0.016)
[`RESULT_V30.md`]. Doubling the training set also made training less stable (three training runs failed a pre-registered
divergence criterion).

### 10. Predicted regions outside the reference are enriched for transcription, but rarely match transcript models

The strict precision counts every prediction outside GENCODE v50 as wrong. We asked, in a pre-registered analysis
with public data (FANTOM5 CAGE, PolyASite 2.0, Snaptron GTEx junctions, CHESS 3 and GTEx long-read transcripts),
how many of these chr20 predictions carry evidence of transcription on the predicted strand, against position-
matched random intervals (100 placements per prediction) [Fig. 5b; `RESULT_BIOVAL.md`].

Of the predictions outside the reference, 59.4% carry such evidence against 29.4% of matched intervals
(P < 0.01 for each of the three models). The cleanest set is predictions with no annotated gene on either strand
within 5 kb: 48.1% against 3.7% (13-fold); swapping the strand lowers support to 32.7%. Predictions made by at
least two of the three models (111 regions) reach 82.9% against 54.5%. About half of the predictions are already
in LncBook 2.1 or RNAcentral.

The support is almost entirely transcription inside the predicted region (spliced RNA-seq reads), not the
predicted structure: 1.1% of predicted introns match an observed junction exactly, no predicted intron chain
matches a long-read transcript, and CAGE peaks lie within 500 bp of 4.7% of predicted starts (0.8% for matched
intervals). The strict precision therefore underestimates how often the model points to transcribed sequence,
but the evidence supports transcribed regions, not new gene models.

### 11. Detection follows CpG-island promoters and long transcription units

Which lncRNA genes are found? In a multivariable analysis of all lncRNA genes on chr18 and chr20, and in a pre-registered
replication on chr9, detection was associated with promoter CpG observed/expected ratio
(standardised coefficients +0.35, +0.34 and +0.32) and with gene span (+0.55, +0.70 and +0.58), but not consistently, after
these, with CAGE support at the start, repeat content or exon number [Fig. 5c; `RESULT_S3.md`,
`RESULT_V20_CHR9.md`]. Detection rates rose from 22–24% in the lowest to 41–42% in the highest CpG quartile and
from 12–13% to 48–54% across span quartiles. In post hoc analyses, lncRNA genes antisense to a coding gene were
harder to detect (−0.26 and −0.29), and lncRNA genes whose start lies in a transposable element (about a quarter)
were harder to detect (detection rate −0.064 on chr20 and −0.067 on chr9), mainly because they lack CpG-island promoters
[`RESULT_S3_TE.md`]. The model thus locates lncRNA-associated regions mainly through gene-like promoters and long
transcription units rather than through exon sequence; sections 12–15 show that the same dependence on CpG-island
promoters governs where the start of an lncRNA gene can be placed.

### 12. Ab initio lncRNA predictions overlap annotated loci but rarely delimit them

Under the locus-level measure used so far (any same-strand exon overlap; Methods), our sequence-only annotator
recovers 30.6% of the 955 GENCODE v50 lncRNA genes on held-out chr20 at a strict precision of 0.42, eleven times
the recall of GENATATOR (2.8%) at the same precision (section 7). This measure does not ask whether a prediction
delimits a gene: one long prediction can touch several genes, and a prediction that covers a small part of a gene
counts as a hit. We therefore also scored every prediction set with a one-to-one span criterion: a prediction and
a reference gene of the same class and strand are matched when each covers at least half of the other's genomic
span, and a maximum bipartite matching makes the pairs unique (S50; Methods).

Under S50 the lncRNA recall of our models falls to **6.4%** on chr20 (three models 4.1%, 10.2%, 5.0%; precision
9.4%) and **6.9%** on chr9 with the calibrated v2.0 models (precision 10.8%) [Fig. 6a; `RESULT_B1.md`,
`B1/B1.json`]. The relative advantage over existing tools remains (GENATATOR 0.6%, SegmentNT 0.8%), but in
absolute terms only about one annotated lncRNA gene in fifteen is delimited by a single prediction. The criterion
is not unattainable: GENCODE v41, scored as a prediction of v50, keeps 39.6% of lncRNA genes under S50 against
47.1% under exon overlap, and our own mRNA predictions keep 55% (GENATATOR mRNA 79%, Tiberius 65%).

The loss comes from the ends. On chr20 the median distance between predicted and annotated lncRNA starts is
7.7–9.5 kb and between ends 9.3–21 kb (mRNA: 1.0–1.9 kb and 1.6–2.4 kb for the same models; GENATATOR mRNA:
79 bp and 21 bp). The models err in different directions: among annotated lncRNA genes hit by a prediction on
chr18, the best-overlapping prediction is too long for 68% with one model and too short for 56% with another
[`RESULT_BND.md`, D1], and on chr20 8–30% of hitting predictions span two or more annotated lncRNA genes.

*We therefore describe the locus-level result as the detection of lncRNA-associated transcribed regions rather
than as the discovery of lncRNA genes, and report both measures throughout.*

### 13. Known ends make lncRNA structure partly recoverable, but only with sub-kilobase precision

Is exon structure the harder problem once the ends are known? Two experiments on the development chromosome
(chr18) say no. First, re-decoding our posteriors with the annotated gene span imposed gives exact lncRNA intron
chains for 0.3% of genes, with the annotated splice sites imposed for 16.8%, and with both for 34.1% (three-model
means) [`RESULT_S1B.md`]. Second, GENATATOR's segmentation model — designed to segment a given transcript
interval — reconstructs the exact intron chain of 37% of 200 multi-exon lncRNA genes when given the annotated
transcript interval and 35.5% when given the gene span [Fig. 6b; `RESULT_S11.md`]. (chr18 is part of GENATATOR's
training data, so this is an upper bound.)

The same model collapses when the interval is imprecise: padding both ends by 1 kb lowers exact chains to 3%,
padding by 5 kb or more to 0%, and trimming 1 kb from each end to 1%. Junction precision falls from 0.49 to 0.12
(+1 kb) and 0.01 (+10 kb). The tolerance of interval segmentation (below 1 kb) is an order of magnitude tighter
than the end error of ab initio lncRNA locus prediction (section 12). Exon structure is therefore not the first
obstacle to end-to-end lncRNA annotation; the ends of the transcription unit are. A recent benchmark of fine-tuned
DNA language models for within-interval segmentation lists gene-boundary identification among the problems that
remain open for all tools (Shmelev et al., 2026), but scores neither transcript starts nor ends; the experiments
here measure how precise the ends must be and how far ab initio predictions are from that precision.

### 14. The posteriors carry no lncRNA end signal, and the cause is not in decoding

We next asked where in our pipeline the end information is lost (chr18, three v1.9 models)
[`RESULT_BND.md`], and repeated every read-out on held-out chr20 under pre-registered replication criteria
[`PLAN_REP20.md`, `RESULT_REP20.md`]. Three explanations that would not require retraining were ruled out on both
chromosomes:

* *Tiling.* Posteriors are computed on 24.6-kb tiles; predicted lncRNA ends fall within 100 bp of a tile seam
  2.1% of the time against 1.0% expected by chance (chr20: 2.3%).
* *Decoding.* Decoding the same posteriors by per-base argmax instead of maximum expected accuracy does not help
  (S50 recall 4.3% against 5.4%; chr20: 5.6% against 6.4%).
* *Unreliable reference ends.* Annotated lncRNA ends are known to be often incomplete and are revised as
  full-length and 5′-end data accumulate (Hon et al., 2017; Lagarde et al., 2017; Kaur et al., 2024). Only 16.8% of chr18 lncRNA genes have a FANTOM5 CAGE peak within 200 bp of their
  outermost start and 32.7% a PolyASite cluster within 100 bp of their end, but genes with both kinds of support
  are delimited only slightly better than genes with neither (S50 6.4% against 5.1%; on chr20, 10.9% of 92 genes
  against 5.7% of 466, below the pre-registered two-fold threshold). Even fully supported reference ends leave
  about nine in ten lncRNA genes undelimited.

What remains is the posterior itself. Taking the summed gene-state posterior of the gene's class and strand, the
annotated start of an mRNA is a step of ~0.4 and the onset of the posterior lies within 1 kb of it for 70% of
detected genes; for lncRNA the step is ~0.05 and the onset lies within 1 kb for 23% (chr20: 25% against 69% for
mRNA). The model knows that a region
belongs to an lncRNA gene but not where the gene begins. Predicted lncRNA starts are nonetheless promoter-like:
8–11% lie within 200 bp of a same-strand CAGE peak, against 0.3–0.6% after shifting them by 50 kb.

### 15. With current sequence models, only CpG-island lncRNA promoters can be located, and lncRNA ends cannot

The boundary component of our model is a single set of start / end logits shared by mRNA and lncRNA genes, trained
on the canonical-transcript ends of both classes. In the training windows, 93% of mRNA start labels but only 31%
of lncRNA start labels have a CAGE peak within 200 bp (polyA support at ends: 89% and 38%); genome-wide, both ends
of the canonical transcript are supported for 82% of mRNA genes and 10.5% of lncRNA genes [`PLAN_END.md` §0].
We tested two remedies with a pre-registered probe on frozen features (chr18; Methods): class-specific end
channels trained on the existing labels (P-ann), and the same channels trained only on evidence-supported ends,
placed at the CAGE or polyA position, with unsupported ends ignored (P-ev). The read-out is the share of reference
starts whose highest-scoring position within ±20 kb lies within 1 kb (5% expected by chance) [`RESULT_END.md`].

For CAGE-supported lncRNA starts, the existing shared logits reach 15.4% (three models). Class-specific channels
reach 40.5% and 21.4% (two seeds), evidence-filtered labels 15.6% and 26.0%; for mRNA starts the probes reach
41–79%. None of the pre-registered outcomes (labels limiting, channel sharing limiting, or not locatable) was met
after the second seed, so we draw no conclusion on the two remedies.

Two patterns, however, held in all seven models (three full models, four probes) on chr18. First, lncRNA starts are
locatable only where the promoter is CpG-rich: across quartiles of promoter CpG observed/expected ratio, the
share within 1 kb stays near chance in the two lowest quartiles (2.6–8.6% in the lowest) and rises to 17–61% in
the highest; in every model the third quartile exceeds the two lowest and the highest exceeds the third [Fig. 6c].
On held-out chr20, read once with the same models, the ordering held in six of seven models (two lowest quartiles
5–10%, highest 16–54%), meeting the pre-registered replication criterion [`RESULT_REP20.md`]. This mirrors what drives detection itself — promoter CpG and transcription-unit length are the
features associated with a lncRNA locus being found (section 11; `RESULT_S3.md`), and lncRNA genes whose start
lies in a transposable element are harder to find mainly because they lack CpG-island promoters
(`RESULT_S3_TE.md`). The dependence on CpG islands is not new in itself: CpG islands have long been the dominant
signal in promoter prediction, and promoters without them were suggested to be unpredictable from local sequence
(Hannenhalli and Levy, 2001). What it implies for lncRNA follows from their origin: most intergenic lncRNAs start
at enhancers rather than at promoters (Hon et al., 2017), and lincRNA core promoters carry simpler motif profiles
than mRNA promoters (Mattioli et al., 2019). Our data quantify this for a genome-scale, sequence-only annotator:
the lncRNA starts that lack a CpG island — most of them — are placed no better than chance.

A dedicated model of transcription initiation shows the same limit. We scored Puffin (Dudnyk et al., 2024), which
predicts CAGE and related initiation signals at base resolution, with the same read-out on its own held-out
chromosomes, chr8 and chr9, under pre-registered criteria [`PLAN_DEDEND.md`, `RESULT_DEDEND.md`]. Puffin placed
71% of CAGE-supported mRNA starts within 1 kb, but only 8.6% of lncRNA starts in the two lowest promoter-CpG
quartiles (7.9% and 9.4%), against 21% in the third and 56% in the highest quartile; for all lncRNA starts it
reached 24% and for CAGE-supported ones 40% [Fig. 6c]. On chr18 and chr20, which were in Puffin's training data,
the gradient was the same (lowest quartile 3–9%, highest 49–51%). The CpG-island dependence is therefore not a
weakness of general-purpose language-model features: a model trained directly on initiation data cannot place
low-CpG lncRNA starts either.

Second, with these features, lncRNA transcript ends are not locatable: 4.5–8.7% of polyA-supported lncRNA ends
on chr18 and 4.7–9.9% on chr20 lie within 1 kb of the top-scoring position, in every model. On chr18 mRNA ends
were equally unlocatable (2.2–12.6%), but on chr20 two probes placed 17–18% of them within 1 kb, so the pre-registered
replication criterion, which required both classes to stay at or below 15%, was not met; for mRNA we therefore
report only that the end signal in the features is weak. The mRNA ends that our full model places correctly
(28–41% of detected genes within 1 kb on chr18, 49–58% on chr20) are thus set mostly by the gene grammar (exon
structure, coding frame). A dedicated model of 3′-end processing does not change this. PolyaID and PolyaStrength (Stroup et al., 2023)
identify cleavage and polyadenylation sites at nucleotide resolution; scored at every position within ±20 kb of
each gene end on the gene's strand (pre-registered: the strongest site among those classified as polyadenylation
sites), they placed 7.1% (chr18) and 11.9% (chr20) of polyA-supported lncRNA ends within 1 kb, and only 12% and
19% of polyA-supported mRNA ends [Fig. 6d; `RESULT_DEDEND.md`]. About 5% of all positions in these windows were
classified as cleavage-competent, so recognising polyadenylation signals does not identify which of them ends a
gene. Both chromosomes may have been in PolyaID's training data, which would favour it.

*Taken together, sections 12–15 identify the localisation of transcription-unit ends as the rate-limiting step of
end-to-end lncRNA annotation: structure becomes partly recoverable once the ends are known to within a kilobase,
but with current sequence models — general DNA language-model features as well as dedicated models of
initiation and polyadenylation — the start can be placed only for CpG-island promoters and the lncRNA end not at
all. Each of these observations except the mRNA end result was replicated on held-out chr20.*

## Discussion

**Detection is not delimitation.** The central observation of this work is a gap of about five-fold between two
measures that the field often treats as one. Under exon overlap, the usual locus-level measure, our annotator
recovers 31% of lncRNA genes on a held-out chromosome; under a one-to-one span criterion, 6–7%. The same criterion
leaves 55% of coding genes and 40% of lncRNA genes when one GENCODE release is scored against the next, so it is
not unattainable. Predictions outside the reference are enriched thirteen-fold for transcription in intergenic
sequence (section 10), so the model is pointing at transcribed sequence; it is the boundaries of the
transcription unit that it gets wrong, by a median of 8–21 kb against 1–2 kb for coding genes. We therefore
describe the result as the detection of lncRNA-associated transcribed regions, and we suggest that ab initio
lncRNA annotators report a span-level criterion alongside exon overlap; under exon overlap alone, a model that
paints long transcribed stretches can look like a gene finder.

**Why the ends.** The decomposition places the failure at the ends rather than in exon structure, and the
pattern of what can and cannot be located matches lncRNA biology. Starts are locatable at CpG-island promoters
and near chance elsewhere. This echoes the long-standing observation that CpG islands dominate sequence-based
promoter prediction and that promoters without them are hard to find from local sequence (Hannenhalli and Levy,
2001), and it matters more for lncRNA than for coding genes because most intergenic lncRNAs initiate at
enhancer-like elements (Hon et al., 2017) whose core promoters have simpler motif content (Mattioli et al., 2019).
lncRNA transcript ends were not locatable with our features, and the end signal for coding genes was weak; the
coding ends our full model places correctly are set mostly by the gene grammar — exon structure and reading
frame — which lncRNA genes lack. lncRNA
transcripts are also less efficiently spliced and, for a subset, terminated in ways less tied to canonical
polyadenylation than mRNAs (Mele et al., 2017; Schlackow et al., 2017) [verify], which would leave even less for a
sequence model to anchor an end on. Detection itself follows the same features: promoter CpG content and the
length of the transcription unit, not exon sequence, are what distinguish the lncRNA genes the model finds
(section 11).

**What would close the gap.** Three directions follow. First, end evidence, where it exists, should enter the
decoder directly: CAGE or long-read 5′ ends and polyA-site data constrain exactly the quantity that sequence
fails to provide, and the grammar can take them as position-dependent bonuses without retraining (Methods,
boundary logits). Second, for genomes without such data, plugging in dedicated sequence models is not enough on its own. A model
of transcription initiation (Puffin; Dudnyk et al., 2024) placed CpG-island lncRNA starts well on its held-out
chromosomes but low-CpG starts at near chance, and a polyadenylation model (PolyaID; Stroup et al., 2023) found
cleavage-competent sites throughout each 40-kb window and picked the gene end for about one lncRNA in ten
(section 15). Initiation models remain useful as start scores for CpG-island genes; for the rest, the missing
information is which of many candidate sites is used, which is a property of the cell rather than of the local
sequence alone, and points to expression data or models conditioned on cell type. Third, the tolerance of
structure recovery to end error is now measured (below 1 kb for exact intron chains, section 13), which gives a
concrete target: an end predictor that reaches sub-kilobase precision for most lncRNA genes would make existing
segmentation models usable end to end.

**Relation to other approaches.** Fine-tuned segmentation models (Shmelev et al., 2026) report that lncRNA
segmentation is harder than coding segmentation and list gene-boundary identification among open problems; our
results quantify that statement for end-to-end annotation and show that, given correct ends to within a kilobase,
such a model already recovers a third of exact lncRNA intron chains. Coverage predictors (AlphaGenome) recover
more lncRNA regions than our annotator when their tracks are segmented, but at a strict precision below 0.1 even
on chromosomes they were trained on. Coding annotators (Tiberius, Helixer) remain far better on coding gene
structure than our joint model, and some of their false positives are lncRNA loci reported as coding genes
(section 2). A practical annotation therefore combines a dedicated coding annotator with a non-coding layer like
ours, and the joint grammar remains useful in that combination: removing the coding pathway lowered lncRNA
performance in all seeds (section 8), and our model reports fewer lncRNA loci as coding than any coding annotator.
The annotation we release follows the results rather than the format: coding genes from Tiberius, lncRNA genes as
loci with calibrated confidence and a start class, and no lncRNA transcript models, because fewer than 2% of the
decoded lncRNA introns matched an annotated intron exactly (Methods).

**Limitations.** (i) We used a single, frozen backbone (GENERanno-1.2b); a backbone fine-tuned end to end or with
a longer receptive field might carry more end information, although the end-localisation probe on frozen features
did not suggest that labels or class sharing were the main limit. (ii) The end-localisation analyses were developed on chr18; the posterior and probe read-outs (sections 14–15)
were then repeated once on held-out chr20 under pre-registered criteria, but the structure-tolerance experiment
(section 13) was not, and the reference-end control on chr20 came close to its threshold (two-fold).
(iii) Models trained with different seeds differ substantially in recall (chr20 L1 0.21–0.40), mostly through
the number of genes reported; strict precision is stable. (iv) Every lncRNA reference used here is incomplete,
and FANTOM5 CAGE, used to identify supported ends, favours broadly and highly expressed transcripts; the
association between CAGE support and promoter CpG was not separated. (v) Evidence for predictions outside the
reference shows transcription in the predicted region, not correct gene models: no predicted intron chain matched
a long-read transcript. (vi) Transfer fails on *C. elegans*, and adding mouse training did not improve transfer
to other species.

**Conclusion.** Sequence-only models can now find a substantial fraction of lncRNA-associated transcription on
unseen chromosomes and species, at a precision that existing tools do not reach. What they cannot yet do is say
where an lncRNA gene begins and ends, and this, rather than exon structure, is what separates detection from
annotation. Making the ends of transcription units predictable — from initiation and termination models, or from
the end evidence that a decoder can take directly — is the step on which ab initio lncRNA annotation now depends.

## Methods

### Evaluation framework

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
true positives under a newer release. Confidence intervals come from bootstrap resamples of genes (1,000 for the benchmark, 2,000 for the held-out reads
of sections 7–15).
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

### Annotator, held-out reads and unified annotation

**Model.** Features: hidden states of layers 17 and 23 of GENERanno-1.2b (CDS-annotation checkpoint), one 4,096-d
vector per 6-mer token, computed per 24,576-bp window and kept frozen. Head: linear projection to 256 dimensions,
bidirectional GRU over tokens, repeated to bases and fused with a two-layer convolution over one-hot DNA, giving
per-base potentials for 17 states and four boundary logits (start / end per strand). Grammar: a first-order
hidden-Markov CRF over the 17 states with motif gates at intron borders, minimum lengths (exon 20, intron 50,
lncRNA 200 nt) and the boundary logits as position-dependent bonuses on gene-entry and -exit transitions
[`PHASE1_PROTOCOL.md` §13–16].

**Training.** 3,000 windows (lncRNA-rich, coding-rich and random) from human chromosomes other than chr18 and chr20;
for v2.0 the 154 windows on chr9 were replaced by windows on other chromosomes before chr9 was read. Labels: GENCODE v50 Ensembl_canonical lncRNA and protein-coding transcripts, repaired to grammar-feasible
paths; unannotated bases within ±500 bp of Snaptron junctions with ≥ 100 samples and ≥ 500 reads are unknown
(marginalised) rather than intergenic. Loss: CRF negative log-likelihood plus binary cross-entropy for the boundary
logits; AdamW, 12 epochs; the model is the uniform average of the weights of epochs 10–12. Three seeds per arm.

**Decoding and calibration.** Exact float64 forward-backward posteriors on overlapping tiles; exon-level maximum
expected accuracy decoding (γ = 3) on the posteriors; lncRNA genes kept when their mean posterior is ≥ 0.6 (v1.9).
v2.0 adds a gene-level logistic calibration fitted on chr18 (features:
posterior confidence, exon count, coding-state posterior mass over the exons, log ORF length, ORF fraction, hexamer
and Fickett scores, log length), with a threshold of 0.36 on the calibrated probability chosen on chr18 by a
pre-registered rule.

**Pre-registration and single reads.** Every decoding and selection rule was frozen in the protocol before the
corresponding test read. Before each test read, the chr18 posteriors were recomputed and required to match the
development posteriors byte for byte (sha256). chr20, T2T chr20, chr9, mouse chr19 and the panel chromosomes were
each read once per model; no parameter was tuned on them.

**Strict span matching** as in sections 12–15 (Methods). For GENATATOR, transcripts of one predicted gene are
merged before span matching; its L1 precision in the text is the published per-transcript value.

**Transcription evidence for predictions outside the reference.** As pre-registered in `PREREG_BIOVAL.md`: a
prediction is supported if the predicted region contains a same-strand spliced junction (Snaptron GTEx, ≥ 5 samples,
≥ 20 reads, not a v50 intron), or overlaps a non-v50 long-read / assembled transcript, or has both a CAGE peak
within 500 bp of its start and a polyA cluster within 500 bp of its end. Controls: 100 placements per prediction at
the same length and strand on the same chromosome, outside v50 genes as the prediction.

**Unified annotation.** For each chromosome we write one GFF3 with three layers [`RESULT_UGFF.md`]. Coding genes are
Tiberius predictions (mammalia_nosofttmasking_v2; CDS only). lncRNA genes are those of the shipped model (seed 42;
chr20 v1.9, chr18 and chr9 v2.0) that pass the reporting threshold and do not overlap a Tiberius CDS on the same
strand; this guardrail, fixed in advance, raised strict precision by 0–0.04 and lowered recall by 0.01–0.05 (three
chromosomes, three seeds), and the values in the text are for the predictions without it. A lncRNA transcript model
is written only when every predicted intron has a mean intron posterior above a threshold chosen on chr18 as the
lowest value giving an exact-intron precision of at least 0.5. No threshold reached it (at 0.95 the precision was
0.04; with all decoded introns written it was 0.003–0.015 on chr18 and below 0.01 on chr20 and chr9), so every
lncRNA gene is reported at locus level, with its decoded exons as an attribute and not as a transcript model. Each
lncRNA gene carries the CpG observed/expected ratio of ±500 bp around its predicted start and a class (CpG-rich:
above the upper quartile of annotated lncRNA starts on chr18, 0.46): on chr18, 73% of CpG-rich predicted starts of
predictions overlapping an annotated lncRNA gene lay within 1 kb of an annotated start of that gene, against 29% of
the others. Small structured RNAs (Rfam, tRNAscan-SE; Phase 0 integration rules) are included for chr20.

### End analyses

**One-to-one span matching (S50 / S25).** Gene span = leftmost to rightmost exon boundary over all transcripts
(reference) or of the prediction (transcripts of one predicted gene merged). Candidate pairs: same class, same
strand (unstranded predictions may pair with either strand), |p ∩ g| / |g| ≥ θ and |p ∩ g| / |p| ≥ θ, θ = 0.5
(S50, primary) or 0.25 (S25). Pairs are made one-to-one by maximum-cardinality bipartite matching (augmenting
paths; checked against exhaustive search on 3,000 random graphs). Recall = matched / reference genes, precision
= matched / predictions. The wording rule (locus discovery / candidate locus discovery / transcribed regions)
was fixed before scoring [`PLAN_B1.md`, `scripts/phase1/b1/b1_strict.py`].

**Interval jitter of segmentation.** 200 chr18 GENCODE v50 lncRNA genes with a multi-exon transcript and a gene
span ≤ 60 kb (numpy rng 0). GENATATOR's segmentation model (genatator-caduceus-ps-multispecies-segmentation,
revision 2716f2b2) was run as in its released pipeline (reverse-complement averaging, 250-kb context, sigmoid,
exon-versus-intron decision, splice filter, rejection of intervals that start or end in intron) on the longest
multi-exon transcript interval, the gene span, the interval padded by 1, 5, 10, 20 kb or trimmed by 1, 5 kb at
both ends. Exact chain = predicted intron chain identical to that of any reference transcript of the gene;
rejected or single-exon outputs count as failures. Stop rule fixed in advance [`PLAN_S11.md`].

**End diagnostics.** Tile seams at k × 20,472 + 2,052 bp (k ≥ 1); error types from the best-overlapping prediction
of each hit gene (both span coverages ≥ 0.5 / prediction too long / too short / shifted); posterior end signal
from the summed posterior of the four gene states of the gene's class and strand (step = mean over 1 kb inside
minus 1 kb outside; onset = nearest crossing of 0.5 within ±50 kb, for genes with mean posterior ≥ 0.5 over the
span) [`PLAN_BND.md`].

**End probes.** Hidden layers of the annotation head (projection, bidirectional GRU, DNA convolution) with eight
end logits ({lncRNA, mRNA} × {start, end} × strand) on frozen GENERanno-1.2b features of the 3,000 v2.0 training
windows; unweighted binary cross-entropy, AdamW (lr 1e-3, weight decay 1e-2, cosine schedule), batch 8, four
epochs, last epoch used (validation windows lie on chr18). P-ann positives: canonical-transcript ends ± 10 bp;
ignored: non-canonical ends ± 50 bp and evidence-masked unannotated sequence. P-ev: canonical ends with a
same-strand CAGE peak (start, ± 200 bp) or PolyASite cluster (end, ± 100 bp) are positives at the evidence
position ± 10 bp; unsupported canonical ends are ignored ± 1 kb. Read-out on chr18 canonical transcripts: Loc1k
= share of ends whose argmax within ± 20 kb on the class- and strand-specific channel lies within 1 kb; strata by
evidence and by promoter CpG o/e (TSS ± 500 bp) quartile. Pre-registered decision rules and the grey-zone rule
(second seed) in [`PLAN_END.md` §4].

**Dedicated end models.** Puffin (github.com/jzhoulab/puffin, commit 9bc81eb; FANTOM_CAGE output, + and − strand
channels; forward pass re-implemented without selene_sdk, weights loaded with all parameters matched) over whole
chromosomes in 200-kb pieces with 325-bp context; PolyaID and PolyaStrength (github.com/zhejilab/PolyaModelsHuman,
commit d70bcd0; published weights, sha256 checked) on 240-nt windows centred on every position within ± 20 kb of
a canonical end, on the gene's strand. Score for an end: Puffin argmax; PolyaID — the maximum PolyaStrength among
positions with PolyaID ≥ 0.75, or the PolyaID maximum if there is none. Truth, evidence, CpG strata and Loc1k as
for the probes; CpG quartiles per chromosome group. Puffin was read on its test chromosomes chr8 and chr9 (primary)
and on chr18 and chr20 (in its training set); PolyaID on chr18 and chr20. Decision rules in [`PLAN_DEDEND.md` §4].

## Figure legends

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

**Supplementary Figure 1 | Treating evidence-supported unannotated sequence as unknown.** Pre-registered test
on frozen features with a minimal decoder, three seeds per arm (dots; short lines, mean): unannotated sequence treated
as background (F), as unknown where pre-2022 RNA-seq junctions support it (U-ev), or everywhere (U-pu). **a**,
Recovery of lncRNA loci first annotated after GENCODE v41 (n = 473, chr20); **b**, AUPRC for exons of those loci
(exploratory, post hoc); **c**, strict precision against v50.

## References (to be completed and formatted)

* Bogard N et al. A deep neural network for predicting and engineering alternative polyadenylation. *Cell* 2019.
* Dudnyk K et al. Sequence basis of transcription initiation in the human genome. *Science* 2024.
* Hannenhalli S, Levy S. Promoter prediction in the human genome. *Bioinformatics* 2001.
* Hon CC et al. An atlas of human long non-coding RNAs with accurate 5′ ends. *Nature* 2017.
* Kaur G, Perteghella T et al. GENCODE: massively expanding the lncRNA catalog through capture long-read RNA
  sequencing. *bioRxiv* 2024, doi:10.1101/2024.10.29.620654.
* Lagarde J et al. High-throughput annotation of full-length long noncoding RNAs with capture long-read
  sequencing. *Nat Genet* 2017;49:1731, doi:10.1038/ng.3988.
* Mattioli K et al. High-throughput functional analysis of lncRNA core promoters elucidates rules governing
  tissue specificity. *Genome Res* 2019.
* Shmelev A et al. GENATATORs: ab initio gene annotation with DNA language models. *bioRxiv* 2026,
  doi:10.64898/2026.06.17.732686.
* Stroup EK et al. Deep learning of human polyadenylation sites at nucleotide resolution. *Nat Commun* 2023.
* Derrien T et al. The GENCODE v7 catalog of human long noncoding RNAs. *Genome Res* 2012.
* Dudnyk K et al. Sequence basis of transcription initiation in the human genome. *Science* 2024.
* Cochran K et al. Dissecting the cis-regulatory syntax of transcription initiation with deep learning. *bioRxiv* 2024.
* Mele M et al. Chromatin environment, transcriptional regulation, and splicing distinguish lincRNAs and mRNAs.
  *Genome Res* 2017 [verify].
* Schlackow M et al. Distinctive patterns of transcription and RNA processing for human lincRNAs. *Mol Cell* 2017
  [verify].
* Tiberius (Gabriel et al., *Bioinformatics* 2024), Helixer (Holst et al., *Nat Methods* 2025), ANNEVO
  (*Nat Methods* 2026), SegmentNT, AlphaGenome, GENATATOR (2025) — citation details to fill from the benchmark draft.
