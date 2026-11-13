# Draft — results sections: sequence-only detection of lncRNA-associated transcribed regions

> Status: first draft, 2026-11-14. Sections 7–11, between the benchmark (`DRAFT_results_benchmark.md`, sections 1–6)
> and end localisation (`DRAFT_results_ends.md`, sections 12–15). Every number points to a result document in
> `docs/phase1_results/` and its JSON. Not yet checked by a second reader; figures planned, not drawn.
> Wording follows the pre-registered rule of `PLAN_B1.md`: locus-level results are reported under both the
> exon-overlap measure (L1) and the one-to-one span measure (S50), and are described as the detection of
> lncRNA-associated transcribed regions, not as the discovery of lncRNA genes.
> Language: English (manuscript); notes for the authors in Chinese at the end.

## Results

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
same precision (gain +0.278, 95% CI 0.256–0.300) [Fig. 6a; `RESULT_P15.md`]. Under the one-to-one span measure the
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
precision [Fig. 6b; `RESULT_BASELINES.md`, `RESULT_S2_AG.md`]:

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
claim competitive coding-gene annotation; in the unified output (section 11, Methods) coding structures are best
taken from a dedicated annotator.

Recall differs between the three models trained with different seeds (chr20 L1 0.214–0.402; S50 0.041–0.102),
mostly because they report different numbers of lncRNA genes (500–891); strict precision is stable (0.41–0.43).
Detection is higher for longer lncRNA genes and genes with more exons (single-exon genes 7.8%, five or more exons
64.0%), and nearly as high for genes first annotated after GENCODE v41 (26.1%) and for genes without detectable
homology to any training window (29.2%) as for all genes [`RESULT_P15.md` §5].

### 8. Each component of the grammar contributes; a richer grammar does not help

Pre-registered ablations, each trained with three seeds and read once on chr20, all lowered the area under the
precision–recall curve for lncRNA genes relative to the same-seed full model [Fig. 6c; `RESULT_P15.md` §7]:
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
recovers 0–2.2% [Fig. 6d; `RESULT_P16.md`]. Strict precision is low on these chromosomes (0.246, 0.065, 0.089 and
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
matched random intervals (100 placements per prediction) [Fig. 6e; `RESULT_BIOVAL.md`].

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
these, with CAGE support at the start, repeat content or exon number [Fig. 6f; `RESULT_S3.md`,
`RESULT_V20_CHR9.md`]. Detection rates rose from 22–24% in the lowest to 41–42% in the highest CpG quartile and
from 12–13% to 48–54% across span quartiles. In post hoc analyses, lncRNA genes antisense to a coding gene were
harder to detect (−0.26 and −0.29), and lncRNA genes whose start lies in a transposable element (about a quarter)
were harder to detect (detection rate −0.064 on chr20 and −0.067 on chr9), mainly because they lack CpG-island promoters
[`RESULT_S3_TE.md`]. The model thus locates lncRNA-associated regions mainly through gene-like promoters and long
transcription units rather than through exon sequence; sections 12–15 show that the same dependence on CpG-island
promoters governs where the start of an lncRNA gene can be placed.

## Methods (additions)

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
v2.0 replaces the 154 training windows on chr9 and adds a gene-level logistic calibration fitted on chr18 (features:
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

## Figure plan (Figure 6)

* a) lncRNA recall and strict precision on chr20, T2T chr20 and chr9, exon overlap and S50 side by side, for our
  three models (points) and GENATATOR.
* b) Same-ruler comparison with SegmentNT, AlphaGenome RNA-seq segmentation, Tiberius, Helixer: lncRNA recall
  against strict precision; inset: lncRNA genes reported as coding.
* c) Ablations: ΔAUPRC per seed, chr18 and chr20.
* d) Cross-species: mouse chr19 and four panel chromosomes, ours against GENATATOR (lncRNA and mRNA).
* e) Transcription evidence for predictions outside the reference, by stratum, against matched controls.
* f) Detection against promoter CpG and gene span quartiles (chr18, chr20, chr9).

---

## 作者备注（中文）

- **措辞**：全节按 B1 事前规则的最低一档写（「lncRNA-associated transcribed regions」），标题和每个位点数字都和 S50 并列。如果你更想保留「locus」一词，可以写成「locus-level detection」，但不能写「discovery of lncRNA genes」。
- **新算的数**：T2T chr20、小鼠 chr19、四个面板物种的 S50 是本次补算的（`scripts/phase1/b1/b1_more.py`，`docs/phase1_results/B1/B1_MORE.json`），只用已读出的预测，没有新读取。原口径数字与已发表的完全一致（核对过）。面板物种的 S50：鸡 0.088、斑马鱼 0.119、果蝇 0.048、线虫 0.007；GENATATOR 依次为 0、0.007、0、0。
- **T2T 上 GENATATOR 的严格 Pr**：按转录本算是 0.259（P1.5 发表值），按基因合并后是 0.277；正文用前者，方法里已注明。
- **会被问到的**：
  1. seed 方差大（Sn 0.21–0.40）。正文已如实报告，未加 seed（按「投稿前不在稳健性上循环」的原则）。
  2. 编码侧结构远弱于 Tiberius。已写明不主张编码注释有竞争力；统一 GFF3 里编码结构建议用专门工具——这一点还需要真正实现（`CONVERGENCE_20261112.md` 第 5 项）。
  3. AlphaGenome 的比较是样本内、且是我们自己的朴素分段，正文已限定。
  4. 第 10 节的证据是「有转录」，不是「结构对」，正文已写明。
  5. v3.0（加小鼠）是负结果，放在正文还是补充材料，待全文篇幅定了再决定。
- **章节编号**：本节为 7–11，端点一节已顺延为 12–15，基准一节为 1–6。图号：本节图 6，端点一节图 7。
