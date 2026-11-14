# Draft — results section: locating the ends of lncRNA transcription units

> Status: first draft, 2026-11-13. Covers the end-localisation line (B1, S11, BND, END). Every number points to
> its result document in `docs/phase1_results/` and its JSON. 2026-11-15: read-through; figure drawn
> (Figure 6). Sections 12–15; they follow the benchmark (`DRAFT_results_benchmark.md`, sections 1–6) and the
> locus-detection results (`DRAFT_results_locus.md`, sections 7–11). Renumbered from 7–10 on 2026-11-14;
> wording and citations revised after the novelty check (`NOVELTY_ends.md`) on the same day.
> Language: English (manuscript); notes for the authors in Chinese at the end.

## Results

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

## Methods (additions)

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

## Figure 6 (drawn)

* a) lncRNA and mRNA recall under exon overlap versus S50 for our models, GENCODE v41, GENATATOR, SegmentNT and
  Tiberius (section 12).
* b) Exact intron chains of GENATATOR segmentation against interval padding / trimming, chr18 (section 13).
* c) lncRNA start localisation by promoter CpG quartile: our seven models on chr18 and chr20, Puffin on chr8 + chr9
  (section 15).
* d) polyA-supported end localisation, lncRNA and mRNA: our seven models and PolyaID (section 15).
* Supplementary: tiling, decoder and reference-support controls (section 14); label support statistics; per-seed
  probe tables.

## References cited in sections 12–15

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

---

## 作者备注（中文）

- **这一节改变了全文的定位**：原 L1 的 0.306 不能再写成「发现了三成 lncRNA 基因位点」。第 12 节末尾的斜体句是事前规则给出的措辞（`PLAN_B1.md` 第 3 节，落在最低一档）。「section L」（v1.9 / v2.0 的位点检出结果）成稿时要同步改口径：两种口径并列，标题和摘要用 S50 能支撑的说法。
- **哪些是检验、哪些是描述**：B1 的措辞规则、S11 的停止规则、END 的判读规则都是事前写定的。BND 是事前写定读数的诊断。END 里的 CpG 分层和终点结论是**描述性**的，7 个模型方向一致，但不是事前判据；正文已写成「held in all seven models」，不要写成检验结论。
- **审稿人可能会问**：
  1. 开发在 chr18 上做（BND、END、S11）。BND 与 END 的读数已于 2026-11-14 在 chr20 上按事前规则复现（`RESULT_REP20.md`）；S11 未复现，且 chr18 对 GENATATOR 是训练染色体（已注明是上界）。
  2. 「CAGE 支持」与「CpG 高」没有拆开；FANTOM5 CAGE 偏向高表达转录本。
  3. END 两个臂的 seed 方差很大（同一臂 mRNA 起点 0.785 对 0.414），所以第 15 节对两种改进办法不下结论。返修时如被要求，可补 seed。
  4. 「终点定位不了」是在这套特征（GENERanno 第 17、23 层，24.6 kb 窗口）和这个头下的结论，不是对序列信息的一般断言。2026-11-14 已按查新把第 15 节标题和结论句收紧为「with (current) DNA language-model features」，并引 APARENT（Bogard 2019）、Stroup 2023，说明专用 polyA 模型未测。
- **与文献的关系**：lncRNA 剪接效率低（Mele 2017）、SPLAIRE（低使用量位点难预测）、GENATATORs（区间内分割）都需要在讨论里引用；本节的新意在于端到端地把瓶颈定位到转录单位的端点，并量化区间分割对端点误差的容忍度。查新已做（`NOVELTY_ends.md`，2026-11-14）：主说法未见先例；CpG 分层写成对 Hannenhalli 2001 的印证并用 Hon 2017（增强子起源）、Mattioli 2019 解释；端点注释不完整（Lagarde 2017、Hon 2017、Kaur 2024）只作背景。禁用说法：「首次发现 lncRNA 端点难注释」「序列不含终点信号」。
- **没写进来的**：S1 / S1R（剪接使用量）那条机制线；它与本节的关系是「端点已知之后，结构恢复受剪接使用量限制」，放在第 13 节之后还是讨论里，待全文结构定下来再决定。
