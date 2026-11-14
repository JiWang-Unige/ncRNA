# Draft — Introduction and Discussion

> Status: first draft, 2026-11-14. Built on `DRAFT_results_benchmark.md` (sections 1–6), `DRAFT_results_locus.md`
> (7–11) and `DRAFT_results_ends.md` (12–15), and on the novelty check `NOVELTY_ends.md`. Every number repeats one
> in the results drafts; no new numbers. References marked [verify] have not been checked against the original
> paper yet. Language: English (manuscript); notes for the authors in Chinese at the end.

## Working title (proposed, replaces the October title)

*Sequence-only annotation finds lncRNA-associated transcribed regions but cannot delimit them: the ends of
transcription units limit ab initio lncRNA annotation*

Alternatives are listed in the author notes.

## Introduction

Most newly sequenced eukaryotic genomes are annotated without matched transcriptomic data for every tissue and
developmental stage, so gene models come in large part from ab initio predictors that read the genome sequence
alone. For protein-coding genes these predictors have improved markedly: deep-learning models combined with
hidden-Markov decoding now recover most human coding genes at gene level and transfer across distant lineages
(Tiberius, Helixer, ANNEVO) [verify exact figures before citing]. Long non-coding RNA (lncRNA) genes are as
numerous as protein-coding genes in the human reference, yet they are absent from the output of these tools.
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

## References to add (beyond the results drafts)

* Derrien T et al. The GENCODE v7 catalog of human long noncoding RNAs. *Genome Res* 2012.
* Dudnyk K et al. Sequence basis of transcription initiation in the human genome. *Science* 2024.
* Cochran K et al. Dissecting the cis-regulatory syntax of transcription initiation with deep learning. *bioRxiv* 2024.
* Mele M et al. Chromatin environment, transcriptional regulation, and splicing distinguish lincRNAs and mRNAs.
  *Genome Res* 2017 [verify].
* Schlackow M et al. Distinctive patterns of transcription and RNA processing for human lincRNAs. *Mol Cell* 2017
  [verify].
* Tiberius (Gabriel et al., *Bioinformatics* 2024), Helixer (Holst et al., *Nat Methods* 2025), ANNEVO
  (*Nat Methods* 2026), SegmentNT, AlphaGenome, GENATATOR (2025) — citation details to fill from the benchmark draft.

---

## 作者备注（中文）

- **标题**：10 月的标题（「Non-coding genes are invisible to ab initio genome annotation…」）以基准为主线；现在主线是「检出 ≠ 界定，端点限速」，所以换了标题。备选：
  1. 「The ends of transcription units limit sequence-only annotation of lncRNA genes」（短，最直接）；
  2. 「Finding but not delimiting: ab initio annotation of lncRNA genes from DNA sequence」；
  3. 保留旧标题作为副标题之一，突出基准。
  选择取决于目标期刊：Nature Methods / Genome Biology 偏方法与基准，可用 2；若偏生物学结论（如 Nature Communications），用 1。
- **叙事顺序**：引言最后一段按「基准 → 方法 → 检出 → 分解」的顺序，对应结果 1–6、7–11、12–15。如果主线换成「端点」，可考虑把第 12 节（S50）提前到位点检出节里，第 13–15 节单独成一节。
- **引言里刻意没写的**：没有说「首次」，没有说「序列不含端点信号」（查新禁用说法）；端点注释不完整只作背景。
- **[verify] 条目**：Tiberius/Helixer/ANNEVO 的具体数字、AlphaGenome 的引用、Mele 2017 与 Schlackow 2017 的具体结论（剪接效率低、终止方式不同）都要读原文确认后再保留。Schlackow 2017 那句若核对不上就删掉，不影响论证。
- **讨论里提出的「下一步实验」**：把 Puffin/ProCapNet 类起始模型、APARENT 类 polyA 模型放进 ±20 kb 搜索窗口，看能否选出真实端点。这是对主说法最直接的补强，也可能被审稿人要求；算力不大（推理为主），是否投稿前做，等你决定。
- **专用端点模型（DEDEND，2026-11-14）**：T0 有效（Puffin mRNA 0.71），T1 = B（低 CpG lnc 起点 0.086），T2 灰区（有 CAGE 的 lnc 起点 0.40），P1 = B（PolyaID lnc 终点 0.07 / 0.12）。讨论第二条方向已据此改写。「which is a property of the cell rather than of the local sequence alone」是解释，不是我们测出来的，审稿人可能挑；如不想冒险可删。
- **局限 (ii)**：2026-11-14 已在 chr20 上复现第 14–15 节（`RESULT_REP20.md`）：R1–R3 成立，R4 只对 lnc 成立（两个探针的 mRNA 终点 0.17–0.18）。正文已按事前规则把「终点对两类都定位不了」改成「lnc 终点定位不了，mRNA 终点信号弱」。第 13 节（S11）未复现。
- **没写进来**：摘要（已另写 `DRAFT_abstract.md`）。统一 GFF3 已实现（`RESULT_UGFF.md`），讨论「Relation to other approaches」末尾加了一句发布形式。
