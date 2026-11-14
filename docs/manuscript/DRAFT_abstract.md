# Draft — Abstract

> Status: first draft, 2026-11-14. Every number is taken from the results drafts (sections 1–15); none is new.
> Length: ~190 words (Nature Communications / Genome Biology limit 150–250). A 150-word version follows.
> Language: English; notes for the authors in Chinese at the end.

## Abstract (~190 words)

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

## Short version (~150 words)

Ab initio gene finders annotate protein-coding genes from sequence alone, but not long non-coding RNA (lncRNA)
genes. We built a joint coding–lncRNA annotator that decodes a gene grammar over frozen DNA language-model
features, and an evaluation framework for an incomplete reference. On held-out chromosomes it detected 30% of
lncRNA genes at a strict precision of 0.42–0.50, ten times the recall of existing tools, and transferred across
vertebrates. Only about 7% of lncRNA genes, however, were delimited by a single prediction, against 55% of coding
genes. The ends of transcription units, not exon structure, account for this gap: given ends within a kilobase,
structure becomes partly recoverable, but lncRNA starts can be located from sequence only at CpG-island promoters,
and lncRNA ends not at all — by general DNA language-model features and by dedicated initiation and
polyadenylation models alike.

---

## 作者备注（中文）

- 摘要的主线是「检出 ≠ 界定，端点限速」，与新标题一致。基准（真值在快速增长、现有工具约 3%）压缩成一句「evaluation framework built for an incomplete and rapidly growing reference」和「ten times the recall of existing tools」。如果目标期刊更看重基准，可以把 v41→v50 只覆盖 49% 这个数字加进去。
- 「30%」是 chr20 的 L1（0.306）与 chr9（0.304）；「0.42–0.50」是 chr20 v1.9 与 chr9 v2.0 的严格精度；「about 7%」是 S50（chr20 6.4%，chr9 6.9%）；「55%」是我们自己的 mRNA S50（chr20）。
- 「also by dedicated models」依据 DEDEND：Puffin 在其留出染色体上低 CpG lnc 起点 0.086，PolyaID lnc 终点 0.07 / 0.12。T2（有 CAGE 的 lnc 起点）是灰区，摘要没有涉及。
- 没有写进摘要的：seed 方差、编码结构弱于 Tiberius、线虫失败、加小鼠训练无益、预测外区域的转录证据（59% vs 29%）。后者如果篇幅允许可加半句，它说明「检出」的区域确实在转录。
- 「transferred to mouse and non-mammalian vertebrates」：小鼠 L1 0.280、鸡 0.318、斑马鱼 0.348；果蝇 0.214 不是脊椎动物，线虫失败，所以只写脊椎动物。
