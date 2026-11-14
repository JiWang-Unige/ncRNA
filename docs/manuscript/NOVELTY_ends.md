# 查新：「端点定位是 lncRNA 从头注释的限速步骤」（2026-11-14）

> 检索：PubMed、Consensus、网页检索、bioRxiv/arXiv 原文。只读了摘要和可获取的正文；结论是"未找到先例"，不是"确认无先例"。
> 对应草稿：`DRAFT_results_locus.md`（§7–11）、`DRAFT_results_ends.md`（§12–15）。

## 结论

拆成四条子论点，逐条看：

| # | 子论点 | 先例 | 判断 |
|---|---|---|---|
| 1 | lncRNA 参考注释的 5′/3′ 端不完整 | Lagarde 2017（CLS）；Hon 2017（FANTOM CAT，CAGE 校正 5′ 端）；Kaur/Perteghella 2024（GENCODE 长读扩容）；Kainth 2023（短读组装末端模糊） | **已知**，只作背景引用 |
| 2 | 纯序列模型：lnc 位点能检出，但界定不了；误差主要在端点（L1 vs S50 约 5 倍差距；端点偏移 8–21 kb vs mRNA 1–2 kb；拼接/解码/真值支持均被排除） | GENATATORs 2026 预印本在局限中写「基因边界（含 UTR）对所有工具仍困难」，但不评 TSS/TES、不拆分检出与界定、不区分 lnc/mRNA 端点 | **未见先例，可作为主新说法**；必须引用 GENATATORs 并说明区别 |
| 3 | 序列只能定位 CpG 岛型 lnc 起点，低 CpG 起点在随机水平 | Hannenhalli 2001：非 CpG 启动子仅靠局部序列可能无法预测（针对一般基因）；Hon 2017：多数基因间 lncRNA 起源于增强子；Mattioli 2019：lincRNA 启动子 motif 更简单、组织特异；Alam 2014：lnc 与编码启动子特征不同；Puffin/ProCapNet/CLIPNET 能从序列解释多数启动子起始 | **机制部分已知，组合新**：「DNA 语言模型特征下 lnc 起点可定位性按 CpG 分层，Q4 0.17–0.61，Q1–Q2 ≈ 随机」未见报道。应写成对 Hannenhalli 的现代、lnc 专属的定量印证，并用 Hon 2017 的增强子来源解释 |
| 4 | 终点（TES）在这套特征下对 lnc 与 mRNA 都定位不了；mRNA 后验里的终点精度来自解码语法 | APARENT（2019）、Stroup 2023（PolyaID 类，核苷酸分辨率 polyA 位点）能在给定窗口内预测切割位点 | **需小心措辞**：只能说「冻结的通用 DNA LM 特征 + 线性/小探针在 ±20 kb 内选不出正确终点」，不能说"序列不含终点信号"。审稿人会拿 APARENT 类模型反驳 |

**总体判断：主新说法（#2）站得住**。它是「从头注释任务」的分解结论：检出 ≈ 可做、界定 ≈ 不可做，瓶颈在端点。已有工作要么是注释/测序侧（#1），要么是启动子/polyA 的给定窗口建模（#3、#4），要么只在局限里泛泛提一句（GENATATORs 2026）。

## 最接近的工作与区别

1. **GENATATORs（Shmelev et al., bioRxiv 2026-06, doi 10.64898/2026.06.17.732686）**
   - 内容：对 DNA LM 做微调分割（exon/intron/CDS/UTR），用零容差区间 F1 和基因级指标；报告 lncRNA 比 mRNA 难；「任何模型只能正确分割约 30–40% 人类基因」；局限中说基因边界（含 UTR）识别仍困难。
   - 区别：不评 TSS/TES；不拆分"位点检出 vs 界定"；没有端点诊断（拼接/解码/真值支持）；没有 CpG 分层。
   - **对我们的影响**：(a) 必须引用；(b) 基准部分目前评的是旧版 GENATATOR 流程，审稿人会问新版。新版主要报告区间内分割，与我们的"端到端位点"任务不同，可在讨论里说明，必要时补跑（另议）。
2. **Hannenhalli & Levy 2001（Bioinformatics）**：CpG 岛是启动子预测的主导信号，非 CpG 启动子仅靠局部序列可能无法预测。我们的 #3 是在 DNA LM 时代、针对 lnc 的定量重现，应正面引用，而不是声称新发现。
3. **Hon et al. 2017（Nature，FANTOM CAT）**：多数基因间 lncRNA 起源于增强子 → 增强子起源通常缺 CpG 岛 → 解释为什么大量 lnc 起点不可定位。这是讨论里的生物学落点。
4. **Puffin（Dudnyk 2024 Science）/ ProCapNet / CLIPNET**：给定候选区域时能在碱基分辨率解释起始。区别：它们预测的是"已知调控元件内的起始分布"，需要全基因组扫描时在大量增强子/启动子中挑出真正的基因起点，这一点它们不评估。讨论里可提出"把这类专用起始模型接入解码"作为前景。
5. **APARENT（Bogard 2019 Cell）/ Stroup 2023（Nat Commun）**：polyA 位点的序列模型。用于 #4 的措辞约束和讨论中的"可接入组件"。

## 对草稿的修改建议

- 引言：把「端点不完整」放进背景（Lagarde 2017、Hon 2017、GENCODE 2024），然后提出本文的问题：「纯序列模型的瓶颈在检出还是界定？」——这是先例没问过的。
- §12–15：CpG 分层段加 Hannenhalli 2001 与 Hon 2017；TES 段改成"在冻结特征和探针下"，并引 APARENT/Stroup，说明专用 polyA 模型是否能解决未测。
- 讨论：GENATATORs 2026 单独一段，说明任务差异（区间内分割 vs 端到端位点+端点）。
- 不能用的说法：「首次发现 lncRNA 端点难以注释」「序列不含 lnc 终点信号」。

## 未覆盖

- 未检索会议论文全文（ICLR/NeurIPS 生物轨道）与中文文献。
- GENATATORs 2026 的表 1/2（各类别数字）未读到，正文只拿到了摘要和部分段落。
- 未检索 AlphaGenome/Borzoi 是否报告过 lnc TSS 定位的分层结果（它们预测 CAGE 覆盖，可能有相关图）。

## 参考

- Lagarde J et al. Nat Genet 2017;49:1731. doi:10.1038/ng.3988
- Hon CC et al. Nature 2017. An atlas of human long non-coding RNAs with accurate 5′ ends.
- Kaur G, Perteghella T et al. bioRxiv 2024. doi:10.1101/2024.10.29.620654
- Carbonell-Sala S et al. Nat Commun 2024;15:5278. doi:10.1038/s41467-024-49523-3
- Kainth AS et al. PLoS Comput Biol 2023.
- Hannenhalli S, Levy S. Bioinformatics 2001. Promoter prediction in the human genome.
- Alam T et al. PLoS ONE 2014.
- Mattioli K et al. Genome Res 2019.
- Dudnyk K et al. Science 2024 (Puffin); Cochran K et al. bioRxiv 2024 (ProCapNet); He AY et al. bioRxiv 2024 (CLIPNET).
- Bogard N et al. Cell 2019 (APARENT); Stroup EK et al. Nat Commun 2023.
- Shmelev A et al. bioRxiv 2026, doi:10.64898/2026.06.17.732686 (GENATATORs).
- Marin FI et al. arXiv 2505.03377 (2025)（编码基因，无 lnc/端点，仅作背景）。
