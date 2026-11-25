# Gemini 第 2 轮互评：要点存档（非逐字原文）

> 页面抽取。回复约 5 万字符，其中大段是 Gemini 泄露出来的英文起草与自检过程；这里只存最终中文内容的要点。原对话：https://gemini.google.com/app/30e30f82d0db1a61
> 过程异常：该对话仍处于 Deep Research 模式，第一次发送后卡在「Next Steps in Verification」9 分钟以上；手动取消生成并刷新后，才显示出完整回复。

## 对 ChatGPT 方向的评价
- **总体判断同意**：缺严格 OOD 基准、条件依赖模型、跨物种深同源、把测序当带噪观测；不需要更大的纯序列 LM。
- **PMID 40205851 的结论要补边界**：RiNALMo / ERNIE-RNA 在部分未见家族上明显好于经典方法（tmRNA F1 0.73 vs 0.38）。真正全线失败（F1 < 0.15）的是超长、样本少的 telomerase 家族；one-hot 只是略胜较弱的 RNAErnie。
- **LncADeep 2.0**：AUROC 0.980 在它的检索片段里属于 LncRFM，不是 LncADeep 2.0【我方按 ChatGPT R1 的出处待核】。「lncRNA 分类已饱和」不准确：GENCODE v46–v47 稳定子集上约 45% 的转录本存在工具间分歧（PMID 42039468）。
  - 【核对：PMID 42039468 是 bioRxiv 预印本，不是同行评审，Gemini 写错。】
- **ChatGPT #1（结构化小 RNA 变异效应）**：强烈质疑。
  - 空间分辨率错配：Gnocchi 按 1 kb 窗口打分，而 miRNA / snoRNA / tRNA 只有 22–200 nt；多数 snoRNA 在宿主基因内含子里，窗口分数会混入邻近外显子、剪接位点的选择压力。
  - 除少数线粒体 tRNA、RNU4ATAC、TERC 外，缺少足量的致病变异真值。
- **ChatGPT #2（条件化体内结构，RASP v2.0）**：高度赞同范式，但不能混合全部 438 个数据集。
  - 其中 81 个是邻近连接、141 个是靶向探测，方法横跨 DMS（只测 A/C）与 icSHAPE（四种碱基），还有 59 套是插补值。
  - 应聚焦单一试剂平台内的跨区室 / 跨细胞比较，例如同一 icSHAPE 平台的 HEK293 体内 / 体外 / 染色质 / 核小体 / 胞质。
- **ChatGPT #3（深同源检索）**：强烈同意，认为检索框定和「剥离邻近蛋白编码基因身份」的消融优于自己原来的方案。
  - 风险：非模式物种没有 CLIP，只能用 RBPNet 预测。
  - 真值不能只靠同线性定义，要用独立证据，例如 Degalez 2026 的 17 组织表达保守性。
- **ChatGPT #4（RNA–RNA 潜变量模型）**：问题真实，但投入产出低。
  - 技术偏好在重复之间本身就可重复，所以「重复间可重复性」作主指标有盲区。
  - 应像 KARR-seq（PMID 38238480）那样，对已知冷冻电镜结构（核糖体、TERC、U3）做精确度 / 召回验证。
- **ChatGPT #5（抗泄漏基准）**：不同意作主线。Brief Bioinform 2025、xvalRNAfolding、RNAGym、NABench、GraphFLA、LncBench 已覆盖，窗口基本关闭。

## 自己漏掉的、因此改变看法的
- **Gnocchi**：只有 3.7% 的 lncRNA 序列受强约束（与整体非编码区的 3.1% 相近），miRNA 为 22.2%。
  - 推论：约 96% 的 lncRNA 序列近乎中性，纯序列 MLM 学通用语法效率很低。
  - 这支持「局部基序 / RBP 锚点检索」和「环境条件建模」。
- **RASP v2.0** 使「序列 + 区室上下文 → 差异反应性」可以直接监督学习，因此下调自己偏理论模拟的 SBI 方向（合并或聚焦）。
- **lncHOME 与 Degalez 2026（NAR，PMID 42483833）**：同源检索竞争比原来估计的更活跃，并提供了量化基线：
  - 13 物种严格同源每物种约 200 个；
  - 人–斑马鱼 570 个 coPARSE-lncRNA。

  它的方向 4 必须采用「检索 + 剥离邻近 PCG」的设计。
- **GPN-Star 已发表于 Nature 2026**：基于多物种比对的约束模型已是强敌，小团队不宜在非编码变异约束上正面竞争。

## 更新后的前 3
1. **区室与 RBP 环境条件化的体内结构反应性 / 差异折叠模型**（融合 ChatGPT #2 与自己的上下文视角）。
   - 数据：RASP v2.0 中同一平台的 HEK293 多区室 icSHAPE，只用非插补的高质量位点。
   - 划分：基因家族 / 相似性聚类留出 × 留一区室。
   - 模型：冻结 RiNALMo / RNA-FM + 轻量交叉注意力条件头。
   - 门槛：Pearson +0.05，差异位点（|Δreactivity| 前 10%）auPRC +0.05。
   - 基线：5-mer / 7-mer 上下文均值反应性表（含体内外偏移校正）；ViennaRNA 未配对概率 + 序列 ridge。
2. **剥离邻近 PCG 身份后的多模态快演化 lncRNA 深同源检索**（融合 ChatGPT #3、自己的 RBP 软对齐和系统发育零模型）。
   - 数据：Degalez 2026 的 13 物种严格同源集（450 多条人类保守 lncRNA）+ lncHOME 人–斑马鱼 570 对；同染色体和非同染色体诱饵；遮蔽邻近 PCG。
   - 独立检验：17 组织表达保守性。
   - 门槛：无 PCG 特征时 Recall@10 比 BLASTN + k-mer 基线高 10 个百分点。
   - 基线：宽松 BLASTN + 6-mer / SEEKR 余弦、LncLOOM 基序重叠。
3. **lncRNA 半衰期与核 / 质分配的统一反式条件建模**（自己原 #2 的精炼）。
   - 假设：核滞留与半衰期是一体两面，可在连续尺度上联合建模。
   - 数据：K562 PRO-seq / RNA-seq 半衰期（PMC7885420）+ lncATLAS 15 细胞系连续 CNRCI + ENCODE RBP 表达。
   - 门槛：未见细胞 × 未见基因的 CNRCI Pearson +0.08（对加性和 k-mer 基线），半衰期 Spearman +0.05（对 LightGBM）。
   - 【注：Gemini 没有提到 PMID 42760417，那篇工作直接显示 lncRNA 半衰期几乎无法从序列预测。】

## 三个开放问题与最便宜的实验
1. **RASP 跨区室差异是真生物学信号还是噪声？**（2–3 天）
   - 做法：下载 HEK293 跨区室 icSHAPE 及其重复（不含插补值），按覆盖度分层；比较「跨区室 Δ 的重复间相关」与「重复内技术差异的幅度」；用 7-mer 上下文均值 + ridge 拟合 Δ。
   - 读法：生物方差不显著大于技术方差，或差异主要由测序深度解释，就降级。
2. **遮蔽 PCG 后，序列侧还有没有同源信号？**（3–4 天）
   - 做法：lncHOME 570 对 + Degalez 严格直系同源，每个查询配 50 条随机非同源 lncRNA（含同染色体背景）；零样本比较 6-mer 余弦、冻结 RNA-FM / RiNALMo 余弦、RBPNet 结合峰相关的 Recall@10。
   - 读法：三者都接近随机，就止损。
3. **深度表征是否超过矩阵分解 + 表格特征？**（2–3 天）
   - 做法：lncATLAS 双线性回归、LightGBM 半衰期，再把冻结的 RNA-LM 嵌入拼进去。
   - 读法：交叉验证增量 < 0.02 Pearson，就不值得堆深度模型。

## 引用核对（Gemini 自报）
- 同行评审：
  - CASP16 核酸（10.1002/prot.70072）；
  - AF3 at CASP16（PMID 40851426）；
  - 定位中间排除（10.1093/nargab/lqae125，PMC11409063）；
  - RNA-LLM 二级结构基准（10.1093/bib/bbaf137，PMID 40205851）；
  - R-scape（10.1371/journal.pcbi.1011262）；
  - Dorado 评估（Genome Biol 2026，PMID 42098827；eLife 2026，110672）；
  - Evo 2（Nature 2026，PMID 41781614）；
  - GPN-Star（Nature 2026，PMID 42717086）；
  - SnoBIRD（PMID 40737086）、KARR-seq（PMID 38238480）、VASA-seq（PMID 35760914）、Saluki（PMID 36419176）、RBPNet（PMID 37542318）。
- 预印本或未验证：
  - GraphFLA（arXiv 2510.24826）、RNAGym、NABench（arXiv 2511.02888）、LncRFM；
  - **arXiv 2610.33960 自认无法核验，不作立项依据**；
  - bpRNA-1m k=0 F1 0.608 来自 BEACON 的 OpenReview 回复，不是同行评审正文。
