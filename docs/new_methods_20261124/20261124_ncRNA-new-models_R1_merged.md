# ncRNA 新模型 / 新方法：第 1 轮合并报告（2026-11-24）

> 来源：ChatGPT 深度研究（`deep_research/20261124_ncRNA-new-models_ChatGPT_R1.md`）、Gemini Deep Research（`deep_research/20261124_ncRNA-new-models_Gemini_R1.md`）。两家只拿到中立提示，没有项目背景。按作者新规则，本轮不做 Claude 调研，本文件由 Claude Code 主持合并。
> 核对：关键论文用 Europe PMC 核对 PMID，仓库用 GitHub API 核对；核不到的标【未核实】。

## 1. 结论先行

两家判断高度一致：
- **不要再训练通用 RNA 基础模型**（RNA-FM、RiNALMo、AIDO.RNA、Evo 2 已覆盖；1–2 张 24 GB GPU 竞争不了）。
- **不要在随机划分上刷分**：家族留出、低同源测试下，RNA 语言模型的优势缩小甚至消失（Brief Bioinform 2025，PMID 40205851；Brief Bioinform 2026 零样本基准，PMID 41789565）。
- RNA 3D 仍依赖模板。CASP16 核酸评估（Proteins 2026，PMID 41165252）的结论是，控制模板后深度学习没有带来实质提升。
- 真正的机会在**改变问题的定义**：
  - 条件依赖：细胞状态、上下文；
  - 系综：不是单一结构；
  - 显式建模实验观测过程：assay 是有偏的观测；
  - 跨物种深同源；
  - 严格的分布外评价。

  把现成大模型当冻结特征抽取器，把有限算力投到这些地方。

## 2. 两家的候选方向

| # | ChatGPT | Gemini |
|---|---|---|
| 1 | 结构感知的 ncRNA 变异效应模型（结构化小 RNA，人群约束 + 共变 + RNA LM，家族留出） | 细胞系条件化的 lncRNA **连续**亚细胞定位（细胞 RBP 表达谱作条件，预测连续 CNRCI，留一细胞系） |
| 2 | 序列 + 细胞上下文的条件化体内结构 / 化学探测模型（RASP v2.0 的 438 个数据集） | lncRNA 序列 → 降解动力学（PRO-seq / 4sU，Saluki 只做 mRNA） |
| 3 | 多模态深同源**检索**（同线性 + 基序 + 结构 + RBP 结合 + 基因组上下文） | 摊销模拟推断（SBI）做化学探测多构象去卷积（对标 DREEM / DRACO） |
| 4 | RNA–RNA 邻近连接数据的潜变量观测模型（显式建模 assay 偏差） | RBP 结合轮廓软对齐的同线性 lncRNA 保守模块发现（LncLOOM 升级 + 树模拟零模型） |
| 5 | 抗泄漏的 RNA 基础模型基准 | — |

**两家共同指向的**（第 2 轮重点）：
- **深同源 / 跨物种功能同源**：ChatGPT #3 ≈ Gemini #4。最近工作是 lncHOME（Nat Genet 2024，PMID 38195860）和 13 物种 lncRNA 直系同源工作流（NAR 2026，PMID 42483833）。
- **结构系综 / 体内结构**：ChatGPT #2 与 Gemini #3 都主张预测系综而不是单一结构，并利用化学探测数据。
- **assay 观测模型**：ChatGPT #4（RNA–RNA）和 Gemini #3（化学探测）都是「潜变量 + 观测过程」思路。

**不值得做**（两家一致）：
- lncRNA vs mRNA 分类器（LncADeep 2.0 AUROC 0.980，Bioinformatics 2026，PMID 41923359）；
- Rfam 随机划分家族分类；
- 从头预训练 RNA 基础模型；
- 单序列 3D 正面竞争 AF3 / Boltz；
- 纯序列二元定位、二元 miRNA 靶标；
- 「contact map 当图片 + U-Net」式纯架构移植。

## 3. 关键事实核对

| 事实 | 出处 | 核对 |
|---|---|---|
| 去掉「中间排除」后 lncRNA 定位准确率掉到 56–59% | NAR Genom Bioinform 2024，PMID 39296930 | 论文存在；具体数字【未逐一核对】 |
| LncRFM（lncRNA 基础模型，跨界编码 F1 0.374） | bioRxiv 2026-11-20，10.64898/2026.11.20.768117 | 存在（4 天前的预印本） |
| lncHOME：同线性 + RBP 结合 → 斑马鱼–人功能保守 lncRNA | Nat Genet 2024，PMID 38195860 | 属实 |
| 13 物种 lncRNA 直系同源工作流 | NAR 2026，PMID 42483833 | 属实；严格标准下每物种约 200 个【数字未核】 |
| Gnocchi：76,156 个基因组，3.7% lncRNA 序列强约束 | Nature 2024，PMID 38057664 | 论文属实；3.7% / 22.2%【未核】 |
| RNA LM 二级结构在低同源下不稳定胜过 one-hot | Brief Bioinform 2025，PMID 40205851 | 属实 |
| RNA LM 零样本基准 | Brief Bioinform 2026，PMID 41789565 | 属实 |
| CASP16 核酸评估：模板依赖 | Proteins 2026，PMID 41165252 | 属实（Gemini 写成 2025） |
| SwitchSeeker：体内 DMS 找结构开关 | Nat Methods 2024，PMID 39014073 | 属实 |
| RASP v2.0：438 个结构探测数据集 | NAR 2025，PMID 39546630 | 论文属实；438【未核】 |
| GPN-Star | **Nature 2026，PMID 42717086** | Gemini 说「预印本」，已过时 |
| Tani 2026：RBP 占有预测 lncRNA 定位 | Int J Mol Sci 2026，PMID 42353306 | 属实 |
| RNAGym | bioRxiv 2025 预印本 | 属实 |
| tRNAscan-SE 最新版 v2.0.13（2026-07-05） | GitHub release | 属实（Gemini 写 2.0.12，过时） |
| 仓库：evo2 4,340★；RNA-FM 396★；RiNALMo 181★（2026-10 有推送）；RhoFold 248★；DRfold2 54★ | GitHub API 2026-11-24 | 属实 |

**各家的问题**：
- **ChatGPT**：导出的 Markdown 里引用只剩 `citeturn…` 占位符，没有链接，需要在第 2 轮请它补 DOI / PMID。「Rfam 15.1 有 4,227 个家族」【未核实】。
- **Gemini**：
  - GPN-Star 发表状态过时；tRNAscan-SE 版本过时；
  - CASP16 年份写成 2025；
  - 引了 ResearchGate、MDPI 综述等弱来源；
  - 「arXiv 2610.33960 miRNA 基准」两家都引了，arXiv 编号【未核实】。

## 4. 主持人判断（Claude Code，结合项目背景；本节不发给两家）

- **Gemini #1（细胞条件化连续定位）与我们的方向 A（`13_RNA_fate/`）几乎重合**。E1 / E1b 已测出：去混杂后序列对 lncRNA 核 / 质比的增益 ΔR² 约 0.11–0.13，与 Tani 2026 是同一类问题。Gemini 补上了 A 缺的东西：
  - 不做中间排除，预测连续值；
  - 用细胞自身的 RBP 表达谱作条件，留一细胞系外推。

  这把「静态序列定位」（已饱和）变成「细胞条件化的定位」，新意比 A 的现状强，而且 A 的数据和代码可以直接复用。第 2 轮应重点让两家评它。
- **Gemini #2（lncRNA 降解动力学）**：方向 A 的检查已指出 Churchman 2024（Mol Cell，PMID 38964322）里的 lncRNA 只有 22 个，真值稀缺。
- **ChatGPT #1（结构化小 RNA 变异效应）**：与已关闭的 sncRNA 副线（`11_sncRNA_tools/`）有交集，那条线判过「tRNA 已解决、Rfam 家族标签少」。真值（致病变异）稀少是主要风险。
- **深同源（两家共同）**：最有「新模型」味道，但金标准极少（lncHOME 的救援实验只有数十个）。
- **SBI 化学探测去卷积（Gemini #3）**：方法新意清楚，真值可以靠合成混合体和公开金标准，但离 ncRNA 生物学问题稍远。
- **抗泄漏基准（ChatGPT #5）**：可以作为任何方向的评价部分，单独成文新意不足。

## 5. 第 2 轮要做的

把对方的前 3–4 名发给每家，要求：
- 逐条评、指出对方的事实错误；
- 更新排序；
- 列出 3 个最能改变方向的开放问题和最便宜的实验；
- ChatGPT 补上关键引用的 DOI / PMID。
