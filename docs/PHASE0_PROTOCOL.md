# Phase 0 协议（v2，2026-10-20）

> **v2（2026-10-20，用户批准 U-ev 为主臂、U-pu 为副臂、F 为对照、3 个 seed；在任何 0.5 结果之前冻结）**：本节取代 §3 中 0.5 的笼统描述。§1、§2、§4 不变。

## 0.5 预注册（v2）

**骨干与特征**
- GENERanno-eukaryote-1.2b-cds-annotator-preview，冻结，取第 23 层隐状态（0.4 中 lnc 外显子的最佳层），以 6-mer token 为单位。
- 用 PCA 降到 512 维：在 200 个随机训练窗口上拟合，所有臂共用，存为 fp16。
- 窗口长 24,576 bp，即 4,096 个 token。

**数据**
- 训练窗口：从 hg38 chr1–17、chr19、chr21、chrX 中抽 2,000 个，一半以 v41 lncRNA 基因为中心，一半均匀随机；N 占比 > 10% 的窗口丢弃。seed 固定为 0。
- 内部验证：chr18 上 300 个窗口，抽法相同。
- 评测：chr20 全长。窗口长 24,576 bp，步长 20,480 bp，只保留中心部分。

**标签**（只用 GENCODE v41，逐 token 取中心碱基）
- 5 类：背景、lnc 外显子、lnc 内含子、mRNA 外显子、mRNA 内含子。不分链。外显子优先于内含子，lnc 优先于 mRNA。

**unknown 区域（U-ev）**
- 证据来源：Snaptron `srav3h`，即 recount3 人类 SRA 剪接点，2021 年发布，早于 v41。
- 保留同时满足以下条件的剪接点：
  - `samples_count ≥ 100` 且 `coverage_sum ≥ 500`
  - motif 为标准型 GT-AG、GC-AG 或 AT-AC（含反链形式 CT-AC、CT-GC、GT-AT）
  - 长度 ≤ 500 kb
  - 剪接点区间与任一 v41 基因（任意类别、任意链；基因范围取外显子最左到最右）没有任何重叠
- unknown 区域 = 每个保留剪接点的 [start − 500, end + 500]，合并，再减去所有 v41 基因范围。

**三个臂**（网络结构、数据和训练步数完全相同）
- **F**：所有位置都按标签监督，未注释区一律当作背景。
- **U-ev**：unknown 区域内的位置，损失改为 −log(p_背景 + p_lnc外显子 + p_lnc内含子)，即在"未转录 / lnc"三者上做边缘化；其余位置与 F 相同。
- **U-pu**：未注释位置视为"未标注"。对 lnc 外显子这一维用 nnPU 风险估计（β = 0），类先验取 π = 训练窗口中 v41 lnc 外显子碱基数 ÷ 未注释碱基数（即假设未标注部分中 lnc 外显子的质量与已标注部分相等）。确定 π 时不参考 v50。其余标签用交叉熵。

**小头与训练**
- 小头：Conv1d(512→256, k=9) → BiGRU(256, 1 层) → Linear(5)，约 1.0M 参数。
- AdamW，学习率 1e-3，6 个 epoch，batch 为 8 个窗口。用 chr18 验证集上 v41 lnc 外显子的逐碱基 F1 选 checkpoint。
- seed 取 {42, 123, 456}。

**解码**（三个臂同一套流程）
1. 逐碱基 lnc 外显子后验（token 后验复制到 6 bp）。
2. 阈值在 chr18 验证集上按 v41 逐碱基 F1 选定。
3. 取连续片段，长度 ≥ 50 bp。
4. 每个片段视为一个单外显子 lnc 预测。

**终点**
- 主终点：chr20 上 Δ 集合（473 个基因）的 L1 Sn，取 3 个 seed 的平均命中率。
- **成立判据**：U-ev − F ≥ +5 个百分点，且按 Δ 基因做配对 bootstrap（1,000 次）的 95% CI 不含 0；**同时** U-ev 在 v50 下的严格 L1 Pr 比 F 下降 ≤ 5 个百分点。
- 次要终点（只报告）：v50 全部 lnc 与 v41 lnc 的 L1 Sn / Pr、mRNA 指标，以及 U-pu 相对 F 的同一组对比。

**运行**
- 特征缓存与训练在 Modal（A10G）上进行。
- 先做一次小规模 smoke（20 个训练窗口、1 个 epoch、只跑 chr20 的 1 Mb），用来验证流程端到端可行；smoke 的数值不读。

# v1.2（2026-10-20）

> **v1.2 变更（2026-10-20，在任何 0.1b / 0.6 结果之前）**：面板物种的 GTF biotype 命名各不相同。
> - 斑马鱼 GRCz11 用旧命名 lincRNA / processed_transcript / antisense。
> - 果蝇、线虫用笼统的 ncRNA。
>
> 因此固定类别映射（`ruler.py` v1.2）：
> - **lnc** = lncRNA、lincRNA、antisense(_RNA)、sense_overlapping、sense_intronic、3prime_overlapping_ncRNA、bidirectional_promoter_lncRNA、macro_lncRNA；另外，ncRNA 与 processed_transcript 中最长剪接转录本 ≥ 200 nt 的基因。
> - **small** = snRNA、snoRNA、miRNA、pre_miRNA、misc_RNA、scaRNA、vault_RNA、Y_RNA、tRNA、piRNA、rRNA、Mt_tRNA、Mt_rRNA、ribozyme、sRNA、scRNA、SRP_RNA、RNase_P_RNA、RNase_MRP_RNA、telomerase_RNA；另外，上述两种笼统类型中 < 200 nt 的基因。
> - **other** = transposable_element（果蝇）、TEC 等；pseudogene 单独计。
>
> 对人类 GENCODE 的 lnc / mrna 结果没有影响（回归测试与 v1 完全一致）。
> 小 RNA 的 Rfam 预测会剔除 lncRNA 家族（218 个）以及 Cis-reg、Intron 类型（Rfam 15.1 `family.txt` 第 19 列）。
> 各物种的真值密度（lnc 基因 / Mb）作为方案 B 的描述性结果报告。

# v1.1（2026-10-20）

> **v1.1 变更（2026-10-20，用户决定）**：
> - 只做动物。
> - RNA 类别以 lncRNA 为主，加入小 ncRNA。
> - 非哺乳动物面板与染色体只按长度选定，选定时尚未看到任何 0.1b 结果（见 §0）。
> - 新增 §3 的 0.6（小 ncRNA）。
> - 判据与阈值没有任何改动。
>
> v1 冻结于 2026-10-20，早于任何模型结果。

> 定位：A+B 合一（用户 2026-10-19 批准）。Phase 0 的目的是用低成本实验回答两个 go/no-go 问题，并冻结评价尺子。
> 代码与输出：Baobab `common/ncrna-annot/`（事实来源，git）→ ont `~/ncrna-annot/`（运行）。本目录 `05_phase0/` 保存协议与结果报告的副本。
> 算力：ont → Modal → Baobab（Baobab 需用户批准）。ont 的 GPU 被其他项目占用时不插队。
> 修改本协议必须新增版本号并写明原因；判据和阈值不能在看到结果后再改。

## 0. 范围（v1.1 确定）

- **物种**：只做动物。GENATATOR 训练用的是 39 个哺乳动物，以下 4 个非哺乳动物都不在其训练集内：

  | 谱系 | 物种 | 组装 | 评测染色体 | 真值 |
  |---|---|---|---|---|
  | 鸟 | 鸡 Gallus gallus | GRCg7b | 5 | Ensembl 当前版本（116） |
  | 鱼 | 斑马鱼 Danio rerio | GRCz11 | 22 | Ensembl 116 |
  | 昆虫 | 果蝇 D. melanogaster | BDGP6.54 | 2L | Ensembl 116（来源是 FlyBase） |
  | 线虫 | C. elegans | WBcel235 | I | Ensembl 116（来源是 WormBase） |

  染色体只按长度选（15–60 Mb，控制 GPU 时间）。
- **RNA 类别**：lncRNA 是主线（方法主张）；小 ncRNA（tRNA、rRNA、snRNA、snoRNA、miRNA 等）作为第二类别，写进统一 GFF3。
- **小 ncRNA 的循环论证问题**：Ensembl/GENCODE 的小 RNA 注释本身就是用 Rfam/Infernal、tRNAscan-SE、miRBase 产生的。用这些工具对它打分是循环的。因此小 ncRNA 部分**不作为方法主张**，只报告工程整合与覆盖情况。若要做方法比较，需要独立证据（例如 small RNA-seq），暂不纳入 Phase 0。

## 1. 数据与划分

- **人类**：hg38（GRCh38）。真值来自 GENCODE，**chr20 是唯一的人类评测染色体**；chr21 只用于 smoke 和调通代码，不报告任何结论性数字。
  - 注意：GENATATOR 的 gene-finding 数据集把 hg38 chr8/20/21 作为验证集（不参与训练，但可能用于选模型），因此 chr20 对它不算完全盲测。在报告里如实说明。
- **时间切分**：
  - `T_old` = GENCODE v41（2022-07）
  - `T_new` = GENCODE v50（最新）
  - `Δ` = v50 中新增、且与 v41 任何外显子都不重叠的 lncRNA 基因
  - 主口径用 comprehensive（全部 isoform）；basic 作为敏感性分析。
- **lncRNA 的定义**：`gene_type == lncRNA`。其余基因按类别归档：protein_coding、small ncRNA（snRNA/snoRNA/miRNA/misc_RNA/scaRNA/vault_RNA/Mt_*/rRNA*）、pseudogene、其他。
- 坐标一律 1-based 闭区间（GFF 约定）。

## 2. 评价尺子（ruler v1，冻结）

预测与真值都转成转录本列表 `(chrom, strand, [exons])`，按类别（lncRNA / mRNA）分别评估。

| 指标 | 定义 |
|---|---|
| **L1 locus** | 预测转录本与真值基因在同链有 ≥1 bp 外显子重叠，即为命中。Sn = 命中的真值基因 / 真值基因数；Pr = 命中任一真值基因的预测转录本 / 预测转录本数 |
| **L2 intron-chain** | 多外显子转录本的内含子链与任一真值 isoform 完全一致。Sn 按真值基因统计（任一 isoform 被命中即算）；Pr 按预测转录本统计 |
| **L3 junction** | 内含子 (start, end, strand) 精确匹配。P/R 在去重后的内含子集合上计算 |
| **L4 exon** | 内部外显子精确匹配，首末外显子只要求内侧边界一致。P/R 在去重后的外显子集合上计算 |
| **L5 single-exon** | 单外显子真值与预测之间，同链 reciprocal overlap ≥ 0.5 即为命中 |
| **L6 base** | 同链外显子碱基 F1 |
| **T type** | 对命中的位点做类别混淆：真值 lncRNA 被预测为 mRNA 的比例，以及反向的比例 |

**真值不确定性口径**（方案 B 的核心）：
- 每个指标分别对 `T_old` 和 `T_new` 计算。
- 预测落在 `T_new` 之外、且不与任何类别的任何基因重叠时，单独记为 `novel`，不算入 Pr 的分母；同时报告"严格 Pr"（novel 计为 FP）和"宽松 Pr"（novel 剔除），两者都报告。
- **不完整度指标** `ΔFP→TP` = 在 `T_old` 下判为 FP、在 `T_new` 下变成 TP 的预测所占比例。

**统计**：按基因做配对 bootstrap（1,000 次）给出 95% CI。

## 3. 实验

### 0.1a GENATATOR 人类 sanity（物种无关，先做）

- 输入 hg38 chr21（smoke）和 chr20。
- 使用官方 pipeline 的默认生物学参数（`config.json`），固定模型版本并记录哈希。float32，单卡 4090。
- smoke 必须端到端跑通（包括 GFF 输出和 summary），并记录时间和显存峰值，据此外推 chr20 的耗时。
- 完成判据：`summary.complete == True`；GFF 合法；能用尺子打分。
- 作为 sanity 的预期：lncRNA L1 Sn 明显 > 0，并与论文报告的量级大体一致（论文写的是 T2T chr20 上的"正确分割基因数"，口径不同，只做量级比较）。

### 0.2 尺子实现与自检

- 用真值对真值自测：`T_new` 对 `T_new` 所有指标都应为 1。
- 用真值对真值的不完整度自测：`T_old` 当作预测、`T_new` 当作真值，得到 L1/L2 的 Sn，这就是"注释本身随时间的召回上限"。
- 再用 GENATATOR 的预测出具第一份正式数字。

### 0.1b GENATATOR 跨谱系（等物种面板确定）

- 每个物种 1–2 条染色体，打分用各物种的 Ensembl/RefSeq lncRNA 注释。
- **跨谱系退化判据**：某物种 lncRNA 的 L1 Sn 或 L2 Sn，低于该物种 mRNA 同指标的 50%，**同时**低于人类 chr20 lncRNA 同指标的 50%，则判为"明显退化"。用 mRNA 作对照，是为了排除"整条染色体都跑坏了"的情况。

### 0.3 编码 SOTA 在 lncRNA 位点的误报审计（等物种面板确定，人类可先做）

- 统计 Tiberius / ANNEVO / Helixer 的预测基因与真值 lncRNA 位点（同链外显子重叠，且不与任何 protein_coding 基因重叠）的重叠比例。
- 只报告描述性统计，不设 gate。

### 0.4 骨干逐层探针（GPU，0.1a 之后）

- 候选骨干：NT-v2-500M、GENERanno-eukaryote（coding 项目用过的版本）、Caduceus-PS（GENATATOR 用的）、PlantCAD2（仅当纳入植物时）。
- 标签：chr21 训练、chr20 评测；lncRNA 的外显子 / 内含子 / 供体 / 受体 / TSS / TES。
- 线性探针，复用 coding 项目的 `layer_probe.py`。报告每层的 AUPRC。
- 只用于选骨干；**不**作为论文主张。

### 0.5 unknown 边缘化最小验证（GPU，go/no-go 核心）

- 选定骨干后冻结，在其上训练小头（同一结构，同一 seed 集合 {42, 123, 456}）。
- 训练标签只用 `T_old`。
  - **F 臂**：未注释区当作背景。
  - **U 臂**：未注释区当作 unknown，做边缘化，复用 M30 的机制。
- 训练区域为 chr1–19 中的若干条染色体（具体在 0.4 后按算力确定，并在跑之前写进 v2 协议）。只用训练集内部验证选 checkpoint；chr20 每臂只评一次。
- **主终点**：chr20 上 `Δ`（v50 新增 lncRNA）的 L1 Sn。
- **成立判据**：U − F ≥ +5 个百分点，且 3 seed 均值的配对 bootstrap 95% CI 不含 0；**同时** U 在 `T_new` 下的严格 Pr 下降不超过 5 个百分点。

### 0.6 小 ncRNA 组件（CPU，与其他实验并行）

- 在 ont 上安装 Infernal + Rfam.cm（当前版本）、tRNAscan-SE 2.x。
- 在人类 chr20 和 4 个面板染色体上运行：cmscan 用 Rfam 的 GA 阈值，排除 lncRNA 家族；tRNAscan-SE 用真核模式。
- 用尺子（small 类：L1 locus、L5 reciprocal overlap）与 Ensembl/GENCODE 的小 RNA 对照。由于存在循环，**只报告覆盖率和运行时间**，作为统一 GFF3 的组件验收，不设 gate。
- 同时统计 GENATATOR 预测的 lnc_RNA 与小 RNA 位点的重叠（小 RNA 是否被误当成 lncRNA 外显子），作为合并规则的依据。

## 4. Gate

| 0.1b 跨谱系明显退化 | 0.5 U 边缘化成立 | 决定 |
|---|---|---|
| 是 | 是 | 全力推进方案 A（跨谱系统一语法 + 部分标签），B 作为评价章节 |
| 是 | 否 | A 的主卖点改为跨谱系迁移本身，unknown 边缘化降为分析 |
| 否 | 是 | 以 B + 方法（真值不完整下的训练）为主线，跨谱系只作展示 |
| 否 | 否 | 以 B（基准 + 评测）为主，重新讨论 |

## 5. 记录

每个实验在 `outputs/PHASE0/<exp>/` 下写 `RESULT.md` 和 `metrics.json`，同时复制到本地 `05_phase0/results/`。失败与"未执行"分开记录。
