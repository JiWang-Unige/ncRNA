# M1 实验 3 计划：AlphaGenome 留出 fold 模型能否定位 lncRNA 起点（冻结 2026-11-17）

> **冻结**：作者 2026-11-17 回复「批准并冻结」，同时批准第 7 节第 2–4 项（下载四个 fold 权重与 Borzoi 区间文件、Colab ≤ 10 CU）。第 2–5 节（fold 归属规则、真值、打分、判读）自此不改。实现细节以 `scripts/m1/ag/`（5e97f4a）为准；如有偏离，记在结果文档里，不改本文件。

> 依据：`docs/review/DECISION_20261117.md` 第 2 节实验 3（判据按 ChatGPT 第二轮修改）。
> 设计与 `PLAN_DEDEND.md`（Puffin，冻结 8deaaff）完全对齐，这样可以在同一批基因上配对比较。
> **状态**：草案。没有下载任何东西，没有产生任何 AlphaGenome 输出。批准事项见第 7 节。
> - 2026-11-17：作者已在 Hugging Face 接受四个 fold 仓库的条款；Baobab 上用作者的登录核对过，四个仓库的文件都能访问（fold-0 修订 4ce9fd40、fold-1 ebbc6ab6、fold-2 8b7e73c4、fold-3 75fa73c8，各约 0.75 GB）。
> - 代码已写好（`scripts/m1/ag/`）：`ag_fold_prep.py` 起点表与 fold 分配，`ag_fold_job.py` Colab 作业，`ag_fold_run.sh` Colab 包装（沿用 S2-AG 的会话、超时、CU 上限与清理），`ag_fold_metrics.py` 打分与判读。尚未运行。
> 只作对照：按 AlphaGenome 条款，输出不用于训练或改进我们的方法。

## 1. 问题

1 Mb 上下文、按细胞类型输出 CAGE / PRO-cap 的最强序列模型，能否在 ±20 kb 内把低 CpG lncRNA 的起点定位到 1 kb 以内？

- **能**：说明序列中确有起点信息，只是我们的骨干和 Puffin 没取到。这对应 DECISION 的结果 A（方向 1）。
- **不能**：「当前序列模型都定不了」这一结论扩展到最强的模型。

只测起点。CAGE / PRO-cap 只反映转录起始，不能用来测终点。

## 2. 模型与 fold

- **模型**：AlphaGenome 1（AG1），Hugging Face `google/alphagenome-fold-0` 到 `fold-3`，每个约 0.75 GB，文件结构与已用过的 all-folds 相同。
  - 官方代码 `alphagenome_research` 的 `create_from_huggingface` 只是换仓库名；我们的 `ag_colab_job.py` 用本地检查点目录加载，换目录即可。
  - 不用 AG2：它去掉了 PRO-cap 输出，而且没有核实是否有留出 fold 的权重。
- **fold 定义**（已在官方代码 `alphagenome/data/fold_intervals.py` 中核对）：沿用 Borzoi 的 8 个序列块（`sequences_human.bed.gz`，Borzoi 提交 5c93582）。每个模型留出 2 块：

  | 模型 | 验证块 | 测试块 |
  |---|---|---|
  | FOLD_0 | fold0 | fold1 |
  | FOLD_1 | fold3 | fold4 |
  | FOLD_2 | fold2 | fold5 |
  | FOLD_3 | fold6 | fold7 |

  all-folds 用全部 8 块训练，所以我们之前的 S2-AG 结果都是样本内的。
- **起点归属哪个模型**：看起点落在哪个 Borzoi 区间、该区间属于哪一块，就用留出这一块的模型。
  - **主分析**：起点在该模型的验证块或测试块内，且离任何训练块区间至少 430 kb。训练窗口是以区间中点为中心的 1 Mb，再加平移增强，标签最远可达约 428 kb。这是保守判据，官方没有规定。
  - **敏感性分析**：只用测试块。验证块可能用于早停。
  - 推理时输入窗口跨到训练块不算泄漏，因为输入序列不是标签。

## 3. 染色体与真值（与 DEDEND 相同）

- **主读数：chr8 + chr9**。这是 Puffin 的测试染色体，可以与 Puffin 在同一批基因上配对比较。
  - chr9 由 Borzoi 的 fold3（FOLD_1 的验证块）和 fold7（FOLD_3 的测试块）覆盖。
  - chr8 属于哪些块，要等下载区间文件后才能算出来。
- **次读数：chr18、chr20**。各块都由某个 fold 模型留出，可以与我们的探针和 Puffin 的样本内结果并排比较。
- **真值**：GENCODE v50 中 Ensembl_canonical 的 lncRNA 和 protein_coding 转录本起点。
  - 有无 CAGE 证据：起点 ±200 bp 内有 FANTOM5 CAGE 峰。
  - CpG 四分位：按每组染色体内的 lnc 起点分别计算。
  - 低 CpG = Q1 + Q2。
  - 以上规则与 `end_metrics.py` 和 DEDEND 相同。

## 4. 打分（事先规定）

- **对每个起点**：在基因所在链上取 [TSS − 20 kb, TSS + 20 kb]。若 argmax 位置离真值起点 ≤ 1 kb，记为命中（Loc1k）。随机期望约 0.05。
- **主得分**：基因所在链上**全部 CAGE 轨道的均值**（每条链 273 条，FANTOM5 的 hCAGE 和 LQhCAGE）。
- **次得分**（只描述）：
  - CAGE 轨道的最大值；
  - PRO-cap 轨道的均值（6 个样本）。
- **对照**（同一批基因，配对比较）：
  - Puffin：直接用 DEDEND 已有的逐基因结果，不重跑；
  - CpG 对照：±20 kb 内以 500 bp 滑窗计 CpG 个数，取 argmax；
  - 随机期望。
- **在虚拟机里只保存小结果**：每个起点的每条 CAGE / PRO-cap 轨道在窗口内的 argmax 偏移和最大值，以及三种汇总得分的 argmax。不保存逐碱基的全部轨道。

## 5. 判读规则（冻结后不改）

| 编号 | 读数（chr8 + chr9，主分析集） | A：序列中有起点信息 | B：最强模型也定不了 | 之间 |
|---|---|---|---|---|
| G0 有效性 | 有 CAGE 的 mRNA 起点 Loc1k | — | — | < 0.60 则只作描述 |
| G1（主） | 低 CpG lnc 起点 Loc1k | ≥ 0.25，且配对差值 AG − Puffin 和 AG − CpG 对照的 95% CI 下限都 > 0 | ≤ 0.12，且 95% CI 上限 < 0.25 | 灰区 |
| G2 | 有 CAGE 的 lnc 起点 Loc1k | ≥ 0.50 | ≤ 0.25 | 灰区 |

- 置信区间：对基因做 2,000 次自助抽样（种子 0）。
- 对照参考值：Puffin 在 chr8 + chr9 上，T1（低 CpG）为 0.086，T2（有 CAGE 的 lnc）为 0.395，mRNA 为 0.707。

## 6. 第二阶段（可选）：细胞类型依赖

要回答：用与 lncRNA 实际表达的细胞类型相匹配的轨道，起点是否更准？

- **只在同一种实验内比较**：都用 hCAGE，比较匹配与不匹配的细胞类型。不能用「最大值对均值」来论证。
- **做法**：
  - 用 FANTOM5 逐样本的 CAGE 峰表达矩阵，按 ontology 把样本聚成 AlphaGenome 的轨道组。
  - 对每个有 CAGE 证据的 lnc 起点，找出它表达的轨道组（匹配）和不表达的轨道组（不匹配）。
  - 分别计算 Loc1k，配对比较。
- **判据**：匹配减不匹配 ≥ 0.10 且 CI 下限 > 0 → 起点定位依赖细胞类型（DECISION 的结果 C 方向）。
- **需要另外批准下载**：FANTOM5 counts 矩阵约 111 MB、sdrf 180 KB、ontology 1.7 MB。
- 第一阶段的数据已经足够做第二阶段（每条轨道的 argmax 都存了），不需要重跑 AlphaGenome。

## 7. 计算与需要作者批准的事

- **硬件**：Colab G4（RTX PRO 6000 96 GB）。上次实测每个 1 Mb 窗口 0.34 s，峰值 31.8 GB。
  - 主分析集需要约 460 个窗口，推理不到 3 分钟；主要时间花在上传权重和加载模型（4 次）。
  - 预计 ≤ 8 CU，上限 10 CU（M1 总上限 15 CU）。
  - 沿用上次的 ssh / scp 流程、硬超时和 watchdog；用完清理会话。
  - ont 的 4090 显存（24 GB）不够跑 1 Mb 窗口。

| # | 需要作者做或批准的事 | 说明 |
|---|---|---|
| 1 | 在 Hugging Face 上**分别接受** `google/alphagenome-fold-0`、`-1`、`-2`、`-3` 的条款 | 四个仓库各自 gated；all-folds 已接受，不代表这四个已接受。只能你本人点 |
| 2 | 下载四个 fold 权重（共约 3 GB） | 像上次一样在 Baobab 登录节点用你的 HF 登录下载，转存到 ont（磁盘剩 91 GB），核对 sha256 |
| 3 | 下载 Borzoi 区间文件 `sequences_human.bed.gz`（约几百 KB） | 来自 GitHub calico/borzoi 的固定提交 |
| 4 | Colab 花费 ≤ 10 CU | — |
| 5 | （可选）第二阶段的 FANTOM5 下载（约 113 MB） | — |

批准第 1–4 项后：先写代码，并在 chr18 的 10 个起点上做烟雾测试（检查加载、链和偏移）；然后冻结本计划；最后一次性读取全部染色体。
