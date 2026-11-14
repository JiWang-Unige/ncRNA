# PLAN DEDEND：专用端点模型能否在 ±20 kb 内选出真实的基因端点（草案，未冻结，2026-11-14）

> 起因：讨论稿把「专用起始 / polyA 模型能否选出真实端点」列为对主说法最直接的下一步实验（`DRAFT_intro_discussion.md`）。END / REP20 只说明**冻结的通用 DNA LM 特征 + 轻探针**不行：lnc 起点只在高 CpG 启动子上可定位，lnc 终点不可定位。
> 要回答：**换成专门为转录起始 / 3′ 端加工训练的序列模型，结论会变吗？** 三种结果都能写：(a) 专用模型能定位低 CpG lnc 起点或 lnc 终点 → 正文结论收窄为「通用 LM 特征缺这一信号」，讨论改成「把专用模型接入解码器是明确的改进路线」；(b) 专用模型也不行 → 「端点限速」从特征问题上升为「仅凭局部序列难以选出基因端点」的更强说法；(c) 灰区 → 描述性报告。
> 性质：只做推理，不训练、不调参；读数与 END 完全相同（Loc1k，±20 kb）。

## 1. 模型

| 用途 | 模型 | 来源 | 训练数据 / 留出 |
|---|---|---|---|
| 起点 | Puffin（Dudnyk et al., *Science* 2024），FANTOM_CAGE 输出，基因所在链，碱基分辨率 | github.com/jzhoulab/puffin（权重 `resources/puffin.pth` 2 MB；仅限学术使用） | FANTOM/ENCODE CAGE、RAMPAGE、GRO/PRO-cap；**测试染色体 chr8、chr9**，验证 chr10，其余训练 |
| 终点 | PolyaID + PolyaStrength（Stroup & Ji, *Nat Commun* 2023） | github.com/zhejilab/PolyaModelsHuman（`PolyaID.h5` 132 MB、`PolyaStrength.h5` 4 MB，Git LFS；未声明许可证） | 3′READS polyA 位点；论文写的是随机划分，**未见按染色体留出**，所以任何染色体都可能在训练集内 |

## 2. 染色体与真值

- **Puffin 主读数：chr8 + chr9**（Puffin 的测试染色体）。chr9 已被我们的 v2.0 检出读过一次，这里不涉及我们的模型，只读 Puffin。
- **Puffin 次读数：chr18、chr20**（在 Puffin 训练集内），用于和 END/REP20 的探针在同一批基因上并排比较；明确标为「样本内，偏向 Puffin」。
- **PolyaID：chr18（开发）+ chr20（读一次）**，与 END/REP20 的探针同一批基因并排。样本内性质如实写出。
- 真值：GENCODE v50 Ensembl_canonical lncRNA / protein_coding 转录本；证据：FANTOM5 CAGE（起点 ±200 bp）、PolyASite 2.0（终点 ±100 bp）；CpG 四分位按每组染色体的 lnc 起点分别算。规则与 `end_metrics.py` 相同。

## 3. 打分（事前规定，不改）

- **Puffin**：对每个基因，在其链上取 [TSS − 20 kb, TSS + 20 kb] 的 FANTOM_CAGE 预测，argmax 位置到真值起点 ≤ 1 kb 为命中（Loc1k）。
- **PolyaID**：对每个位置，以它为中心的 240 nt 窗口（基因所在链）给出 PolyaID 概率与 PolyaStrength。主得分：PolyaID ≥ 0.75 的位置中 PolyaStrength 最大者；若窗口内没有 PolyaID ≥ 0.75 的位置，取 PolyaID 概率最大者。命中定义同上（到真值终点 ≤ 1 kb）。次得分：PolyaID 概率（并列时取 PolyaStrength 更大者）。
- 随机期望 ≈ 0.05。

## 4. 判读规则（冻结后不改）

| 编号 | 读数 | 结论 A（专用模型能做到） | 结论 B（专用模型也做不到） | 之间 |
|---|---|---|---|---|
| T0 有效性 | Puffin，有 CAGE 的 mRNA 起点 Loc1k（chr8+9） | — | — | < 0.60 则 Puffin 的读数视为无效，只作描述 |
| T1 | Puffin，低 CpG（Q1+Q2）lnc 起点 Loc1k（chr8+9） | ≥ 0.25 | ≤ 0.12 | 灰区 |
| T2 | Puffin，有 CAGE 的 lnc 起点 Loc1k（chr8+9） | ≥ 0.50 | ≤ 0.25 | 灰区 |
| P0 有效性 | PolyaID，有 polyA 的 mRNA 终点 Loc1k（chr18、chr20 各自） | — | — | 只描述，不设门槛（mRNA 终点同样有多个 PAS） |
| P1 | PolyaID，有 polyA 的 lnc 终点 Loc1k，chr18 与 chr20 **都**满足 | ≥ 0.30 | ≤ 0.15 | 灰区 |

正文写法：
- T1 = A → 删去「只有 CpG 岛型启动子可定位」中的「仅凭序列」含义，改为「通用 LM 特征下」；讨论把 Puffin 列为可接入的组件。T1 = B → 第 15 节加一句「即使专为转录起始训练的模型，在留出染色体上也只能定位 x% 的低 CpG lnc 起点」。
- P1 = A → 「lnc 终点不可定位」限定为「在通用 LM 特征下」，讨论把 PolyaID 列为组件。P1 = B → 加一句「专用 polyA 模型在 ±20 kb 内同样选不出 lnc 终点」，并说明这可能反映的是同一区域内多个可用 PAS，而非序列不含 polyA 信号。
- 灰区：描述性报告，不加实验。

## 5. 计算与下载（需用户同意）

- 下载：上表两个仓库（共约 3 MB 代码）+ PolyaID / PolyaStrength 权重（136 MB）；Puffin 的 Python 依赖（selene_sdk、torch_fftconv、docopt，装进新的 `~/envs/puffin`，复用已有 torch 版本）；PolyaID 的新环境 `~/envs/polya`（TensorFlow 2.19，conda-forge + pip，约 2–3 GB）。不改动任何已有环境（TE 的环境尤其不碰）。
- GPU：Puffin 全染色体前向（chr8、9、18、20，两条链）预计 < 1 GPU 小时；PolyaID 只算端点 ±20 kb 的并集（每条染色体约 4–6 千万个窗口），预计每条染色体 1–2 GPU 小时。全部经 `gpu_gate.sh`，给 TE 让路。磁盘：输出只保留每个端点窗口的 argmax 位置与得分，< 1 GB。
- 产出：`RESULT_DEDEND.md`，`docs/phase1_results/DEDEND/*.json`；草稿第 15 节与讨论按第 4 节更新。
