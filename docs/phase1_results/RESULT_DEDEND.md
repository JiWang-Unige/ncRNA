# RESULT DEDEND：专用端点模型能否在 ±20 kb 内选出真实的基因端点（2026-11-14）

> 计划：`PLAN_DEDEND.md`（冻结 8deaaff，用户「冻结并下载」）。代码 `scripts/phase1/dedend/`（`puffin_run.py`、`polya_run.py`、`dedend_metrics.py`、`dedend_queue.sh`）。原始输出 `docs/phase1_results/DEDEND/DEDEND.json`；ont `outputs/PHASE1/DEDEND/`。
> 模型：Puffin（jzhoulab/puffin 9bc81eb，权重 2 MB，14 个参数张量全部匹配）；PolyaID + PolyaStrength（zhejilab/PolyaModelsHuman d70bcd0，权重 sha256 与 LFS 记录一致）。新环境 `~/envs/puffin`（复用 generanno 的 torch）、`~/envs/polya`（TF 2.19，Python 3.12）。
> 运行：ont 4090，GPU 闸门。Puffin 四条染色体约 6 分钟；PolyaID chr18 + chr20 约 20 分钟（每条链 1,600–2,400 万个位置）。
> 实现偏离（不影响设计）：Puffin 分块从 2 Mb 改为 200 kb（FFT 卷积在 2 Mb 时超出 7 GB 显存上限，第一次运行失败）。冒烟：Puffin 在 chr18 40 个有 CAGE 的 mRNA 起点上 Loc1k 0.85（中位 55 bp）；PolyaID 在多数 mRNA 终点附近 100 bp 内给出 > 0.9 的概率。

## 结论（按冻结规则）

| 编号 | 读数 | 值 | 规则 | 结果 |
|---|---|---|---|---|
| T0 | Puffin，有 CAGE 的 mRNA 起点 Loc1k（chr8+9，n = 1,316） | **0.707** | ≥ 0.60 才有效 | 有效 |
| T1 | Puffin，低 CpG（Q1+Q2）lnc 起点 Loc1k（chr8+9） | **0.086** | A ≥ 0.25；B ≤ 0.12 | **B：专用起始模型也定位不了低 CpG lnc 起点** |
| T2 | Puffin，有 CAGE 的 lnc 起点 Loc1k（chr8+9，n = 708） | 0.395 | A ≥ 0.50；B ≤ 0.25 | 灰区 |
| P0 | PolyaID，有 polyA 的 mRNA 终点 Loc1k | chr18 0.122；chr20 0.188 | 只描述 | — |
| P1 | PolyaID，有 polyA 的 lnc 终点 Loc1k | chr18 **0.071**；chr20 **0.119** | A 两条都 ≥ 0.30；B 两条都 ≤ 0.15 | **B：专用 polyA 模型也选不出 lnc 终点** |

## 细节

### Puffin（起点，Loc1k）

| 组 | lnc 全部 | lnc 有 CAGE | lnc 无 CAGE | CpG Q1 | Q2 | Q3 | Q4 | mRNA 全部 | mRNA 有 CAGE |
|---|---|---|---|---|---|---|---|---|---|
| chr8+9（Puffin 测试集） | 0.237（3,031） | 0.395（708） | 0.189 | 0.079 | 0.094 | 0.214 | **0.561** | 0.652 | 0.707 |
| chr18（Puffin 训练集） | 0.207（929） | 0.335（173） | 0.177 | 0.030 | 0.099 | 0.190 | 0.509 | 0.742 | 0.754 |
| chr20（Puffin 训练集） | 0.213（955） | 0.391（233） | 0.155 | 0.088 | 0.100 | 0.168 | 0.494 | 0.679 | 0.708 |

- 与我们的探针（chr18，有 CAGE 的 lnc 起点 0.16–0.41；CpG Q4 0.17–0.61；mRNA 0.41–0.79）在同一范围。专用模型在 mRNA 上更稳定，但 **CpG 梯度几乎一模一样**：Q1–Q2 接近随机，Q4 约一半。
- lnc 起点的中位偏差仍约 7 kb。

### PolyaID（终点，Loc1k）

| 染色体 | 得分 | lnc 有 polyA | lnc 全部 | mRNA 有 polyA | mRNA 全部 |
|---|---|---|---|---|---|
| chr18 | 主（PolyaID ≥ 0.75 中 PolyaStrength 最大） | 0.071（311） | 0.062 | 0.122（230） | 0.112 |
| chr18 | 次（PolyaID 最大，PolyaStrength 打破并列） | 0.116 | 0.086 | 0.187 | 0.165 |
| chr20 | 主 | 0.119（362） | 0.084 | 0.188（490） | 0.172 |
| chr20 | 次 | 0.152 | 0.108 | 0.304 | 0.288 |

- 窗口内约 5% 的位置 PolyaID ≥ 0.75（例如 chr20 负链 131 万 / 2,365 万）：能加 polyA 的位点到处都是，选出「基因用的那个」才是难点。
- mRNA 终点也只有 0.12–0.19（主得分），说明我们整模型里 mRNA 终点的精度（检出基因中 28–58% 在 1 kb 内）来自基因语法，而不是 polyA 信号。
- 次得分略高（chr20 mRNA 0.30），但判读按事前规定用主得分；次得分下 lnc 也只有 0.12–0.15。

## 对论文的影响（已改）

- 第 15 节标题改为「With current sequence models …」。加入 Puffin（留出染色体）与 PolyaID 的结果；终点段原先的「专用模型未测」改为实测结果；结论斜体句扩展为「通用 LM 特征以及专用起始 / polyA 模型」。方法加「Dedicated end models」，图 7 加 d。
- 引言：分解一句加上「与专用起始和 polyA 模型一样」。
- 讨论「What would close the gap」第二条：从「最直接的下一步实验」改为实测结论——直接接入专用模型不够；起始模型可作为 CpG 岛基因的起点得分；其余起点 / 终点的缺失信息在于「许多候选位点中用哪一个」，指向表达数据或按细胞类型条件化的模型。（这一句是解释，不是测出来的，已在作者备注中标出。）

## 局限

- PolyaID 未按染色体留出，chr18 / chr20 可能在其训练集中（偏向 PolyaID，结论方向不受影响）。
- 只测了每类一个专用模型；ProCapNet、CLIPNET、APARENT2 未测。
- Puffin 用的是 FANTOM_CAGE 输出；其他输出（RAMPAGE、PRO-cap）未测。
