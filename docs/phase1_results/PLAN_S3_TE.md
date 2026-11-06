# S3-TE 计划：TE 来源 lncRNA 的分层检出（定稿，2026-11-06，计算之前）

> 来源：用户经看板维护会话转来（「直接传回主对话把第一点做了」）。
> 只作描述：不重训，不新读模型，只用已有的读取结果。chr20 / chr9 不用于任何调参。
> 外部下载（UCSC hg38 `rmsk.txt.gz`，155.6 MB）走决定卡 `s3-te-rmsk`，批准后才下载与计算。
> 背景：S3 只把「软屏蔽碱基占比」当作一个位点特征检验过（−0.02，不相关）。软屏蔽包括简单重复和低复杂度序列，不等于 TE；TE 来源的 lnc 从未作为单独的层报告。

## 1. 数据

- **TE 区间**：UCSC hg38 `database/rmsk.txt.gz`（RepeatMasker）。
  - 只取 repClass ∈ {LINE, SINE, LTR, DNA, Retroposon, RC}，不分链。
  - 不算 TE 的类：带「?」的不确定类、Simple_repeat、Low_complexity、Satellite、各类 RNA、Unknown。
- **真值**：GENCODE v50 的 lnc 基因（ruler 判为 lnc 的）。
  - chr20：955 个，用于 v1.9 的读取；
  - chr9：1,425 个，用于 v2.0 的读取。
- **预测**：已有读取，不重读：
  - chr20：v1.9 主模型 g1b-42 / 123 / 456-ke（`outputs/PHASE1/P15/chr20`），τ = 0.6；
  - chr9：v2.0 的三个模型（`outputs/PHASE1/V20/read/chr9`），校准概率 ≥ τ_v20 = 0.36，与 `v20_score.py` 的报出规则相同。

## 2. 分层（计算前写定）

对每个真值 lnc 基因：

- **外显子 TE 占比 f**：把所有转录本的外显子合并，计算其中落在 TE 区间内的碱基比例。
- **TSS 在 TE 内**：Ensembl_canonical 转录本的 TSS 碱基落在 TE 区间内。没有 canonical 标记时，取该基因最长的转录本。

| 分层 | 定义 |
|---|---|
| 主分层（三层） | TE 少：f < 0.1；部分：0.1 ≤ f < 0.5；TE 多：f ≥ 0.5 |
| 次分层（两层） | TSS 在 TE 内 / 不在 |
| 「TE 来源」 | f ≥ 0.5，或 TSS 在 TE 内 |

## 3. 读数

每条染色体：

1. 各层的真值基因数；各层的 lnc L1 Sn（3 seed 平均）与漏检率（1 − Sn）。
   - 95% CI 按真值基因重抽 2,000 次，`default_rng(0)`，统计量为 3 seed 平均。
2. 「TE 多 − TE 少」与「TSS 在 TE 内 − 不在」的 Sn 差，以及差值的 CI（同一次重抽）。
3. 漏检基因中「TE 来源」的占比（三个 run 合并计算），与全部真值 lnc 中该占比的比值（富集）。
4. 补充：各层的外显子数、成熟长度中位数。用来判断差别是否只是长度造成的。

## 4. 判读（写定）

- **「TE 来源的 lnc 更难检出」**：chr20 与 chr9 上，「TE 多 − TE 少」的 Sn 差都 < 0，且 CI 上界都 < 0。
- 只在一条染色体上成立，记为「不一致」。
- 两条都不成立，记为「没有证据表明 TE 来源影响检出」。
- 次分层（TSS 在 TE 内）同样判读，单独报告。
- 无论结果如何都写进论文。如果「更难」成立，作为位点特征分析（S3）的补充：模型靠启动子 CpG 与长转录单元找位点，而 TE 来源的启动子不同于 CpG 岛启动子。

## 5. 产出与算力

- 代码 `scripts/phase1/s3/s3_te.py`，结果 `docs/phase1_results/RESULT_S3_TE.md` 与 `S3TE/`。
- 只用 ont CPU（nice 19，≤ 8 线程），约数分钟。不占 GPU，不影响正在跑的 v3.0 训练队列。
