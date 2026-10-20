# RESULT 0.3（部分）：编码 SOTA 在 lncRNA 位点上的表现——人类 chr20、Helixer（2026-10-20）

> 状态：**部分完成**。只有 Helixer，而且只覆盖人类。
> - Tiberius 和 ANNEVO 在 hg38 chr20 上的现成预测不存在：旧项目只留下日志。需要补跑，属于 GPU 任务，排在 0.1a 之后，计划走 Kaggle 或 ont。
> - 面板物种的审计等 0.1b 完成后再做。

## 输入

- 预测：Helixer 0.3.6，模型 `vertebrate_v0.3_m_0080`，hg38 chr20。
  - 复用旧项目的输出：`coding_rna_ab_initio/runs/helixer_Homo_sapiens_7945934/`（2026-05-15）。
  - 已复制到 `ncrna-annot/data/baselines/human/`。
- 真值：GENCODE v50 chr20。打分用 ruler v1.2，按 mRNA 类计算。

## 结果

| 指标 | 数值 |
|---|---|
| Helixer 预测的 mRNA | 531 |
| 编码基因 L1 Sn / 严格 Pr | 0.920 [0.896, 0.942] / 0.940 [0.919, 0.959] |
| 编码基因 L2 intron-chain Sn / Pr | 0.245 / 0.257 |
| **被任一 Helixer 预测触及的 lncRNA 基因** | **31 / 955（3.2%）** |
| 其中，该预测不碰任何编码基因（等于把 lncRNA 报成 mRNA） | 14 个 lncRNA 基因 |
| Helixer 位点级假阳性（不命中任何编码基因的预测） | 32 个：**18 个只落在 lncRNA 上（56%）**，8 个落在假基因上，6 个落在完全未注释区 |

## 解读

1. 对这一个代表性 SOTA 编码注释器来说，**约 97% 的 lncRNA 位点完全不可见**。这为"编码 SOTA 不关注非编码基因"提供了定量依据，比"从未关注"这种措辞更站得住。
2. 它犯错时，**一半以上的位点级假阳性是把已注释的 lncRNA 报成了编码基因**。这与 ANNEVO 自报的"18.3% 假阳性与已注释非编码转录本重叠"方向一致（口径不同，不能直接比较）。这正是"编码 ⊗ lncRNA 乘积语法，让两类在配分函数里竞争"的动机。
3. 局限：
   - 只有一个工具和一条染色体，样本小（18 / 32）。
   - Helixer 的 vertebrate 模型并非最新版本。
   - 结论需要 Tiberius 和 ANNEVO 来补强。
