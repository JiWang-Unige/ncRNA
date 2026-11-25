# PROTOCOL_DEC1：Decima 在 lncRNA 上有没有「盲区」（2026-11-25）

作者 2026-11-25「做 Decima 诊断」。本协议在计算任何分层统计量之前冻结。此前只看了元数据：obs / var 列名、gene_type 计数、fold / dataset 计数。查看表头时意外打印出 3 个蛋白编码基因的 `pearson` 值（STRADA、ETV4、USP25），不涉及 lncRNA，也不涉及任何分层。

## 1. 数据
- Decima 补充文件 1（Zenodo 15092691，`Supplementary_file_1.h5ad`，658.6 MB，md5 0391874d…，与 Zenodo 一致），存于 ont `~/ncrna-annot/data/ctx18/decima/`。
- 8,856 个伪 bulk × 18,457 个基因；X 为实测（log1p），`layers/preds` 为 Decima 预测；var 有 `gene_type`、`dataset`（train / val / test）、`mean_counts`、`gene_length`、`frac_nan`、`pearson`。
- 基因类型：protein_coding 16,344，lncRNA 2,084，其余 29 个不用。

## 2. 指标
- 主指标：每个基因跨伪 bulk 的 Pearson（预测 vs 实测），即「细胞特异性」能力。用 X 与 preds 重算，并与 var `pearson` 核对（差异中位数应 < 0.01，否则报告并以重算值为准）。
- 主分析只用 **test 基因**（Decima 训练时未见），val 作复核。
- 协变量（都由实测计算或取自 var，与预测无关）：log 平均表达 `mean_counts`、实测跨伪 bulk 的标准差、log 基因长度、`frac_nan`。

## 3. 假设与门槛
- **H1 原始差距**（描述）：test 基因上 protein_coding 与 lncRNA 的 Pearson 中位数之差，基因 bootstrap 2000 次 95% CI。
- **H2 lncRNA 特有差距（主门槛）**：在 test 基因上做 OLS：Pearson ~ is_lnc + 自然样条（log mean_counts，4 自由度）+ 自然样条（实测 SD，4 自由度）+ log 长度 + frac_nan。
  - is_lnc 系数 ≤ −0.05，且基因 bootstrap 95% CI 上限 < 0 → **通过**：Decima 在 lncRNA 上有不能被表达量和变异度解释的盲区。
  - 系数在 (−0.02, 0] 或 CI 跨 0 → **不过**：lncRNA 预测差主要是因为表达低、变异小，不是 lncRNA 特有。
  - −0.05 < 系数 ≤ −0.02 且 CI < 0 → 灰区，只报告。
  - 复核：同样模型在 val 基因上方向一致（不进门槛）。
  - 稳健性（不进门槛）：按 mean_counts 十分位内配对比较中位差。
- **H3 覆盖缺口**（描述）：Decima 纳入的 lncRNA 占 GENCODE lncRNA 的比例 vs 蛋白编码的比例。分母用 TS total-RNA var 的 biotype 计数（lncRNA 18,039，protein_coding 20,017），是近似值。
- **探索（不进门槛）**：lncRNA 按与最近蛋白编码基因的关系分组（重叠且反向链 = 反义；10 kb 内 = 邻近；其余 = 基因间），只在 Decima 自己的 18,457 个基因坐标内计算，报告各组中位 Pearson。

## 4. 决策
- H2 通过 → 写下一步计划：用 TS total-RNA 检验盲区是否在非 poly(A) / 核内 lncRNA 上更大，并评估在开源 Decima 上微调的可行性（需作者批准算力）。
- H2 不过 → 「Decima 的 lncRNA 盲区」作为切入点不成立（差距只是表达量问题），报告作者，不再追。
- 灰区 → 只报告，交作者决定。
- 结果只读一次；脚本运行前冻结（sha256 见 `PROTOCOL_DEC1.sha256`）。

## 偏离 1（2026-11-25，读结果之前）
- 首次运行在第一步崩溃：`a.X` 以稀疏矩阵存储，`np.asarray` 未转成稠密数组（IndexError，未产生任何统计量）。
- 修复：稀疏矩阵先 `.toarray()`。门槛、指标、协变量都不变；重新冻结。

## 偏离 2（2026-11-25，读结果之前）
- 第二次运行仍崩溃：补充文件 1 **没有实测矩阵 X**，只有 `layers/preds`。另下载了 HF `Genentech/decima-data` 的 `metadata.h5ad`（3.27 GB，sha256 b6d1d743…，与 HF etag 一致），它同样只有预测层（preds、v1_rep0–3），没有实测值。公开文件里拿不到实测表达。
- 因此：
  - 主指标改为 var 中作者给出的每基因 `pearson`（Decima 自己算的预测 vs 实测、跨伪 bulk 的相关），不再重算。
  - 协变量「实测 SD」无法计算，主模型改为：r ~ is_lnc + cr(log mean_counts, 4) + cr(n_tracks, 4) + log 长度 + frac_nan（n_tracks 取自 var，是该基因有数据的轨道数）。
  - 敏感性（不进门槛）：再加 cr(预测 SD, 4)。预测 SD 来自模型本身，若 Decima 对 lncRNA 低估变化，会把真实盲区吸收掉，所以只偏向「不过」，不会制造假阳性。
- 门槛数值和决策规则不变；重新冻结。
