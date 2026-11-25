# B0-lite 协议：区室间 RNA 结构反应性变化是否可预测（v1，2026-11-25，下载完成、计算之前冻结）

> 来源：`16_new_methods_20261124/20261125_ncRNA-new-models_SYNTHESIS.md` 第 4 节；作者 2026-11-25：「Run pilot」。
> 状态：作者 2026-11-25 选择「D0 + B0-lite」：只用 108 MB 处理后文件，**不做 G1**。本文件在任何标签或模型结果计算之前冻结（sha256 见 `PROTOCOL_B0.md.sha256`）。若 B0-lite 通过 G2 + G3，再申请 B0-full 补 G1。

## 1. 数据核查结果（2026-11-25，只查元数据和一个 1.7 MB 文件的前几行）
- GSE117840（Sun et al., NSMB 2019，PMID 30886404）处理后文件共 108 MB：
  - 每个细胞（HEK293、mES）× 区室（染色质 ch / 核质 np / 胞质 cy）× 条件（in vivo / in vitro）各**一个已合并的反应性文件**（icSHAPE 分数，每个核苷酸一个值）；
  - HEK293 染色质 in vivo 只覆盖 4,686 个转录本。
- **处理后文件没有分重复的反应性。**原始数据（SRA PRJNA483400）每个区室有 2 个 DMSO（D）、2 个 in vivo NAI-N₃（N）、2 个 in vitro（T）文库：
  - HEK293 全部 18 个文库约 280 GB；
  - 只取 HEK293 核质与胞质的 in vivo（D + N，8 个文库）约 94 GB，共约 26 亿条读段，需要在 ont 上用 icSHAPE-pipe 重新处理（CPU 预计 1–2 天）。
- 所以 G1（用独立重复构造的差异谱是否可重复）**只能靠原始数据**；只用处理后文件做不了。

## 2. 两个可选方案
- **方案 B0-lite（只下 108 MB 处理后文件）**：
  - 能做：G0（规模）、G2 / G3（可预测性与方法增益，以合并值为真值）。
  - 做不了：G1。替代的弱检查是「in vitro 差异 vs in vivo 差异」和「人 HEK293 与小鼠 mES 同源转录本上区室差异的方向一致性」，但都不是重复。
  - 风险：差异里混有测量噪声时，G2 / G3 的「可预测」可能只是学到了覆盖度和噪声。
- **方案 B0-full（再下约 94 GB 原始数据并重新处理）**：G0–G3 全部可做，且能用独立重复构造输入和评价目标，避免模型学到参考值里的噪声。

## 3. 任务、基线、门槛（两个方案共用，冻结时按选定方案删改）
- 主任务：HEK293 核质 → 胞质（in vivo）；次任务：染色质 → 核质、in vitro → in vivo；mES 作独立数据集（G4）。
- 标签：两个区室都实际观测到的核苷酸（非 NULL），Δr = r_target − r_ref。
- 划分：按基因留出（5 折，按基因 ID 分组；去除与训练集共享 ≥ 50 nt 完全相同片段的转录本）。
- 子任务：B0-S（序列 + 目标区室身份）；B0-E（序列 + 参考区室实测值 + 目标区室身份）。
- 基线：零变化；区室均值偏移；只用参考反应性回归；LightGBM（参考反应性、局部 5-mer、GC、ViennaRNA 未配对概率、转录本类型、合法协变量）；小型 1D CNN（不用 RNA LM）。
- 门槛：
  - G0：≥ 300 个合格 lncRNA 基因且 ≥ 30,000 个两区室共同观测核苷酸；不足则改为通用 RNA 方法、单独报告 ncRNA 子集。
  - G1（仅 full）：独立重复构造的差异谱 Pearson ≥ 0.40；不过 → STOP。
  - G2：留出基因上，最佳模型差异 Pearson ≥ 0.30 且 ΔR² > 0。
  - G3：比最强简单基线 +0.05，基因块 bootstrap（2,000 次）95% CI 下限 > 0。
  - G4：mES 上方向一致（不进主判定，只限定结论范围）。
- 不使用 RASP v2.0 的插补值或 StructureImpute / PrismNet 的预训练权重（它们在这些数据上训练过）。

## 4. 资源
- B0-lite：下载 108 MB，ont CPU 半天，GPU ≤ 2 小时（小 CNN）。
- B0-full：下载约 94 GB（ENA），ont CPU 1–2 天重处理，再加 B0-lite 的计算。

## 5. B0-lite 冻结细节（2026-11-25 补充，优先于上文的通用描述）
- 文件：GSE117840 处理后 icSHAPE 分数（`data/context_pilot/icshape/`，12 个文件）。主任务用 293_np_vivo → 293_cy_vivo；次任务 293_ch_vivo → 293_np_vivo、293_cy_vitro → 293_cy_vivo；G4 用 mes_np_vivo → mes_cy_vivo。
- 转录本序列：GENCODE v41（ont 已有）按去版本号 ENST 取剪接后序列；只保留 v41 长度与 .out 文件长度列相同的转录本。转录本类型取 v41 transcript_type。
- 核苷酸：两个区室都非 NULL；去掉转录本首尾各 30 nt（icSHAPE 末端不可靠）。标签 Δr = r_target − r_ref。
- 划分：按 gene_id 的 5 折 GroupKFold（种子 0）。
- 输入：
  - B0-E：参考区室反应性；±12 nt 序列 one-hot；±25 nt GC；RNAplfold（W = 150，L = 100，u = 1）未配对概率；转录本类型；log 转录本长度；参考区室 RPKM（.out 第 3 列）。**不使用目标区室的任何量**（包括 RPKM）。
  - B0-S：同上但去掉参考区室反应性和参考 RPKM。
- 模型：
  - 简单基线：零变化；区室均值偏移（训练集 Δ 均值）；只用参考反应性的线性回归（仅 B0-E）。
  - 学习模型：HistGradientBoostingRegressor（max_iter = 300，learning_rate = 0.05，种子 0，不调参）；小型 1D CNN（输入 ±25 nt one-hot + 参考反应性轨道，3 层卷积，ont 4090，≤ 2 GPU 小时）。
- 评估：留出折的核苷酸级 Δr Pearson 与 R²（1 − SSE/SS）。全体与 lncRNA 子集分开报告。
- 门槛（B0-lite 版）：
  - G0：两区室共同观测的 lncRNA 基因 ≥ 300 且核苷酸 ≥ 30,000；不足则主结论改为通用 RNA，lncRNA 只作描述。
  - G2：主任务上最佳学习模型 Δr Pearson ≥ 0.30 且 R² > 0（全体）。
  - G3：最佳学习模型比最佳简单基线 Pearson 高 ≥ 0.05，基因块 bootstrap（2,000 次）95% CI 下限 > 0。
  - G4：mES 主任务上 G3 的方向一致（只限定结论范围）。
- 判定：G2 + G3 过 → 申请 B0-full 补 G1；G2 或 G3 不过 → B 方向停止（不下载原始数据）。
- 环境：特征和 HistGBR 用 ont `~/envs/rnafate`（sklearn 1.3.2）；RNAplfold 用 `~/envs/xaudit/bin/RNAplfold`；CNN 用 `~/envs/te_benchmark`（torch 2.5.1）。nice 19，≤ 16 核。
