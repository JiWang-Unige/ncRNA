# PROTOCOL_DEC2：用 Tabula Sapiens total-RNA 的真实测量检验 Decima 在 lncRNA 上的盲区（2026-11-25）

作者 2026-11-25「做 DEC2」。本协议在计算任何 Decima 与 TS 之间的比较之前冻结。此前只看了元数据：TS 73 个上下文的名字、Decima 健康轨道的 organ / cell_type 名字，以及 Decima 包的源码。

## 1. 数据
- **实测**：D1 的 `pseudobulk.npz`（TS total-RNA，73 个上下文，Cells 与 Nuclei 两种 CPM；单一供体）。只用 Cells 的 log1p(CPM) 作实测；Nuclei 只用来算核富集。
- **Decima 预测**：
  - 已在 Decima 基因表内的基因：直接用补充文件 1 的 `layers/preds`（作者给出的预测）。
  - 不在表内的基因：用开源权重 rep0（HF `Genentech/decima-model`，`rep0.safetensors`）自己预测。窗口按 Decima 的约定构建：524,288 bp，TSS 位于第 163,840 位，负链取反向互补，基因掩码从 TSS 覆盖到基因末端（截到窗口末端）。坐标取 GENCODE v50 基因层面的起止。
- **对应关系**（写死在脚本 `ORG`、`CT` 两个字典里）：
  - TS 组织 → Decima organ：Bladder→urinary，Large_intestine / Small_intestine→gut，Lung→lung + airway，Mucosa / Salivary_gland / Tongue→mouth，Muscle→muscle，Prostate→prostate，Vasculature→vasculature，Testis 不对应。
  - TS 细胞类型 → Decima cell_type 按 CL 名字逐一对应，见脚本；无对应的（如 Pericyte、Spermatid）去掉。
  - 只用 Decima 的 healthy 轨道。一个 TS 上下文对应多条轨道时，取这些轨道预测的平均。
  - 对应上的上下文少于 20 个 → 停，报告。

## 2. 指标
- 每个基因：跨对应上下文的 Pearson（Decima 预测 vs TS 实测 log1p CPM），即「跨细胞类型的特异性」。
- 基因纳入：TS biotype 为 lncRNA 或 protein_coding；在 ≥ 10 个对应上下文中 CPM_cell ≥ 1；实测 SD > 0。
- 协变量（都来自 TS 实测）：log 平均 CPM、实测 SD、log 基因长度。
- 核富集 μ_g：73 个上下文里 y = log2((CPM_nuc+1)/(CPM_cell+1)) 在观测条件（(CPM_nuc+CPM_cell)/2 ≥ 2）下的平均，至少 10 个上下文。

## 3. 假设与门槛（基因 bootstrap 2000 次，95% CI）
- **S0 管线核对（前提）**：对 GPU 集合里已在表内的基因，自己用 rep0 预测的 8,856 条轨道与作者给出的预测之间的 Pearson，中位数 ≥ 0.95。不过 → 停，查管线，不读 H3。
- **H1 total-RNA 下的 lncRNA 差距**：Decima 留出基因（test + val）上，r ~ is_lnc + cr(log 平均, 4) + cr(SD, 4) + log 长度。is_lnc 系数 ≤ −0.05 且 CI 上限 < 0 → 过；−0.05 < 系数 ≤ −0.02 且 CI < 0 → 灰区；否则不过。train 基因另报，不进门槛。
- **H2 核富集越强越差**：表内所有 lncRNA（含 train，加 C(dataset) 协变量）上，r ~ μ_g + 同样协变量。μ_g 系数 × μ_g 四分位距 ≤ −0.05 且 CI 上限 < 0 → 过；同一模型在蛋白编码上另报作对照。
- **H3 表外 lncRNA**：GPU 集合里，表外 lncRNA（随机抽 ≤ 3,000 个，seed 0）与表内留出 lncRNA（全部纳入的）比较，两者都用 rep0 自己预测。r ~ is_out + 同样协变量。系数 ≤ −0.05 且 CI < 0 → 过。
  - GPU 集合另加 500 个表内留出蛋白编码基因（seed 0），用于 S0 和报告。
- 门槛判定只看上述系数；中位数只报告。

## 4. 决策
- H1、H2、H3 任一通过 → Decima 在 total-RNA / lncRNA 上的盲区成立，写下一步计划（在开源 Decima 上加入 total-RNA 数据微调；算力和数据需作者批准）。
- 全部不过（含灰区）→ Decima 盲区这个切入点关闭，报告作者。
- 已知局限：单一供体；TS 是 Smart-seq 类全长 total-RNA，Decima 训练多为 10x poly(A)，平台差异本身就是我们要测的东西，但也会带来与 lncRNA 无关的整体偏差（由 is_lnc 对照蛋白编码来扣除）。
- 结果只读一次；脚本运行前冻结（sha256 见 `PROTOCOL_DEC2.sha256`）。
