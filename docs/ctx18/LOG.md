# 18 上下文重新考虑（2026-11-25）

主题短语：`ncRNA-context-models`

## 起因（作者 2026-11-25 原话要点）
- DNA-only 预测 lncRNA 改为「序列 + 细胞类型 / 环境」一起输入；
- 方向 D：参考 AlphaGenome 2（解决不同细胞类型下的表达），和 ChatGPT 商议；
- 方向 B：公开数据应不止 Sun 2019 一项研究，和 ChatGPT 梳理。

## AlphaGenome 2 要点（Claude 读博客与预印本前半部分）
- Gayoso … Avsec，bioRxiv 10.64898/2026.11.11.765857（2026-11-17）；博客 2026-11-11。
- 输入 1 Mb DNA + 条件（本体字典，或 scVI 式编码器把几十个细胞的单细胞谱映射成条件嵌入）+ 物种；输出链特异 CAGE / RNA-seq 覆盖、剪接、基因 TPM（GENCODE v46，含 lncRNA 基因，摘要来源）、染色质、甲基化、接触图。
- 1.07 亿单细胞谱训练编码器；TenK10K sc-eQTL 在 90% 精度下召回 30%（AG1 16%）。
- 非商业 API 可用；权重是否公开未核实。未见专门的 lncRNA 评估。

## 对话
- ChatGPT 深度研究（ncRNA项目内）：https://chatgpt.com/g/g-p-6b00508701cc8191a5d833658a0c32e4-ncrnaxiang-mu/c/6b06a62b-5f54-83ed-9431-5bea9bab5e7d
- 提示：`PROMPT_ChatGPT_EN.txt`（带项目背景，非中立轮；作者这次只要求 ChatGPT）。

## R1 结果（2026-11-25）
- 报告完成（7 分钟，31 引用）。窗格隐藏 → 导出失败，存要点：`deep_research/20261125_ncRNA-context-models_ChatGPT_R1_points.md`；原文在对话第二条回复。
- 主持人判断与核对：`20261125_ncRNA-context-models_Claude_R1.md`。
- 排序：D 改为「测量方式感知的上下文定位交互模型」第 1；条件化先验第 3；B 第 6。
- 下一步待作者：B 零成本审计（5 MB）→ D 重新预研（TS total-RNA 21 GB）→ AG2 先验测试（需作者注册 API）。

## B 审计与 D1（2026-11-25）
- B：GSE236224 只有 262 个基因、4 个 lncRNA → B 作为 lncRNA 方法停止（`B_AUDIT_GSE236224.md`）。
- D1（TS total-RNA，冻结 d99d867）：G0 过（0.456）、G1 过（57.6%，但蛋白编码对照 68.9%）、**G2 不过**（最佳中位 r 0.150，灰区；R² 为负）→ 按协议停 D（`RESULT_D1.md`）。

## Decima 诊断 DEC1–DEC3（2026-11-25）
- DEC1（Decima 自身数据，94e6666）：lncRNA 差距完全由表达低解释（校正后 +0.070）；Decima 只覆盖 11.6% 的 lncRNA（`RESULT_DEC1.md`）。
- DEC2（TS total-RNA）：H1 过（校正后 −0.108）；H2 不过（核富集越强越好）；S0 偏离 1 后过（0.992）；H3 不过（−0.011，事后放行）→ 形式 blind_spot，仅 H1 支撑（`RESULT_DEC2.md`）。
- DEC3（ENTEx total vs GTEx poly(A)，冻结 4391be8）：R1、R2 都不过 → not_replicated；但蛋白编码中位 r ≈ 0，检验本身无分辨力，留出 lncRNA 只有 25 个（`RESULT_DEC3.md`）。
- 综合：Decima 盲区只在单一供体 TS 上见到，机制假设不成立，第二套数据无法复核 → 按协议不立项，交作者确认关闭。
