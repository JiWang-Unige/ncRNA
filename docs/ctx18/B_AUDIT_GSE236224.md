# B 零成本审计：GSE236224（2026-11-25）

- 文件：`GSE236224_Normalized_reactivity.csv.gz`（5.1 MB），存在 ont `~/ncrna-annot/data/ctx18/b_audit/`。
- 内容：
  - 共 169,073 行，四个阶段：H9 54,053 行、D0 46,473、D1 37,595、D7 30,952；
  - 只有 **262 个基因**，其中蛋白编码 241 个、未匹配 16 个、**lncRNA 4 个**（DANCR、GAS5、SNHG5、MALAT1）；
  - 两个以上阶段都有的基因里，lncRNA 只有 1 个。
- 原始数据（9 个样本，12.5 GB，hESC / NPC 的 DMS 与对照）每个样本约 1.2 GB，测序深度浅；即使重新处理，lncRNA 覆盖也不可能到 100 个。
- 结论：按 ChatGPT R1 与主持人约定的门槛（≥ 100 个有重复的 lncRNA），**方向 B 作为 lncRNA 方法停止**。
  - Sun 2019 只有 16–39 个，这里 4 个；巨噬细胞 smartSHAPE 要 267 GB 且没有对照，不再追。
  - 低丰度 lncRNA 在全转录组结构探测里覆盖不足，是这类数据的普遍局限（RASP v2.0 需要靠插补，PMID 39546630）。
