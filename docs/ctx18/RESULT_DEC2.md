# RESULT_DEC2：Tabula Sapiens total-RNA 上的 Decima lncRNA 盲区（2026-11-25，只读一次）

- 协议 `PROTOCOL_DEC2.md`（含偏离 1），脚本 `scripts/dec2.py`。输出 `DEC2.json`（终版）、`DEC2_firstread.json`（首读，S0 按原设定未过、H3 未算）。
- 对应上的上下文与基因数见 `DEC2.json`（`n_ctx`、`n_genes_ok`）。

## 结果
| 项 | 结果 | 判定 |
|---|---|---|
| **H1** 留出基因：中位 r lncRNA / 蛋白编码 | 0.073（n 390）/ 0.242（n 3,324） | |
| H1 校正后 is_lnc 系数 | **−0.108**（CI −0.129 至 −0.085）；train −0.091 | **过** |
| **H2** 核富集每 IQR 效应（lncRNA） | **+0.050**（CI +0.036 至 +0.064）；蛋白编码 +0.079 | 不过（方向相反：越核富集越好） |
| **S0**（偏离 1 后）自预测 rep0 vs 公开 v1_rep0 | 中位 0.992（n 882）；vs 集成 0.947 | 过 |
| **H3** 表外 vs 表内留出 lncRNA（均为 rep0） | 中位 r 0.031（n 2,780）vs 0.056（n 386）；校正后 −0.011（CI −0.032 至 +0.010） | 不过（**事后放行**） |

- 重跑后 H1、H2 与首读完全一致（脚本自检 `H1 same True H2 same True`）。
- 决策（按协议「任一通过」）：**blind_spot**，仅由 H1 支撑。

## 读法
- 在 TS total-RNA 上，Decima 对 lncRNA 的跨细胞类型特异性预测明显差于蛋白编码，且校正表达、变异度、长度后仍差约 0.11。这与 DEC1（Decima 自己数据上校正后不差）形成对比。
- 但 H2 不支持「核内 RNA 是原因」：核富集越强，预测反而越好（蛋白编码同向且更强）。
- H3：表外 lncRNA 并不比表内 lncRNA 更差（差距由表达低解释）。
- 局限：单一供体；Smart-seq 全长平台 vs Decima 的 10x poly(A)；H1 差距可能来自平台 / 供体，而非 total-RNA 本身 → 交 DEC3 在第二套数据上复核。
