# RESULT_DEC3：第二套 total-RNA 复核（ENTEx total vs GTEx poly(A)，2026-11-25，只读一次）

- 协议 `PROTOCOL_DEC3.md`（687a4014…），脚本 `scripts/dec3.py`（0a0784e8…），冻结提交 4391be8。输出 `DEC3.json`。
- 15 个 organ 全部对上（每 organ 28–1,141 条 Decima 轨道）。

## 结果
| 项 | 蛋白编码 | lncRNA |
|---|---|---|
| 留出基因数（test + val） | 2,365 | **25** |
| 中位 r_tot（Decima vs ENTEx total） | −0.029 | −0.093 |
| 中位 r_pA（Decima vs GTEx poly(A)） | −0.045 | −0.145 |
| train 中位 r_tot / r_pA | −0.023 / +0.044 | −0.009 / +0.056 |

| 门槛 | 系数（95% CI） | 判定 |
|---|---|---|
| R1 total 下 is_lnc | −0.011（−0.121 至 +0.110） | 不过 |
| R1 poly(A) 另报 | −0.042（−0.167 至 +0.085） | — |
| R2 (r_tot − r_pA) 下 is_lnc | +0.059（−0.053 至 +0.170） | 不过 |

- **按协议决策：not_replicated。**

## 读法（如实说明：这个检验本身没有分辨力）
- 蛋白编码基因的中位 r 在 total 和 poly(A) 下都接近 0（−0.03 / −0.04）。也就是说，在「organ 内所有轨道预测取平均 vs bulk 组织 TPM」这个比较里，Decima 连蛋白编码的跨器官变化都没有对上。协议没有事先设阳性对照，这是设计缺口。
- 可能原因（未检验，不再追）：Decima 轨道按细胞类型伪 bulk，organ 内取平均后被各 study 的细胞组成主导，不等于 bulk 组织；15 个 organ 的跨器官相关自由度也小。
- 在两种测量下都满足表达门槛、又在 Decima 留出集里的 lncRNA 只有 25 个，CI 宽到 ±0.11。
- 所以 DEC3 **既没有复现、也没有否定** DEC2-H1。形式结论按协议记为 not_replicated；实际含义是「第二套数据无法检验」。
- 综合 DEC1–DEC3：Decima 的 lncRNA 差距只在单一供体的 TS total-RNA 上看到（H1 −0.108），机制假设（核富集 H2、表外基因 H3）都不成立，第二套数据无法复核。按协议，Decima 盲区不立项。
