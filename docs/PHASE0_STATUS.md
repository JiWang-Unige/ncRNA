# Phase 0 进度总览（2026-10-20）

| 实验 | 状态 | 一句话结论 | 报告 |
|---|---|---|---|
| 0.2 评价尺子 | ✅ | ruler v1.3（GENCODE / Ensembl / RefSeq GFF3、小 RNA、无链模式）经合成与真实数据自检。GENCODE lncRNA 从 v41 到 v50 翻倍（17,960 → 34,677），v41 只覆盖 v50 中 49% 的位点；mRNA 为 99.8% | `results/RESULT_0.2_ruler.md` |
| 面板真值概况 | ✅ | 各谱系 lnc : mRNA 真值比例相差 30 倍以上（0.05–1.74）；斑马鱼 GRCz11 仍用旧命名 | `results/RESULT_truth_density.md` |
| 0.1a GENATATOR 人类 | ✅ | 在它自己的 T2T 基准上，lncRNA 转录本召回 < 1%（官方 metric，k = 250），位点召回 3.5%；hg38 chr20 + GENCODE v50 下位点召回 0.028、mRNA 内含子链 0.708；瓶颈在区间发现 | `results/RESULT_0.1a_0.3_T2T_benchmark.md` |
| 0.3 编码 SOTA 审计 | ✅（人类） | Tiberius、ANNEVO、Helixer、AUGUSTUS 都不输出 lncRNA；会把 3–52 个 lnc 位点报成编码基因；Helixer 位点级假阳性中 56% 是 lncRNA | 同上，加 `RESULT_0.3_fp_audit_partial.md` |
| 0.1b 跨谱系 | ✅ 5 个物种全部完成 | lnc 在各谱系都约等于 0；mRNA 迁移在远缘谱系崩溃（线虫 0.006）；原 gate 没有区分度，已如实记录 | `results/RESULT_0.1b_cross_lineage.md` |
| 0.4 骨干探针 | ✅ | 选定 GENERanno-1.2b-cds（与 coding 项目、Carbon-A 同支）。"CDS 监督会压制非编码"的假设被推翻；lnc 外显子的序列信号弱（AUPRC 是先验的 2.4 倍），必须靠结构约束 | `results/RESULT_0.4_probe.md` |
| 0.6 小 ncRNA | ✅ | 经典结构 RNA 基本全召回；缺口在 miRNA、piRNA，以及 Alu/SRP 伪命中，已整理成整合规则 | `results/RESULT_0.6_smallrna.md` |
| 0.5 unknown 边缘化 | ✅ **gate 未通过** | 主臂 U-ev 比 F 高 +1.7 个百分点，CI 为 [+0.4, +3.2]，显著但未达到预注册的 +5。终点在最小解码下饱和（Sn 约 0.8，Pr 约 0.14）。探索性 AUPRC 显示两种 unknown 处理都稳定但小幅改善了新增 lncRNA 的排序（相对 +9–11%） | `results/RESULT_0.5_unknown.md` |

## 算力与花费（截至 2026-10-20）

- **Modal**（`ji-wang-unige`），估算约 14 美元；11:17 workspace 因额度用尽被停用（与 TE 项目共用）：
  - GENATATOR：chr21 一次，加面板 5 次
  - 小 RNA 扫描：约 42 CPU-h
  - 4 个骨干探针
- **Kaggle**（`jiunige`）：2 次 T4 smoke。
- **ont**：CPU 用于安装、数据准备和打分；GPU 只用于 0.5 的小头训练，共 57 分钟，是在 GPU 空闲后才开始的。
- **Baobab**：0 credits（只做了文件存储和同步）。

## 对论文主线的影响

- **方案 B 的骨架已基本成形**：
  - 真值在快速扩张（lncRNA 位点 4 年翻倍）
  - 各谱系真值差异巨大
  - 所有现有工具端到端几乎恢复不了 lncRNA，且编码工具会误报 lncRNA
  - 跨谱系迁移崩溃
- **方案 A 的门槛被重新定义**：只要在人类等内分布物种上把 lnc 召回做到明显不为 0（例如位点召回 > 20%）且 precision 可接受，就是领域首个可用结果；跨谱系是加分项，不再是前提。
- 0.5 是 A 是否采用"部分标签 / unknown"作为核心卖点的判定实验。

## Phase 0 gate 判读（按协议 §4，如实记录）

| 0.1b 跨谱系明显退化 | 0.5 U 边缘化成立 | 协议给出的决定 |
|---|---|---|
| 形式上成立（4 个物种中 3 个），但没有区分度：人类自身 lnc 召回就只有 0.028 | **不成立**（+1.7 个百分点，低于 +5） | "A 的主卖点改为跨谱系迁移本身，unknown 边缘化降为分析" |

补充说明：0.1a、0.3 揭示的"端到端 lncRNA 结构注释在所有工具上都接近 0"比跨谱系更根本。因此建议的 Phase 1 主线是：**编码 ⊗ lncRNA 结构化语法解码器 + 方案 B 基准**，跨谱系与 unknown 处理作为两个分析或组件。具体方向待与用户讨论。
