# RESULT 0.1a + 0.3：现有 ab initio 工具在 lncRNA 上的表现（2026-10-20）

> 这是 Phase 0 迄今最重要的结果。
> - 第一部分用的是 GENATATOR 团队**自己的**真值和**自己公开的**各工具预测，我们没有自己跑任何一个工具。
> - 第二部分用我们在 Modal 上跑出的 GENATATOR chr21 smoke，与第一部分互相印证。

## 1. GENATATOR leaderboard 面板：T2T chr20（NC_060944.1）

- **真值**：leaderboard 的 `chr20.gff`（RefSeq），只含 mRNA 和 lnc_RNA 基因，共 434 个 lncRNA 基因、546 个 mRNA 基因。
  - 来源：HF Space `AIRI-Institute/genatator-ab-initio-annotation-leaderboard`，revision `d32d4dc`。
- **预测**：GitHub `alexeyshmelev/genatator-ab-initio-leaderboard-predictions`，commit `f26384d`，共 13 份。
- **打分**：ruler v1.3。SegmentNT 与 NTv3 的输出没有链信息，因此用无链模式打分，这与 leaderboard 对它们的 `use_strand=False` 处理一致。

| 模型 | 有链 | lnc L1 Sn（位点） | lnc L1 Pr | lnc L2 Sn（内含子链） | 被任一预测触及的 lnc 基因 | 被报成 mRNA 的 lnc 基因 | mRNA L1 Sn / Pr | mRNA L2 Sn / Pr |
|---|---|---|---|---|---|---|---|---|
| AUGUSTUS | ✔ | 0（不输出 lnc） | – | 0 | 62 | **52** | 0.905 / 0.546 | 0.139 / 0.080 |
| Helixer | ✔ | 0 | – | 0 | 12 | 4 | 0.883 / 0.890 | 0.245 / 0.240 |
| Tiberius | ✔ | 0 | – | 0 | 9 | 3 | 0.932 / 0.919 | 0.452 / 0.468 |
| ANNEVO | ✔ | 0 | – | 0 | 15 | 8 | 0.963 / 0.899 | 0.325 / 0.304 |
| **GENATATOR** | ✔ | **0.035** | 0.259 | **0.007** | 42 | 22 | 0.932 / 0.947 | **0.648** / 0.362 |
| NTv3 100M / 650M（32 kb） | ✘ | 0.512 / 0.465 | 0.135 / 0.234 | 0.000 / 0.005 | 252 / 228 | – | 0.965 / 0.971（Pr 0.43 / 0.63） | 0.04 / 0.09 |
| NTv3 100M / 650M（1 Mb） | ✘ | 0.514 / 0.447 | 0.144 / 0.200 | 0.000 / 0.005 | 252 / 221 | – | 0.960 / 0.965 | 0.04 / 0.12 |
| SegmentNT（4 个版本） | ✘ | 0.071–0.136 | 0.09–0.18 | 0 | 93–126 | – | 0.94–0.96 | ≈0.05 |

**用他们官方的 metric 引擎在其主口径（k = 250 bp）下核对 GENATATOR**：

| | interval-level P / R | segmentation-level P / R |
|---|---|---|
| lnc_RNA | 0.049 / **0.009** | 0.037 / **0.007** |
| mRNA | 0.425 / 0.707 | 0.272 / 0.496 |

在 k = 0 时，lnc_RNA 的 P 和 R 都是 0。官方引擎解析不了 ANNEVO、Helixer、Tiberius 的 GFF，但这三个工具本来就不输出 lncRNA，不影响结论。

## 2. 我们在 Modal 上跑的 GENATATOR（hg38 chr21，smoke，按协议不作为结论）

- 运行：76 分钟，A10G，fp32，显存峰值 22 GB（靠 pipeline 自动缩小 batch）。
- 对上游做了一处兼容修复：`dtype` 被吞掉的 bug，见 `genatator_cli.py`；生物学参数没有任何改动。
- 输出：408 个 mRNA、46 个 lnc_RNA。5,939 个候选区间中保留 782 个，其中 727 个被判为 mRNA。
- 对 GENCODE v50 打分：
  - lnc L1 Sn 0.029，L2 0.007
  - mRNA L1 Sn 0.819 / Pr 0.934，L2 0.562
- 590 个 lncRNA 基因中只有 33 个被任一预测触及：**瓶颈在区间发现阶段，而不是类型分类**。
- 用官方 T2T 预测作参照，lnc : mRNA 的输出比例（46 : 408 对 81 : 1089）同一量级。说明我们这边的运行没有跑偏，lncRNA 输出稀少是已发布 pipeline 本身的特性。

## 3. 结论

1. **端到端地从基因组中找出 lncRNA 基因结构，现有方法基本没有解决。**
   - 唯一号称做这件事的完整流程（GENATATOR），在它自己的基准、自己的口径下，lncRNA 转录本召回 < 1%，位点召回 3.5%。
   - 所有蛋白编码 SOTA 都不输出 lncRNA，而且会把 3–52 个 lncRNA 位点报成编码基因。
   - NTv3 能标出约一半 lncRNA 所在的区域，但给不出任何正确的内含子链，precision 也很低。
2. 这把最初的立论（"SOTA 不关注非编码基因"）以可量化、可复核的形式重新立住了。论文措辞改为："**现有 ab initio 注释器，包括唯一明确面向 lncRNA 的 GENATATOR，在端到端设定下几乎无法恢复 lncRNA 基因结构**"，而不是"从未有人做过"。
3. GENATATOR 论文中 lncRNA 的正面结果，应当来自"已知转录本区间内的分割"任务，也就是我们旧项目追的那个基准；而不是全基因组的 gene finding。这一点需要在论文中引用原文核实后再表述。
4. **对协议 gate 的影响（如实记录，不改判据）**：
   - 0.1b 的"跨谱系明显退化"判据，以人类 lnc 指标为参照；人类本身已经接近 0，所以这个判据失去区分度。
   - 0.1b 继续运行，作为方案 B 基准的一部分：它能给出各谱系上 mRNA 的迁移表现，以及 lnc 的绝对水平。
   - gate 的最终判读，留到 0.1b 和 0.5 都完成后与用户讨论。
5. **对方案 A 的含义**：主张不再需要依赖"跨谱系退化"。方法只要能在人类等内分布物种上，把 lncRNA 的位点召回或内含子链召回做到明显不为 0（例如 > 20%），而 precision 又可以接受，就已经是领域内的首个可用结果。0.5（unknown 边缘化）现在是 A 的核心 go/no-go。
