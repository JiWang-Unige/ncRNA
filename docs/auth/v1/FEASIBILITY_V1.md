# AUTH v1 第一步：跨细胞系可行性（只查元数据，2026-11-24）

> 来源：`PLAN_AUTH.md` 第 6 节决策表，「G1 与 G2 都通过 → 写 v1：跨细胞系可行性（只查元数据）」。作者 2026-11-24 选「查跨细胞系数据」。
> 范围：只查元数据，**没有下载任何数据文件**。查了 ENCODE 门户（实验、文件、长读长 QC）、GEO（E-utilities）、FANTOM5 的 hg38 / mm10 重处理目录、LRGASP 文档和论文（PMC11543605）。
> 原始查询结果在 `meta/`，汇总在 `V1_CANDIDATES.json`。

## 1. 需要什么

要在另一个细胞系上重复 AUTH，需要同一种细胞同时有：

1. **起点证据**：带链信息的 5′ 帽测序（CAGE 一类），最好有重复；
2. **终点证据**：3′ 端测序（QuantSeq / PAS-seq 一类）；
3. **至少 2 个长读长重复**，有 TALON 定量和 GTF，用同一规则建真值（两个重复都观察到的 lnc）和「不确定」口径；最好还有第 3 个重复，用来跑 StringTie / IsoQuant 作参照。

WTC11 三样都来自 LRGASP 同一批细胞：ssCAGE 6 个样本（GSE185917），QuantSeq 3 个重复，PacBio 3 个重复。

## 2. 结论

- **没有第二个细胞系像 WTC11 那样三样齐全、来自同一批细胞。** LRGASP 只公开了 WTC11 的 ssCAGE；H1、H1 定形内胚层（H1-DE）、小鼠 ES 有同一批的 QuantSeq 和 PacBio，但没有同一批的 CAGE。LRGASP 论文提到小鼠 ES 也用了同样本 CAGE，但 GEO 和 ENCODE 里都找不到，可能只在 Synapse（需要登录）。
- **能拼出来的组合有 4 个**，都需要从别处补起点证据，或者换起点测序类型。最可行的是 **GM12878** 和 **H1**，各有一个明显短板（见第 3 节）。
- 3′ 端和长读长都来自 Mortazavi 实验室的 ENCODE4 流程，格式与 WTC11 相同，可以直接复用代码。CAGE 要换来源（FANTOM5），峰要自己重新定义，这是迁移时最大的改动。

## 3. 候选

长读长深度用 ENCODE QC 里的全长非嵌合读段数（FLNC，百万）。

| 细胞系 | 起点证据 | 终点证据 | 长读长（PacBio） | 同一批细胞？ | 主要问题 |
|---|---|---|---|---|---|
| WTC11（已用） | LRGASP ssCAGE，6 个样本 | QuantSeq 3 个重复（每个 8.7 MB） | 3 个重复：2.62 / 1.71 / 3.09 | 三样都是 | — |
| **GM12878** | FANTOM5 hCAGE 3 个重复（hg38 CTSS，每个 22–25 MB，较深）；另有 ENCODE RAMPAGE 2 个重复、ENCODE2 RIKEN CAGE | QuantSeq **只有 1 个重复**，位点文件 1.7 MB（约为 WTC11 每个重复的 1/5） | Sequel II 2 个重复：2.14 / 2.54；另有 Sequel I 4 个重复：0.69 / 0.79 / 0.68 / 0.75 | 否（CAGE 约 2011 年，长读长 2019–2020 年） | 终点证据浅、只有 1 个重复；和 WTC11 细胞类型不同（淋巴母细胞） |
| **H1 hESC** | FANTOM5 hCAGE 3 个重复（CD34+ 分化时间序列第 0 天，每个 6.1 / 8.1 / 19 MB，较浅）；ENCODE2 RIKEN CAGE（2011 年，整细胞 polyA+ 2 个重复） | LRGASP QuantSeq 3 个重复（每个 6.9 MB） | LRGASP 3 个重复：**0.67 / 0.25 / 1.34**（偏浅、不均） | QuantSeq 和 PacBio 是；CAGE 不是 | 长读长浅：真值只能用 rep1 ∩ rep3，留给参照的 rep2 只有 0.25 M；CAGE 浅 |
| H1-DE | **没有 CAGE**，只有 PRO-cap 2 个重复（ENCODE4，2025 年，新生转录起点） | LRGASP QuantSeq 3 个重复（每个 8.1 MB） | LRGASP 3 个重复：1.59 / 1.22 / 1.55（均衡） | QuantSeq 和 PacBio 是 | 起点测序类型不同：PRO-cap 测新生 RNA，含大量不稳定转录，不能直接当 CAGE 用 |
| 小鼠 ES（F121-9） | 同一株没有公开 CAGE；FANTOM5 有别的 ES 株（ES-OS25 3 个重复、ES-46C 2 个重复） | LRGASP QuantSeq 3 个重复（mm10） | LRGASP 3 个重复：4.06 / 2.53 / 2.14（最深） | QuantSeq 和 PacBio 是 | 跨物种；CAGE 来自别的 ES 株；要另下 mm10 基因组、GENCODE 小鼠注释和重复序列注释 |
| K562 / HepG2 | FANTOM5 3 个重复、ENCODE RAMPAGE / CAGE | ENCODE 里**没有** 3′ 端测序；GEO 只有别的实验室的零散数据（如 HepG2 3′READS） | K562 Sequel II 2 个重复：2.86 / 2.98 | 否 | 缺 3′ 端证据，要从原始读段自己处理；不建议 |

其他已排除的：
- ENCODE 的 PAS-seq 在人和小鼠里只有上面 5 个细胞系，另有人角质形成细胞（没有长读长）和小鼠海马、大脑皮层（PAS-seq 各 1 个重复，有长读长，但没有 CAGE，又是跨物种组织）；
- 有长读长和 PRO-cap 的组织（心、肝、肺等）都没有 3′ 端测序；
- GEO 里的 ssCAGE、nanoCAGE、SLIC-CAGE 系列（人和小鼠共 30 多个）没有一个与上述长读长细胞系对应。

## 4. 迁移时要先解决的方法问题

这些不影响「能不能做」，但会影响 v1 计划怎么写。

1. **CAGE 峰的定义**：WTC11 直接用 LRGASP 给的 count.bed 峰。FANTOM5 只给单碱基 CTSS 和一套全 FANTOM5 的 DPI 峰（`hg38_fair+new_CAGE_peaks_phase1and2.bed.gz`，5.2 MB）。「峰宽」的系数在模型里排第二（−1.358，仅次于跨度内 QuantSeq 位点数 +1.366），「主峰占比」也在模型里（−0.429）；峰的定义一变，这两个特征的分布就变。
2. **计数的量纲**：「ssCAGE 合计读段」「QuantSeq 信号」是绝对数，「支持样本数」「QuantSeq 支持重复数」依赖样本数（6 个对 3 个，3 个对 1 个）。直接把 WTC11 模型搬过去，需要先换成 CPM、比例或分位数。
3. **两种问法要分开**：
   - 零样本迁移：WTC11 上训练好的模型原样用到新细胞系，问的是「学到的东西能不能搬」；
   - 重新学习：在新细胞系的训练染色体上重新训练，问的是「这种真伪信号是不是普遍可学」。
   两者都只能在测试染色体（chr2、6、10、16、22）上评估：训练染色体上的 lnc 位点在两个细胞系里大量相同，用了会泄漏。
4. **真值深度**：「错」的定义是与真值重复里任何转录本都不重叠。真值重复越浅，被算成「错」的真实低表达转录本越多，精度绝对值越低。H1 的真值重复（0.67 M 和 1.34 M）明显浅于 WTC11（1.71 M 和 3.09 M）；GM12878 的 Sequel II 两个重复（2.14 M 和 2.54 M）与 WTC11 相当。
5. **长读长参照**：GM12878 可以把 Sequel I 的 4 个重复合起来（约 2.9 M）当参照，深度和 WTC11 rep1（2.62 M）相当，但平台不同；H1 的 rep2 只有 0.25 M，拿它做参照对长读长不公平。

## 5. 建议

先不写完整的 v1 计划，再做一步只下载小文件的描述性检查（不训练、不打分、不设门槛），用数字决定选哪个细胞系：

- **GM12878 和 H1 各下小文件**：TALON 定量与 GTF、QuantSeq 位点、FANTOM5 CTSS 和 DPI 峰集；不下 BAM。
- 在 ont 上算：
  1. 测试染色体上的真值位点数，按 W24 规则；
  2. 证据覆盖：真值起点 100 bp 内有 CAGE 峰、终点 100 bp 内有 QuantSeq 位点的比例，与 WTC11 对照；
  3. 各特征（CPM 化后）的分布和 WTC11 的差别。
- 这一步不碰测试集上的预测，所以不影响之后冻结的迁移测试。
- 选定细胞系后再写 v1 计划（零样本迁移为主，重新学习为辅），门槛和读取方式照 AUTH 的做法冻结。

下载量（都放 ont，不花钱；ENCODE 和 FANTOM5 都是公开数据）：

| 细胞系 | 文件 | 大小 |
|---|---|---|
| GM12878 | Sequel II 2 个重复的 TALON GTF + 定量（ENCFF871ODI、ENCFF993PKH、ENCFF079JHG、ENCFF386OOG）；QuantSeq 位点 2 个版本（ENCFF259EZO、ENCFF448CNL）；FANTOM5 CTSS 3 个；RAMPAGE 起点 3 个 | 约 0.25 GB |
| H1 | 3 个重复的 TALON GTF + 定量（ENCFF721HVH、ENCFF814SOM、ENCFF236CGS、ENCFF580WQM、ENCFF563JZI、ENCFF014CCY）；QuantSeq 位点 3 个（ENCFF970YYF、ENCFF138ZWV、ENCFF215NMJ）；FANTOM5 CTSS 3 个 | 约 0.25 GB |
| 共用 | FANTOM5 hg38 DPI 峰集 | 5.2 MB |

合计约 0.5 GB，下载加计算预计 1–2 小时，只用 ont CPU。之后如果进入正式测试，还要下参照用的 BAM（GM12878 Sequel I 4 个，0.47 GB；或 H1 全部 3 个，0.66 GB）。

**不建议**：K562 / HepG2（缺 3′ 端证据）；小鼠 ES 暂缓（跨物种、要另下基因组和注释；如果作者想做跨物种，可以单独再议）。H1-DE 可以作为「换成 PRO-cap 起点」的压力测试，但不适合当主测试。

## 6. 需要作者决定

1. 是否批准下载第 5 节的小文件（约 0.5 GB，ENCODE 和 FANTOM5 公开数据，放 ont），做描述性检查。
2. 候选范围：GM12878 和 H1 都查（建议），还是只查一个。
3. 是否把 H1-DE（PRO-cap）留作之后的压力测试。
