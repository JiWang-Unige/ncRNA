# RESULT_UGFF — 选择性统一 GFF3（2026-11-14）

> 计划：`PLAN_UGFF.md`（合并与打分前写定）。脚本：`scripts/phase1/ugff/{tiberius_queue.sh, ugff.py, ugff_queue.sh}`。
> 数据：`docs/phase1_results/UGFF/RULES.json`（chr18 开发）、`UGFF/<chrom>/SUMMARY.json`。
> 产物：ont `outputs/PHASE1/UGFF/<chrom>/unified.gff3`。sha256 前 12 位：chr18 6a6d9665fe70、chr20 9d8815def8ae、chr9 6a181f1ae4f1。
> 没有重新读取任何模型。chr20 用的是 v1.9 已读出的预测，chr9 用的是 v2.0 已读出的预测。

## 1. 结论

1. **三条染色体的统一 GFF3 都已生成，格式检查 0 个问题**：ID 唯一，每个 Parent 都存在，子特征都落在父特征范围内，按起点排序。重新读入后的条数与写出的一致。

   | 染色体 | 编码（Tiberius）转录本 | lncRNA 基因 | 小 RNA | 行数 |
   |---|---|---|---|---|
   | chr18（开发） | 313 | 442 | —（未跑） | 7,408 |
   | chr20 | 565 | 422 | 570 | 13,098 |
   | chr9 | 896 | 811 | —（未跑） | 19,027 |

2. **内含子阈值为 none：所有 lncRNA 都只以位点级输出（`structure=unresolved`）。**
   - **强制写出的内含子几乎都是错的。** chr18 上，按解码结果写成转录本时，与 v50 lncRNA 内含子精确一致的比例，三个 run 分别只有 0.26%、0.48% 和 1.5%。
   - **提高阈值也救不回来。** 阈值升到 0.95，三 run 平均精确率也只有 4.1%，离 0.5 的门槛很远。
   - **完整一致的内含子链几乎为零。** 三个 run 合计只有 1 条（seed 456）。
   - **chr20 和 chr9 一样。** 强制输出时内含子精确率为 0.97% 和 0.34%，完整一致的内含子链为 0。
   - **这正是正文的结论**：lncRNA 能检出，但结构界定不出。软件的输出方式也照此处理，不伪造转录本模型，解码出的外显子只放在 `decoded_exons=` 属性里。

3. **编码护栏用敏感度换精度。** 下表是三个 run 的范围。

   | 染色体 | 剔除的 lnc 预测 | 其中只命中 mRNA | 其中命中 lnc | lnc L1 Sn 变化 | 严格 Pr 变化 |
   |---|---|---|---|---|---|
   | chr18 | 31–58 | 10–35 | 14–22 | −0.012 至 −0.027 | 0 至 +0.027 |
   | chr20 | 79–130 | 54–98 | 16–25 | −0.011 至 −0.030 | +0.028 至 +0.043 |
   | chr9 | 84–161 | 37–86 | 37–61 | −0.019 至 −0.046 | +0.009 至 +0.028 |

   - 被剔除的大多是 lncRNA 外显子落在 Tiberius CDS 上的情况，以只命中 mRNA 的为主。
   - 但 chr9 上被剔除的预测里有 37–61 个命中了真实 lncRNA，可能是与编码基因重叠的 lncRNA，或者是 Tiberius 的假阳性。
   - 护栏是预先固定的，没有调；正文里的 Sn 和 Pr 仍按未加护栏的预测报告。

4. **起点的 CpG 分类在开发集上有区分度（chr18，描述性）。**
   - **切点**：取 chr18 上 v50 lncRNA 基因起点 ±500 bp 的 CpG o/e，上四分位数是 0.459。
   - **两类的起点准确度差很多**：在命中 v50 lncRNA 的预测中，CpG 富集类的预测起点有 73%（157/215）落在所命中基因某个转录本的 5′ 端 1 kb 以内；其他类只有 29%（117/408）。
   - **与正文一致**：只有 CpG 岛启动子的起点能从序列定位。所以 `start_class=CpG_rich` 可以当作「起点相对可信」的标记。
   - **解读要注意**：这里的标准比较宽松（命中基因的任一转录本 5′ 端都算），所以 29% 不能和正文 Loc1k 的随机水平 0.05 直接比较。
   - **各染色体的 CpG 富集比例**：chr18 30%、chr20 37%、chr9 35%。

5. **编码层的质量就是 Tiberius 本身的水平。**

   | 染色体 | mRNA L1 Sn | mRNA L1 严格 Pr | 内含子链 Sn | 内含子链 Pr |
   |---|---|---|---|---|
   | chr18 | 0.963 | 0.821 | 0.500 | 0.469 |
   | chr20 | 0.964 | 0.917 | 0.557 | 0.585 |
   | chr9 | 0.959 | 0.818 | 0.529 | 0.512 |

   chr20 的数与 RESULT_BASELINES 一致。chr18 和 chr9 的 Tiberius 是这次新跑的，用的是发布时的设置，什么都没调。

## 2. 输出格式（每个 lncRNA 基因一行 `lncRNA_gene`）

`ID, biotype=lncRNA, confidence, support=k/3, structure=unresolved, n_decoded_exons, min_intron_conf,
promoter_cpg_oe, start_class=CpG_rich|other, end_precision=low, model=ncrna_annot_v1.9|v2.0, decoded_exons=a-b,...`

- **编码层**：`gene` → `mRNA` → `exon` 和 `CDS`，来源列为 Tiberius；只有 CDS，没有 UTR。
- **小 RNA 层（只有 chr20）**：`ncRNA_gene` → `ncRNA` → `exon`。与 lncRNA 同链重叠的写 `Note=overlaps_lnc:<ID>`。

## 3. 局限

- **小 RNA 层只有 chr20 有。** chr18 和 chr9 没跑过 Rfam/tRNAscan（CPU，约数小时），需要时再补。
- **版本不统一。** chr20 是 v1.9，chr18 和 chr9 是 v2.0；要统一，得在 chr20 上重新读一次 v2.0，按规则这次没有做。
- **只发 seed 42 的预测**，另外两个 seed 只用来算 `support`。
- **lncRNA 没有转录本模型**，下游工具若要求每个 gene 都有 transcript 子特征，需要另行处理。这是有意的选择。

## 4. 算力

ont：Tiberius chr18 约 9 分钟，chr9 约 14 分钟（GPU，经 gpu_gate）；合并与打分每条染色体 6 秒到 1 分钟（CPU，nice 19）。
第一次排队时路径里的 `~` 在引号中没有展开，三个 CPU 步骤返回 127，什么都没算；改用 `$HOME` 后重跑成功。
