# 全稿通读记录（2026-11-15）

> 全稿：`MANUSCRIPT_DRAFT.md`，由 `scripts/manuscript/assemble.py` 从各节草稿拼成。以后只改各节草稿，再重新拼装。
> 本次没有新增任何数字，只核对、对齐和改措辞。

## 1. 已经改了的

**图号重排**：
- 原图 5（unknown 处理）→ 补充图 1。
- 原图 6 → 图 5（位点检出），原图 7 → 图 6（端点）。
- 两个作图脚本的输出名和 `fig_data*.json` 的键名都已同步，重新生成后数值不变（只有键名的 diff）。
- `LEGENDS.md` 的标题和编号已对应修改。

**正文中的图引用逐条对齐了实际画出来的面板**：
- 第 2 节：补上图 2b（官方指标）和图 2c（lncRNA 被报成编码）的引用；图 2c 加上 GENATATOR 的 22。
- 第 5 节：补上图 4。
- 第 6 节：补上补充图 1。
- 第 7 节：SegmentNT / AlphaGenome / Tiberius / Helixer 的同尺比较原来引「图 6b」，实际图里没有这一面板，改为补充表 1。
- 第 8 节：消融原来引「图 6c」，图里没有（之前已决定放补充表），改为补充表 2。
- 第 9–11 节：改为图 5a / 5b / 5c。
- 第 15 节：Puffin 原来引「图 7d」，实际画在 6c 里，已改正；PolyaID 改为 6d。

**数字和措辞的不一致**：
1. 第 2 节「46 lncRNA loci reported only as mRNA」写错了。46 是「只落在 lncRNA 位点上的 mRNA 预测数」；按第 7 节的口径（被预测 mRNA 覆盖、且该预测不命中任何参考 mRNA 的 lncRNA 基因）应为 34（`RESULT_BASELINES.md` 表中的「34 / 46」）。现在两个数都写出来了，与第 7 节一致。
2. 第 6 节：AUPRC 的提升是事后探索（`RESULT_0.5_unknown.md` 中「探索性分析（事后进行，不参与 gate）」），原文与预注册终点写在同一句里，现已注明「in a post hoc analysis」，与补充图 1 的图注一致。
3. 引言原写「lncRNA genes are as numerous as protein-coding genes」。v50 中 lncRNA 34,677、编码 20,033，已改为「more numerous」并给出数字。
4. 方法：自助法次数原只写 1,000，实际 Phase 1 读取用的是 2,000，现已分别注明。
5. 第 3 节：补上 GENATATOR 在人类 chr20 上的 mRNA 召回 0.97 作对照，与图 3 一致。
6. 方法「Decoding and calibration」与「Training」两段重复写了「换掉 chr9 上 154 个窗口」，删掉了后一处。

**状态行**：
- 各节草稿顶部「figures planned, not drawn」已改为图已画。
- 端点一节作者备注第 1 条（「如被要求可在 chr20 上复现」）已更新为「已复现，S11 未复现」。

**核对过、无需改动的**：
- 摘要、引言、讨论里的所有数字都与结果节一致：30.6% / 30.4%、0.42–0.50、6.4% / 6.9%、55%、39.6%、13 倍、8–21 kb、PolyaID 约 1/10。
- 图 1 的染色体范围是 chrX 0.394 到 chr16 0.580；chr9 为 0.396，与 chrX 几乎相同，正文写 chrX 没错。
- 图 3 的 lncRNA 预测数：非哺乳动物 0–13 与正文一致（人类 83 不在此范围内）。

## 2. 需要你决定的

1. **篇幅**：正文约 6,050 词（引言 700、结果 4,260、讨论 1,090），方法约 1,630 词，图注约 640 词。Nature Communications 的正文一般不超过 5,000 词，要压掉约 1,000 词；Genome Biology 没有硬性上限。可删的候选：

   | 候选 | 约省词数 | 代价 |
   |---|---|---|
   | 第 4 节（小 RNA）移入补充材料 | 100 | 基本没有：它本来就不是方法主张 |
   | 第 6 节并入方法或补充图 1 的图注 | 90 | 第 7 节需要一句话交代 unknown 处理 |
   | 第 12 节开头重复第 7 节的数字 | 80 | 无 |
   | 第 15 节（924 词，最长）的 Puffin / PolyaID 细节移入方法 | 200–250 | 结论句保留 |
   | 第 8 节消融只保留结论和一个数字，其余放补充表 2 | 100 | 无 |
   | 第 9 节 v3.0（加小鼠训练）的负结果移入补充材料 | 120 | 作者备注里原本就把它列为待定 |
   | 第 14 节三个排除项写得更简洁 | 100 | 无 |

   合计约 800–850 词，再加上讨论部分压缩 150–200 词，即可降到约 5,000 词。

2. **期刊**：上面这些删法取决于投哪本期刊（标题的选择也取决于此，见引言草稿的作者备注）。

## 3. 还没解决的（不需要新算力）

- **[verify] 引用 5 处**：Tiberius / Helixer / ANNEVO 的具体数字、AlphaGenome、Mele 2017、Schlackow 2017、ANNEVO 的 18.3%。需要逐条读原文核对。
- **参考文献不全**：目前只列了端点一节和引言新增的条目。还缺：GENCODE（v41 / v50）、Snaptron / recount3、FANTOM5、PolyASite 2.0、CHESS 3、LncBook、RNAcentral、Infernal / Rfam、tRNAscan-SE、GENERanno、Nucleotide Transformer、Caduceus、SegmentNT、AlphaGenome、AUGUSTUS、Vipsania、GENATATOR（2025）、Helixer。
- **「lncRNA 被报成编码」有三套数字**：分别来自不同染色体和不同参考（T2T chr20 RefSeq：Tiberius 3 / Helixer 4；hg38 chr20 v50：15 / 14 / GENATATOR 34；第 2 节 Helixer 的「32 个假阳性中 18 个落在 lncRNA」也是 hg38）。正文已分别标明出处，但审稿人可能会觉得乱。建议在第 2 节只保留 T2T 的数字，hg38 的放进补充表 1。
- **补充材料还没成形**：补充表 1（同尺比较）、补充表 2（消融），以及第 14 节的对照、标签支持统计、逐 seed 探针表（端点一节图计划里列为补充材料）都只有出处，还没做成表。
- **第 13 节（S11）没有在 chr20 上复现**，已写在局限 (ii) 里。
- **统一 GFF3 的模型版本不一致**（chr20 用 v1.9），方法里已写明。
