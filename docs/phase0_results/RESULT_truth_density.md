# 面板物种真值概况（2026-10-20，只描述注释，不涉及模型）

真值来源：Ensembl 116（果蝇、线虫为 Metazoa 63 同期文件），人类为 GENCODE v50。类别按协议 v1.2 映射（ruler v1.2）。每个物种都做了真值对真值自检，所有类别的 L1–L6 = 1.0。

| 物种 / 染色体 | 长度 | lnc 基因 | lnc / Mb | mRNA 基因 | small | lnc : mRNA | 原始 biotype 说明 |
|---|---|---|---|---|---|---|---|
| 人 chr20 | 64.4 Mb | 955 | 14.8 | 548 | 196 | 1.74 | GENCODE lncRNA |
| 鸡 5 | 59.5 Mb | 694 | 11.7 | 991 | 58 | 0.70 | 标准 lncRNA |
| 斑马鱼 22 | 39.1 Mb | 135 | 3.5 | 1,007 | 76（含 55 个 rRNA） | 0.13 | 旧命名：lincRNA 35、processed_transcript 89、antisense 11 |
| 果蝇 2L | 23.5 Mb | 619 | 26.3 | 2,657 | 283 | 0.23 | FlyBase 笼统 ncRNA 804 个，按 200 nt 拆分；另有 953 个 transposable_element 归为 other |
| 线虫 I | 15.1 Mb | 150 | 9.9 | 2,878 | 988 | 0.05 | WormBase ncRNA 809 个，按 200 nt 拆分；加 lincRNA 23、antisense 14 |

（人类 small 数取 GENCODE chr20 的 misc_RNA、miRNA、snRNA、snoRNA、rRNA、scaRNA 合计。）

## 含义

- 各谱系 lncRNA 真值的密度和"相对编码基因的比例"差异超过 30 倍（lnc : mRNA 从 0.05 到 1.74）。这既反映生物学差异，也反映注释投入的差异：斑马鱼 GRCz11 的基因集多年未更新。
- 由此带来的后果：
  - 在斑马鱼、线虫上，"precision" 主要衡量的是注释缺口，而不是模型错误。按协议，这些物种以 Sn 为主，严格 / 宽松两种 precision 都报告。
  - 0.1b 的退化判据以同物种的 mRNA 作对照，正是为了把"注释差"和"模型差"分开。
- 这一张表本身就是方案 B 的第二个动机：跨谱系的 lncRNA 真值不只是不完整，而且不完整的程度在谱系之间差别极大。
