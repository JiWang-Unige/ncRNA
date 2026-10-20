# RESULT 0.6：小 ncRNA 组件（Infernal/Rfam + tRNAscan-SE，2026-10-20）

> 按协议 v1.2，这一部分只是描述性的组件验收，不设 gate，也不构成方法主张：真值本身就是由同类工具产生的，存在循环。

## 设置

- **工具**
  - Infernal 1.1.5 cmscan，参数 `--rfam --cut_ga --nohmmonly --oclan --clanin`，Rfam 15.1
  - tRNAscan-SE 2.0.12，`-E` 真核模式
- **运行**：Modal CPU（`ji-wang-unige`），按 1 Mb 分块并行（重叠 10 kb）。
  - 代码：`scripts/phase0/smallrna_modal.py`、`score_smallrna_panel.py`
  - 输出：Baobab `outputs/PHASE0/0.6_smallrna/`
- **预测过滤**：剔除 clan 内被压掉的重叠命中（`=`）、Rfam lncRNA 家族（218 个），以及 Cis-reg、Intron 类型。
- **打分**：ruler v1.2 的 small 类，L1 指同链外显子重叠即为命中。

## 结果

| 染色体 | 真值 small | 预测数 | L1 Sn [95% CI] | L1 Pr（严格） | 新位点预测 | cmscan CPU-h（每 Mb） |
|---|---|---|---|---|---|---|
| 人 chr20 | 196 | 570 | 0.781 | 0.268 | 397 | 17.6（0.27） |
| 鸡 5 | 58 | 53 | 0.379 [0.26, 0.50] | 0.415 | 31 | 9.1（0.15） |
| 斑马鱼 22 | 76 | 842 | 0.553 | 0.050 | 798 | 7.8（0.20） |
| 果蝇 2L | 283 | 133 | 0.527 | 0.992 | 1 | 3.9（0.16） |
| 线虫 I | 988 | 194 | 0.116 | 0.933 | 4 | 3.7（0.24） |

### 按真值 biotype 的召回（命中 / 总数）

| 染色体 | 召回情况 |
|---|---|
| 人 | snRNA 45/46，snoRNA 24/28，rRNA 5/5，misc_RNA 52/68，**miRNA 26/48** |
| 鸡 | snoRNA 10/10，snRNA 2/2，misc_RNA 4/4，**miRNA 4/40** |
| 斑马鱼 | miRNA 10/10，snRNA 2/2，snoRNA 3/7，rRNA 27/55 |
| 果蝇 | tRNA 41/41，snRNA 15/16，**snoRNA 6/41**，笼统短 ncRNA 87/185 |
| 线虫 | tRNA 66/66，snRNA 16/16，snoRNA 17/74，miRNA 8/33，**piRNA 0/98，笼统短 ncRNA 4/696** |

总花费：5 条染色体约 42 CPU-h，Modal 估计约 3 美元。

## 解读

1. **经典结构化 RNA（tRNA、snRNA、rRNA、脊椎动物 snoRNA）几乎全部召回**。这部分可以直接整合进统一的 GFF3。
2. **召回的缺口集中在三类**：
   - miRNA：鸡只有 4/40；人类、线虫约一半。Rfam 的 GA 阈值对 miRBase 发夹覆盖有限，需要接 miRNA 专用工具或 miRBase 比对。
   - 线虫 piRNA 和 WormBase 笼统的短 ncRNA：没有结构模型可用，属于从序列很难注释的类别。
   - 果蝇 snoRNA（6/41）。
3. **低 precision 主要来自真值缺口和重复序列，而不是工具误判**：
   - 人类 397 个新位点中，357 个是 Metazoa_SRP 命中，即由 7SL 衍生的 Alu。整合规则需要过滤重复序列内的 SRP / U6 假阳性：参照 Rfam 自己的做法，要求全长匹配，并排除与 RepeatMasker 注释重叠的命中。
   - 斑马鱼 798 个新位点主要是 tRNA（377 个）和 5S rRNA 簇，而 Ensembl GTF 根本不收录 tRNA 基因。
4. **对统一 GFF3 的整合规则（待写入软件设计）**：
   - Rfam 与 tRNAscan 的命中直接作为 small 层输出。
   - SRP / U6 等易受重复序列干扰的家族，需要额外过滤。
   - miRNA 另接专用组件。
   - 与 lncRNA 外显子重叠的小 RNA（例如 snoRNA 宿主基因）要保留两者，并在 GFF3 中用 `Parent` / `Note` 标注关系。
