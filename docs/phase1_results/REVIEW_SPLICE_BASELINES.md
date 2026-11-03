# 剪接信号进入解码 / 基线可行性 / 转录证据：桌面调研（2026-11-03）

> 用户 2026-11-03 决定卡 `p15-next` 第 ① 条（调研 Tiberius、Helixer、SegmentNT 怎么把剪接信号接进解码，写成对比）和第 ④ 条（SegmentNT / Helixer 加进 chr20 对照的可行性）、第 ③ 条（公开转录证据数据源）。桌面调研，只读网页与论文，未下载或安装任何东西；每条事实附来源，查不到的标「未核实」，运行时间均为估算。

范围：Part A 剪接信号如何进入基因结构解码；Part B 在 hg38 chr20 上加入 SegmentNT / Helixer（及 Tiberius）作为基线的可行性；Part C 验证 chr20 上"GENCODE v50 之外"的新 lncRNA 候选所需的公开转录证据。
约定：每条事实后附来源 URL；拿不到一手确认的标「未核实」；我自己的推算标「估算」。所有网页只读，未下载大文件、未安装任何东西（只对文件做了 HTTP HEAD 取大小，并对 Snaptron 做了 3 个 1 Mb 的小查询估体量）。
部分页面（PMC、Nature）被反爬挡住，改用 Europe PMC 的 fullTextXML 接口读同一篇论文正文（如 https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12615259/fullTextXML），引用时仍给出 PMC 文章 URL。

---

## 0. 速览（最重要的 10 条）

1. **SegmentNT（human 版）chr20 = 官方测试集，不泄漏**：训练 = 除 chr20/21 外所有染色体，chr22 = 验证，chr20+chr21 = 测试，标签为 GENCODE V44（去掉 level 3）(https://www.biorxiv.org/content/10.1101/2024.03.14.584712v2.full ；https://pmc.ncbi.nlm.nih.gov/articles/PMC12615259/)。**但我们的 dev 集 chr18 在它的训练集里**，所以 SegmentNT 只能在 chr20 上比。multi-species 版的人类染色体划分论文未写（未核实）。
2. **Helixer 公开 vertebrate 模型在 human 上泄漏**：Tiberius 官方模型配置写明 "Delphinapterus leucas and Homo sapiens were used during training or validation of Helixer" (https://github.com/Gaius-Augustus/Tiberius/blob/main/model_cfg/mammalia_nosofttmasking_v2.yaml)；Tiberius 论文也写测试的三个物种都在 Helixer 训练的 315 个脊椎动物基因组里 (https://pmc.ncbi.nlm.nih.gov/articles/PMC11645249/)。Helixer 的 Zenodo v2 新增了 `mammal_v0.3_m_0100.h5`（36.3 MB），论文称 mammal 模型与 Tiberius 使用同一批物种、human 为测试物种 → 很可能不泄漏，但物种表（Suppl. Table 21）我看不到，**标未核实** (https://zenodo.org/records/17850139 ；https://pmc.ncbi.nlm.nih.gov/articles/PMC13076211/)。
3. **Tiberius：human 完整留出**（论文 test = human/cow/beluga；当前 v2 配置的 29 个训练物种里没有 human）→ 无泄漏 (https://pmc.ncbi.nlm.nih.gov/articles/PMC11645249/ ；上面的 yaml)。
4. **三者都不是 lncRNA 注释器**：Helixer / Tiberius 只预测编码基因（最长 CDS 转录本），不预测 lncRNA。SegmentNT 有 "lncRNAs" 类，但 3 kb 模型 lncRNA 的测试 MCC < 0.1，30 kb 模型的 lncRNA 数值正文未给（只在 Fig. 2a）(https://pmc.ncbi.nlm.nih.gov/articles/PMC12615259/)。
5. **许可**：SegmentNT 代码 + 权重均 CC BY-NC-SA 4.0 (https://github.com/instadeepai/nucleotide-transformer/blob/main/LICENSE.md ；https://huggingface.co/InstaDeepAI/segment_nt)；Helixer 代码 GPL-3.0、权重 CC BY 4.0 (https://github.com/usadellab/Helixer ；https://zenodo.org/records/17850139)；Tiberius 代码 MIT（仓库 LICENSE）但网页服务写 Apache-2.0，权重许可未写（未核实）(https://github.com/Gaius-Augustus/Tiberius/blob/main/LICENSE ；https://bioinf.uni-greifswald.de/tiberius/)。
6. **算力都能在一张 RTX 4090 上完成**：下载量 SegmentNT 2.24 GB、Helixer ~36 MB、Tiberius ~30 MB；chr20（64,444,167 bp，https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes）估算：SegmentNT 约 5–45 min（取决于是否重叠窗口/反向互补平均），Helixer 约 10–15 min，Tiberius 约 2–10 min（均为估算，见 B 节）。
7. **SpliceAI 自身对我们有泄漏**：SpliceAI 训练染色体 = 2,4,6,8,10–22,X,Y，测试 = 1,3,5,7,9（OpenSpliceAI 论文转述）(https://elifesciences.org/articles/107454)。**chr18（dev）和 chr20（test）都在 SpliceAI 训练集里**（针对 mRNA 位点）。SpliceAI 真正"没见过"的 chr1,3,5,7,9 恰好在我们的训练集里，可用来校准权重。
8. **许可陷阱**：SpliceAI 源码是 PolyForm Strict 1.0.0（只允许非商业使用，不允许修改、做衍生作品或分发）；模型 CC BY-NC 4.0 (https://github.com/Illumina/SpliceAI ；https://raw.githubusercontent.com/polyformproject/polyform-licenses/1.0.0/PolyForm-Strict-1.0.0.md)。
9. **最接近"SpliceAI 类分数 + 动态规划解码"的现成做法是 minisplice**（Heng Li 2025）：把 GT/AG 的校准 log-odds 直接加进 minimap2/miniprot 的 DP (https://arxiv.org/html/2506.12986v2)。我没有找到发表的"SpliceAI + HMM/CRF 解码 lncRNA 结构"的方法（多轮搜索，未发现）。
10. **最小转录证据下载集**（约 120–180 MB + Snaptron API 流式查询）：FANTOM5 hg38 CAGE peaks（5.5 MB）、PolyASite 2.0 GRCh38.96（50.6 MB）、Snaptron gtexv2 的 chr20 剪接位点（API）、可选 GTEx long-read FLAIR GTF（14 MB）、CHESS 3.1.3 GRCh38 GTF（20 MB）、LncBook v2.1 GTF（28 MB）。详见 C 节。

---

## A. 剪接信号如何进入解码

### A1. Tiberius（Gabriel et al. 2024, Bioinformatics）
来源：论文 https://pmc.ncbi.nlm.nih.gov/articles/PMC11645249/ ；代码 https://github.com/Gaius-Augustus/Tiberius
- **信号从哪来**：没有专门的剪接 head。CNN + 双向 LSTM 输出 15 个逐碱基类别，类别与 HMM 的 15 个状态一一对应：intergenic 1、intron 3、exon（按阅读框）3、acceptor-site（ASS，每个阅读框一个）3、donor-site（DSS）3、start/stop 2。论文正文对"11 个 exon 状态"与图注不一致，以图 2 为准。
- **结合方式**：可微 HMM 层，**端到端**（训练时梯度穿过 HMM 回到网络）。HMM 只有 24 个未训练的参数，转移参数按 intergenic/intron/exon 的经验平均长度设定，emission 只有一个平滑参数。"HMM 设计为强制规范剪接位点模式、起止密码子"，但 GT/AG 是硬约束还是加权、细节在 Supplement，正文未说（未核实）；正文还说"允许的剪接模式及其相对权重原则上可由用户指定甚至学习"。
- **学什么/固定什么**：网络学；HMM 转移与允许的模式固定。
- **效果**：端到端 HMM 微调使 F1 提高 gene 级 +2.6、exon 级 +1.1 个百分点（相对不带 HMM 训练）。人类（当前 v2 mammalia 无软屏蔽模型，gffcompare）exon F1 87.5、gene F1 65.9；同表 Helixer 为 74.0 / 26.1 (https://github.com/Gaius-Augustus/Tiberius/blob/main/model_cfg/mammalia_nosofttmasking_v2.yaml)。论文正文对 exon/gene 匹配判据在 Supplement（正文未定义"边界是否须完全一致"，未核实）。
- **lncRNA**：不处理。训练标签只用每个基因"最长 CDS 的转录本"，且只能训练一基因一转录本 (PMC11645249；README)；2026 预印本说评估排除 UTR 因为 Tiberius 不预测 UTR (https://www.biorxiv.org/content/10.64898/2026.04.24.720536v1.full)。
- 推理：约 500 kb 的 tile 独立处理，tile 边界内的基因再用约 1 Mb 做二次 pass (PMC11645249)。

### A2. Helixer + HelixerPost（Holst et al.，Nature Methods 2025/26）
来源：论文 https://pmc.ncbi.nlm.nih.gov/articles/PMC13076211/ ；HelixerPost https://github.com/usadellab/HelixerPost ；Helixer https://github.com/usadellab/Helixer
- **信号从哪来**：CNN+BiLSTM 逐碱基输出 4 个 genic 类（intergenic、UTR、CDS、intron）+ 4 个相位类（none、0、1、2）。**没有剪接位点类**。
- **结合方式**：**后处理（post-hoc）HMM**，网络训练时不知道 HMM。HelixerPost 先用滑窗（示例 100 bp）的 genic 均值定候选区（edge 阈值 0.1 开/关区，peak 阈值 0.8 才保留），再在每个候选区内找"罚分最小"的基因模型：罚分 = 与 Helixer 预测的不一致 + 违反生物学知识的状态转移。论文 Fig. 3 图注：intron 状态按"内含子起始/延续、**splice motif**、外层状态"拆成 **60 个子状态**，所以剪接 motif 是在 HMM 里通过读基因组序列加入的，不是网络输出。
- **学什么/固定什么**：网络学；HMM 罚分按"成熟的生物学知识"手工编码（正文未明说是否训练，未核实；正文把它和"带训练的状态/转移似然"的传统 HMM 对比）。
- **效果**：HelixerPost 相对网络直出（HelixerBW）只有小幅提升，如脊椎动物 base-wise phase F1 0.8641 → 0.8829（43/45 测试物种略升）；脊椎动物测试组均值 exon F1 0.7405、intron F1 0.6912、**intron-chain F1 仅 0.1740**、transcript F1 0.1977。human（经 Tiberius 的 yaml 基准）exon F1 74.0、gene F1 26.1。
- **lncRNA**：不涉及（只给"主要基因模型=最长蛋白异构体"，非编码 RNA 仅作为先前工作被引用）。
- 注意：human 在 Helixer vertebrate 模型的训练/验证集中（见 B 节）。

### A3. SegmentNT（InstaDeep；de Almeida et al., Nature Methods 2025）
来源：https://pmc.ncbi.nlm.nih.gov/articles/PMC12615259/ ；预印本 v2 https://www.biorxiv.org/content/10.1101/2024.03.14.584712v2.full
- **信号从哪来**：专用输出通道。14 类逐核苷酸概率（含 splice donor、splice acceptor、exon、intron、lncRNA 等），剪接位点标签用 HISAT2 的 extract_splice_sites.py 从 GENCODE V44 抽取（level-3 转录本剔除）。
- **结合方式**：**没有解码器**。作者明确写 "We do not add further constraints on predictions such as the fact that one nucleotide belongs only to one element"；评估只用 0.5 阈值的核苷酸级 F1/MCC + 区域级 SOV；与 AUGUSTUS 比较时也是这些指标，**没有 exon-level / intron-chain 精确匹配**。
- **学什么/固定什么**：全部端到端学习，无外部结构先验。
- **效果（剪接位点）**：在 SpliceAI 的 mRNA 测试集上 auPRC：acceptor 0.93 vs SpliceAI 0.96 / Pangolin 0.94，donor 0.93 vs 0.94 / 0.94；全染色体测试集 MCC 高于二者。预印本 v2 给出的 MCC 为 acceptor 0.75 vs 0.48、donor 0.76 vs 0.42；某聚合页对后续版本转述为 SpliceAI 0.67/0.59 — 版本间数字不一致，**以最终版为准（未核实）**。
- **strand**：标签对正负链基因合并（"including protein-coding genes and lncRNAs in both positive and negative orientation"），**输出不带链信息**。
- **lncRNA**：有 "lncRNAs" 类，但是 locus 类；exon / intron 类不区分 mRNA 与 lncRNA。3 kb 模型 lncRNA 与 CTCF 位点是最难的，测试 MCC < 0.1；10/30 kb 增益"主要来自蛋白编码基因、lncRNA 和内含子"，但 30 kb 的 lncRNA 数值只在 Fig. 2a（我没取到）。multi-species 版只有 7 个基因元素，不含 lncRNA。
- **已知偏差**：同一核苷酸的概率依赖它在输入窗口里的位置（首/中/尾），且有 24 nt 周期振荡（一篇独立预印本，用 30 kb 模型）；窗口中心最可靠，边缘方差大 (https://www.biorxiv.org/content/10.1101/2025.04.09.647946.full)。

### A4. AUGUSTUS / BRAKER / GeneMark（外部证据以"奖励/强制"进入）
- **AUGUSTUS hints**：外部信息以 GFF hint 提供，可局部强制或建议结构；extrinsic.cfg 里的 bonus 与 malus "是乘到基因结构后验概率上的因子"，bonus=1 / malus=1 即关闭。hint 类型包括 ass / dss（剪接受体/供体）、intron / intronpart、exon / exonpart、CDS、UTR、tss / tts 等；精确 intron hint 只有完全吻合才得 bonus，intronpart 对区间内每个内含子碱基给 bonus，exonpart 的 malus 按碱基累乘因此应接近 1；dss/intron 等 hint 还会让 AUGUSTUS 预测 GC 供体 (https://sources.debian.org/src/augustus/3.5.0%2Bdfsg-2/docs/tutorial/data/extrinsic.bug.cfg)。权重是配置文件里按来源设定的常数，**不学习**。本质上就是"在 log 后验上加 log(bonus) 的外部势"。
- **GeneMark-ETP / BRAKER3**：RNA-seq 与蛋白证据先产生高置信基因，再把筛选后的 hints "plugged in … into the Viterbi algorithm"（强制还是软奖励正文未区分）；第 3 类 RNA-seq 内含子 hint **只有在 ab initio 预测也支持时才采用**，文中说这样可过滤来自"表达的非编码 RNA"的假阳性内含子 hint (https://pmc.ncbi.nlm.nih.gov/articles/PMC11216313)。→ 对编码基因注释器而言，lncRNA 的剪接证据是"噪声"；这对我们是反向的好信号：证据源不会天然偏向编码。BRAKER3 = GeneMark-ETP + AUGUSTUS + TSEBRA 合并 (https://pmc.ncbi.nlm.nih.gov/articles/PMC10312602)。
- lncRNA：这些工具都不建模（GeneMark-ETP 仅把 ncRNA 内含子当过滤对象）。

### A5. 2023–2026 其他相关方法
- **GeneDecoder（Marin, Pultz, Boomsma, arXiv 2025）**：潜变量线性链 CRF，观测标签 E/I/D/A/NC（外显子、内含子、供体、受体、非编码），隐状态加链与密码子位置共 25 个，转移限制为生物学合法、正负链共享转移权重；Viterbi 精确解码 + forward-backward 后验；特征网络（空洞 CNN + BiLSTM）与 CRF 端到端训练；未使用 GT-AG 显式规则，也不处理 lncRNA。消融（拟南芥子集，MCC/F1）：LSTM 单独 0.79/0.75 → LSTM+CRF 0.86/0.84；文中称没有 CRF 时剪接位点常被标错 (https://arxiv.org/html/2505.03377)。
- **minisplice（Heng Li, arXiv 2025）**：7,026 参数的 1D CNN 对基因组中每个 GT / AG 估计剪接概率；原始得分分 50 个 bin，每个 bin 的经验概率相对 null 换算成 2·log2 的 log-odds，作为 d(i)、a(i) 项进入 minimap2/miniprot 的 DP 递推（来自论文公式 2，"加/减"的方向是我读公式的推断）；只用 GT/AG；奇数染色体 80% 基因训练、偶数染色体全部 GT/AG 测试；效果：斑马鱼蛋白比对到人类时，未注释剪接位点占比 14.01% → 4.37%，直接 RNA-seq 的 minimap2 junction 错误率 1.4% → 1.0% (https://arxiv.org/html/2506.12986v2 ；代码 https://github.com/lh3/minisplice)。
- **Splam（Chao et al., Genome Biology 2024）**：剪接位点深度模型，用于过滤比对/转录本中的坏 junction，不是结构解码器 (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03379-4)。
- **OpenSpliceAI（eLife 2025）**：PyTorch 版 SpliceAI，RefSeq MANE v1.3 / GRCh38 重训；human 版沿用 SpliceAI 的染色体划分；默认 biotype = protein-coding，未评估 lncRNA (https://elifesciences.org/articles/107454)；GPL-3.0 代码，权重许可页面未说明（未核实）(https://github.com/Kuanhao-Chao/OpenSpliceAI)。
- **GENATATORs（bioRxiv 2026-06）**：DNA 语言模型微调（Caduceus PS、GENA-LM large，250 kb 上下文，多物种）。逐核苷酸 5 类（exon、intron、CDS、5'UTR、3'UTR），**不使用 HMM/CRF/Viterbi**，也未提 SpliceAI；有一个"splice site filter"（排除两侧二核苷酸在真实数据里从不出现的外显子）、反向互补平均（Caduceus 约 1.5×）、CDS 启发式。声称各 DNA LM 都能检出 lncRNA 基因（Tiberius 不能）。**人类划分：训练 = 除 8、20、21 外的染色体，验证 = chr8、chr21，测试 = chr20（T2T 组装）**。与本项目的设定（chr20 测试、含 lncRNA）高度重合，建议在论文里并列讨论 (https://www.biorxiv.org/content/10.64898/2026.06.17.732686v1.full)；其 lncRNA 外显子/内含子/locus 的具体数字未取到（表格未随文本给出）。
- **ANNEVO（Nature Methods 2026）**：MoE 基因组语言模型，纯 ab initio；"Gene structure decoding" 在 Supplementary Method 3，我没取到剪接/Viterbi 细节（未核实）(https://www.nature.com/nature-index/article/10.1038/s41592-026-03036-7)；Tiberius 配置称它训练时没见过 human（同上 yaml）。
- **NTv3 annotate（InstaDeep 2025–26）**：有独立解码器与过滤（CDS<201 bp、含内部终止密码子的基因被丢弃），只输出 gene/CDS/intron，不含 lncRNA、UTR、exon 行；mammals 模型在 28 个哺乳动物上训练；1 Mb 窗口峰值显存 25.3 GB，文档明说"24 GB 的卡跑不了"；HF 权重需申请（gated）、许可页未写 (https://github.com/instadeepai/nucleotide-transformer/blob/main/docs/ntv3_genome_annotation.md)。→ 4090 上不可行，不列入基线。
- **没找到的**：2023–2026 年把 SpliceAI 类分数接入 HMM/CRF 去解码 gene/lncRNA 结构的已发表工作。最近的是 minisplice（DP 中的 log-odds）与 GeneDecoder（CRF 但自带特征网络）。Illumina 页面提到 "SpliceAI 2"（预测位点使用频率、junction 连接与全长异构体），属商业许可的线索，未核实 (https://help.connected.illumina.com/annotation/data-sources/splice-ai-2.md)。

### A6. lncRNA 剪接信号的特殊性（影响 SpliceAI 能否直接复用）
- lincRNA 剪接效率低于 mRNA，U2AF65 结合较弱、剪接相关 motif 较弱，polypyrimidine tract 的嘧啶更少（Melé et al. 2017, Genome Research）(https://genome.cshlp.org/content/27/1/27.abstract)。→ 在 mRNA 上训练的 SpliceAI 对 lnc 位点的分数分布可能整体偏低且校准不同（推断，未核实；建议先在训练染色体上直接量化）。
- SpliceAI 随包提供的是 GENCODE V24 canonical 注释 (https://github.com/Illumina/SpliceAI)；其训练转录本类型（principal transcripts of protein-coding genes）出自 Cell 论文，本次未抓到原文（未核实）。

### A7. SpliceAI 对我们 chr18 实验的直接含义
- 训练染色体 2,4,6,8,10–22,X,Y；测试 1,3,5,7,9（"the training set is defined by chromosomes 2, 4, 6, …, 22, X, and Y … testing set … 1, 3, 5, 7, and 9"，OpenSpliceAI 论文引述）(https://elifesciences.org/articles/107454)。OSAIMANE 沿用同一划分。
- 后果：(i) chr18/chr20 上 mRNA 位点的 SpliceAI 分数偏乐观（模型见过这些位点的标签）；(ii) 在我们的训练染色体里，只有 chr1,3,5,7,9 的分数是"诚实的"；(iii) 如果让 CRF 在全部训练染色体上端到端学权重，会学到对 SpliceAI 过度信任的系数，换到 chr18/chr20 就失配。
- 许可：源码 PolyForm Strict 1.0.0 ——只许非商业使用，**不许分发或"修改/做衍生作品"**；模型 CC BY-NC 4.0；预计算分数仅限学术/非营利 (https://github.com/Illumina/SpliceAI ；PolyForm 条款摘要 https://raw.githubusercontent.com/polyformproject/polyform-licenses/1.0.0/PolyForm-Strict-1.0.0.md)。仓库已于 2026-04-20 归档（只读）。

### A8. 对照表

| 方法 | 剪接信号来源 | 与结构解码器的结合 | 学习 vs 固定 | 精确度（来源见上） | lncRNA |
|---|---|---|---|---|---|
| Tiberius | 网络逐碱基 15 类（含 ASS/DSS 状态），无独立 head | 可微 HMM 层，端到端 | 网络学；HMM 转移/模式固定（24 参数） | human exon F1 87.5 / gene F1 65.9；端到端 HMM 比不带 HMM 训练 +2.6 gene / +1.1 exon F1 | 否（只编码，最长 CDS 转录本） |
| Helixer + HelixerPost | 网络 4 类+相位；剪接 motif 由 HMM 读序列加入（intron 拆 60 子状态） | 后处理 HMM（Viterbi 罚分最小化） | 网络学；HMM 罚分手工 | 脊椎动物 exon F1 0.74、intron F1 0.69、intron-chain F1 0.17；human gene F1 26.1 | 否 |
| SegmentNT | 专用 donor/acceptor/exon/intron 通道 | 无解码器，0.5 阈值 | 全端到端 | 只有核苷酸级 F1/MCC 与 SOV；无 exon/intron-chain 精确匹配 | 有 lncRNA 类，3 kb MCC<0.1；不分 lnc/mRNA 的 exon/intron |
| AUGUSTUS hints | 外部 ass/dss/intron/exonpart 提示 | 乘到结构后验的 bonus/malus | 全固定（配置） | —（依赖证据质量） | 否 |
| GeneMark-ETP | RNA-seq 内含子 + 蛋白 hint | 强制进 Viterbi；ab initio 须佐证 | 基因模型自训练；hint 规则固定 | —（BRAKER3 human gene F1 72.2） | 把 ncRNA 内含子当噪声过滤 |
| GeneDecoder | 网络直接输出 D/A 标签 | 线性链 CRF（Viterbi + FB），端到端 | 全学（转移受限于合法图） | 拟南芥子集 MCC 0.79→0.86（加 CRF） | 未提 |
| minisplice | 小 CNN 对每个 GT/AG 打分 | log-odds 加入对齐 DP | 分箱校准，DP 参数固定 | miniprot 未注释 junction 14.01%→4.37% | 未提 |
| GENATATORs | 逐核苷酸 5 类，无专门剪接 head | 无 HMM/CRF；splice-site filter | 全学 | 未取到 | 有（chr20 T2T 测试） |

### A9. 对我们设计的 3–5 条建议
1. **势的权重要"校准后再加"**：不要把 SpliceAI 原始概率当 log-potential。按 minisplice 的做法把分数分箱，换算成相对 null 的 log-odds（或等渗校准的 logit），再乘以每个边界状态各自的标量 λ（donor/acceptor × mRNA/lncRNA）。AUGUSTUS 的 bonus/malus 也是同构（对 log 后验加常数），但它是手设的，我们应该用数据拟合。
2. **拟合 λ 的染色体要避开 SpliceAI 见过的**：只在 chr1,3,5,7,9 上拟合/校准（SpliceAI 测试染色体，同时属于我们的训练集），再在 chr18 上做一次性确认，最后 chr20 只读一次。否则学到的系数对 chr18/chr20 偏乐观（见 A7）。
3. **plug-in 先行，端到端后置**：我们的 CRF 是 log-linear，加特征很自然；Tiberius 的端到端 +2.6 gene F1 说明联合训练有价值，但 SpliceAI 在大多数训练染色体上已拟合训练位点，联合训练会过度信任它。建议先 post-hoc 拟合少量标量（几个 λ + 偏置）；若要端到端，改用只在我们训练染色体上训练的小型剪接 CNN（minisplice 量级 7k 参数），这样没有泄漏，也绕开 PolyForm Strict 的"禁止修改/分发"。
4. **lncRNA 与 mRNA 分开校准，并先做诊断**：lncRNA 剪接信号弱（Melé 2017），mRNA 训练的分数对 lnc 位点很可能偏低。先在训练染色体上看 lnc 真位点的 SpliceAI 分数分布与 AUPRC，再决定是否对 lnc 用不同的 λ/偏置；这也是可写进论文的 negative/positive 控制。
5. **预期收益在"边界精确度"而不是 locus 检出**：Helixer 的 intron-chain F1 只有 0.17、SegmentNT 无解码，说明"逐碱基分类 + 简单后处理"拿不到精确 intron chain；Tiberius/Helixer/minisplice 的共同点是 motif 门控 + 序列级剪接信号。我们已有 GT/GC-AG gating，SpliceAI 主要应改善候选 donor/acceptor 的排序（exact ±0），评估指标要用 junction / intron-chain 而非 locus 级。

---

## B. 在 hg38 chr20 上加入 SegmentNT / Helixer / Tiberius 基线

chr20 = 64,444,167 bp (https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes)；约占 hg38 的 2.1%。

### B1. 汇总表

| 项目 | SegmentNT | Helixer | Tiberius |
|---|---|---|---|
| 代码许可 | CC BY-NC-SA 4.0 (https://github.com/instadeepai/nucleotide-transformer/blob/main/LICENSE.md) | GPL-3.0；HelixerPost 也是 GPL-3.0 (https://github.com/usadellab/Helixer ；https://github.com/usadellab/HelixerPost) | MIT（仓库 LICENSE）(https://github.com/Gaius-Augustus/Tiberius/blob/main/LICENSE)；网页服务页写 Apache 2.0（冲突）(https://bioinf.uni-greifswald.de/tiberius/) |
| 权重许可 | CC BY-NC-SA 4.0 (https://huggingface.co/InstaDeepAI/segment_nt) | CC BY 4.0 (https://zenodo.org/records/17850139) | 未写（未核实） |
| 训练数据 | 人类 GENCODE V44（去 level 3）+ ENCODE SCREEN 调控元件 | 公开基因组 + 参考注释，v0.3 脊椎动物模型 | 37 个哺乳动物组装（论文）；当前 v2 为 29 个训练物种，RefSeq 注释，每基因取最长 CDS 转录本 |
| 划分 | 人类：chr20+chr21 测试，chr22 验证，其余训练 | 按物种划分；human 在 vertebrate 模型的训练或验证集 | 按物种划分；human/cow/beluga 全程留出测试 |
| **chr20 泄漏？** | human 版：**不泄漏**（标签层面）；multi-species 版：未说明（未核实）；注意 chr18 在其训练集 | vertebrate 模型：**泄漏**（编码基因）；mammal 模型：推测不泄漏（未核实） | **不泄漏** |
| 预测 lncRNA？ | 有 "lncRNAs" 类（locus 级），3 kb MCC<0.1；multi-species 版没有 | 否（编码，最长蛋白异构体） | 否 |
| 模型大小/变体 | segment_nt（30 kb，562M 参数）、segment_nt_multi_species；预印本另有 3 kb/10 kb、SegmentEnformer/Borzoi（HF 上只放了前两个） | 4 个谱系（fungi/land_plant/vertebrate/invertebrate），另有 Zenodo v2 的 mammal 模型 | 9 个谱系配置（mammalia ×3、vertebrates 等） |
| 权重文件大小 | pytorch_model.bin 2,237.5 MB；JAX pytree 2,250 MB（各一份）(https://huggingface.co/api/models/InstaDeepAI/segment_nt?blobs=true) | vertebrate_v0.3_m_0080.h5 36.3 MB；mammal_v0.3_m_0100.h5 36.3 MB (https://zenodo.org/records/17850139) | tiberius_nosm_weights_v2.tar.gz 29.8 MB；tiberius_weights_v2.tar.gz 29.9 MB；vertebrates 36.4 MB（HTTP HEAD 实测，https://bioinf.uni-greifswald.de/bioinf/tiberius/models/） |
| 上下文长度 | 30 kb 训练（5,001 个 6-mer token），可零样本到 50 kb（需改 rescaling_factor）；输入不能含 N；DNA token 数需被 4 整除 (https://huggingface.co/InstaDeepAI/segment_nt) | 脊椎动物推荐 subsequence 213,840 bp (https://github.com/usadellab/Helixer) | 默认 seq_len 400,050；推理约 500 kb tile |
| 框架 | JAX/Haiku（GitHub）+ PyTorch（HF，`trust_remote_code`，需从源码装 transformers）；无官方 Docker 说明 | TensorFlow；Docker/Singularity 见 github.com/gglyptodon/helixer-docker | TensorFlow ≥2.17、Python ≥3.12；Singularity（自动拉取）、Docker `gaiusaugustus/tiberius` |
| GPU 显存 | 官方未给（估算：fp32 权重 2.2 GB，30 kb 窗口 batch 1–4 远低于 24 GB） | 11 GB 够（GTX 1080 Ti），8 GB 也应可跑 | ≥8 GB；README 示例 RTX 3090 (25 GB) batch 8 |
| 输出 | 每窗口 14 通道逐核苷酸概率；**无链信息**；不加互斥约束 | GFF3 | GTF/GFF3（也可输出 CDS/蛋白 FASTA） |

### B2. SegmentNT 细节
- **14 类**：protein-coding gene、lncRNA、5'UTR、3'UTR、exon、intron、splice acceptor、splice donor、polyA signal、组织不变/特异启动子、组织不变/特异增强子、CTCF 结合位点 (https://huggingface.co/InstaDeepAI/segment_nt ；https://www.biorxiv.org/content/10.1101/2024.03.14.584712v2.full)。基因元素来自 **GENCODE V44**（我们的是 v50；GENCODE 2025 论文称其 lncRNA 目录已"大幅扩充并重构"，见 https://repositori.upf.edu/items/8000aec6-2d52-4c4e-9f7d-dcb718e6daa3 ；建议评估时记录 v44/v50 标签差异）。
- **泄漏评估**：human 版的监督标签不含 chr20（测试）。两点保留：(i) 骨干 NT-v2-500m-multi-species 在 850 个 NCBI 基因组（每属 1 个，另"加选若干被深入研究的基因组"）上自监督预训练，是否含 GRCh38 未核实（只涉及序列，不涉及标签）(https://huggingface.co/InstaDeepAI/nucleotide-transformer-v2-500m-multi-species ；https://huggingface.co/datasets/InstaDeepAI/multi_species_genomes)；(ii) multi-species 版由 human SegmentNT-30kb 在 human+mouse/chicken/fly/zebrafish/worm 上微调，"每个物种分训练/验证/测试染色体"，但人类具体染色体未写。**建议只用 human 版 `InstaDeepAI/segment_nt` 作为 SegmentNT 基线**。
- **lncRNA 的可用性**：locus 级可用（lncRNA 通道），但 exon/intron 通道是 mRNA+lnc 合并；要得到"lnc 的外显子-内含子结构"只能取 lncRNA 概率高的区间内的 exon/intron 概率，且 SegmentNT 对 lncRNA 的核心判别力弱（3 kb 模型 MCC<0.1）。
- **窗口与速度（估算）**：
  - 官方说法：预印本 v2 称 SegmentNT-10kb 在单张 A100 上 0.16 ms 处理一条 10 kb 序列；InstaDeep 博客称 30 kb 输入 0.009 s 出 420,000 个预测，均未说明 batch/精度 (https://www.biorxiv.org/content/10.1101/2024.03.14.584712v2.full ；https://instadeep.com/2025/10/genome-annotation-with-segmentnt/)。
  - **这些数字与 FLOP 量级不符**：编码器约 5 亿参数 × 5,001 token × 2 ≈ 5 TFLOP/窗口（不含注意力与 U-Net），A100 BF16 峰值约 312 TFLOPS，物理下限约十几毫秒，0.009 s 超过峰值，故应视为摊销/特殊条件下的数字，**不要拿来规划**。
  - 我的估算（4090，BF16 约 40% 利用率）：每个 30 kb 窗口约 0.1–0.3 s。chr20 不重叠 2,148 窗口 ≈ 4–11 min；步长 15 kb（重叠 50%，取中间部分以避开边缘偏差）4,296 窗口 ≈ 9–22 min；再做正反互补平均 ≈ 18–45 min。无需切分显存，batch 4 以内即可。以上均为估算，未实测。
- **作者自己怎么"解码"**：没有。论文说该模型"允许直接预测潜在转录本异构体"，评估用 0.5 阈值 + 核苷酸级 F1/MCC + SOV，且对 AUGUSTUS 的比较没有 exon-level 精确匹配 (PMC12615259)。
- **要在 locus / intron-chain / junction 尺子下比较需要的解码（建议两档）**：
  1. 弱基线（"作者式"）：窗口拼接（步长 15 kb 取中央）→ 阈值 0.5 → exon/intron 连通域；donor/acceptor 取局部极大且 >τ；用 GT..AG vs CT..AC 判链；在 lncRNA 概率高的区域内组装 intron chain。
  2. 强基线（公平对照）："SegmentNT 概率作为 emission + 我们的 CRF/Viterbi/MEA"。因为 SegmentNT 输出无链信息、不保证 exon/intron 互斥，必须有外部解码器才能产生 gene model。
  - 注意事项：位置依赖偏差与 24 nt 周期（见 A3）；输入不能含 N（需对 N 区段切断）；只跑正向序列即可得到两个方向的元素，但是否做反向互补平均需自行决定。

### B3. Helixer 细节
- **模型**：vertebrate_v0.3_m_0080（README 表中的脊椎动物默认）；Zenodo v2（2025-12-07）新增 `mammal_v0.3_m_0100.h5` (https://zenodo.org/records/17850139)。README 未列 mammal 谱系，需用 `--model-filepath`（我据 README 推断，未实测）。
- **human 状态**：Tiberius 官方 yaml：Homo sapiens 在 Helixer 训练或验证中；Tiberius 论文："测试的三个物种都在 Helixer 训练用的 315 个脊椎动物基因组中" (https://github.com/Gaius-Augustus/Tiberius/blob/main/model_cfg/mammalia_nosofttmasking_v2.yaml ；https://pmc.ncbi.nlm.nih.gov/articles/PMC11645249/)。Helixer 自己的论文说 mammal 模型"与 Tiberius 用同一批物种"，human 为测试物种，但物种表在 Suppl. Table 21（不可见）(https://pmc.ncbi.nlm.nih.gov/articles/PMC13076211/)。→ 若用 vertebrate 模型，**chr20 编码基因泄漏**，结果只能作为"乐观上界"；mammal 模型应优先，但需向作者/补充表确认。
- **只是编码基因**：不预测 lncRNA；对 chr20 的比较应限定在 mRNA 位点（可作为"编码 vs lnc"对照）。
- **运行**：人类全基因组单线程 8.5 h (Nature Methods)；Tiberius 基准表里 Helixer 人类 10:56 (h:mm)（硬件未写）→ chr20 约 10–15 min（估算：按 2.1% 线性折算，含 CPU 端 HelixerPost）。GPU 11 GB 即可。输出 GFF3。

### B4. Tiberius 细节
- **human 完整留出**：论文测试物种 human/cow/beluga，验证 leopard/rat，且训练/验证里没有 Hominidae、Ruminantia、Cetacea (PMC11645249)；当前 v2 mammalia 配置的训练物种列表里没有 human，并写明 "Test species … not seen during training" (上面 yaml)。→ 作为 chr20 基线无泄漏。
- **只编码**：不预测 UTR 与 lncRNA；每基因最长 CDS 转录本；可作为 "lncRNA 只能被编码模型错判为 intergenic" 的负对照。
- **运行**：论文基准 A100 80 GB + 48 线程，人类约 1.5–2 h；README：RTX 3090 (25 GB) batch 8、RTX 2070 (8 GB) batch 2 → 4090 属同一档；chr20 约 2–10 min（估算）。
- **备注**：mammalia 有 softmasking / nosoftmasking / clamsa 三个配置；README 建议默认用不需要软屏蔽的模型 (https://github.com/Gaius-Augustus/Tiberius)。

### B5. 其他可选基线（简述）
- NTv3 annotate：需 ≥32 GB 显存，4090 不可行，仅编码，权重 gated（见 A5）。
- GENATATORs：与本项目最像（含 lncRNA、chr20 测试），但代码/权重是否公开、是否用 hg38 未核实；建议作为 related work 并核对其 chr20（T2T）评测口径。
- ANNEVO：仅编码，human 留出（上面 yaml）。

### B6. 需要用户决定的事项
1. **许可是否可接受（Nature 系列论文）**：SegmentNT 为 CC BY-NC-SA（非商业；ShareAlike 只在分发"改编材料"时触发，推理输出是否算改编材料我未能确认）；对论文基线用途通常没问题，但不要再分发其改动过的权重。Helixer GPL-3.0 / CC BY 4.0 无顾虑。Tiberius 权重许可未写、代码许可文档互相冲突，建议向作者确认。SpliceAI 源码 PolyForm Strict（不许修改/分发）：只能原样调用或自己写推理代码加载 CC BY-NC 权重；更稳妥的是自训小剪接 CNN。
2. **泄漏处理**：(a) SegmentNT 仅在 chr20 比较，chr18 数字不可用；(b) Helixer 用 vertebrate 模型须标注"human 在其训练/验证集，乐观"，或改用 mammal 模型并核实；(c) SpliceAI/OpenSpliceAI 在 chr18/chr20 的 mRNA 位点均已见过，需在论文里声明并用 A9-2 的校准办法；(d) GENCODE v44（SegmentNT 标签）vs v50（我们的标签）的差异需说明，或把评估限制在 v44 已有的基因。
3. **算力**：三者单卡 4090 均可在 1 小时内完成 chr20（估算）；是否做 SegmentNT 反向互补平均/重叠窗口。
4. **下载与大小**：SegmentNT HF `segment_nt` pytorch_model.bin 2.24 GB（可选 hg38_data.bed.gz 33.7 MB 标签文件）；Helixer mammal + vertebrate 权重各 36.3 MB；Tiberius 权重约 30 MB；容器镜像大小未核实。（另需安装 transformers 源码版 / TF 环境，本次未安装。）
5. **SegmentNT 用哪一档解码**：只报"作者式阈值解码"，还是同时报"SegmentNT 概率 + 我们的解码器"。后者更公平，也可作为"解码器贡献"的消融。
6. **SegmentNT 变体**：human 版（建议）还是 multi-species（划分未知，不含 lncRNA）。
7. **是否在论文中引用 GENATATORs/NTv3/ANNEVO 作为 related work 而不跑**。

---

## C. 验证 chr20 上新 lncRNA 候选的公开转录证据

候选 = 预测的 lnc locus 且不在 GENCODE v50。所有来源均为 hg38/GRCh38（除另标）。

| 来源 | 是什么 | 下载 URL 与文件 | 大小 | 许可 | 能否便宜取 chr20 | 循环性风险 |
|---|---|---|---|---|---|---|
| FANTOM5 CAGE peaks | hg38 的 CAGE 峰（phase 1+2），BED9，thickStart = 代表性 TSS | https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz（README 推荐此文件；另有 liftover 版 5.3 MB） | 5.5 MB（HEAD 实测 5,500,816 B） | CC BY 4.0 (https://fantom.gsc.riken.jp/5/data/) | 整文件即可；非 tabix，但 5.5 MB，zcat+awk 取 chr20 | 低：原始数据早于 v50；但 GENCODE 一直用 CAGE/polyA 佐证 TSS/TES（2012 年论文称 35% 的 TSS 有 CAGE 支持、62% 蛋白编码基因有 polyA 标注）(https://genome.cshlp.org/content/22/9/1760)，所以"已被 GENCODE 注释"的基因更易有 CAGE；对"v50 之外的候选"反而可作独立证据。我们的模型不吃 CAGE |
| PolyASite 2.0 | 3' 端测序整合的 polyA 位点图谱（人/鼠/线虫） | https://polyasite.unibas.ch/download/atlas/2.0/GRCh38.96/atlas.clusters.2.0.GRCh38.96.tsv.gz | 50.6 MB（实测 50,621,563 B；服务器 TLS 证书已过期，HTTP 抓取需忽略证书） | 数据许可未核实（论文本身 CC BY-NC 的说法来自二手页） | 整文件取回后 awk 取染色体 "20"（列：chrom, chromStart, chromEnd, name, score, strand, rep, frac_samples, nr_prots, annotation, gene_name, gene_id, repSite_signals + 各样本 TPM；染色体名无 "chr" 前缀，注释基于 Ensembl 96）(https://pubmed.ncbi.nlm.nih.gov/31617559/) | 低-中：与 GENCODE 的 TES 佐证来源有重叠，但与我们的标签/模型无关 |
| PolyA_DB v4.1 | 3'READS+ 与长读长数据的 polyA 位点库（人约 140 万 PAS） | https://exon.apps.wistar.org/polya_db/v4/download/4.1/HumanPas.zip（hg38） | 55.8 MB | 页面未写（未核实） | 整文件 | 与 PolyASite 重叠（独立性有限），作为二选一/交叉验证 |
| Snaptron（recount3 剪接位点） | 数十万份公共 RNA-seq 的 junction 数据库，hg38；`gtexv2` ≈3300 万 junction / 约 1.9 万 GTEx 样本；`srav3h` ≈2.28 亿 junction / 约 31.6 万 SRA 样本 (http://snaptron.cs.jhu.edu/data.html) | REST：`http://snaptron.cs.jhu.edu/gtexv2/snaptron?regions=chr20:START-END&rfilter=samples_count>:10`（我实测通；返回 TSV：chromosome、start、end、length、strand、annotated、left_motif、right_motif、samples_count、coverage_sum 等）(http://snaptron.cs.jhu.edu/wsi.html) | 无需整库下载；整库 bgz：gtexv2 junctions 21 GB（tbi 4.9 MB）、srav3h 127 GB (http://snaptron.cs.jhu.edu/data/gtexv2/)。chr20 的 API 体量按 3 个 1 Mb 窗口（samples_count≥10：1,002–2,926 行、5–24 MB/Mb）外推约 0.3–1.5 GB 文本（估算，样本列占大头；收紧到 samples_count≥20 可明显缩小） | 条款见仓库 license（未核实） | 原生支持区域查询（tabix/Lucene），**最便宜** | 中：junction 的 "annotated" 标记所参照的 GENCODE 版本未核实（数据目录里有 gencode.v25/v26 注释文件）；GENCODE 注释本身会参考 RNA-seq，但 GTEx/SRA 读段对我们模型是外部数据 |
| recount3 | Snaptron 的底层汇总（GENCODE v26/v29、RefSeq、FANTOM6_cat 注释） | https://rna.recount.bio/docs/ ；Bioconductor `recount3` / `snapcount` | 按研究分文件，整体很大 | 未核实 | 通过 snapcount 查询（本质走 Snaptron） | 同 Snaptron；与 Snaptron 二选一 |
| intropolis | 21,504 份 SRA 样本的 junction（hg19，Rail-RNA） | https://github.com/nellore/intropolis | hg19 TSV 6.6 GB（gz）；hg38 liftover 约 6.87 GB；已被 Snaptron 取代 | 未核实 | 否（需整文件 + liftover） | 不推荐 |
| GTEx long-read（v9） | ONT 长读长 88 个样本，FLAIR 合并后的转录本（经 TransDecoder 过滤） | https://storage.googleapis.com/adult-gtex/long-read-data/v9/long-read-RNA-seq/flair_filter_transcripts.gtf.gz（README：https://storage.googleapis.com/adult-gtex/long-read-data/v9/long-read-RNA-seq/README.txt） | 14.0 MB（GCS 清单） | GTEx 开放访问（具体条款未核实） | 整文件 | 低-中；**TransDecoder 过滤偏向有 ORF 的转录本，对 lncRNA 可能漏**（README 原文 "filtering using transdecoder"）；论文 https://www.nature.com/articles/s41586-022-05035-y |
| ENCODE4 long-read | PacBio/ONT 长读长，TALON 以 GENCODE V29 为种子给出每个重复的 GTF | ENCODE 门户，GRCh38 "transcriptome annotations" 172 个 GTF（API 查询统计） | 合计约 10.0 GB，中位约 59 MB/个 | ENCODE 数据使用政策（未核实） | 否（整文件，gz 无索引） | 低-中；不建议放进最小集 |
| LRGASP | 长读长方法比较联盟数据 | 数据存放位置我没确认（论文 https://pmc.ncbi.nlm.nih.gov/articles/PMC11543605） | 未核实 | 未核实 | 未核实 | 未核实 |
| CHESS 3.1.3 | 基于约 1 万个 GTEx RNA-seq 的人类基因目录，含 lncRNA 基因 | https://github.com/chess-genome/chess/releases/download/v.3.1.3/chess3.1.3.GRCh38.gtf.gz | 20.0 MB（GitHub release API） | 未核实 | 整文件 | 中-高：CHESS 3 含全部 MANE、并为多数 RefSeq/GENCODE 基因保留至少一个转录本，故与 GENCODE 大量重叠；只看"CHESS 有、v50 无"的子集才有意义 (https://doi.org/10.1186/s13059-023-03088-4) |
| LncBook 2.1 | 人 lncRNA 目录，52.6 万转录本（v2.0 为 9.5 万基因/32.4 万转录本） | https://ngdc.cncb.ac.cn/lncbook/downloads → lncRNA_LncBookv2.1_GRCh38.gtf.gz | 28 MB（v2.0 15.2 MB） | 学术免费，商用需授权 (https://ngdc.cncb.ac.cn/lncbook/help) | 整文件 | 高：整合 GENCODE v33、RefLnc、CHESS v2.2、FANTOM-CAT 等；用作"别处已知"而非独立证据 |
| RNAcentral | 非编码 RNA 综合库 | https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/genome_coordinates/bed/homo_sapiens.GRCh38.bed.gz（Release 28，2026-10） | 59.2 MB（GFF3 版 153.8 MB） | 未核实 | 整文件 | 高（聚合多个库，成员构成未核实） |
| NONCODE | lncRNA 目录 | http://www.noncode.org/download.php | 未核实（本机访问返回 403/空）；v5 的 hg38 GTF 命名为 NONCODEv5_human_hg38_lncRNA.gtf（二手资料）；2016 版人类 167,150 转录本/101,700 基因 (https://pmc.ncbi.nlm.nih.gov/articles/PMC4702886) | 未核实 | 未核实 | 高 |

### C2. 最小集建议（体量最小、循环性最低）
1. **FANTOM5 hg38 CAGE peaks**（5.5 MB）→ TSS 支持。
2. **PolyASite 2.0 GRCh38.96**（50.6 MB，取 chr "20"）→ TES 支持；PolyA_DB v4.1（55.8 MB）作为可选交叉。
3. **Snaptron `gtexv2` 的 chr20 junction**（API，按窗口分批；用 samples_count 过滤控制体量）→ 内含子支持，且自带 strand、left/right motif。
4. 可选：GTEx long-read FLAIR GTF（14 MB）和 CHESS 3.1.3 GTF（20 MB）→ "整条 intron chain 精确匹配"（注意 FLAIR 经 TransDecoder 过滤，对 lnc 偏保守）。
5. "别处已知"标签（不作为独立证据）：LncBook v2.1 GTF（28 MB）、RNAcentral GRCh38 BED（59 MB）。
合计固定文件约 120–180 MB（不含 Snaptron 流式查询）。

### C3. 预注册的支持标准草案（可调，需在看 chr20 之前冻结）
- **候选定义**：预测 lnc locus 与 GENCODE v50 任何基因无同链重叠（反义重叠单独标记）；两端各留 5 kb 缓冲，避免借用邻近已注释基因的 TSS/TES。
- **S1（TSS）**：同链 FANTOM5 CAGE peak 的代表 TSS（BED thickStart）落在预测 TSS ±100 bp 内（±50 / ±300 作敏感性分析；主终点只选 ±100）。
- **S2（junction）**：每个预测内含子的供体/受体坐标与 Snaptron `gtexv2` 中同链 junction 精确匹配，samples_count ≥ 5（敏感性 ≥10）且 coverage_sum ≥ 20；多外显子候选要求"全部内含子"或"≥50% 内含子"（二选一，预先写死，另一个作敏感性）。
- **S3（TES）**：同链 PolyASite 2.0 cluster 的代表位点落在预测 TES 的 −50 ~ +200 bp 内（敏感性：要求 frac_samples ≥ 0.05 或存在 PAS 信号）。
- **S4（可选，整链）**：预测 intron chain 与 GTEx FLAIR 或 CHESS 3.1.3 中非 v50 转录本完全一致。
- **综合判据**：多外显子候选 = S2 且（S1 或 S3）；单外显子候选 = S1 且 S3。主终点 = 综合支持比例；各单项作次要终点。
- **匹配零模型**：在 chr20 上从 v50 基因 ±5 kb 之外、且不与我们预测 locus 重叠的区域，抽取与候选长度分布一致（同时匹配外显子数/内含子长度分布可通过随机放置同一批预测结构模板实现）的"洗牌基因间 locus"，N ≥ 1,000 次，随机链向；比较支持比例并给经验 p；另做"链翻转"零模型检验 CAGE/polyA 的链特异性，以及"随机 GT…AG 对（内含子长度分布一致）"零模型检验 S2。
- **阳性/阴性对照**：同一套标准应用于 chr20 上已注释的 v50 lncRNA（阳性对照）和同长度的 v50 基因间区（阴性对照），用来标定"可预期的支持率"。因为 v50 之外的候选已经被 GENCODE 注释流程过滤过（GENCODE 本身就用 CAGE/polyA/长读长），候选的支持率应显著低于已注释 lncRNA，这是选择效应，不是方法缺陷。
- **循环性提示**：我们的模型/标签 (GENCODE v50) 从未使用上述证据；但 GENCODE、CHESS、LncBook、FANTOM-CAT 之间互相引用（LncBook 汇入 GENCODE v33、CHESS v2.2、FANTOM-CAT），所以"目录类"来源只能用于标注"别处已知"，不能算作独立证实。

---

## 未核实 / 来源冲突清单
- SegmentNT multi-species 版的人类训练/验证/测试染色体（论文只给人类版 chr20/21/22）。
- NT-v2-500m-multi-species 预训练的 850 个基因组是否含 GRCh38。
- SegmentNT 30 kb 模型 lncRNA 的 MCC 数值（Fig. 2a）；SpliceAI 对比 MCC 在 v2 与后续版本数字不一致。
- Helixer mammal 模型的物种表（Suppl. Table 21）与 human 是否留出；Helixer 的 HelixerPost 参数是否训练。
- Tiberius 权重许可；Tiberius 代码许可 MIT（仓库）vs Apache 2.0（网页）；GT/AG 在 HMM 里是硬约束还是加权。
- ANNEVO 的剪接/Viterbi 细节；GENATATORs 的代码/权重公开与 lncRNA 数字。
- SpliceAI 训练转录本是否仅为 protein-coding principal transcripts（Cell 论文原文）。
- OpenSpliceAI 权重许可；Snaptron / CHESS / RNAcentral / PolyASite / PolyA_DB / ENCODE long-read 的数据许可条款；LRGASP 数据位置；NONCODE 现行下载文件（本机被 403）；LNCipedia 下载链接（尝试的 URL 返回 404，未纳入）。
- 所有 chr20 运行时间均为估算，未实测。

## 关键来源索引
- Tiberius 论文 https://pmc.ncbi.nlm.nih.gov/articles/PMC11645249/ ；仓库 https://github.com/Gaius-Augustus/Tiberius ；v2 配置 https://github.com/Gaius-Augustus/Tiberius/blob/main/model_cfg/mammalia_nosofttmasking_v2.yaml ；多谱系预印本 https://www.biorxiv.org/content/10.64898/2026.04.24.720536v1.full
- Helixer 论文 https://pmc.ncbi.nlm.nih.gov/articles/PMC13076211/ ；权重 https://zenodo.org/records/17850139 ；https://zenodo.org/records/10836346
- SegmentNT 论文 https://pmc.ncbi.nlm.nih.gov/articles/PMC12615259/ ；预印本 https://www.biorxiv.org/content/10.1101/2024.03.14.584712v2.full ；HF https://huggingface.co/InstaDeepAI/segment_nt ；https://huggingface.co/InstaDeepAI/segment_nt_multi_species ；GitHub https://github.com/instadeepai/nucleotide-transformer ；偏差预印本 https://www.biorxiv.org/content/10.1101/2025.04.09.647946.full
- GeneDecoder https://arxiv.org/html/2505.03377 ；minisplice https://arxiv.org/html/2506.12986v2 ；GENATATORs https://www.biorxiv.org/content/10.64898/2026.06.17.732686v1.full ；OpenSpliceAI https://elifesciences.org/articles/107454
- AUGUSTUS extrinsic.cfg https://sources.debian.org/src/augustus/3.5.0%2Bdfsg-2/docs/tutorial/data/extrinsic.bug.cfg ；GeneMark-ETP https://pmc.ncbi.nlm.nih.gov/articles/PMC11216313
- Snaptron http://snaptron.cs.jhu.edu/data.html ；FANTOM5 https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks/ ；CHESS https://github.com/chess-genome/chess/releases ；LncBook https://ngdc.cncb.ac.cn/lncbook/downloads
