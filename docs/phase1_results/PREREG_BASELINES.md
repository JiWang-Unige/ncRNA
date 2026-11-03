# chr20 对照扩充：Tiberius 与 SegmentNT 的运行参数与解码规则（读 chr20 前写定，2026-11-03）

> 来源：用户 2026-11-03 下午在对话中回复「按照推荐推进」，即决定卡 `baselines-chr20` 第一项。
> 本文件提交后才下载模型、读 chr20；下面的参数此后不改。每个对照在 hg38 chr20 上只跑一次，用与 P1.5 相同的尺子（ruler v1.3，`p15_score.score_preds`，bootstrap 2,000）打分。
> 两个对照都没有在 chr18 或其他人类染色体上调过任何东西，这里也不调。chr9 不碰。

## 1. Tiberius（编码基因对照，人类在其训练中完整留出）

- 代码：GitHub `Gaius-Augustus/Tiberius` 克隆时的 commit（记进结果）；装在 ont 新建的 uv 环境 `~/envs/tiberius`（Python 3.12，TensorFlow 按其 requirements）。
- 模型：`mammalia_nosofttmasking_v2`，即不需要软屏蔽的哺乳动物模型（README 推荐的默认）；权重文件 sha256 记进结果。
- 输入：`data/human/chr20.fa`（hg38 chr20，md5 与 Phase 0 下载时的 MD5SUMS 一致；大小写原样，模型不看软屏蔽）。
- 参数：其余全用默认；`--batch_size 2`，以把显存压在约 8 GB 以内（只影响速度，不改预测）。
- 输出转换：GFF3 中每个转录本的外显子作为一个 mRNA 预测（Tiberius 只输出 CDS，没有 UTR），全部保留，不设阈值。
- 读数：mRNA L1–L4；两个护栏口径（lnc 基因被报成 mRNA、mRNA 预测只命中 lnc），对应看板 L4.b3「编码工具把 lnc 误报成编码」。Tiberius 不预测 lnc，lnc 指标为 0，不报。

## 2. SegmentNT（唯一带 lncRNA 类、在 chr20 上不泄漏的深度学习注释器）

- 模型：Hugging Face `InstaDeepAI/segment_nt`（人类版，30 kb），revision `1048ad869036bb8e50a15d82b323ebeb52e8ab92`，只用 `pytorch_model.bin` 及其代码文件。按 PyTorch 版加载（`trust_remote_code`）；ont 上先用 generanno 环境（只读，transformers 4.49）试；装不上时另建 uv 环境，版本记进结果。
- 输入切分：
  - chr20 按 N 切成不含 N 的片段，序列转大写。
  - 长度 ≥ 30,000 bp 的片段用 30,000 bp 窗口（5,000 个 6-mer token，默认 rescaling 2.44140625 正好对应）、步长 15,000 bp 平铺，最后一个窗口对齐片段末端。每个碱基取它离窗口中心最近的那个窗口的预测。
  - 短于 30,000 bp 的片段不跑，输出记 0，总长度写进结果。
  - 只跑正向序列，不做反向互补平均（作者的用法）。
- 概率：每个位置每类 softmax 后取「是」的概率；14 类中用 `protein_coding_gene`、`lncRNA`、`exon`。
- 作者式解码（主读数，阈值 0.5）：
  1. lnc 基因：`lncRNA` 概率 ≥ 0.5 的连通区间为一个基因跨度；跨度内 `exon` 概率 ≥ 0.5 的连通区间为外显子；跨度内没有外显子则整个跨度作单外显子。成熟长度 < 200 nt 的去掉（与我们的过滤相同）。
  2. mRNA 基因：`protein_coding_gene` 概率 ≥ 0.5 的连通区间为跨度，外显子同上，不设长度过滤。
  3. 链：数相邻外显子之间每个内含子的两端二核苷酸，GT…AG（及 GC…AG）记 +，CT…AC（及 CT…GC）记 −，取多数。没有内含子或票数相同的基因，两条链各报一份（主读数）；敏感性读数只保留能定链的基因。
  4. 不做平滑、不合并短间隔、不用剪接位点通道。
- 敏感性：阈值 0.3 与 0.7 重做 1–3 步；链规则的敏感性见第 3 步。
- 读数：lnc 与 mRNA 的 L1 Sn / 严格 Pr / 宽松 Pr、L2–L4、两个护栏口径。另报 SegmentNT 标签版本说明：其训练标签为 GENCODE v44，我们的真值是 v50。
- 不做「SegmentNT 概率 + 我们的解码器」：那需要在 chr18 上调参，而 chr18 在 SegmentNT 的训练集里。

## 3. Helixer

- 不重跑。沿用 Phase 0 在 hg38 chr20 上的 Helixer 输出，标注「人类在其训练 / 验证集里，结果偏乐观」。

## 4. 算力与顺序

- ont，经 `gpu_gate.sh`；排在自训剪接网络训练之后。估计 Tiberius 约 10–30 分钟、SegmentNT 约 20–60 分钟 GPU（未实测）。
- 下载：Tiberius 代码与权重约 30 MB、TensorFlow 环境约 3 GB；SegmentNT 权重 2.24 GB。下载与环境安装在本文件提交之后进行。
