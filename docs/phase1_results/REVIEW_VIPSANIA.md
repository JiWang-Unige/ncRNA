# Vipsania（无监督深度基因预测）系统梳理：能给 lnc 借什么（2026-11-03）

> 起因：用户 11-03 晚在对话中提出，ANNEVO、Carbon-A 等编码基因方法已经参照过、在 lnc 上没有好效果，Vipsania 还没有系统梳理，建议梳理参照。
>
> 范围：方法、结果、与我们已参照方法的异同、对 lnc 的逐条评估、建议的下一步。
>
> 只读：读了 GitHub 仓库的 README、docs/、configs/ 和 vipsania/model/ 下的源码（raw 文本），以及 bioRxiv 预印本全文（网页）。没有安装或运行任何东西，也没有下载模型。论文中的数字经网页摘要工具读取，**进论文前要对原文逐条核对**；从代码推断的内容标「推断」。
>
> 来源：仓库 https://github.com/Gaius-Augustus/Vipsania （MIT，2026-07 建，最近推送 10-21）；预印本 Krieg, Becker, Saenko, Diehl, Stanke，"Vipsania: Unsupervised Deep Gene Finding"，bioRxiv 2026-08-30，https://www.biorxiv.org/content/10.64898/2026.08.26.747235v1 。与 Tiberius 同组（Greifswald，Stanke）。

---

## 0. 速览

1. **它是什么**：只预测**蛋白编码基因**的 ab initio 基因预测器。训练时**从不看任何注释**，只做掩码碱基预测（遮 5% 的碱基，让模型从上下文猜）。网络中间插一层可微 HMM；基因结构作为 HMM 的隐状态「涌现」出来，因为它有助于猜对被遮的碱基。
2. **人类是它的测试物种，不在训练集里**。我们用到的人、小鼠、鸡、斑马鱼、果蝇、线虫都是测试物种。训练集里离人最近的是倭黑猩猩，而且只用序列。
   - 人类 locus F1：预训练 0.412，目标基因组微调后 0.470；exon F1 0.706 / 0.719（仓库 `docs/test_species.tsv`，Vipsania 1.0.0）。
   - 论文自述：脊椎动物是 17 个类群里唯一不领先的，原因是基因和内含子长。
3. **只建模编码基因**：
   - 18 态 HMM：基因间 1、内含子 6（每个相位两个）、编码外显子 3、受体位点 3、供体位点 3、起始 1、终止 1。没有 UTR 态，也没有非编码外显子态。
   - 评估只看蛋白编码基因的 CDS。论文没有提到 lncRNA、ncRNA 或假基因。
4. **结构能「涌现」，靠的是密码子语法**：
   - HMM 的发射被手工设定为：只允许规范剪接位点和起止密码子，外显子内不许出现读框内终止密码子。
   - 论文的解释是，在依赖基因结构的位点上（外显子边界、第三位密码子），被遮的碱基更好猜。
   - lnc 外显子没有这种可被利用的语法。我的推断：Vipsania 这类无监督结构学习不能直接解决 lnc 外显子结构（我们的 L2 ≈ 0）。这与两次剪接诊断的结论一致：lnc 剪接位点在序列上难以区分，SpliceAI 的 AP 是 0.140，自训网络 0.161。
5. **最值得借的一点：不依赖我们外显子结构的「编码证据」**：
   - v2.0 草案 §3.3 的编码潜能特征（最长 ORF、六聚体分、Fickett）都在**我们自己预测的转录本**上计算。可 v1.9 的 lnc 剪接点只有 0.4–1.5% 是对的，L2 = 0，等于在几乎全错的外显子链上算 ORF。
   - 像 Vipsania 这样的基因预测器会在同一位点上**独立**找出最优的编码基因结构，从而给出「这里像编码基因」的证据，不需要我们的结构是对的。
   - 这正对着 v2.0 要解决的精度问题：lnc 预测落在编码基因上。chr20 护栏 10–12 个；P1.6 线虫 seed 42 的 2,167 个 lnc 预测里有 1,197 个落在编码基因上。
6. **但人类上 Tiberius 更强，而且 ont 已经装好**：
   - chr20 上 Tiberius 的 mRNA L1 Sn 0.964 / 严格 Pr 0.917、L2 0.557，人类在它的训练中完整留出。
   - Vipsania 的独特优势（不要注释，每个类群都有模型）在非模式物种上更有用。
   - 所以要验证第 5 点，诊断应同时比较两者，否则审稿人会问「为什么不用 Tiberius」。
7. **次要可借**：
   - 作为 chr20 上「没见过任何标签」的编码对照，F5g 护栏面板多一行。
   - 写作上可借两种分析：按与训练集的距离画斜率；用训练集重叠作对照。论文报告 Tiberius 在自己训练过的真菌物种上高 3.8 个百分点，在没训练过的物种上低 5.9 个百分点。这支持我们把 Helixer 的 chr20 结果标为偏乐观。
8. **不建议现在借**：
   - 把 HMM 放进网络中间、用掩码预测去学 lnc 结构：要重新设计并预训练，按第 4 点预计学不到 lnc 信号。
   - 把 Vipsania 的中间表示作为我们模型头部的输入：要重训，留到 v2.x 以后。
   - 在目标基因组上做无监督微调：我们的骨干有 1.2B 参数，代价高，暂缓（发散）。

---

## 1. 方法要点

### 1.1 训练目标

- 遮罩式语言模型（BERT 式交叉熵）。窗口 T = 20 kb，最多遮 5% 的碱基，N 不遮。重复区（软屏蔽小写）的损失权重为 1e-4（`configs/base_*.json`）。
- 「spliced loss」：用 HMM 估出的外显子边界概率，经一个矩形核平滑后，加重边界附近猜错碱基的损失。
  - 实现在 `vipsania/model/loss.py`：`SplicedCrossentropyAdapter` 按后验质量给每个位置的样本权重加项。
  - 权重 γ 前 100 epoch 为 0，随后 5 epoch 内升到 0.05。
  - 论文说**脊椎动物模型不用 spliced loss**，改用一个辅助损失：按转移矩阵算出内含子与基因间区的期望长度，惩罚两者差异过大；并训练 3000 epoch。
  - 作者认为这不违反「无监督」，因为边界概率来自只靠碱基预测训练出来的 HMM。

### 1.2 架构

- 主体由双向 LRU（线性循环单元）、SwiGLU 前馈层和 LayerNorm 组成，单碱基分辨率，没有卷积编码器。
  - 25M 模型：16 层，宽 320。仓库还提供 10M 模型（10 层，宽 256）和 2M 模型。
- HMM 层放在网络中间（25M 模型在第 8 层后）。
  - 前后向算法向量化、可微。
  - 状态后验经一个宽 9 的卷积（readout）加回残差流，后面的层可以利用它（`vipsania/model/hmm.py`、`base.py`）。
  - 正反两条链用两份参数共享的 HMM 同时建模。

### 1.3 HMM 的状态与参数

- 18 态 = Tiberius 的 15 态，外加每个相位两个内含子态，以便同时学短内含子和长内含子。
  - 具体排列是我按 `loss.py` 的下标推断的：`p[..., -8:]` 为 8 个位点态，即起始、3 个受体、3 个供体、终止；`p[..., -11:]` 在此基础上再加 3 个外显子态；`p[..., 1:]` 是除基因间外的全部状态。
- 发射以三联体编码，每个状态看到左右各两个碱基。
  - 三联体发射矩阵**手工设定且固定**：强制规范剪接位点和起止密码子，外显子内禁止读框内终止密码子。
  - 学习的部分是一个从残差流到发射字母表的 softmax 映射。
- 转移矩阵、发射字母表映射和初始分布都参与训练，且不做权重衰减。
- 配置示例（`configs/base_10M.json` 的 `hmm` 段）：初始外显子长 1,000，内含子长 [100, 10,000]，基因间长 10,000；非编码态共享发射；读框之间共享。

### 1.4 目标基因组微调与推理

- `vipsania annotate <模型> genome.fa --finetune`：在要注释的基因组上继续做 10 个 epoch 的无监督训练。论文给的时间是 A100 上约 160 分钟。超过 1 Gbp 的基因组只用重复比例低于 25% 的序列。
- 推理：在 200 kb 的 tile 上做 Viterbi，tile 边界只允许落在基因间区或内含子里，接缝处重新预测。果蝇（144 Mb）在 A100 上不到 6 分钟。
- 输出：GFF3、GTF 或 GenePred，可以同时写出蛋白和 CDS 序列。
- `vipsania embed`：导出任意层的残差流（HDF5）。25M 模型默认取第 7 层、第 2 条（从 0 计），也就是 HMM 层的输出。
- Python 接口里 `toggle_inference(mode=POSTERIOR)` 可以直接取 HMM 的状态后验（`base.py`；CLI 没有直接暴露，属推断）。
- 依赖：Python 3.12、TensorFlow（带 CUDA）。pip 安装约 5 GB，模型约 100 MB，首次使用时自动下载到 `~/.cache/vipsania/models`。

### 1.5 训练规模

- 共 17 个类群，1,559 个基因组，约 910 Gbp。每个类群一个模型，脊椎动物模型用了 200 个物种。
- 25M 模型训练 10 天、1,000 epoch（论文）。

---

## 2. 结果中与我们相关的部分

| 物种（Vertebrata 模型 etb1go6q） | 模式 | exon F1 | locus Sn | locus Pr | locus F1 |
|---|---|---|---|---|---|
| 人（GCF_000001405.40） | 预训练 | 0.706 | 0.407 | 0.418 | 0.412 |
| 人 | 微调 | 0.719 | 0.450 | 0.492 | 0.470 |
| 小鼠 | 预训练 / 微调 | 0.733 / 0.749 | | | 0.366 / 0.439 |
| 鸡 | 预训练 / 微调 | 0.804 / 0.822 | | | 0.438 / 0.482 |
| 斑马鱼 | 预训练 / 微调 | 0.833 / 0.851 | | | 0.452 / 0.528 |

来源：仓库 `docs/test_species.tsv`。参考注释的版本没有写明，从 GCF 编号看应是 RefSeq（推断）。另外两个测试物种：果蝇（Insecta 模型）locus F1 为 0.672 / 0.779，线虫（Nematoda 模型）为 0.719 / 0.752。

论文的比较（只读到了正文摘要，**需核对原文**）：

- 「与对手持平或更好」的物种数：对 Tiberius 38/60，对 ANNEVO 82/92，对 Helixer 91/92，对 GeneMark 123/126。17 个类群的平均 BUSCO：Vipsania 86.2%，Tiberius 80.4%，ANNEVO 79.8%，Helixer 75.7%。
- 与训练集的序列一致性每下降 10 个百分点，Helixer、ANNEVO、Tiberius 的 F1 分别降约 3、4、5 个点；Vipsania 的斜率显著更平。
- 脊椎动物是唯一不领先的类群；Tiberius 在重复序列少的基因组上更好。
- 没有做「去掉 HMM 层」的消融，所以 HMM 层是否必要没有直接证据。

与我们 chr20 的数字**不能直接比**：

- 尺子不同。他们用全基因组 RefSeq CDS，我们用 GENCODE v50 和 P1.5 尺子。
- Tiberius 的人类 gene F1 65.9（Tiberius 仓库 mammalia v2 配置）同样不能和 Vipsania 的 locus F1 0.41–0.47 直接比。不过两者都在人类上测过，量级上 Tiberius 更强，与论文自述「脊椎动物不领先」一致。

---

## 3. 与我们已参照方法的对比

| 方法 | 监督信号 | 结构怎么进模型 | 有无非编码外显子 / UTR / lnc 态 | 人类是否泄漏 | 我们已做的 |
|---|---|---|---|---|---|
| Tiberius | 注释（每基因最长 CDS 转录本） | 末端可微 HMM，端到端训练（比不带 HMM 的 gene F1 高 2.6） | 无，15 态只有编码 | 人类完整留出 | chr20 对照（RESULT_BASELINES） |
| Helixer | 注释 | 网络不知道 HMM，后处理 HMM（HelixerPost） | 有 UTR 类，无 lnc | 公开 vertebrate 模型训练用过人类 | chr20 重打分（偏乐观） |
| ANNEVO | 注释（MoE 基因组语言模型） | 解码细节在补充材料，未核实 | 只做编码 | 人类留出（其 yaml） | Phase 0 T2T 评测；FP 审计引用其 18.3% |
| SegmentNT | 注释（GENCODE V44 分割） | 没有解码器，阈值 0.5 | 有 lncRNA 类（我们在 chr20 上核查的逐碱基 AUROC 只有 0.649） | chr20 是其测试集 | chr20 对照 |
| Carbon-A（HF Bio，2026-10） | 注释（由 GENERanno-eukaryote-1.2b-cds 注释器衍生，与我们的骨干同支） | 官方解码器未发布，社区复现用 Viterbi | 只做 CDS，无 UTR | 未核实 | Phase 0 选骨干时参照；文献核查（NCRNA_ANNOTATION_LANDSCAPE_2026） |
| GENATATOR | 注释 | 生成式分割 | 有 lnc_RNA 类，但 chr20 lnc Sn 只有 0.028 | 未核实 | 主对照 |
| **Vipsania** | **无（只用序列）** | **HMM 在网络中间，后验回流残差流，用掩码预测训练** | **无**：18 态只有编码，评估只看 CDS | **无标签，人类不在训练集** | 本文档 |

一句话：截至 2026 年，监督的（Tiberius、Helixer、ANNEVO、Carbon-A）和无监督的（Vipsania）深度基因预测器都只建模编码基因；唯一带 lnc 类的 SegmentNT 在 chr20 上 lnc 位点 Sn 0.108、严格 Pr 0.113，远低于 v1.9（0.306 / 0.422）。Vipsania 没有改变「lnc 结构没人做好」这个判断，反而从侧面支持了它：结构之所以能在不看标签的情况下学出来，是因为编码基因有密码子语法。

---

## 4. 对 lnc 的逐条评估

| # | 可借的点 | 对应我们的问题 | 评估 | 建议 |
|---|---|---|---|---|
| A | 用它的编码基因预测和 HMM 后验作为**基因层面的编码证据**，进 v2.0 校准 | 精度：lnc 预测落在编码基因上；草案的 ORF 特征算在错误的外显子链上 | 不依赖我们的结构；无标签，chr9 不会泄漏。要和 Tiberius（已装）比，人类上 Tiberius 可能更好 | **推荐做 chr18 小诊断（§5）**，按判据决定是否改 §3.3 |
| B | 作为 chr20「无标签」编码对照 | F5g 护栏：多少 lnc 被报成编码 | 用同一把尺子，几分钟 GPU。能回答「监督工具很少把 lnc 报成编码，是从标签里学来的，还是序列语法本身就能做到」 | 可与 A 共用一次安装，描述性 |
| C | 用它的受体 / 供体位点后验做剪接信号 | 结构：lnc 剪接位点区分 | 受读框约束，lnc 内含子大多不在 ORF 里。预计 AP 不会高于 SpliceAI（0.140） | 只在 A 的同一次运行里顺带报告 AP，不进任何方案 |
| D | HMM 放进网络中间、用掩码预测训练（结构从序列中涌现） | 结构：lnc L2 ≈ 0 | 涌现靠密码子语法，lnc 没有；我们的 AB6 G2-PH（监督的 HSMM 训练）在 chr20 已失败（lnc Sn 0.036） | 不做 |
| E | 它的中间表示（HMM 层输出，320 维）作为我们模型头部的额外输入 | lnc 外显子信号弱（Phase 0.4：AUPRC 只是先验的 2.4 倍） | 要重训，改模型结构 | 暂缓（发散），v2.0 之后再考虑探针 |
| F | 在目标基因组上做无监督微调 | 跨物种：P1.6 线虫 0.073 | 我们的骨干 1.2B 参数，代价高；P1.6 已冻结完成 | 暂缓（发散） |
| G | 按与训练集的距离画斜率；用训练集重叠作对照 | 写作 | 我们 P1.6 只有 4 个物种，画不了斜率；重叠对照可以引用来支持 Helixer 偏乐观 | 写进论文的讨论和局限 |

关于 A 的补充：线虫上它的意义可能更大（Nematoda 模型的线虫 locus F1 0.72–0.75；我们线虫 lnc 严格 Pr 只有 0.008，一半以上的预测落在编码基因上）。但 P1.6 已冻结，v2.0 也只针对人类 chr9。如果 A 在人类上有效，跨物种的应用要另起协议。

---

## 5. 建议的 chr18 小诊断（待决定卡 `vipsania-ref` 批准；判据在运行前写定）

> **11-03 晚已运行**（用户选选项 1）。健全性通过（Vipsania mRNA L1 Sn 0.865）。V 和 T 都不有效（ΔAUPRC −0.013 / +0.002），v2.0 按原草案。位点 AP 0.037。详见 `RESULT_VIPSANIA.md`。

**目的**：判断「不依赖我们外显子结构的编码证据」能不能帮 v2.0 校准，区分 v1.9 的 lnc 真阳性和假阳性。只用 chr18（开发染色体），chr9 不碰，chr20 / T2T 不用于任何选择。

**输入**：
- v1.9 三个主模型（g1b-42 / 123 / 456-ke，第 10–12 epoch 权重平均）在 chr18 上的 MEAX lnc 预测（γ = 3，τ 过滤之前；已有）。
- 编码证据来源两个，都在 hg38 chr18 上跑一次，两条链：
  - **V（Vipsania）**：Vertebrata 模型 etb1go6q，预训练，不做目标基因组微调。输入原样使用我们的 hg38（软屏蔽，sha c1dd8706）。输出 GFF3；若 Python 接口可行，另存编码态后验（3 个外显子态加 8 个位点态之和，每链一条，float16）。
  - **T（Tiberius）**：ont 上已装（b3ea48f2），输出 GFF3。

**单位与目标**：每个 lnc 预测为一行。目标 y 为严格 L1 命中，即与 v50 lnc 基因同链外显子重叠，与 DRAFT_V20 §3.4 相同。

**特征**：
- 基础集 B：都不需要训练表。
  - v1.9 置信度（lnc 态后验均值）；
  - 外显子上 mRNA 态后验均值（若 chr18 后验 dump 可用）；
  - log 成熟长度、外显子数；
  - 预测转录本上的最长 ORF 长度和占比、Fickett 分。
  - 六聚体分要用训练染色体统计频率表，属于 v2.0 开发，本诊断不算，作为局限写明。
- V：
  - 预测位点跨度内、同链被 Vipsania CDS 覆盖的碱基比例；
  - 同链是否有 Vipsania 基因重叠（0/1）；
  - 若有后验：预测外显子上编码态后验的均值。
- T：同链被 Tiberius CDS 覆盖的比例，以及是否重叠（0/1）。

**评估**：
- chr18 按 v1.9 开发时的前后两半切分（< 40 Mb 和 ≥ 40 Mb）。
- 在一半上拟合 logistic（L2，C = 1，特征标准化），在另一半上算 AUPRC，两个方向都做，每个 seed 分别做。
- 比较 B、B + V、B + T、B + V + T。

**判据（写定）**：
- 某个来源 X「有效」的条件：B + X 相对 B 的折外 ΔAUPRC，3 个 seed 平均 ≥ 0.02，且 6 个（seed × 半边）中至少 5 个 > 0。
- V 和 T 都有效时，取平均 Δ 大的那个；两者相差 < 0.01 时取 T（不新增依赖）。B + V + T 相对单个来源再提升 ≥ 0.01 时才两者都用。
- 有效 → 修改 DRAFT_V20 §3.3，把该来源加为编码潜能特征，在 v20-freeze 卡上说明后再请你冻结。
- 都无效 → v2.0 按原草案，诊断结果写进论文。

**健全性门槛（先看，不过则先查用法）**：Vipsania 在 chr18 上的 mRNA L1 Sn ≥ 0.5（同一把尺子）。它在人类上的 locus F1 约 0.41，位点层面的 Sn 应明显更高。

**只报告、不进任何方案**：
- V 的受体 / 供体位点后验在 lnc 真值剪接位点上的 AP，与 SpliceAI 的 0.140、自训网络的 0.161 用同一函数 `own_splice.site_ap` 比较。
- V 在 chr18 上把多少 v50 lnc 基因报成编码（护栏口径 1 / 2）。

**资源**：
- 下载：pip 包 vipsania 及其依赖（TensorFlow + CUDA 约 5 GB），装进 ont 上独立的 Python 3.12 环境（ont 有 /usr/bin/python3.12 和 uv）；Vertebrata 模型约 100 MB。驱动 595.58 满足要求。
- GPU：TensorFlow 默认会占满显存。要用包装脚本先开 `TF_FORCE_GPU_ALLOW_GROWTH`，并把单进程显存上限设在约 7 GB（30%），走 gpu_gate.sh，给 TE 让路。若在上限内跑不动（减小 `-B` 和 `-T` 后仍不行）就停下来报告。
- Tiberius 在 chr20 上用过约 8.7 GB，超过上限。chr18 改用 batch 1；仍超就只做 V，并如实写明。
- 预计 GPU 合计约 1 小时，CPU 分析半天，墙钟约 1 天。

**若同时批准 chr20 对照（B 项）**：用同一安装在 hg38 chr20 上跑 Vipsania（预训练，与 Tiberius 一样不微调），用 P1.5 尺子和 RESULT_BASELINES 的口径打分（lnc 与 mRNA 的 L1 / L2、两种护栏），只描述，进 F5g。

---

## 6. 论文里怎么引用（草句）

- Related work：「Recent deep gene finders—supervised (Tiberius, Helixer, ANNEVO, Carbon-A) and, most recently, unsupervised (Vipsania, trained by masked-nucleotide prediction with an in-network HMM)—model protein-coding genes only; their state spaces contain no non-coding exon or UTR states, and their evaluations score CDS only.」
- Discussion：「Vipsania shows that coding gene structure can emerge from sequence alone, because the codon grammar makes nucleotides at structure-dependent positions predictable. lncRNA exons lack such a grammar; consistent with this, the splice sites of lncRNAs were poorly separable from motif-matched decoys by two sequence-only splice models (AP 0.14–0.16), and no method recovered lncRNA intron chains (L2 ≈ 0).」
- 局限 / 基线：引用其「训练集重叠」分析，支持把 Helixer（训练含人类）的 chr20 结果标为偏乐观。

---

## 7. 未核实与局限

- 论文正文经网页摘要工具读取，具体数字（胜出物种数、BUSCO、斜率、3.8 / 5.9）进论文前需对照原文和补充表。论文正文没给脊椎动物与 Tiberius 的逐物种数字。
- 18 态的具体排列是从代码下标和「Tiberius 15 态 + 每相位两个内含子态」推断的，论文图 2 没有读到文字版。
- 能否直接取 HMM 后验（Python 接口）、在 7 GB 显存上限内能否跑 200 kb 的 tile，都要装上后才知道。
- 测试物种的参考注释版本没有写明（推断为 RefSeq），它的 locus F1 与我们 GENCODE v50 + P1.5 尺子的数字不可比。
- 结论 4 和表 4 的 C、D 是推断（lnc 外显子缺乏可被掩码预测利用的语法），没有直接实验。如果做 §5 的诊断，C 项的 AP 能给一个间接检验。
