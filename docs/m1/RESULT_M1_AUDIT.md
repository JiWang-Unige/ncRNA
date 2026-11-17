# M1 实验 1：GENERanno-cds 骨干训练数据审计（2026-11-17）

> 依据：`docs/review/DECISION_20261117.md` 第 2 节实验 1。只查资料，不用 GPU。
> 要回答的唯一问题：现有 chr20 / chr9 的结果，能否称为对骨干的**严格监督样本外测试**？

## 结论：否

我们用的骨干权重就是 **Carbon-A-1.2B**。它的 CDS 监督训练数据包含人类 GRCh38 的全部染色体，没有任何染色体留出。所以 chr20 / chr9 只对我们的轻量头和语法是留出的，对骨干不是。

## 证据（2026-11-17 核对）

1. **权重完全相同**。
   - ont 上我们实际加载的 `GENERanno-eukaryote-1.2b-cds-annotator-preview`（HF 快照 1d0ef45），`model.safetensors` 的 sha256 是 `5e0a58184ac9…f71c`（4,614,654,504 字节）。
   - `HuggingFaceBio/Carbon-A-1.2B` 的 `provenance.json` 写明：来源仓库就是这个 preview（修订 b0483c2），`model.safetensors` 原样保留，sha256 与字节数与上面完全一致。
   - 因此「preview 就是 Carbon-A」不是推断，而是逐字节一致。
2. **CDS 监督数据包含人类**。
   - Carbon-A 模型卡写明：任务是逐碱基、分链的 CDS 预测；训练数据为 2,055 个 RefSeq GCF 组装（2,047 个分类单元），其中哺乳动物 242 个。
   - 数据桶 `HuggingFaceBio/Carbon-A-training-data` 的说明写明：98,304 bp 的滑动窗口，不含 CDS 的窗口也保留（所以非编码背景同样参与训练）。全文没有任何按染色体留出的描述。
   - 评测集 `eval/seen` 收录人类 GCF_000001405.40（GRCh38.p14），技术报告（bioRxiv 10.64898/2026.10.20.760448）对 seen 的定义是「from the training snapshot」，并明确说人类在其中。
3. **我们用过的其他物种也都在训练集里**。`eval/seen` 里还有：
   - 小鼠 GCF_000001635.27
   - 鸡 GCF_016699485.2
   - 果蝇 GCF_000001215.4
   - 线虫 GCF_000002985.6
   - 非洲爪蟾 GCF_000004195.4
   - 斑马鱼不在 seen 列表中，但「其他脊椎动物」有 444 个组装，很可能也在训练集里（未逐一核实）。

## 三种情况

| 情况 | 结论 | 说明 |
|---|---|---|
| 预训练见过人类序列 | 是 | 自监督，本身可以接受 |
| CDS 监督见过 chr20 / chr9 的标签 | **是** | 人类全基因组都在 CDS 训练快照里，没有染色体留出 |
| 见过同源序列 | 是 | 242 个哺乳动物组装，灵长类和其他哺乳动物的 CDS 都在内 |

## 对已归档稿件各项说法的影响

| 说法 | 影响 |
|---|---|
| 「chr20 / chr9 一次性留出测试」 | 只对头和语法成立。论文里必须写明：骨干见过这两条染色体的 CDS 标签 |
| mRNA 的检出和结构数字（chr20 mRNA L1 0.86 等） | 对骨干是样本内的，可能偏高。统一 GFF3 的编码层本来就用 Tiberius，主线不受影响 |
| lncRNA 检出 | 骨干从没见过 lncRNA 标签（只有 CDS）。但它见过这些序列「哪里是 CDS、哪里不是」，可能间接帮助区分 lnc 与 mRNA。影响大小未知 |
| 跨物种「零样本」（小鼠 28%，鸡、斑马鱼、果蝇、线虫） | 只对头是零样本。骨干见过这些基因组的 CDS 标签，「跨物种迁移」要改写 |
| 端点诊断（「后验里没有 lnc 端点信号」） | 结论方向不变：见过更多信息的骨干都定位不了，说明问题不在泄漏 |

## 下一步：是否需要「自监督骨干」对照

DECISION 原定：无法核实时加一项自监督骨干对照。现在结论是「否」，对照是否还有必要取决于实验 2：

- **实验 2 判为「撤回」**（简单检测器在两条测试染色体上都达到我们的 80% 以上）：复杂模型对 lnc 检出的增益主张本来就撤回了，骨干是否泄漏不再影响任何 lnc 主张。**不做对照**，只在文中如实交代骨干的训练数据。
- **实验 2 判为「保留」或灰区**（尤其是低 CpG 子集上我们更好）：需要分清这部分优势是骨干的真实能力，还是骨干见过测试染色体带来的。两个候选：

| 方案 | 做法 | 成本 | 局限 |
|---|---|---|---|
| A（推荐先做） | **骨干没见过的组装**：Carbon-A 的 `eval/unseen` 是 2026 年才进入 RefSeq 的 14 个组装，其中哺乳动物有长颈鹿（T2T）、牦牛、欧洲棕背䶄（*Myodes glareolus*），都带 NCBI Gnomon 的 lnc_RNA 注释。取其中一条染色体，跑我们的模型（只推理）和简单检测器，用 RefSeq lnc_RNA 作真值 | 约 1 GPU 小时；需下载一条染色体的 FASTA 和 GFF（几十 MB，需批准） | 组装没见过，但近缘物种（牛、小鼠 / 大鼠）见过；真值是 Gnomon 而不是 GENCODE |
| B | **换自监督骨干重训头**：`GenerTeam/GENERator-v2-eukaryote-1.2b-base`，6-mer 分词和 2048 维与现骨干相同，头不用改；技术报告说 Carbon-A 由 GENERator-1.2B 派生（是否正是这个检查点未核实） | 3 个 seed 约 12 GPU 小时；ont 只剩 91 GB 磁盘，存不下整套特征缓存，要改成在线计算特征 | 它是单向（因果）注意力，与现骨干不完全可比；它的预训练同样含人类序列（自监督，可以接受） |

## 可以直接写进论文的说法

- 骨干（GENERanno-eukaryote-1.2b-cds-annotator-preview，与 Carbon-A-1.2B 权重相同）在含人类 GRCh38 的 RefSeq CDS 标签上做过监督训练。
- 我们的头和语法只在 chr1–17、19、21、22、X 的窗口上训练，chr18 开发，chr20 / chr9 测试；这个划分只保证头部层面干净。
- lncRNA 标签从未用于骨干训练。

## 来源

- Carbon-A-1.2B 模型卡与 `provenance.json`：https://huggingface.co/HuggingFaceBio/Carbon-A-1.2B
- 训练数据桶说明与 `eval/seen`、`eval/unseen` 列表：https://huggingface.co/buckets/HuggingFaceBio/Carbon-A-training-data
- 技术报告：https://www.biorxiv.org/content/10.64898/2026.10.20.760448v1.full
- preview 模型卡（只写了 "More technical details are coming soon"）：https://huggingface.co/GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview
- 组装与物种的对应：NCBI Datasets API（`/genome/accession/<GCF>/dataset_report`）
- 查过但没有找到的：GENERanno-eukaryote-1.2b-base（仓库不存在）；GENERanno bioRxiv v3（只讲原核）；GitHub GenerTeam/GENERanno（没有训练数据说明）；preview 的 HF 讨论区和 GitHub issues（都是空的）。
