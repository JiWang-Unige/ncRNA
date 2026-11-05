# S2 计划：AlphaGenome 对照（定稿，2026-11-05，任何 AlphaGenome 输出产生之前）

> 决定卡 `alphagenome-route`：用户 11-05 选第一项，同意用 Colab CU，上限 20 CU。权重许可：用户的 Hugging Face 账号 JiUNIGE2024 已接受非商业条款。
> 只作对照。按条款，AlphaGenome 的输出不用于训练或改进我们的方法。只上传公开的 hg38 序列，不涉及任何封存数据。

## 1. 模型与运行

- **权重**：`google/alphagenome-all-folds`，Hugging Face 修订版本 a8f293a，701 MB。
  - 在 Baobab 登录节点用用户的登录状态下载，转存到 ont 的 `data/alphagenome/hf_all_folds`。12 个文件的 sha256 见 `hf_all_folds.sha256`，在 ont 上核对一致。
- **代码**：官方 `alphagenome_research`（JAX，Apache-2.0）。在 Colab 虚拟机里从 GitHub 安装，提交号记入运行日志。
  - 加载用 `dna_model.create(<本地检查点>)`，预测用 `predict_sequence(..., organism=HOMO_SAPIENS)`。
- **硬件**：Colab H100，拿不到时用 A100。
  - 先做烟雾测试：一个 1 Mb 窗口，实测显存、单窗口耗时和 CU 速率，然后才正式跑。
  - 全程 CU 上限 20，超过即停。用完按惯例清理会话（`colab_cleanup.py --orphans`），并确认 `colab usage` 显示 0 个分配。
- **输入**：hg38 的 chr20、chr18、chr7、chr9。
  - 虚拟机从 UCSC 下载各染色体 FASTA。大写序列的 md5 必须等于 ont 上本项目 hg38（sha256 c1dd8706）对应染色体的 md5，否则停下。
- **窗口**：1,048,576 bp，步长 524,288 bp。每个窗口只保留中间的 524,288 bp；染色体两端的窗口保留到端点，超出部分用 N 补齐。
- **保存的输出**（float16，按碱基）：
  - 剪接位点：`SPLICE_SITES` 的 4 条轨道，按名字与链对应到我们的列（受体 +、供体 +、受体 −、供体 −）。
  - 转录：`RNA_SEQ` 中带链的轨道，按链取所有组织的最大值，每条链一条。不带链的轨道不用。
- **训练重叠**：all_folds 用全部人类染色体（含 chr9、chr20）和 GTEx / ENCODE 的 RNA-seq 训练过。所以下面所有读数对 AlphaGenome 都是**样本内**的、偏乐观，论文中如此标注。

## 2. 分析 A：剪接位点并入 S1（描述性）

- 用 `s1_analyze.py` 原样运行于 chr20 / chr18 / chr7 / chr9，模型名为 `alphagenome`。
- 偏移：在 chr20 上按 v50 mRNA 位点 AP，在 −2..+2 中为每列各选一次，固定后用到其他染色体。规则与 Pangolin 相同。
- 报告：
  - lnc / mRNA 的 AUROC；
  - 使用量匹配后的差距；
  - 各使用量层的 AUROC；
  - 文献口径的 top-k。
- 不设门槛。要回答的问题是：用 RNA-seq 训练的模型，是否也在低使用量位点上变弱，lnc 与 mRNA 在匹配后是否仍然差别很小。

## 3. 分析 B：lnc 位点检出对照（描述性）

- **片段**：每条链上，转录信号 ≥ 阈值 t 的连续碱基为一个片段。间隔 ≤ 50 bp 的片段合并，长度 < 200 bp 的丢弃。每个片段当作一个单外显子的 lnc 预测。
  - AlphaGenome 不区分编码和非编码，所以命中 mRNA 的片段按严格 Pr 的规则算错。
- **阈值**：只在 chr18 上选。
  - 候选为 chr18 全部碱基转录信号的 50%、60%、70%、80%、85%、90%、95%、97.5%、99%、99.5% 分位数。
  - 选法与我们的 τ 规则相同：在 lnc 严格 Pr ≥ 0.43 的阈值中取 lnc L1 Sn 最大者；都达不到时取严格 Pr 最高者，并如实报告。
  - 选定后在 chr20 和 chr9 上各读一次。
- **尺子**：P1.5，`p15_score.score_preds`，bootstrap 2,000。
- **对比**：chr20 对 v1.9（Sn 0.306 / 严格 Pr 0.422），chr9 对 v2.0（Sn 0.304 / 严格 Pr 0.498）。只作描述，并注明 AlphaGenome 是样本内的。

## 4. 产出

- `results/RESULT_S2_AG.md`；S1 / S1R / chr9 的跨模型表各加一行；图 5g 加一行。
- 运行记录：虚拟机类型、实测 CU、显存峰值、单窗口耗时、代码提交号、权重 sha256。
