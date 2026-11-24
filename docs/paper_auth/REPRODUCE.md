# 复现说明（草案，2026-11-24）

> 用途：论文的 Data / Code availability 和补充表 S2 的底稿。投稿时把仓库打一个版本标签并在 Zenodo 归档（需作者本人用自己的账户操作），把 DOI 填进正文。

## 1. 数据（全部公开）

| 用途 | 来源 | 编号 / 文件 |
|---|---|---|
| WTC11 起点证据 | GEO GSE185917（LRGASP ssCAGE） | GSM5625727–GSM5625732 的 `*.count.bed.gz` 与 `*.ctss.bed.gz` |
| WTC11 终点证据 | ENCODE ENCSR322MWL（QuantSeq） | ENCFF603KSV、ENCFF937BYX、ENCFF004EPK（polyA sites） |
| WTC11 真值与「不确定」 | ENCODE ENCSR507JOF（LRGASP PacBio） | rep2：ENCFF527KWD（定量）+ ENCFF798DWV（GTF）；rep3：ENCFF486EKE + ENCFF422CGK |
| WTC11 长读长参照、E1 分子连接 | 同上 | rep1：ENCFF373TKM（BAM） |
| H1 起点证据 | FANTOM5 hg38 重处理 CTSS | CNhs14067、CNhs14068、CNhs13964（H1 day00）；峰集 `hg38_fair+new_CAGE_peaks_phase1and2.bed.gz` |
| H1 终点证据 | ENCODE ENCSR813AOB | ENCFF970YYF、ENCFF138ZWV、ENCFF215NMJ |
| H1 真值 | ENCODE ENCSR271KEJ | rep1：ENCFF580WQM + ENCFF721HVH；rep3：ENCFF014CCY + ENCFF236CGS |
| 编码掩码 | GENCODE | v41 comprehensive annotation（protein_coding 基因） |
| 基因组 | UCSC hg38 | 软屏蔽 `hg38.fa.gz` |

下载脚本和 md5：`scripts/auth/v1_fetch.sh`（H1 / GM12878 / FANTOM5），`docs/m1/W24/ACCESSIONS.json`（WTC11）。

## 2. 代码与运行顺序（`scripts/auth/`）

| 步骤 | 脚本 | 产出 |
|---|---|---|
| 1 真值 | `auth_truth.py` | WTC11 真值位点 |
| 2 区域开发 | `auth_dev.py` | 冻结的 AUTH 模型（`MODEL.json`） |
| 3 区域测试 | `auth_score.py --mode test`（长读长参照：`auth_queue.sh`） | AUTH 结果 |
| 4 跨细胞系检查 | `v1_check.py`、`v1_check_supp.py` | GM12878 / H1 可行性 |
| 5 迁移 | `v1_dev.py` → `v1_score.py --mode freeze / test` | v1 结果 |
| 6 配对 | `e1_labels.py` → `e1_dev.py` → `e1_score.py --mode freeze / test` | E1 结果 |
| 7 事后分析 | `e1_posthoc.py`、`e1_posthoc_g10ci.py`、`g1_region.py`、`paper_extra.py` | 补充分析 |

共用库：`auth_lib.py`、`v1_lib.py`，以及 `scripts/m1/w24/`（W24 证据与真值规则）、`scripts/m1/simple_det.py`、`scripts/phase1/b1/b1_strict.py`、`scripts/phase1/s3/s3_features.py`。

冻结记录：AUTH、v1、E1 各自的 `FROZEN`（文件 sha256），在 `docs/auth/RESULT/`、`docs/auth/v1/RESULT/`、`docs/auth/e1/RESULT/`。

## 3. 环境

- Python 3.10.19，numpy 2.2.6，scipy 1.15.3（不依赖 GPU，也不需要深度学习框架）。
- samtools 1.25；StringTie 2.1.7；IsoQuant（`auth_queue.sh` 调用，版本号待核对）。
- 全部计算在一台 32 核 Linux 服务器上用 CPU 完成。主要步骤每步几分钟到半小时。

## 4. 要随论文一起归档的中间文件（Zenodo）

| 内容 | 大小 |
|---|---|
| 冻结模型：AUTH `MODEL.json`、v1 `MODEL_V1.json`、E1 `MODEL_E1.json` | < 1 MB |
| 候选单位与特征（`outputs/AUTH/pairs/`） | 约 210 MB |
| E1 分子连接标签（训练 + 测试） | 约 66 MB |
| 全部结果 JSON 与日志 | < 10 MB |

## 5. 还没做的

- 仓库加许可证（建议 MIT）和面向读者的 README；
- 投稿时打版本标签，并在 Zenodo 归档拿 DOI（需作者本人操作）；
- 把服务器上的路径（`~/ncrna-annot`）在 README 里说明为可配置的根目录（脚本已支持环境变量 `NCRNA_ROOT`）。
