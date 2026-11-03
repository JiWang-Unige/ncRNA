# 预注册：chr20 上 GENCODE v50 之外的候选 lnc 有没有转录证据（冻结版）

> 状态：**冻结，2026-11-03**。用户 11-03 下午在对话中回复「按照推荐推进」，即决定卡 `bioval-data` 第一项（下载最小集，冻结后分析）。本文件提交后（commit 记在协议记录与看板）才下载证据数据。
> 下载后只计算文件大小与 sha256，写进 `docs/phase1_results/BIOVAL/MANIFEST.json` 并提交，之后才打开文件内容。此后判据不改。
> 草案（9abecd2）到本版的改动都在看任何证据之前：写死下载地址与 Snaptron 查询；补坐标约定的核对方法、阳性参照的口径、取不到数据时的处理，以及分析代码与输出位置。
> 来源：用户 2026-11-03 在决定卡 `p15-next` 上的决定，第 ③ 条：「判据预先写定，结果不回流到模型选择」。

## 1. 目的与边界

- 描述性分析：P1.5 冻结模型在 chr20 上报出、但 GENCODE v50 没有的 lnc 预测，有多大比例带有公开数据里的转录证据，并和位置匹配的随机对照比较。
- **只用已读出的预测**：`outputs/PHASE1/P15/chr20/<run>/pred.gff3`（configs/p15_v19.json，lnc 置信度 ≥ τ* = 0.6，长度 ≥ 200 nt）。不重新解码 chr20，不改 γ / τ / 任何设定。
- **不回流**：结果不用于选模型、选存档、选阈值或改解码，也不用于 v2.0 开发。v2.0 的测试数据另行锁定（决定卡 `v20-testdata`）。

## 2. 候选集

- **候选**（每个 seed 各自一套）：lnc 预测里所有外显子都不与同链任何 v50 基因（任何类别）的外显子重叠的，即 ruler 的 “novel” 预测。seed 42 有 180 个；另两个 seed 的数目在定稿前不看，以免影响判据。
- **分层**：
  - S_inter：基因区外。预测基因体 ±5 kb 内两条链都没有 v50 基因，以免借用邻近已注释基因的 TSS / TES。
  - S_other：其余候选，即落在 v50 基因体内（内含子里或反义）或其 5 kb 侧翼内。
- **主分析**：按 seed 分别计算，报 3 个 seed 的平均和逐 seed 数字。
- **次分析**：一致候选。至少 2 个 seed 在同链有外显子重叠的候选，合并为一个位点。

## 3. 证据与判据（看证据前写定）

数据源（下载地址写死；版本以文件 sha256 为准，写进 MANIFEST.json）：
- S1：FANTOM5 hg38 CAGE peaks（phase 1+2）：`https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz`（HEAD：5,500,816 B，2017-10-27）
- S3：PolyASite 2.0：`https://polyasite.unibas.ch/download/atlas/2.0/GRCh38.96/atlas.clusters.2.0.GRCh38.96.tsv.gz`（HEAD：50,621,563 B，2020-04-22；染色体名无 chr 前缀）
- S2：Snaptron `gtexv2` 的 chr20 剪接 junction。按 1 Mb 窗口查询 `https://snaptron.cs.jhu.edu/gtexv2/snaptron?regions=chr20:<s>-<e>&rfilter=samples_count>:5&fields=snaptron_id,chromosome,start,end,length,strand,annotated,left_motif,right_motif,samples_count,coverage_sum`，窗口覆盖整条 chr20（不按候选取区间，避免按结果取数），按 snaptron_id 去重。
- S4：GTEx v9 长读长 FLAIR 转录本 `https://storage.googleapis.com/adult-gtex/long-read-data/v9/long-read-RNA-seq/flair_filter_transcripts.gtf.gz`（HEAD：14,022,086 B，2023-12-08）；CHESS 3.1.3 `https://github.com/chess-genome/chess/releases/download/v.3.1.3/chess3.1.3.GRCh38.gtf.gz`（HEAD：19,980,963 B，2024-07-31）
- 「别处已知」标签，不算独立证据：LncBook 2.1 `https://ngdc.cncb.ac.cn/lncbook/files/lncRNA_LncBookv2.1_GRCh38.gtf.gz`；RNAcentral 当前版 `https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/genome_coordinates/bed/homo_sapiens.GRCh38.bed.gz`（版本号取同目录 release_notes.txt，记进 MANIFEST）。
- 某个来源下载失败（重试 3 次）时：该项记为「不可得」，其余照常；主终点中缺的那项按不支持计，并在结果里写明。

| 编号 | 证据 | 支持条件（主分析） | 敏感性分析 |
|---|---|---|---|
| S1 | CAGE 转录起点 | 同链 CAGE 峰的代表 TSS（BED thickStart）落在预测 TSS ±500 bp 内；TSS 取转录方向的 5′ 端 | ±100、±1,000 bp |
| S2 | RNA-seq 剪接 | 至少一个同链 junction 两端都在预测基因跨度内，samples_count ≥ 5 且 coverage_sum ≥ 20，且不是 v50 已注释内含子 | samples_count ≥ 10；S2-exact：至少一个预测内含子与 junction 坐标完全一致 |
| S3 | polyA 位点 | 同链 PolyASite 簇的代表位点落在预测 TES ±500 bp 内 | TES −50 至 +200 bp |
| S4 | 长读长 / 转录组组装 | GTEx FLAIR 或 CHESS 3.1.3 中不属于 v50 的同链转录本，有外显子与预测外显子重叠 ≥ 二者中较短者的 50% | 只用 GTEx FLAIR；整条内含子链完全一致 |

- **主终点**：「有转录证据」= S2 或 S4 或（S1 且 S3）的候选比例。理由：剪接读段或完整转录本直接说明该处有加工过的转录；CAGE 与 polyA 同时出现，才算界定了一个转录本的两端。只有 CAGE 一端的证据较弱（例如增强子 RNA、双向启动子），只作为次终点报告。
- **次终点**：S1–S4 各自的比例、S2-exact，以及「别处已知」的比例。

## 4. 随机对照与统计

- **位置对照**：每个候选在 chr20 上随机放 100 次。外显子–内含子结构不变，链随机。放置均匀，排除 N 区，也不与任何预测位点重叠。
  - S_inter 的对照只放在与所有 v50 基因相距 > 5 kb 的区域。
  - S_other 的对照放在同类区域，即 v50 基因体内或其 5 kb 侧翼内。
  - 随机数 seed 0。
- **链翻转对照**：原位点、反向链，用来检验 S1–S3 的链特异性。
- **统计**：报告候选支持率、对照支持率、二者之比，以及单侧经验 p 值。p 值 = 100 次对照（每次给全部候选各放一个）中支持率 ≥ 观察值的比例。
- **阳性参照**：同一套判据用于 chr20 上 v50 的 lnc 基因，以及 τ* 之后命中 v50 lnc 的预测（真阳性）。
  - GENCODE 自己就用 CAGE、polyA 和长读长做佐证，所以 v50 之外的候选支持率预计低于已注释 lnc。这是选择效应，判读时说明。

## 4b. 坐标约定与核对（只用 v50，不看候选）

- 预测与 v50 都用 ruler 的 1-based 闭区间坐标；内含子 (s, e) 为第一个和最后一个内含子碱基。
- Snaptron 的 start / end 按同一约定（第一个和最后一个内含子碱基）。核对：chr20 上 v50 所有转录本的内含子中，能在 Snaptron 结果里找到坐标完全一致 junction 的比例，在偏移 −1、0、+1 下各算一次；取比例最高的偏移，并要求它比另外两个高出一倍以上，否则停下来报告。这一步只用 v50 内含子。
- CAGE：BED thickStart 为 0-based，代表 TSS = thickStart + 1。PolyASite：代表位点列 `rep`（0-based）+ 1。S1 / S3 窗口为 ±500 bp，1 bp 的偏差不影响判定。
- 候选的 TSS = 转录方向的 5′ 端（+ 链取最小坐标，− 链取最大坐标），TES = 3′ 端。

## 4c. 阳性参照的口径

- 对象：chr20 上 v50 的 lnc 基因（所有转录本合并成一个基因跨度与外显子集合），以及 τ* 之后命中 v50 lnc 的预测（真阳性）。
- 判据与候选相同，只去掉两个「不属于 v50」的排除条件：S2 允许 junction 是 v50 已注释内含子，S4 允许转录本与 v50 相同。否则阳性参照会按定义为 0。

## 4d. 代码与输出

- 下载：`scripts/phase1/bioval_fetch.sh`（ont，nice 19；只下载与计算 sha256，不读内容）。
- 分析：`scripts/phase1/bioval.py`（CPU，nice 19）。输出 `outputs/PHASE1/BIOVAL/RESULT.json`，汇总写进 `docs/phase1_results/RESULT_BIOVAL.md`。
- 随机对照的放置：先选位置（均匀、排除 N 区、不与任何 seed 的任何预测位点重叠、满足分层的区域条件），整套外显子结构平移过去，链随机；同一候选 100 次，numpy `default_rng(0)`，候选按 seed、染色体坐标排序后依次放置。

## 5. 不允许的事

- 看到证据后再改候选定义、窗口、阈值、对照规则或主终点。
- 用结果挑模型、挑 seed、挑 τ，或改任何解码设定。
- 只报支持率高的分层或 seed：所有分层和 seed 一律报告。
