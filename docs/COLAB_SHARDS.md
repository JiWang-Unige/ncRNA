# Colab 分片路线（P1.1 正式特征，2026-10-22 起）

用户决定（看板 colab-topup）：先用现有 Colab 余额（约 36 CU）把分片流程验证通，再决定是否充 200 CU。

## 为什么要分片
正式特征（GENERanno-1.2b-cds 第 17+23 层，4,096 维，fp16）全量约 290 GB，放不进 Colab 磁盘；Colab VM 又是临时的。
所以把缓存切成自描述的分片，在 VM 里算一片、下载一片、删一片，最终在 ont 的 `data/p1_feat/` 合并成与
`extract_p1.py` 完全相同的 memmap 布局。VM 里不放任何密钥：下载走 `colab download`。

## 脚本
- `scripts/phase1/extract_p1_shard.py`：在任何 GPU 机器上算一片（train 窗口段或染色体 tile 段），写 `<stem>.npy` + `<stem>.json`
  （sha256、形状、tile 索引、GPU、计算精度、每窗口秒数、可选与参考分片的余弦）。stem 带全局索引。
- `scripts/phase1/assemble_p1_shards.py`：ont 上校验 sha256、写进 memmap、维护 `*.cover.json`；`--ref-dir` 做跨机器余弦检查
  （协议 v1.4 规则：中位余弦 ≥ 0.99，否则停下来查精度，不做判读）。
- `scripts/phase1/colab_p1_shards.sh`：ont 上的驱动：开 VM → 上传 src + 所需 codes 片段 → 后台跑 → 轮询 → 下载 → 校验 → 释放 VM（trap 保证清理）。
- `tests/test_shards.py`：CPU 桩测试（tile 行号、全局索引、sha256 拒绝、cover 文件）。

## 验证步骤（约 0.3–0.5 CU）
1. `bash scripts/phase1/colab_p1_shards.sh train 0:20 T4` → `data/p1_shards/colab/train_00000_00020.*`
2. `bash scripts/phase1/colab_p1_shards.sh chrom 0:10 T4 1800 chr18` → `chr18_t00000_00010.*`
3. 参考分片：ont 4090（bf16，GPU 空闲时，约 1 分钟）`python3 scripts/phase1/extract_p1_shard.py --kind train --codes data/p1/train_codes.npy --range 0:20 --out data/p1_shards/ont`
4. `python3 scripts/phase1/assemble_p1_shards.py --shards data/p1_shards/colab --ref-dir data/p1_shards/ont --check-only`
   通过判据：sha256 一致；中位余弦 ≥ 0.99；下载吞吐可接受（全量 290 GB 在 ≤ 10 MB/s 时要 8 小时以上，需改走分片即训即删或 HF 中转）。
5. 记录每窗口秒数 → 外推 3,000 窗口 + 三条染色体的 CU 消耗，写进看板 L3.features，再决定充值。

## 全量时的存储策略（验证后定）
训练集 3,000 窗口约 100 GB 可一次放进 Colab 付费运行时磁盘（约 235 GB）；chr18 / chr20 / T2T chr20 逐条处理（提取 → 评估 → 删除）。
若训练也在 Colab 做，分片需能在新 VM 里重新装载：候选是 Google Drive（CLI 是否能挂载待核）或 HF 私有数据集中转。
