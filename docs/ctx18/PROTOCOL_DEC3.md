# PROTOCOL_DEC3：第二套 total-RNA 复核（ENTEx total vs GTEx poly(A)，2026-11-25）

作者 2026-11-25「先找第二套 total-RNA 数据复核」。本协议在计算任何 Decima 与 ENTEx / GTEx 的比较之前冻结。此前只看了：ENTEx 文件清单与组织名、GTEx 列名、Decima ENCODE 轨道的组织名。

## 1. 数据
- **total-RNA**：ENTEx（ENCODE，4 位供体）total RNA-seq 的基因定量文件，GRCh38 / GENCODE V29，共 97 个文件（91 个实验），md5 全部与 ENCODE 一致，存于 ont `~/ncrna-annot/data/ctx18/entex/gq/`。用 TPM。
- **poly(A)**：GTEx v8 各组织中位 TPM（`GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz`，6.95 MB）。ENTEx 供体本身是 GTEx 供体，但这里用的是 GTEx 全体中位数，不是同一批样本。
- **Decima**：补充文件 1 的 `layers/preds`（作者集成预测），只用表内基因；按 organ 取所有轨道（不限疾病状态，脑除外）的平均预测。
- 已知：Decima 训练轨道里有 ENCODE 对 ENTEx 组织做的单细胞核数据（10x，poly(A) 捕获），所以 Decima 见过这些组织，但没见过 bulk total-RNA。

## 2. 对应（写死在脚本 `MAP` 里）
- 单位是 Decima 的 organ。ENTEx 组织与 GTEx 组织分别并入 organ，组织内先对供体 / 重复取平均，再在 organ 内对组织取平均。
- 映射：adipose（皮下脂肪、网膜脂肪垫 / GTEx 两种脂肪）；adrenal；breast；gut（横结肠、乙状结肠、胃、Peyer 斑、食管三处 / GTEx 结肠两处、胃、回肠末端、食管三处）；heart（左心室、右心耳 / GTEx 两处）；liver；lung；muscle（腓肠肌 / 骨骼肌）；ovary；pancreas；prostate；skin（两处 / 两处）；spleen；uterus；vasculature（胸主动脉、升主动脉 / 主动脉）。
- 甲状腺、胫神经、睾丸、阴道在 Decima 里无对应 organ，去掉。
- 少于 10 个 organ → 停。

## 3. 指标
- 每个基因：跨 organ 的 Pearson，Decima 预测 vs log1p(TPM)；total 和 poly(A) 各算一个：r_tot、r_pA。
- 基因纳入：biotype（GENCODE v50）为 lncRNA 或 protein_coding；表内；两种测量下都在 ≥ 5 个 organ 中 TPM ≥ 1，且 SD > 0。
- 协变量：log 平均 TPM、SD、log 基因长度。R1 用 total 的；R2 用两种测量各自的。

## 4. 门槛（基因 bootstrap 2000 次，95% CI）
- **R1 复核 DEC2-H1**：Decima 留出基因（test + val）上，r_tot ~ is_lnc + cr(log 平均, 4) + cr(SD, 4) + log 长度。系数 ≤ −0.05 且 CI 上限 < 0 → 过。同一模型用 r_pA 另报。
- **R2 是否与测量方式有关**：留出基因上，(r_tot − r_pA) ~ is_lnc + 两种测量的 cr(log 平均, 3) + cr(SD, 3) + log 长度。系数 ≤ −0.03 且 CI 上限 < 0 → 过。
- train 基因另报，不进门槛。

## 5. 决策
- R1 不过 → TS 上的 lncRNA 差距没有在第二套数据上复现，Decima 盲区不立项，报告作者。
- R1 过、R2 过 → 差距与 total-RNA 测量有关，复现成立，写微调计划交作者批准。
- R1 过、R2 不过 → 差距在外部 bulk 数据上普遍存在（poly(A) 也有），不是 total-RNA 特有。盲区成立但性质改为「外部泛化」，交作者决定。
- 已知局限：organ 层面比较掺入细胞组成差异；GTEx 中位数来自不同供体；V29 与 v50 的基因 ID 用去版本号匹配。
- 结果只读一次；脚本运行前冻结（sha256 见 `PROTOCOL_DEC3.sha256`）。
