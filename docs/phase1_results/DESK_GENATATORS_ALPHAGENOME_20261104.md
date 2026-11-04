# Desk check: GENATATORs appendix and AlphaGenome access (2026-11-04)

> Written by a subagent from the GENATATORs main text + supplement PDFs and AlphaGenome docs / model cards / terms (read only; no terms accepted, nothing downloaded). UNVERIFIED items still need checking. Used by SCREEN_SUBJOURNAL.md.


Read directly: GENATATORs main PDF (https://www.biorxiv.org/content/10.64898/2026.06.17.732686v1.full.pdf, 18 pp), Supplementary PDF (https://www.biorxiv.org/content/biorxiv/early/2026/06/21/2026.06.17.732686/DC1/embed/media-1.pdf, App. A-L) and Fig. 2 image (.../F2.large.jpg). Numbers come from raw text, not summaries.

## Q1 GENATATORs

**(a) lncRNA accuracy.** No "intron chain" or separate transcript-level metric exists. Two exact-match metrics are used:
- Gene-level (Sec. 3.1, App. I.6): a gene counts only if ALL predicted exon intervals equal those of ANY annotated isoform. This is stricter than an intron chain (terminal exon ends must match too).
- Interval-level (Table A15): exact-boundary exon P/R/F1.

Human T2T chr20 (NC_060944.1), lncRNA, ground truth 434 genes (Fig. 2A, Table A16, Table A15):

| Model | Gene-level exact exon structure | Interval exon P / R / F1 |
|---|---|---|
| Caduceus-PS GENATATOR (250 kb) | 111 (25.6%) | 0.62 / 0.52 / 0.56 |
| GENA-large GENATATOR (1 Mb) | 95 (21.9%) | 0.53 / 0.51 / 0.52 |
| AlphaGenome (1 Mb) | 15 (3.5%) | 0.34 / 0.15 / 0.20 |
| Tiberius | 0 | 0.22 / 0.01 / 0.02 |
| AUGUSTUS / Helixer | 0 / 0 | 0.08 / 0.02 / 0.03; 0.19 / 0.02 / 0.03 |
| SegmentNT, NTv3 | 0-4 (read off Fig. 2A) | R 0.05-0.14 |

Other species: Table 2 / Table A16 / Fig. 4C (e.g. A. thaliana 310 of 513; Pan troglodytes 74 of 284, Caduceus-PS).

**There is no locus-level number.** Inputs are the annotated transcript sequences (App. I.2, I.4, I.6), so gene boundaries are given. This is not whole-chromosome locus discovery. Whether flanking context was included is not stated (UNVERIFIED). A reviewer could read their 25.6% exact-structure rate as "GENATATOR recovers lncRNA exon structure, ours recovers none", even though the tasks differ.

**(b) Evaluation set.** T2T-CHM13 chr20 only (Sec. 4 "Input data"). The annotation is NCBI RefSeq on GCF_009914755.1 (App. C), not GENCODE (no GENCODE version is mentioned). Only mRNA and lncRNA genes are kept; training uses the longest transcript per gene, up to 250 kb.
- Training = all human chromosomes except 8, 20, 21. **chr9 is in training; chr20 is held out.**
- Caveat: chr8/20/21 were the validation set for checkpoint selection, and chr20 was also the final metric set ("We did not use a separate test set", App. C).
- Multispecies variants add 38 mammals (Table A6, includes mouse). Orthology leakage is not discussed (UNVERIFIED).
- Weights and data: https://huggingface.co/collections/AIRI-Institute/genatator. Model licenses were not checked (UNVERIFIED).

**(c) AlphaGenome and Tiberius.** Both are compared on chr20 (table above). AlphaGenome used 1 Mb windows, the splice-site track, and acceptor-donor pairing at threshold 0.7 (App. I.2). That crude exon calling may understate AlphaGenome (UNVERIFIED). The checkpoint used, and whether chr20 was in its training data, are UNVERIFIED. Sec. 4.4 says GENATATORs identify more than 3x as many correctly segmented genes as Tiberius or AlphaGenome overall. Tiberius is 0/434 on lncRNA.

**(d) "First" claims.** There is no "first end-to-end lncRNA structure" claim. The only "first" is test-time RC augmentation (Sec. 4.2). Closest wording to avoid contradicting:
- Abstract: "segmentation of previously intractable non-coding transcripts and untranslated regions".
- Conclusions: lncRNA "remain inaccessible to current state-of-the-art tools such as Tiberius".

The ICLR 2025 workshop version (https://airi.net/articles/genatator-de-novo-gene-annotation-with-dna-language-model) already claimed non-coding segmentation. So we must NOT claim "first ab initio lncRNA exon-structure annotation". Locus discovery is the only defensible novelty, with caveats.

## Q2 AlphaGenome

**(a) Weights.** Public but gated: Hugging Face (google/alphagenome-all-folds, fold-0..3) or Kaggle. You must log in and accept the non-commercial Model Terms; approval is automatic (https://github.com/google-deepmind/alphagenome_research README; HF model card). Code is Apache-2.0; weights are license "other" (non-commercial). Nature "Code availability" confirms the weights are in alphagenome_research.
- Size (HF API file listing): all-folds is 734.8 MB; fold-0 is 746.4 MB. All 5 models are about 3.7 GB (folds 1-3 not checked, UNVERIFIED).
- Terms (read via third-party mirror scancode-licensedb, key alphagenome-model-tou-2025): no commercial use, and no training ML models on the Assets (distillation counts as a "derivative"). **Using AlphaGenome outputs as training features or pseudo-labels for our method would be prohibited. Using it only as an evaluation baseline appears fine (UNVERIFIED legally).**

**(b) Local inference.** JAX. The README says "at least an NVIDIA H100" and gives no VRAM. Community benchmark (Tu, Buendia, Kundaje, Mostafavi, https://zenodo.org/records/20272910, not peer-reviewed): 1 Mb needs about 35-41 GB peak (H200 34.6 GB at 198 ms; A40 35.8 GB at 1.33 s). A 24 GB RTX 6000 (Turing; JAX bf16 unsupported) ran only the PyTorch path, reaching 524 kb at 12.5 GB and not 1 Mb. An RTX 4090 was not tested (UNVERIFIED), but 24 GB is below about 35 GB, so 1 Mb is very unlikely; 524 kb or less is plausible. Supported lengths: 16 kb, ~100 kb, ~500 kb, 1 Mb (docs FAQ). Nature: under 1 s on an H100.

**(c) Hosted API.** Needs an API key (alphagenome.google/api, Google sign-in). Free for non-commercial use. Limited to individuals and non-commercial organisations (https://deepmind.google.com/science/alphagenome/terms, modified 2026-09-08). Rates "vary based on demand"; intended for small-to-medium analyses, not more than ~1M predictions (README). Barred: HIPAA data, commercial use, and training variant-effect-type models on outputs. Nothing found against sending public reference intervals (UNVERIFIED: Google APIs ToS not read). Interval queries use hg38; for T2T chr20 we would send raw sequence (docs FAQ says arbitrary sequences of supported lengths are accepted). Chr20 is about 64 windows of 1 Mb.

**(d) Output tracks** (https://google-deepmind-alphagenome.readthedocs-hosted.com/exploring_model_metadata.html, Table 1, human):
- SPLICE_SITES: 4 tracks (donor+, acceptor+, donor-, acceptor-), probabilities, 1 bp, stranded, tissue-agnostic.
- SPLICE_SITE_USAGE: 734 tracks, 282 biosamples, per strand.
- SPLICE_JUNCTIONS: 734 rows, ~367 distinct tracks (strand belongs to the junction); up to 512 donors x 512 acceptors per strand.
- RNA_SEQ: 667 tracks, 285 biosamples, 1 bp, some stranded (+/-/.).
- Biosample ontology is UBERON/CL/CLO/EFO (ENCODE and GTEx sourced).

**(e) lncRNA evaluation.** The authors report none: the Nature page text has no "lncRNA" (checked on nature.com), and the bioRxiv v2 text has one eQTL-feature mention (last ~42k characters unread). The only quantitative evidence found is GENATATORs (15/434 exact; interval F1 0.20). The Nature Methods 2026 cuTAR paper (https://www.nature.com/articles/s41592-026-03071-4, Extended Data Fig. 2) shows splice sites qualitatively only. MIMIC (arXiv 2604.24506) compares junction prediction against AlphaGenome, but the full text was too large to fetch (UNVERIFIED whether it covers lncRNA).
