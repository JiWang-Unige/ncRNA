# Context-conditioned ncRNA modelling after AlphaGenome 2: what is still genuinely open?

## Scope, evidence status, and repository access

This assessment is current to **25 November 2026**. I treated the three project results in your prompt as **author-reported primary project evidence**, not independently reproduced results: the failure of the DNA-only prior to add loci beyond WTC11 ssCAGE+QuantSeq; the D0 localization interaction signal but failed condition transfer; and the B0-lite icSHAPE compartment-change gain with the ±12-nt ablation concern. Those results are strong enough to influence project selection, but because I could not inspect the underlying outputs in the named project folders, I do not treat the numerical values as externally verified.

I could inspect the public-facing `JiWang-Unige/ncRNA` repository only to a limited extent through web indexing. The accessible repository describes an ab initio lncRNA/small-ncRNA annotation project, but searches did **not** expose `08_M1`, `16_new_methods_20261124`, or `17_context_pilot`; a direct container clone was also unavailable in this environment. The most likely explanations are that those folders are on a non-indexed branch, private/unpushed, or simply newer than the indexed snapshot. I therefore could not audit scripts, split definitions, intermediate tables, or raw output for D0/B0-lite. This is an important limitation: recommendations below about your pilots rely on the experimental summaries you supplied rather than code-level verification.

The external landscape has changed materially in 2026. AlphaGenome 2, Decima, HELIX, Otari, IsoPrism and new single-cell localization resources mean that several ideas that would have been genuinely novel in 2024–2025 are now either occupied or need a sharper formulation. AlphaGenome 2 itself is a bioRxiv preprint, DOI **10.64898/2026.11.11.765857**, rather than a peer-reviewed paper as of the cutoff. Its primary paper states that it takes up to 1 Mb of human or mouse DNA plus a tissue/cell/state condition and predicts multimodal regulatory outputs, and that its expression-derived conditioning can generalize to unseen human single-cell datasets without fine-tuning. citeturn809455search7turn403801search6 The public API launched on 11 November 2026 for non-commercial research; its terms explicitly say that direct model outputs and Atlas scores **must not be used to train other machine-learning models**. Live inference is intended for thousands rather than millions of predictions, and expression conditioning accepts up to 256 cells per call. citeturn403801search0turn403801search1turn403801search2turn403801search6

I attempted to inspect the AG2 paper as a PDF, including PDF-page rendering as required for figure/table verification, but the web renderer could not render the roughly 22-MB PDF. I therefore relied on the paper's parsed full text plus the official API documentation and repository. This is sufficient for its model/output definitions and evaluation claims, but not for pixel-level inspection of supplemental figures.

The critical overall conclusion is:

> **Do not build “AlphaGenome for lncRNAs.” That idea has effectively arrived. Do not build another generic context-conditioned isoform or lncRNA-localization classifier either. The strongest remaining opportunity is an assay-aware model of the *context-dependent ncRNA fate interaction*—especially nuclear retention/localization—where the model is explicitly trained to predict the gene×context deviation, not the gene's static localization propensity. A second, riskier opportunity is de-novo context-conditioned transcript architecture for previously unannotated ncRNAs. The current lncRNA structure-change direction should not be a primary project unless it passes a stringent replicate/noise gate first.**

## Context-conditioned transcription and transcript architecture

### What AlphaGenome 2 actually closes

AG2 makes the basic Q1 premise—“DNA plus cellular context rather than DNA alone”—no longer novel. It is trained on human and mouse single-cell atlases plus bulk functional genomics, uses a 1-Mb sequence and condition representation, and predicts RNA-seq/CAGE-like coverage, gene expression, splicing, chromatin accessibility, methylation and related modalities. The official documentation describes stranded RNA outputs, meaning overlapping antisense transcription is at least representable at coverage level. citeturn809455search7turn403801search7 Google reports training on more than 100 million cell profiles across 68 tissues. citeturn809455search0turn403801search5

This is not an isolated result. **Borzoi** already predicts tissue/cell-specific RNA-seq coverage directly from DNA and resolves transcription, splicing and polyadenylation effects; DOI **10.1038/s41588-024-02053-6**, Nature Genetics 2025. citeturn338414search0turn338414search2 **Enformer Celltyping** conditions long-range sequence models on accessibility from a new cell type; DOI **10.1038/s41467-024-54441-5**, PMID **39550354**. citeturn338414search3turn338414search4 **scooby**, DOI **10.1101/2024.09.19.613754**, predicts scRNA-seq coverage and scATAC-seq at single-cell resolution using a cellular embedding. citeturn338414search5turn338414search6 **Decima**, DOI **10.1038/s41592-026-03102-0**, was trained on more than 22 million cells and explicitly predicts cell-type- and disease-state-specific expression from surrounding DNA. citeturn338414search7

Therefore the broad idea “condition sequence models on cell state to rescue context-specific lncRNA prediction” is **already done at the conceptual-method level**. A paper whose main contribution is applying AG2/Decima/scooby to lncRNAs would be a benchmark/resource paper, not a new modelling method.

### The lncRNA-specific benchmark gap is real, but it is not enough by itself

I found no dedicated AG2 publication that systematically stratifies predictive performance into protein-coding versus lncRNA, low-expression versus high-expression, ubiquitous versus cell-type-restricted, antisense versus intergenic, and annotated versus novel transcription. AG2's gene-expression head is annotation-indexed rather than an unrestricted transcript-discovery head; its documentation distinguishes gene-level expression outputs from genomic coverage and splice-junction outputs. citeturn403801search6turn403801search7 This creates a meaningful evaluation gap because a model can have excellent overall gene-expression accuracy while performing badly on sparse, highly specific lncRNAs.

That benchmark is worth doing, but **a benchmark alone is not the method paper you want**. It becomes method-worthy only if it discovers a systematic failure mode and you introduce a correction that generalizes across datasets. Examples would be calibrated low-expression detection, a model of missing-not-at-random assay sensitivity, or a latent transcript-existence model integrating orthogonal evidence.

There is also a serious circularity issue with AG2 expression conditioning. If the input cell-state representation is built from the same scRNA-seq profile whose lncRNA expression you then ask the model to predict, target-gene counts can leak target information into the condition embedding. AG2's documentation describes conditioning from raw UMI gene counts rather than an explicitly target-gene-masked encoder. citeturn403801search6 Therefore a legitimate lncRNA benchmark must include **target-gene-masked conditioning**, preferably also a stronger test using only protein-coding/RBP/context-marker genes to define context. Otherwise high performance on lncRNA expression can be partly tautological.

A decisive pilot is cheap. Select a post-training or otherwise clearly non-overlapping single-cell/total-RNA dataset; create 20–50 biological contexts and 500–2,000 lncRNAs stratified by expression, specificity, antisense/intergenic status and annotation confidence. Compare AG2 reference-condition predictions, expression-conditioned predictions with all genes, and expression-conditioned predictions with each target lncRNA masked. The gate I would use is **context-delta Spearman/Pearson ≥0.25 for lncRNAs after target masking and a statistically significant gain over gene-mean/context-mean baselines in held-out contexts**. If the apparently large advantage collapses when the target gene is removed from the conditioning profile, that is a useful negative result but not a new modelling direction.

### Context-conditioned full transcript architecture is partly occupied

The seemingly strongest Q1 gap—full TSS–splice-chain–TES structure rather than marginal coverage—is narrower than it first appears. In 2026, several models moved directly into context-specific isoforms.

**HELIX**, DOI **10.1038/s43588-026-00988-w**, predicts tissue- and condition-specific splicing and full-transcript isoform usage by integrating pre-mRNA sequence with RBP expression and can transfer to cell-type-specific isoforms. citeturn282178search1 **Otari**, DOI **10.1016/j.xgen.2025.101126**, is trained on long-read transcriptomes from 30 tissues/brain regions and predicts full-length isoform abundance across both protein-coding genes and lncRNAs. citeturn282178search0turn282178search7 **IsoPrism**, DOI **10.1093/bib/bbag621**, PMID **42969855**, uses transcript structure plus expression context to transfer long-read-derived isoform proportions into short-read single-cell/spatial data. citeturn282178search2turn282178search4 **Mach-1**, DOI **10.1101/2024.08.26.609813**, learns long-context pre-mRNA representations and predicts isoform structure and abundance from sequence. citeturn282178search3

So **“predict context-specific isoforms” is already crowded**.

The remaining defensible niche is substantially harder and more specific:

**DNA + cellular context → an unordered set of de-novo transcript objects, including previously unannotated ncRNAs, each with TSS, complete splice chain, TES, abundance and uncertainty.**

That is different from allocating abundance among known transcript candidates. It also addresses a genuine identifiability problem: marginal RNA-seq/CAGE coverage and local splice-junction probabilities do not uniquely determine which TSS belongs to which TES when multiple isoforms overlap. AG2 can predict splice junctions and stranded coverage, but those marginals are not equivalent to full-molecule pairing. This problem requires full-length/long-range supervision—long-read RNA sequencing, paired end evidence, or explicit latent set inference—not merely another coverage decoder. AG2 could be an external comparator or frozen inference prior, but **its current terms prohibit using its outputs to train your model**. citeturn403801search1turn403801search2

On one RTX 4090, I would not begin by training a long-context transformer. A minimal pilot should instead ask whether the problem contains recoverable signal. Freeze a public sequence encoder or use compact local sequence encoders; represent candidate TSS, splice junction and TES nodes as a graph; condition with a low-dimensional cell-state embedding; train a small assignment/set-prediction head against long-read transcript structures. Hold out **both contexts and transcript structures**, and score exact intron-chain recovery, TSS/TES tolerance, transcript-level F1 and novel-locus precision. **Go** only if context improves exact full-transcript recovery by at least ~5 absolute F1 points over sequence+gene-mean and known-isoform baselines, with the gain persisting for previously unseen transcript structures and lncRNAs. **Kill** if nearly all gain is simply choosing among GENCODE transcripts or if de-novo transcript F1 remains close to zero.

### Evidence fusion is more plausible than retraining AG2

Your previous result—that a context-free DNA prior added zero loci beyond ssCAGE+QuantSeq—is exactly the right warning. A context-conditioned prior is worth one **small falsification experiment**, not a new project commitment.

The clean method question is not “does AG2 predict lncRNAs?” but:

**Can a context-conditioned sequence prior improve posterior recovery of transcript structures when direct evidence is sparse, without overriding direct evidence when it is strong?**

Because of AG2's terms, use it only as a frozen external prior or comparator, not as teacher labels for training. citeturn403801search1turn403801search2 A Bayesian/log-odds evidence-fusion pilot can be completely non-neural: endpoint likelihoods from CAGE/3′-end data, splice evidence from RNA-seq, and an independently calibrated AG2 prior from CAGE/RNA-seq/splicing tracks. Evaluate on held-out loci and, crucially, on contexts absent from the evidence used to calibrate the method.

My gate would be **≥10% relative improvement in precision at fixed recall or ≥0.05 absolute transcript-level recall at the same false-discovery rate**, plus a positive gain specifically among weak-evidence loci. If the conditioned prior again contributes essentially zero unique correct structures, stop immediately. Given your previous result, I assign this a higher probability of producing a clean negative answer than of becoming the main paper.

### Mouse and small ncRNAs

Mouse by itself is **not a novelty claim**: AG2 is explicitly trained on human and mouse. citeturn809455search7 A narrower gap exists because AG2's user-supplied expression-conditioning workflow is documented specifically for human scRNA-seq, making arbitrary unseen-state conditioning less symmetric across species. citeturn403801search6 But “port context conditioning to mouse” is an engineering extension, not a compelling ncRNA method unless tied to a biological capability unavailable in human.

Small ncRNAs are a more genuine gap. AG2's main gene-expression outputs are oriented toward annotated gene-expression targets, while many small ncRNAs require dedicated library chemistries, mature-end processing and multimapping-aware quantification. The fact that Tabula Sapiens total-RNA data detect tRNA, snoRNA and snRNA compartmentalization illustrates how different this observation model is from conventional poly(A) scRNA-seq. citeturn117681search3turn117681search4 A context-conditioned small-ncRNA model could be novel, but it is effectively a new project with assay-specific modelling and should not be chosen merely because AG2 has less coverage there.

**Q1 judgment:** run the AG2 lncRNA benchmark as a **diagnostic experiment**, not as the paper. Keep de-novo context-conditioned transcript architecture as a high-risk second-line method. Do one short conditioned-prior fusion test because it directly follows your M1 result. Do not invest in generic context-conditioned expression prediction, generic context-specific isoform usage, or “mouse AG2 for lncRNAs.”

## Context-dependent localisation and RNA fate

### The data landscape is much better than D0 suggested

The original ENCODE/lncATLAS setting is genuinely small. lncATLAS is based on nuclear and cytoplasmic poly(A)+ RNA-seq from **15 human cell lines**, usually with two independent biological replicates per compartment; DOI **10.1261/rna.060814.117**, PMID **28386015**. citeturn844982search0turn844982search1 It is valuable for clean compartment ratios but cannot support a rich cell-state representation by itself.

Other direct localization datasets are deep but similarly narrow in biological context:

| Resource | Direct measurement | Contexts / replication | ncRNA relevance | Assessment |
|---|---|---|---|---|
| **lncATLAS / ENCODE**; DOI 10.1261/rna.060814.117, PMID 28386015 | Nuclear vs cytoplasmic poly(A)+ fractionation | 15 human cell lines; generally 2 biological replicates/fraction | Explicit lncRNA resource | Clean labels, too few contexts for modern condition-transfer learning. citeturn844982search0turn844982search1 |
| **CeFra-seq**, ENCODE **ENCSR885DVH**, GEO **GSE88622**; DOI 10.1016/j.ymeth.2017.05.017 | Nuclear, cytosolic, endomembrane, insoluble, extracellular fractions | K562; isogenic experiment | Transcriptome-wide, including ncRNAs if expressed | Fine localization resolution, effectively one cell context. citeturn844982search2turn844982search3 |
| **APEX-seq**, GEO **GSE116008**, BioProject **PRJNA476786**; Cell 2019 | Proximity-labelled RNA in 9 subcellular locales | HEK293T; 2 labelled + 2 control biological replicates per construct; several drug/time perturbations | Diverse RNA classes and isoforms | Excellent spatial labels and perturbations, but one engineered cellular background. citeturn844982search4turn844982search5 |
| **Halo-seq**, GEO **GSE172281**, BioProject **PRJNA722745**; DOI 10.1093/nar/gkab1185, PMID 34875090 | Chromatin-proximal, nucleolar, cytoplasmic proximity labelling | One human system plus perturbation including leptomycin B | Transcriptome-wide | Useful direct calibration for export/retention, not broad context transfer. citeturn844982search6turn844982search7turn844982search8 |
| **MERFISH spatial transcriptome**, DOI 10.1073/pnas.1912459116, PMID 31501331 | Direct molecule coordinates relative to nucleus/ER | Cultured cells, ~10,050 targeted genes; 9,050 robustly detected | Includes lincRNAs and nuclear-enriched transcripts | Strong orthogonal validation, not a many-context training set. citeturn761738search3 |
| **scRLP**, DOI 10.64898/2026.02.04.703776, PMID 42051298; reported BioProject-like accession **PRJNA1312488** | Compartment-targeted RNA editing at single-cell level | Cytoplasm/ER/plasma-membrane localization in individual cells | Full-length poly(A) RNA, thousands of RNAs/cell | Very important conceptually: cell-to-cell spatial heterogeneity exists. Still preprint and engineered-system data. The manuscript oddly calls PRJNA1312488 a GEO accession; treat the repository mapping as needing manual verification. citeturn761738search0turn761738search1 |

The genuinely important 2026 development is **Tabula Sapiens total-RNA**. The preprint, DOI **10.64898/2026.03.23.713770**, PMID **41929159**, generated single-cell and single-nucleus total-RNA data across **22 tissues**, with **27 whole-cell and 21 nuclear samples**, and **11 tissues measured by both modalities**; 105,471 cells and 121,398 nuclei passed QC and 79 cell types were annotated. The authors explicitly compare matched whole-cell and nuclear measurements to infer cell-type-dependent nuclear compartmentalization of ncRNAs. Gene counts and metadata are available at figshare DOI **10.6084/m9.figshare.31827436**. citeturn117681search3turn117681search4

This is much closer to what D needed: dozens of biologically distinct tissue×cell-type contexts, total-RNA capture of noncoding classes, and a paired whole-cell/nucleus design.

But it is **not direct fractionation**. A whole-cell versus nucleus ratio conflates true localization with modality-specific capture, cell dissociation, nuclear isolation, transcript abundance, intronic/pre-mRNA content, RNA integrity and potentially different cell-composition sampling. The dataset's own design is strongest as a broad *weak-label* source, not a gold-standard localization assay. citeturn117681search4 Moreover, the initial total-RNA expansion is essentially one donor, so “79 cell types” is not equivalent to 79 independent biological replicates. Tabula Sapiens 2.0, DOI **10.1016/j.cell.2026.08.010**, PMID **42685689**, supplies 1.1 million cells from 24 donors and 28 tissues for context representations and donor robustness, but it is conventional single-cell transcriptomics rather than paired localization truth. citeturn761738search2turn761738search5

### A naive conditioned localization model is already done

A major novelty warning is **BiGM-lncLoc**, DOI **10.1007/s12539-024-00679-y**, PMID **39724386**, published in 2025. It explicitly formulates **cell-specific lncRNA subcellular-localization prediction** as a multi-graph meta-learning problem, combining nucleotide-sequence/localization data with cell-line expression. It reports 82.4–94.7% accuracy after removing overlapping samples between cell lines. citeturn117681search0

Therefore the proposed idea

> sequence + cellular expression embedding → lncRNA nuclear/cytoplasmic localization

is **not novel enough** by itself.

Most other recent predictors remain static localization classifiers rather than genuine context-transfer models, so there is still room beyond BiGM-lncLoc; but the method needs to target the biological quantity your D0 established: the **gene × context interaction**.

### The method I would pursue: assay-aware latent context-specific RNA localization

The stronger formulation is:

**Learn a latent, continuous localization/retention propensity \(L_{g,c}\) for RNA \(g\) in biological context \(c\), while explicitly modelling each assay as a noisy observation operator of that latent state. Predict only the context-specific deviation \(L_{g,c}-\bar L_g\) from RNA sequence/architecture plus an independently defined cellular-state embedding.**

This differs from BiGM-lncLoc in four substantive ways.

First, the target is a **continuous context interaction**, not a categorical compartment label. Your D0 already says most easy predictive signal is likely the gene mean; a credible method must not obtain good scores merely by learning that MALAT1 is nuclear and GAPDH mRNA is cytoplasmic.

Second, it can integrate heterogeneous evidence without pretending it is equivalent: ENCODE fractionation observes physical abundance ratios; APEX/Halo/scRLP observe proximity/editing processes; MERFISH observes molecule coordinates; scRNA-versus-snRNA is a modality-dependent proxy. Each gets a separate likelihood/calibration layer.

Third, context embeddings can come from large single-cell atlases such as Tabula Sapiens 2.0, while localization supervision comes from the much smaller direct/proxy datasets. This is more data-efficient than demanding the same assay in hundreds of contexts. citeturn761738search2turn761738search5

Fourth, target leakage can be controlled. The context encoder should exclude the target RNA itself and ideally all ncRNA abundance features, using protein-coding trans regulators, RBP expression and/or low-dimensional cell-state factors. That directly tests whether trans state actually explains the interaction rather than handing the answer to the model.

This is also more biologically defensible than mixing “localization”, “stability”, “unspliced fraction” and “nuclear retention” as though they were interchangeable. Spliced/unspliced ratios are influenced by transcription and kinetics; metabolic labeling measures synthesis/decay; whole-cell/nucleus ratios are localization proxies; APEX/Halo are proximity assays. These should be auxiliary observations of related latent RNA-fate variables, not pooled labels.

### A decisive D pilot

Use the matched Tabula Sapiens total-RNA dataset first because it gives breadth without requiring terabytes of raw processing. citeturn117681search4 Restrict to tissue×cell-type combinations represented in both cell and nucleus with adequate pseudobulk depth. For each RNA \(g\) and context \(c\), estimate an error-aware nuclear-enrichment statistic with bootstrap uncertainty. Remove the gene mean to obtain the interaction target.

Construct context features from protein-coding expression only. Start with a deliberately small model: ridge regression, factorized bilinear model, and a two-tower MLP where one tower encodes RNA sequence/architecture and the other encodes context. There is no reason to start with a foundation model.

The split must be **leave-tissue-out and, where possible, leave-donor/assay-out**, not random gene×context entries. Also include three controls: gene-mean-only; nearest-context transfer; and shuffled-context labels within tissue.

My go/no-go gates are:

| Test | Go threshold | Kill signal |
|---|---:|---:|
| Held-out-context prediction of interaction residual | median Pearson/Spearman **≥0.20–0.25** | <0.10 |
| Variance explained beyond gene mean | **interaction \(R^2>0.05\)** with bootstrap CI >0 | ≤0 |
| Context contribution | condition-aware model beats sequence/gene-only model with gene-block CI excluding 0 | gain collapses to zero |
| Donor/tissue transfer | effect retains sign and meaningful magnitude out of tissue | only within-tissue performance |
| Direct-assay external validation | significant concordance with ENCODE/APEX/Halo/MERFISH subset | no cross-assay concordance |
| Negative control | shuffled/null context performance ≈0 | material predictive signal on null |

The **key success criterion is not a high global correlation**. It is recovering the residual context-dependent movement after subtracting each RNA's mean localization.

On your hardware this is highly feasible. The first pilot is matrix-level statistics plus compact neural models and should fit comfortably within 24 GB GPU memory; most work is CPU preprocessing. I would budget tens of GB rather than hundreds of GB if working from processed Tabula Sapiens matrices. That storage estimate is my engineering estimate; I did not verify the full current figshare package byte size.

If this pilot reaches interaction \(r\sim0.2–0.3\) out-of-tissue, D becomes much stronger than your original 13-cell-line experiment. If it repeats D0 at \(r<0.1\), **kill D rather than adding model complexity**. At that point the likely explanation is that context-specific localization depends on variables not recoverable from steady-state transcriptomic context alone.

## Context-dependent RNA structure

### The public-data situation is broader than GSE117840, but not as favorable for lncRNAs as it first appears

The best broad catalogue is **RASP v2.0**, DOI **10.1093/nar/gkae1117**, PMID **39546630**. It contains 438 structure-probing datasets from 85 articles: 216 transcriptome-wide probing datasets, 141 target-specific sets and 81 RNA–RNA interaction datasets across 24 species and 39 methods. Importantly, RASP imputes missing structural scores for 59 low-coverage datasets. Those imputed values are useful for exploration but **must not be used as ground truth** in your proposed method because they would introduce model-generated labels into model evaluation. citeturn310658search0turn310658search1turn310658search3 RMDB is complementary, curated mainly at the individual construct/chemical-mapping level; DOI **10.1093/bioinformatics/bts554**, PMID **22976082**. citeturn310658search4turn310658search5turn310658search6

I separate candidate datasets below into **qualifying**, **partial**, and **not suitable for the specific lncRNA-delta task**. “Coverage unknown” means I did not find a primary-source count and will not infer one.

| Dataset | Contexts / method | Biological replication | lncRNA evidence | Raw-data / status | Verdict |
|---|---|---|---|---|---|
| **Sun et al. 2019 icSHAPE**, GEO **GSE117840**, BioProject **PRJNA483400**, DOI **10.1038/s41594-019-0200-7**, PMID **30886404** | HEK293 and mES; chromatin, nucleoplasm, cytoplasm; in vivo/in vitro NAI-N3 | Study reports replicated probing; public processed tracks merge replicates in ways that complicate independent-delta construction | Your project finds only 16–39 jointly covered lncRNA genes in relevant contrasts | BioProject-scale raw data are very large; project metadata report hundreds of gigabases. Your ~94-GB estimate appears to refer to a relevant subset rather than the entire deposited project. | **Qualifying biologically, poor lncRNA sample size and expensive.** |
| **smartSHAPE**, GEO **GSE155961**, SRA **SRP276858**, BioProject **PRJNA656214**; DOI **10.1016/j.fmre.2021.12.007**, PMID **38933905** | Two mouse colonic macrophage subtypes under inflammatory setting; NAI-N3; HEK293 input benchmarking | HEK293 benchmark contains replicated libraries; macrophage replicate structure is not sufficiently explicit in the GEO summary for me to certify independent biological pairs without inspecting run-level metadata | lncRNA count not reported | 36 GEO samples; raw SRA available; processed macrophage shape files 23–34 MB each | **Promising partial candidate; replication and lncRNA coverage need run-level audit before use.** citeturn117681search1turn117681search2 |
| **Human neurogenesis icSHAPE**, GEO **GSE236224**, BioProject **PRJNA989350**; DOI **10.1016/j.molcel.2021.09.027**, PMID **34655516** | Paper analyzes hESC→NPC→early neuron→neuron structural remodeling | Public GEO currently exposes clear replicated hESC/NPC control/DMS samples; the deposited sample list does not transparently provide all four reported stages | lncRNA coverage not reported; study emphasis is primarily regulatory mRNA structure | 9 GEO samples, raw SRA available | **Useful replicated biological transition, but only partial public representation and unknown lncRNA yield.** citeturn257340search3turn257340search4 |
| **PrismNet structuromes**, GEO **GSE145805**, BioProject **PRJNA608297**, DOI **10.1038/s41422-021-00476-y** | HEK293, HeLa, K562, HepG2, H9, mES; icSHAPE | Public design largely gives one DMSO and one NAI sample per cell line rather than independent biological replicate pairs | lncRNA count not reported | Very deep libraries, raw SRA available | **Excellent context breadth, fails your replicate criterion.** citeturn257340search2 |
| **PARIS**, GEO **GSE74353**, BioProject **PRJNA299809**; DOI **10.1016/j.cell.2016.04.028**, PMID **27180905** | Duplex maps in HEK293T, HeLa, mES; psoralen crosslink/ligation | HEK293T and mES have multiple libraries; associated icSHAPE controls include explicit replicates | Strong lncRNA exemplar: XIST; broad lncRNA count not reported | 14 GEO samples, raw SRA available | **Qualifies for duplex/interactome modelling, not for nucleotide-reactivity delta.** citeturn257340search1turn257340search6turn257340search7 |
| **SPLASH**, DOI **10.1016/j.molcel.2016.04.028**, PMID **27184079** | Human LCL, hESC, retinoic-acid differentiation, HeLa; psoralen duplex mapping | Published resource includes replicated libraries for several states; hESC→RA gives a true biological transition | Long-range interactions are enriched in ncRNAs, but transcriptome-wide lncRNA count is not stated | Public raw/processed interaction data available through study resources | **One of the better context-change datasets if you change the target from reactivity to RNA duplex architecture.** citeturn257340search0turn257340search5 |
| **RIC-seq** | RBP-mediated RNA proximity/interaction, not chemical accessibility | Dataset-dependent | may capture ncRNAs | varies | **Not interchangeable with SHAPE/DMS. Treat as interaction topology, not “structure reactivity.”** |
| **nanoSHAPE / nanopore-based probing / SMS candidates** | direct-read chemical or structural signals | study-dependent | I did not verify a multi-context, replicated, lncRNA-rich public dataset satisfying all three criteria | unverified | **No qualifying dataset I can responsibly recommend from this search.** |

Two methodological distinctions matter. PARIS/SPLASH/RIC-seq measure **duplexes or proximity contacts**, whereas icSHAPE/SHAPE-MaP/DMS-MaPseq measure **chemical accessibility/reactivity**. They are not interchangeable labels. A model predicting context-specific duplex formation could actually be more biologically interpretable than B0's local reactivity delta because long-range pair formation is a discrete molecular event, but those assays are sparse and capture-biased. citeturn257340search0turn257340search6

Second, high coverage across protein-coding transcripts does not imply useful lncRNA coverage. Your 16–39-lncRNA observation in GSE117840 is exactly the failure mode expected for low-abundance transcripts in transcriptome-wide structure probing. RASP's decision to impute low-coverage measurements underscores how pervasive missingness is, but those imputations cannot solve the independent-label problem for a methods paper. citeturn310658search0

### B0's ±12-nt result is not yet evidence for biological context remodelling

The local-sequence ablation is a major red flag rather than a minor caveat. Chemical probing and reverse-transcription pipelines have local base- and sequence-dependent efficiencies. A predictor can therefore improve \(\Delta r\) by learning assay-specific expected measurement error around each nucleotide, especially when the reference reactivity itself is noisy. In that scenario, the model is performing **heteroscedastic denoising**, not learning context-dependent folding.

The current comparison “reference reactivity alone” versus “reference reactivity + local sequence” is therefore too weak. The correct null is a **replicate-aware chemistry/noise shrinkage model** that already knows nucleotide identity, local k-mer, read depth, RT-stop/mutation propensity and reference uncertainty.

The decisive experiment must remove shared noise between feature and label. Suppose reference context has replicates \(A_1,A_2\) and target has \(B_1,B_2\). One cross-fit can use \(A_1\) as the reference-reactivity feature and define the label from independently measured \(B_2-A_2\). Rotate all replicate assignments. Never use the same reference measurement simultaneously as both an input and part of the delta label.

Then run **within-condition pseudo-deltas**, such as \(A_2-A_1\), through exactly the same model. If local sequence predicts those pseudo-changes, the sequence gain is measurement bias. Coverage and expression distributions should be matched or downsampled between real and pseudo contrasts.

I would require all of the following before committing to B:

| Replicate/noise test | Required result |
|---|---:|
| Cross-fitted replicate delta reliability | ICC/Spearman roughly **≥0.5** for the evaluated nucleotide/window set |
| Local-sequence gain on same-context pseudo-deltas | **≤0.01 Pearson** beyond a chemistry-aware shrinkage baseline |
| Gain on true context contrast | **≥0.05 Pearson** beyond that baseline and at least ~3× the pseudo-delta gain |
| ncRNA support | **≥100 lncRNA genes**, preferably ≥100 independently covered informative nt/gene, with ≥2 biological replicates/context |
| Gene-level uncertainty | gene-block bootstrap CI excludes zero; position-block sensitivity also required |
| External validation | same directional structural change for a nontrivial subset in an independent study or different chemistry |
| Positive controls | known remodeled elements recovered; negative housekeeping controls remain stable |

Cross-chemistry validation should not compare raw SHAPE and DMS numbers one-to-one: SHAPE-type reagents and DMS report different chemistry. Compare window-level ranks, accessibility/pairedness direction, or validated structural transitions instead.

If GSE155961 or GSE236224 cannot supply ≥100 well-covered lncRNAs with true biological replicate pairs, **I would stop trying to make B a lncRNA method**. You could instead pivot B toward all transcripts, or toward duplex topology, but then it ceases to be the ncRNA-specific project initially envisioned.

**Q3 judgment:** with currently verifiable public data, lncRNA context-induced reactivity prediction is **not yet a viable primary direction**. It is worth one replicate audit because a clean pass would be interesting. But it is currently the direction most likely to consume hundreds of gigabytes and weeks of preprocessing only to reveal that the +0.11 is chemistry/noise shrinkage.

## Ranked project decision

The ranges below are **subjective conditional estimates**, not measured probabilities. They assume competent execution, a genuinely held-out validation dataset, no hidden leakage, and a target journal whose scope fits computational RNA/genomics methods. They estimate the chance that the direction could support a credible JCR-Q1 methods paper **if the stated pilot gate passes**; before the pilot, the effective probability is lower.

| Rank | Direction | Novelty against 2025–2026 closest work | 4090 feasibility | Decisive pilot / kill criterion | Conditional Q1-paper potential |
|---|---|---|---|---|---:|
| **1** | **Assay-aware, context-conditioned ncRNA localization / nuclear-retention interaction model** | Generic context-specific localization is already occupied by BiGM-lncLoc, DOI **10.1007/s12539-024-00679-y**. Novelty comes from continuous gene×context interaction, broad tissue/cell states, heterogeneous-assay observation model and held-out-context transfer. citeturn117681search0 | **High**. Processed matrices + compact two-tower/hierarchical models fit easily. | Tabula Sapiens total-RNA weak labels; protein-coding-only context embedding; leave-tissue-out. **Go:** interaction r≥0.20–0.25 and R²>0.05 with external direct-assay concordance. **Kill:** r<0.10 / R²≤0. | **45–65% if gate passes**; roughly **25–40% before pilot**. |
| **2** | **De-novo context-conditioned ncRNA transcript-object model** | HELIX, Otari and IsoPrism already cover context-dependent isoform usage; Mach-1 covers transcript architecture. Novelty survives only if the model discovers **unannotated full transcript objects** rather than choosing among known isoforms. DOI **10.1038/s43588-026-00988-w**, **10.1016/j.xgen.2025.101126**, **10.1093/bib/bbag621**, **10.1101/2024.08.26.609813**. citeturn282178search0turn282178search1turn282178search2turn282178search3 | **Moderate** for graph/set decoder with frozen embeddings; low for training a foundation model. | Long-read-supervised graph pilot. **Go:** ≥5 absolute F1 improvement on unseen lncRNA transcript structures and meaningful novel-locus precision. **Kill:** gain exists only for catalogued isoform allocation. | **35–55% if strong de-novo gate passes**; ~20–35% before pilot. |
| **3** | **Context-conditioned prior + cheap endpoint evidence** | Not intrinsically novel; novelty would be calibrated probabilistic fusion for weak-evidence/novel ncRNA structures. AG2 cannot legally be used as teacher output for ML training under current API terms. citeturn403801search1turn403801search2 | **Very high**; no major training required. | Re-run the M1 logic with frozen conditioned predictions. **Go:** ≥10% relative precision improvement at fixed recall or ≥0.05 absolute recall at fixed FDR, specifically in weak-evidence loci. **Kill:** zero/near-zero unique correct loci, as before. | **25–45% if surprisingly positive**; because of your previous negative result, prior chance is modest. |
| **4** | **AG2-class lncRNA failure benchmark + correction** | No dedicated lncRNA-focused AG2 benchmark found; but benchmarking alone is not a new model. AG2 DOI **10.64898/2026.11.11.765857**. citeturn809455search7 | **Very high**, API-scale experiment. | Target-masked conditioning; stratify low expression, specificity, antisense, lncRNA. **Go as method only if** a reproducible failure mode supports a correction with external generalization. | **20–40% with a real correction; 10–20% as benchmark only**. |
| **5** | **Context-dependent small-ncRNA model** | More genuinely open because small-RNA protocols/processing differ from AG2-style gene expression. Tabula Sapiens total-RNA demonstrates strong cell-type specificity of ncRNA classes. citeturn117681search3turn117681search4 | **Moderate**; data harmonization dominates compute. | Choose one class first—tRNA, snoRNA or miRNA—not “all small ncRNAs.” Need ≥20–30 contexts and an orthogonal held-out assay. | **30–50%**, but effectively a new project and slower than D. |
| **6** | **Context-induced lncRNA chemical-structure change** | Biological question is still interesting; model space less crowded than transcription. Data quality, not model novelty, is the bottleneck. RASP v2.0 confirms many datasets exist but does not solve lncRNA replication/coverage. citeturn310658search0 | **Low–moderate** because raw-read processing/storage dominates; multi-study audit can reach several hundred GB. | Cross-fitted replicate experiment above. **Kill immediately** if null pseudo-delta gain >0.01, true gain <0.05 above chemistry baseline, or <100 lncRNAs survive. | **15–30% if the gate passes**; currently I would assign **<15–20% pre-gate**. |
| **7** | **Mouse-only context-conditioned lncRNA prediction** | AG2 already supports human and mouse sequence modelling; species change alone is not innovation. citeturn809455search7 | High | Only pursue if tied to a capability unavailable in human, not as a port. | **<15%** as a standalone method. |

The option I would **not** rank separately is a generic “RNA-fate foundation model.” That phrase is attractive but currently too broad. Localization, nuclear retention, stability, translation and splicing have different observation processes, and large models such as HELIX already occupy pieces of that space. citeturn282178search1 A multi-task RNA-fate model becomes compelling only after one fate dimension—most plausibly localization—shows transferable gene×context signal. Starting with all fates would add missing-data and batch-identifiability problems before proving the central hypothesis.

The highest-value sequence of work is therefore:

1. **Run the D re-pilot first**, using Tabula Sapiens matched sc/sn total RNA as broad weak localization supervision. Do not reproduce D0 with more cell types; change the statistical target to the interaction residual and the modelling framework to an assay-aware latent localization variable. The input context must exclude the target ncRNA. Tabula Sapiens gives 79 annotated cell types across 22 tissues in the total-RNA resource and is sufficiently different from the 13-cell-line ENCODE matrix to make this a real falsification test rather than incremental tweaking. citeturn117681search4

2. In parallel, run a **small AG2 audit**, not a development project: 20–50 contexts, a few hundred representative lncRNAs, target-masked versus unmasked conditioning, and an explicit coding-vs-lncRNA performance decomposition. This should require thousands—not millions—of predictions, which matches the API's intended scale. citeturn403801search2turn403801search4 The purpose is to determine whether there is a concrete ncRNA failure mode worth exploiting.

3. Run the **conditioned-prior fusion test** on your existing annotation framework as a one- or two-day statistical experiment once AG2 predictions are available. Because your context-free prior already failed, the burden of proof should be high. Do not redesign the pipeline unless the prior contributes unique correct loci.

4. Audit **GSE155961 and GSE236224 at the run level** before downloading a large collection of structure data. Count GENCODE lncRNAs satisfying a pre-specified per-replicate coverage criterion. If neither reaches roughly 100 replicated lncRNAs, stop B immediately. smartSHAPE has biologically interesting macrophage contexts but the accessible GEO summary does not let me certify the relevant biological-replicate layout; GSE236224 has explicit replicated hESC/NPC samples but only nine deposited samples despite the paper discussing four differentiation stages. citeturn117681search1turn257340search3

5. Only if D fails and the transcript-architecture pilot shows unusually strong signal should you commit to **de-novo context-conditioned transcript-object prediction**. The novelty claim must explicitly exclude known-transcript selection and should be evaluated on unannotated long-read-supported lncRNAs. Otherwise HELIX/Otari/IsoPrism reviewers will correctly say the central idea is already present. citeturn282178search0turn282178search1turn282178search2

The most consequential change from the earlier project framing is that **the field no longer needs another model showing that sequence plus cell state predicts transcription or isoform usage**. AG2, Decima, scooby, HELIX and Otari have moved that frontier. citeturn809455search7turn338414search5turn338414search7turn282178search0turn282178search1 What remains underdeveloped is modelling ncRNA properties that arise **after transcription** and separating their stable RNA-intrinsic component from a transferable cellular-context interaction.

On the evidence available today, **context-dependent localization/nuclear retention is the best bet** because your own D0 has already established that the interaction variance is biologically nontrivial, while the new 2026 total-RNA and single-cell localization resources materially improve the data regime. The correct response to D0 is not “try a bigger neural net”; it is “increase genuine context diversity, remove the gene main effect, define context without target leakage, and model assay-specific measurement.” If that still fails out-of-tissue, the result is informative enough to kill the direction with confidence.

By contrast, **B currently has the wrong risk profile** for the main thesis project: its pilot signal is encouraging numerically but exactly aligned with the local chemistry/noise confound that a sequence window would exploit, and public lncRNA coverage remains thin despite a large apparent universe of structure-probing studies. RASP v2.0 makes discovery easier, but its 438 datasets should not be confused with 438 replicated lncRNA-rich biological contrasts. citeturn310658search0turn310658search3 Until the cross-replicate pseudo-delta test says otherwise, I would treat B0's +0.11 as **unresolved measurement signal rather than evidence that sequence predicts context-induced RNA refolding**.

Finally, the most defensible potential paper-level contribution is not “predict where an lncRNA is” but something closer to:

**A transferable, assay-aware model of context-dependent RNA localization that learns which RNA sequence/architecture features interact with cellular state, separates constitutive localization from state-dependent relocalization, and transfers to unseen tissues/cell states and measurement technologies.**

That is both narrower and more novel than a generic RNA-fate foundation model, and it directly attacks the unexplained interaction component exposed by D0 rather than optimizing a task already solved mostly by RNA identity.