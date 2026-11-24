# LITERATURE_P0a: positioning of "evidence-constrained transcript identification with ambiguity certificates"

Date of review: 2026-11-24. Reading only; nothing was posted or sent.

## 0. Method, sources, and verification status

- Full text read (open access via the PMC article pages or Europe PMC): Hiller 2009 (PMC3167695), Ferrer-Bonsoms 2022 (PMC8896638), Khan 2022 (PMC9807076), Dias SafeMFD 2023 (PMC10628435), Scallop2 (PMC9879047), Bookend (PMC9245221), IsoTools (PMC10287928), Dias MFD-ILP 2022 (PMC9700332), Zheng 2022 (PMC8892959), Chen/Shao catfish-LP preprint (PMC12713669, partial grep only).
- Abstract/metadata only (PubMed record via the PubMed tool): Williams 2023 (PMID 35104222), Caceres 2022 (PMID 34847041), Scallop (PMID 29131147), TIF-seq (PMID 23615609, protocol PMID 24967623).
- Web pages read: flowpaths documentation (algbio.github.io/flowpaths), the GitHub READMEs of algbio/mfd-safety and algbio/flowpaths, arXiv:2511.19153 abstract, bioRxiv Shannon (10.1101/039230) and bioRxiv saturating subflows (10.64898/2025.12.11.693570) pages.
- Europe PMC REST was intermittently unavailable (HTTP 500/503) during the session; PMC article pages were used instead. No claim below comes from memory, but items marked UNVERIFIED were not confirmed against the primary text.
- Section references are given as the section names in the papers where I could identify them; the PMC text extraction lost some section numbering, so some are by heading description.

### Identifier verification (all confirmed in PubMed metadata unless stated)

| # | Paper | PMID | DOI | Verdict |
|---|---|---|---|---|
| 1 | Hiller, Jiang, Xu, Wong. Identifiability of isoform deconvolution from junction arrays and RNA-Seq. Bioinformatics 2009;25(23):3056-9 | 19762346 | 10.1093/bioinformatics/btp544 | confirmed |
| 2 | Ferrer-Bonsoms, Morales, Afshar, Wong, Rubio. On the identifiability of the isoform deconvolution problem... Bioinformatics 2022;38(6):1491-1496 | 34978563 | 10.1093/bioinformatics/btab873 | confirmed |
| 3a | Khan, Kortelainen, Caceres, Williams, Tomescu. Improving RNA Assembly via Safety and Completeness in Flow Decompositions. J Comput Biol 2022;29(12):1270-1287 (extended version of the RECOMB 2022 paper "Safety and Completeness in Flow Decompositions for RNA Assembly", arXiv 2201.10372, found by search, not read) | 36288562 | 10.1089/cmb.2022.0261 | confirmed |
| 3b | Caceres, Mumey, Husic, Rizzi, Cairo, Sahlin, Tomescu. Safety in Multi-Assembly via Paths Appearing in All Path Covers of a DAG. IEEE/ACM TCBB 2022;19(6):3673-3684 (the "SafePC" paper) | 34847041 | 10.1109/TCBB.2021.3131203 | confirmed |
| 4 | Dias, Caceres, Williams, Mumey, Tomescu. A safety framework for flow decomposition problems via integer linear programming ("SafeMFD"). Bioinformatics 2023;39(11):btad640 | 37862229 | 10.1093/bioinformatics/btad640 | confirmed |
| 5a | Zhang, Shi, Shao. Accurate assembly of multi-end RNA-seq data with Scallop2. Nat Comput Sci 2022;2(3):148-152 | 36713932 | 10.1038/s43588-022-00216-1 | confirmed |
| 5b | Shao, Kingsford. Accurate assembly of transcripts through phase-preserving graph decomposition (Scallop). Nat Biotechnol 2017;35(12):1167-1169 | 29131147 | 10.1038/nbt.4020 | confirmed (abstract only) |
| 6 | Schon, Lutzmayer, Hofmann, Nodine. Bookend: precise transcript reconstruction with end-guided assembly. Genome Biol 2022;23(1):143 | 35768836 | 10.1186/s13059-022-02700-3 | confirmed |
| 7 | Lienhard et al. IsoTools: a flexible workflow for long-read transcriptome sequencing analysis. Bioinformatics 2023;39(6):btad364 | 37267159 | 10.1093/bioinformatics/btad364 | confirmed |

---

## 1. Requested papers

### 1.1 Hiller et al. 2009 (PMID 19762346)

- Problem: given a known set of K isoforms, when are isoform abundances uniquely determined from exon/junction array data or RNA-seq feature counts (Abstract; Results).
- Observation model: linear model y_ij = sum_k omega_ik * delta_kj * phi_j, where j indexes "unique probe classes" (k-mers / read classes with identical isoform membership), delta_kj is 0/1 membership, phi_j the (probe/sampling) rate. For RNA-seq the sampling rates are treated as known (Methods, RNA-seq paragraph: "Model (1) translates ... to RNA-Seq", reads in class j Poisson-like).
- Solution space / what is computed: a theorem giving sufficient (and, for the unknown-affinity case, necessary) conditions for identifiability: with known rates, one needs K independent probe classes (Delta_1 invertible); with unknown affinities, a bipartite graph G(D) built from D = Delta_1^-1 Delta_2 must be connected (Theorem, Methods). Identifiability is a yes/no property of a GIVEN candidate isoform set (the RefSeq isoforms of a gene), not de novo. R script for checking supplied (http://biogibbs.stanford.edu/~djhiller/nonid_test.R, as stated in the text; link not tested).
- Ends: not modeled as observations. Candidate isoforms already contain their ends; Fig. 2c states alternative TSS / alternative polyadenylation are not identifiable on arrays but are identifiable in RNA-seq if the sampling rates are known.
- Alternative solutions / counterexamples: Fig. 1 gives a hand-built example showing that two different abundance vectors yield the same data; no systematic counterexample output.
- Evidence-increment guidance: qualitative only. Discussion: 3% of genes remain non-identifiable in RNA-seq; in 90% of those, two alternative events are separated by a constitutive exon, so "junction arrays and single end RNA-Seq ... only able to assess local properties"; they suggest paired-end data may help (speculation, not tested) or quantifying events rather than isoforms.
- Parsimony: none; the isoform set is assumed given (they advise small candidate sets because identifiability collapses as the number of transcripts grows, Fig. 3 and Discussion).
- Complexity/scalability: linear algebra / graph connectivity on a per-gene matrix; trivial.
- Software: R script (see above).
- Data: 2256 alternatively spliced human RefSeq genes (hg18), 1118 mouse; Human Exon array, HJAY array, simulated array, RNA-seq. Result: RNA-seq 97% identifiable, exon array 26%, HJAY 69%.

### 1.2 Ferrer-Bonsoms et al. 2022 (PMID 34978563)

- Problem: for a gene with a given transcript set, is the structure and concentration of its transcripts uniquely recoverable from paired-end reads of a given read length and fragment-length distribution? Used to choose the fragment length (Abstract).
- Observation model (Methods, "Model" paragraphs): reads/read pairs are mapped to "bins" (sorted subexon sequences that a single read can span; the bin graph depends on read length) and to pairs of bins. Counts are proportional to the summed concentration of transcripts containing the bin (pair), times an effective length m that depends on read length, fragment-length distribution and distance between bins in each candidate transcript. Matrix form y = M theta, M of size (#bin pairs) x (#paths).
- Solution space / computation: a gene is identifiable iff rank(M) equals the number of start-to-end paths of the bin graph (Results/Methods, whole-genome analysis). So the candidate set is ALL start-to-end paths of the splicing graph, not an annotated list; identification is of structure and concentration jointly, and the test is a rank test (a yes/no per gene). Two libraries are handled by stacking M1 and alpha*M2 (Methods); rank can only increase.
- Ends: explicitly assumed known beforehand ("The transcript start and end sites are known beforehand. There can be several start and end sites for each gene", Hypotheses list); the intro suggests CAGE / polyA data for this. Every start-to-end combination is a candidate path, so this is the closest formal treatment of the TSS-TES pairing issue, but at gene level and with ends given, not as separate observations that may or may not be paired.
- Alternative solutions / counterexamples: none output; rank deficiency is reported, not an alternative decomposition. A supplement figure (Fig. S1) shows concentration confidence intervals widening (not read in detail).
- Evidence-increment guidance: YES at the library-design level: read length, fragment-length choice, and combining a short- and a long-fragment library (Tables 1-2, Fig. 3). Optimal average fragment length 400-600 nt for coding genes, 150-200 nt for lncRNAs; combining two libraries of very different fragment length increases identifiable genes. Not tied to a specific ambiguity inside a gene.
- Parsimony: none (no sparsity assumption; full-rank requirement over all paths). Note that this is far stricter than parsimony-based identifiability; only <30% of complex genes (>4 paths) are identifiable even under the best setting (Results, Table 1).
- Complexity: rank of an (#bin pairs x #paths) matrix per gene; path count can be exponential (they use 11 fragment-length points to speed up the convolution).
- Software: GitHub JFerrer-B/transcriptome-identifiability (Abstract).
- Data: GENCODE 38 (60,619 genes), simulated read/fragment length grids; no real-read experiments for identifiability.

### 1.3 Safe&Comp: safety in flow decompositions (Khan et al. 2022; PMID 36288562; companion Caceres et al. 2022, PMID 34847041)

- Problem: given a DAG with edge flows (conservation holds), report subpaths that appear in EVERY flow decomposition ("safe" paths); report all maximal ones ("complete") (Abstract; Section 2).
- Observation model: the splice graph is the flow; edge weights are read coverages assumed to satisfy flow conservation exactly (error-free flow). No separate read-pair/long-read observations; Caceres 2022 considers "constrained path covers" for path covers rather than flows.
- Solution space: ALL flow decompositions of ANY number of paths (Khan Section 2: "in every possible flow decomposition"), not only minimum ones (SafeMFD paper states this contrast explicitly, Introduction).
- Computed: local characterization of safe paths by "excess flow" (Theorem 1: P is w-safe iff excess flow >= w); O(|P|) test; all maximal safe paths in O(mn) (compact output) or O(mn^2) (explicit) (Theorems 1-2, Section 4). A whole s-t path is a special case of a safe path. Verification of a GIVEN path and de novo enumeration are both supported. Output = safe subpaths (sub-transcript contigs), not transcripts with abundances; the problem formulation of Zheng et al. 2022 (AND-Quant), cited there, gives a quadratic algorithm for "does a path appear in all decompositions".
- Ends: only a unique source s and sink t; no TSS/TES evidence modeled. The assumption is that every s-t path is a feasible transcript.
- Counterexamples: the proof of Theorem 1 constructs a decomposition in which an unsafe path carries only excess-flow weight, but no counterexample is an output of the software (not reported).
- Evidence-increment / design: none.
- Parsimony: none (all FDs); Caceres 2022 uses the path-cover length parameter l = k+1.
- Complexity: polynomial.
- Software: https://github.com/algbio/flow-decomposition-safety (listed in the paper), C++.
- Data: simulated "perfect" splice graphs: Catfish (Shao & Kingsford 2017) human/mouse/zebrafish and Reference-Sim (Ensembl GRCh104 transcripts with simulated abundances). Result: ~50% more coverage than extended unitigs at 100% precision; the paper notes Williams et al. found 20% of human genes have multiple minimum FD solutions on error-free data.

### 1.4 SafeMFD (Dias et al. 2023; PMID 37862229)

- Problem: safe paths for MINIMUM flow decompositions (NP-hard) (Abstract).
- Model: exact-flow DAG; ILP model M of k solution paths (binary edge-per-path variables); can add subpath constraints ("every path in R is a subpath of some P_i", constraints 4a-4b, Section 2.1.3) which can model long-read/paired-end evidence.
- Computed: a "safety test" by ILP maximizing the number of tested paths avoided (GroupTest, Algorithm 1); an outer two-pointer/trimming algorithm finds all maximal safe paths. Both verification of a given path and enumeration of maximal safe paths. The test returns the set N of unsafe paths (those avoided in some MFD). Note: the ILP solution that avoids the path is an MFD that acts as a witness, but the paper's output is the set of safe paths; I found no statement that the witness decomposition is returned or presented (UNVERIFIED for the software; README says output is only "maximal safe paths", GitHub algbio/mfd-safety).
- Hardness: computing all safe paths for MFDs is intractable (Turing reduction from UNIQUE 3SAT, in the Supplementary material, per Introduction).
- Ends: only source/sink.
- Parsimony: exactly minimum k. The paper states explicitly (Results/Discussion) that using "at most k' >= k paths" "is still easily achievable with our framework by just changing the ILP blackbox", but leaves it as future work; this is the "minimum + delta" idea of the planned method.
- Evidence increments: none; but "Table 2" shows that adding subpath constraints from half of the transcripts (prefix with three nontrivial junctions) raises coverage ~3% and precision ~1% (evidence-gain measured empirically, not predicted).
- Scalability: 99.8% of >27,000 non-trivial Catfish graphs in 1.5 h total; 2 min budget per graph then fallback to SafeFlow.
- Software: github.com/algbio/mfd-safety (Python, NetworkX + Gurobi; GPL-3.0 per README).
- Data: Catfish, RefSim, StringTie-derived graphs from SRR307903; all "perfect" superpositions of transcripts. Result: up to 90% of full RNA transcripts recovered at 99% precision. Precision drops for graphs from many transcripts because the true set is then less likely to be an MFD (Results).

### 1.5 Scallop (2017) and Scallop2 (2022)

- Problem: reference-based assembly of a parsimonious transcript set from a weighted splice graph preserving long-range phasing paths (Scallop abstract: "preserves long-range phasing paths ... parsimonious set ... minimizing coverage deviation"); Scallop2 extends to multi-end reads (Smart-seq3 barcodes, paired-end) (Abstract).
- Observation model (Scallop2 Methods, "Constructing the splice graph and multi-end phasing paths"): weighted splice graph G=(V,E,w) plus a set C of multi-end phasing paths (read groups sharing a barcode or a pair).
- Computed: heuristic, de novo. Problem 1: compute s-t paths P covering all multi-end phasing paths with d(P,f) and |P| as small as possible. Steps: (1) bridge each multi-end phasing path into a single-end path by the maximum-bottleneck path in the DAG (dynamic programming), (2) refine splice graph by identifying false starting/ending vertices from reads that fail to bridge (a pseudo-edge of weight 0.5 is added; paths using false ends are labelled "transcript fragment"), (3) phase-preserving decomposition integrating multiple decompositions.
- Ends: modeled as source/sink-adjacent vertices; false starts/ends are inferred and downgraded (Methods "Determining false starting/ending vertices"). No separate TSS/TES evidence stream.
- Alternative solutions: not reported; one solution only. Bridging picks one path where several bridging paths exist ("there might be multiple bridging paths", Methods).
- Guidance on extra data: none.
- Parsimony: yes (minimize |P| and coverage deviation) by heuristics, no guarantees.
- Scalability: reasonably fast; slower than StringTie2 (Supplementary Tables 1-2, not read).
- Software: github.com/Shao-Group/scallop2 (and scallop2-test).
- Data: Smart-seq3 HEK293T (192 cells) and mouse fibroblast (369 cells), Illumina ENCODE10; gffcompare accuracy. Scallop2 gives 3.6% more matching transcripts than StringTie2 in HEK293T.

### 1.6 Bookend (Schon et al. 2022; PMID 35768836)

- Problem: assemble full-length transcripts using 5' and 3' end information in reads (end-labeled reads from Smart-seq2/3 template-switch and polyA reads, CAGE/nanoPARE, long reads) (Abstract; Results "A framework for end-guided transcript assembly").
- Model: reads are aligned and stored as end-labeled reads (ELR) recording 5' labels, splice donors/acceptors, gaps, 3' labels; Start/End tag clusters become TSS/PAS boundaries; each locus is partitioned into "frags"; a Membership Matrix (reads x frags) records overlap/exclusion; reads related by "extends/excludes/contained" relations become an overlap graph with global source/sink (Methods, Assembly algorithms).
- Computed: heuristic de novo assembly: collapse linear chains, resolve containment, then greedy paths from the heaviest element maximizing a score (weight, similarity of coverage distribution, variance) subject to keeping source and sink reachable. Also condensed assembly / meta-assembly over single cells.
- Ends: YES, central. Separate end datasets (e.g. nanoPARE TSS reads, 3P-seq) can be combined with internal reads in "hybrid assembly"; the Discussion mentions augmenting unlabeled RNA-seq with tissue-matched FANTOM5 CAGE. However, connection between a start and an end is only through reads that span, the greedy path choice otherwise decides: no representation that a TSS-TES pairing is undetermined.
- Alternative solutions / uncertainty / guidance: none. Output is one transcript set.
- Parsimony: implicit via greedy heaviest-path peeling.
- Software: pip install bookend-rna; github.com/Gregor-Mendel-Institute/bookend (MIT); Zenodo 10.5281/zenodo.6486387.
- Data: Arabidopsis floral bud (Smart-seq2, PacBio FLNC, nanoPARE), mouse ES cells (Smart-seq2 single cells, SIRVs), validated with SLIC-CAGE and 3P-Seq. Reported 2979 novel isoforms vs 8886 (StringTie2), 17,400 (Scallop), 18,934 (Cufflinks) (Results).

### 1.7 IsoTools (Lienhard et al. 2023; PMID 37267159)

- Problem: long-read transcriptome analysis workflow (transcript reconstruction by grouping reads with identical intron chains, quality control, ASE detection via segment-graph bubbles, differential splicing tests, novel transcript classification into SQANTI-style subcategories) (Abstract; Methods).
- Saturation analysis (Methods 2.3.1): a negative-binomial CDF model of the probability of observing a transcript at least r times, with parameters n (sequencing depth), r (minimum coverage) and p (proportional to the transcript's cellular concentration); given two parameters the model suggests the third, e.g. the required depth. This answers "how deep to sequence for a transcript of concentration c to be DETECTED", i.e. a sampling-depth criterion (Results: at 1 TPM, >70% detection; at >=2 TPM near 100% for 2-4 M reads).
- Observations: full-length reads; the segment graph (nodes disjoint exonic segments, edges adjacency in >=1 transcript). Not a constraint solver; no solution space of transcript sets, no ambiguity states; TSS/PAS usage is handled as alternative start/poly-A sites in differential tests, and CAGE peaks are used for validation.
- Ends: not modeled in assembly except through reads; no TSS-TES pairing uncertainty.
- Alternative solutions / design: none beyond depth; the saturation model is a read-count detection model, not an identifiability model.
- Software: GitHub and PyPI (isotools.readthedocs.io).
- Data: human hepatocytes +/- valproic acid, 2,615,181 and 4,200,885 aligned HiFi reads; short-read RNA-seq and ENCODE HepG2 CAGE peaks as comparison.

---

## 2. Additional closer work found (2018-2026)

Ranked roughly by closeness to the planned method.

### A. Zheng, Ma, Kingsford 2022: "Deriving ranges of optimal estimated transcript expression due to nonidentifiability", J Comput Biol 29(2):121-139 (DOI 10.1089/cmb.2021.0444; PMC8892959; PMID not looked up, listed as UNVERIFIED; bioRxiv 10.1101/2019.12.13.875625). CLOSEST.
- Problem: when optimal abundance estimates are non-unique, give for each transcript the minimum and maximum abundance over ALL optimal estimates (Abstract; Intro). The paper says earlier works (Lacroix 2008; Hiller 2009; Roberts & Pachter 2013) only list/flag non-identifiable transcripts but do not give value ranges (Intro).
- Model: (i) reference-transcript quantification by LP over sufficient statistics (junction or phasing-path totals; Section 2.2), (ii) graph quantification: splice graph S-T paths all allowed; ranges from flow decompositions by a max-flow "introspection" algorithm (AND-Quant / OR-Quant, Section 2.4, Theorems 3-4, "Constructing G+ and G- takes linear time"), plus hybrid between the two (Section 2.5). Extends to phasing paths of multi-junction fragments (building on Ma, Zheng, Kingsford 2020/2021, "Exact transcript quantification over splice graphs", WABI 2020 / Algorithms Mol Biol 2021, DOI 10.1186/s13015-021-00184-7, abstract only read).
- Solution space: all flow decompositions of any size (no parsimony); abundance ranges (lower/upper) per transcript, and for sets of edges; polynomial.
- Ends: source S and sink T only (S = start, T = termination); no TSS/TES evidence or pairing.
- Alternative solution/counterexample: not output; ranges only. A lower bound > 0 implies the transcript is in every decomposition (a "NECESSARY" label in the planned method's terms) and an upper bound of 0 means incompatible, but the paper reports ranges, not labels (my reading).
- What-to-sequence guidance: none; the Discussion only says that full-length transcript sequencing may close the gap.
- Data: 16 Human Body Map samples (35-50% of transcripts have non-unique optimal estimates; 20-47% have inconclusive ranking among sibling isoforms), MCF10 and CD8 T cell DE reliability.
- Code: not verified (UNVERIFIED; the Kingsford lab page for the paper exists).
- Implication: abundance lower/upper bounds under non-identifiability ARE already published, but without parsimony and without end evidence.

### B. Williams, Tomescu, Mumey 2023: "Flow Decomposition With Subpath Constraints", IEEE/ACM TCBB 20(1):360-370 (PMID 35104222; DOI 10.1109/TCBB.2022.3147697; WABI 2021 version LIPIcs.WABI.2021.16). Abstract only read.
- Problem: flow decomposition where each given subpath (from long reads/paired-end) must be contained in at least one decomposition path; first formalization with heuristic and FPT algorithms; subpath constraints increase correct ground-truth recovery by 13% (exact minimal) / 30% (heuristic) on instances with larger solution sets. This is the constraint model "every observed long-read subpath must be covered by some transcript". Single solution output, no ambiguity labelling. Ends: source/sink only.

### C. Dias, Williams, Mumey, Tomescu 2022: "Efficient Minimum Flow Decomposition via Integer Linear Programming", J Comput Biol 29(11):1252-1267 (PMC9700332, full text read).
- Quadratic-size ILP for MFD, plus MFD with subpath constraints (MFDSC), inexact flows (intervals; MIFD), and imperfect flow with total-error bound or squared error (Section 3.1-3.4). Gives the exact solver substrate for "coverage with error tolerance + subpath constraints + min or min+delta paths" (ILP is for a given k or upper bound k, Discussion). Each instance solved in under 13 s. Code github.com/algbio/MFD-ILP (stated in text). Not an ambiguity method: one solution; a web summary claimed that it can report all optima, but I could not find that sentence in the full text (UNVERIFIED, treat as not supported).

### D. Flowpaths library (algbio; pip install flowpaths) and related: Sena & Tomescu, "Fast and Flexible Flow Decompositions in General Graphs via Dominators", arXiv:2511.19153 (Nov 2025; a JCB 2026 version appears in PubMed as PMID 42758647, not read); ESA 2025 (DOI 10.4230/LIPIcs.ESA.2025.55) and SEA 2024 (DOI 10.4230/LIPIcs.SEA.2024.14) safe-sequence optimizations (cited by the flowpaths README, not read).
- Library of MILP models: Minimum Flow Decomposition (DAG and cyclic), k-FD, k-Least Absolute Errors, k-Minimum Path Error, k-Minimum Discordant Nodes, Minimum Path Cover, with subpath constraints (DAGs) / subset constraints (cyclic), and "Additional start/end nodes" (flowpaths docs and README). "Safe" sequences are used only as solver speedups preserving optimality (README: "safety optimizations ... preserve global optimality"). Docs mention no enumeration of optimal solutions, no uniqueness test, no ranges/uncertainty (checked on the documentation home page and README: "Not mentioned"). Implication: the planned method can be built on flowpaths, and no existing feature duplicates it.

### E. Chen, Talesara, Thakkar, Shao 2025 (bioRxiv 10.64898/2025.12.11.693570): "Minimum flow decomposition guided by saturating subflows" (catfish-LP; github Shao-Group/catfish-LP). LP-based identification of "safe merge operations" / saturating subflows in the graph, then a greedy decomposition. It reasons about which decomposition equations are always valid ("superficial equations" not supported by any minimum decomposition), i.e. a form of safety reasoning inside a heuristic; no output of alternatives or guidance. Only partially read.

### F. Kannan, Hui, Mazooji, Pachter, Tse: "Shannon: An Information-Optimal de Novo RNA-Seq Assembler", bioRxiv 10.1101/039230 (Feb 2016); preprint only in the pages I could see (journal status UNVERIFIED).
- Gives a necessary condition for unique reconstruction of a transcriptome from reads of length L by any algorithm (error-free, infinite coverage) and a sufficient condition for its algorithm (generic abundances), in terms of repeats / bottlenecks of the k-mer graph; e.g. for the HESC transcriptome "no more than 95% of the transcripts can be reconstructed at read length 100". So it is an information-theoretic identifiability result with read length as the experimental parameter. No TSS/TES, no per-transcript label, no counterexamples. Predates 2018 but is the closest "optimal experimental parameter" analysis besides Ferrer-Bonsoms.

### G. Pelechano, Wei, Steinmetz 2013 (Nature 497:127-131; PMID 23615609; protocol Nat Protoc 2014, PMID 24967623) TIF-seq. Experimental evidence only: by jointly determining both ends of millions of RNA molecules it gives direct 5'-3' pairing (">26 major transcript isoforms per protein-coding gene in yeast"). It defines the assay that would resolve TSS-TES ambiguity; it is not a computational identifiability tool. Not found in the literature: any method that computes which genes' TSS-TES pairing is unresolved from separately measured end sets (UNVERIFIED negative; absence of evidence from the searches described in section 5).

### Not read or not found but cited in the above papers (UNVERIFIED as to content)
- Bernard et al. 2014 (Bioinformatics; PMID 24813214 matched by title search, not opened); Lacroix et al. 2008 (WABI); Salzman, Jiang, Wong 2011 (Stat Sci; sufficient statistics; PMID 24307754 matched by search but not opened); Roberts & Pachter 2013 (identifiability flag; only via Zheng's description); Ray et al. 2020 lrCaptureSeq (cited by Ferrer-Bonsoms).

---

## 3. Gap table

Rows: the three functions of the planned method. Columns: papers. Y = yes, P = partial, N = no. Short justification and citation in each cell.

Legend of columns: H09 = Hiller 2009 (PMID 19762346); F22 = Ferrer-Bonsoms 2022 (34978563); K22 = Khan et al. Safe&Comp (36288562); D23 = SafeMFD (37862229); S2 = Scallop/Scallop2 (29131147, 36713932); B22 = Bookend (35768836); I23 = IsoTools (37267159); Z22 = Zheng-Ma-Kingsford (DOI 10.1089/cmb.2021.0444); W23 = Williams subpath constraints (35104222); Dd22 = Dias MFD-ILP (PMC9700332); Sh16 = Shannon preprint.

### Function 1: transcript-level three-state labels (NECESSARY / POSSIBLE / INCOMPATIBLE) with abundance lower/upper bounds

| Paper | Verdict | Justification |
|---|---|---|
| H09 | P | Yes/no identifiability of a gene's given isoform set, not per-transcript labels or bounds; no ends (Theorem, Methods). |
| F22 | P | Gene-level rank test over ALL start-end paths with ends given; no per-transcript state, no bounds (Methods; Table 1). |
| K22 | P | NECESSARY only, for subpaths or whole s-t paths, over all FDs of any size (Theorem 1); no POSSIBLE/INCOMPATIBLE, no abundance bounds in the output (excess flow gives a guaranteed minimum flow on a path: "w-safe" for w up to f_P, a lower bound; upper bounds not reported). |
| D23 | P | NECESSARY only, under minimum parsimony, with optional subpath constraints (Section 2.1.3); `k' >= k` extension named as easy future work (Section 4 text); no POSSIBLE/INCOMPATIBLE or bounds. |
| S2 | N | One heuristic solution, no states (Scallop2 Problem 1). |
| B22 | N | One solution via greedy paths. |
| I23 | N | No identification states; saturation is detection probability (Methods 2.3.1). |
| Z22 | Y/P | Per-transcript min/max abundance over all optimal estimates (polynomial, max-flow), including phasing paths; implicitly gives necessary (min>0) and incompatible (max=0) but no parsimony bound, no end evidence, no explicit labels; whether stated as labels UNVERIFIED (not seen). |
| W23 | N | Constraint handling for one solution. |
| Dd22 | N | One solution under exact/inexact/imperfect flow; supports the constraint types needed. |
| Sh16 | N | Reconstructability bound only, no per-transcript labels. |

### Function 2: counterexamples (an alternative feasible solution without the transcript) and certificates (observations that force it)

| Paper | Verdict | Justification |
|---|---|---|
| H09 | P | A toy non-identifiable example (Fig. 1); no algorithmic output. |
| F22 | N | Rank deficiency reported; no alternative solution. |
| K22 | P | Proof of Theorem 1 constructs a decomposition with only excess-flow weight on P (existence of counterexample), the certificate is the excess flow value; software outputs safe paths only. |
| D23 | P | The ILP solution avoiding the path is a counterexample MFD, but only the safe set is documented as output (README). Certificate: none beyond "all ILPs infeasible". |
| S2 | N | Not applicable. |
| B22 | N | Not applicable. |
| I23 | N | Not applicable. |
| Z22 | P | Proofs reconstruct the decomposition attaining each bound, but the tool reports ranges (the reconstructed flow is part of the proof of Theorems 3-4, Appendix); a witness output is not claimed. |
| W23 | N | None. |
| Dd22 | N | None (one optimal solution; claim "reports all optima" not found, UNVERIFIED). |
| Sh16 | N | None. |

### Function 3: missing-evidence diagnosis (which new observation type resolves the ambiguity) and prediction of the most useful extra sequencing

| Paper | Verdict | Justification |
|---|---|---|
| H09 | P | Qualitative: add paired-end reads or quantify events instead of isoforms (Discussion), and adding junction probes (Fig. 1b) fixes the toy case; no algorithm. |
| F22 | Y/P | Strongest existing: quantifies how read length, fragment length and a second library change gene-level identifiability (Tables 1-2, Fig. 3); not per ambiguity, not for ends, not for long-read subpaths, candidate set is all start-end paths. |
| K22 | N | None. |
| D23 | P | Empirical gain from subpath constraints (Table 2) but no prediction of which constraint to add. |
| S2 | N | None (but it shows phasing information from multi-end reads improves accuracy). |
| B22 | N | None; says utility limited by availability of end-labeled data. |
| I23 | P | Depth needed to detect a transcript of given concentration (saturation model); not about structural ambiguity. |
| Z22 | N | None (only remarks that full-length sequencing would help). |
| W23 | P | Shows empirically that extra subpath information helps; no guidance on which. |
| Dd22 | N | None. |
| Sh16 | P | Read-length thresholds for reconstructability (information-theoretic), no per-ambiguity guidance. |

---

## 4. Frank verdict

Is there a defensible novelty gap? Yes, but it is narrower than the plan suggests, and function 1 is largely covered if described in isolation. The defensible gap is the COMBINATION: (i) an evidence-only splice graph with separate start/end sets and observation constraints; (ii) a tri-state classification over a parsimony band (min .. min+delta) including POSSIBLE and INCOMPATIBLE; (iii) witnesses (certificate or counterexample); and (iv) a diagnosis of which new observation type would change a label, with an experiment-ranking application. I found no paper doing (iv) at the level of individual ambiguities, and none doing (i)-(iii) together with end evidence.

What is already done (cannot be claimed as new):
- NECESSARY for subpaths and whole s-t paths over all flow decompositions (Khan 2022) and over minimum flow decompositions (SafeMFD 2023), with subpath constraints already supported in SafeMFD (Section 2.1.3, README).
- Per-transcript abundance lower and upper bounds under non-identifiability (Zheng-Ma-Kingsford 2022), polynomially via max-flow, including multi-junction phasing paths; the transcript-level "range" idea is theirs. This is the strongest overlap: a reviewer will say "NECESSARY is lower bound > 0 in Zheng; INCOMPATIBLE is upper bound = 0."
- The statement that "up to k' >= k paths" can replace exactly-k in the safety framework is made by SafeMFD authors (as future work); so "minimum + delta" is not by itself novel, only its evaluation and the POSSIBLE/INCOMPATIBLE side are.
- Exact ILP models with subpath constraints, inexact flows and error tolerance (Dias 2022; flowpaths library) exist; the solver layer is not novel.
- Library-design identifiability (read/fragment length, two-library combination) at the gene level (Ferrer-Bonsoms 2022), and read-length reconstructability bounds (Shannon 2016).
- End-aware assembly (Bookend 2022, Scallop2 2022 false-start/end handling) exists, but returns one solution.

Most novel function: 3 (missing-evidence diagnosis and experiment ranking), then 2 (explicit counterexamples and certificates as outputs), then 1 (labels). The TSS-TES-pairing question for separately measured CAGE and QuantSeq ends is, in what I read, handled only implicitly: Ferrer-Bonsoms assume ends known and all start-end pairs candidates (so pairing is in the rank test as unresolved when rank deficient), Zheng et al. have only S and T, Bookend has end tags but picks one greedy pairing. None report "TSS-TES pairing is undetermined for transcripts X, Y" nor name the TIF-seq-type observation needed (Pelechano 2013 gives the assay only).

Precise differences the paper should claim:
1. Input is evidence, not an annotation and not a given flow: graph built only from observed starts, ends, segments, junctions; start/end evidence restricts path endpoints (a constraint family absent in K22, D23, Z22, where all s-t paths are allowed).
2. Solution space is a parsimony band with explicit delta, plus optional coverage tolerance, and three states including INCOMPATIBLE (absent from all feasible solutions), which none of the safety papers report; abundance bounds are computed over this band (Zheng computes bounds over all FDs, i.e. delta = infinity).
3. Every label comes with a witness: forcing observation set (certificate) or an alternative feasible solution (counterexample). Existing work outputs the safe set or ranges only.
4. Missing-evidence diagnosis: for each ambiguous transcript, the minimal new observation types (a spanning read across two branches, a TSS-TES-linked molecule, a long read over k events, more molecules for abundance) that change its state, and a ranking of candidate sequencing additions by expected ambiguity reduction; to be validated by simulation and by real held-out evidence (e.g., TIF-seq-like or long-read data), ideally showing predictions confirmed.
5. Quantify robustness of labels to delta; SafeMFD itself reports that precision of minimum-based safety falls on complex graphs because the truth is less likely to be an MFD (Results), so delta-aware labels address a documented weakness.

Work that could make the plan redundant: no single paper does so. Closest threats are (a) Zheng-Ma-Kingsford 2022 (ranges) combined with Khan 2022 / Dias 2023 (necessary paths with subpath constraints): together they cover function 1 and part of 2 for the no-end, all-FD or exact-minimum cases; (b) a SafeMFD follow-up that implements "at most k'" and returns counterexamples, which the authors say is straightforward (UNVERIFIED whether this exists, since a later ESA 2025 / JCB 2026 line of the Tomescu group was not fully read); (c) Ferrer-Bonsoms-type design analysis extended to per-gene end pairing. To stay safe, the paper should (1) lead with function 3 and the evidence-driven diagnosis, (2) state clearly that safe paths and abundance ranges are special cases (delta = infinity, no end constraints), and (3) show that on end-separated evidence the labels differ materially from those under the unconstrained all-s-t-path model.

## 5. Search coverage and limits

- Queries run (Europe PMC, PubMed, web): safety / safe paths in flow decompositions, subpath constraints, flowpaths, AND-Quant, nonidentifiability ranges, information-theoretic transcript reconstruction, identifiability with long reads, experimental design for transcript assembly, TIF-seq pairing, enumeration of all optimal flow decompositions.
- Not covered / UNVERIFIED: whether any 2026 paper (e.g. from the Tomescu, Kingsford, Shao groups) implements counterexamples or "at most k'" safety; full text of Williams 2023 and Caceres 2022; the ESA 2025 and SEA 2024 papers; the JCB 2026 dominator paper; long-read assemblers (IsoQuant, Bambu, StringTie2 long-read mode) were not examined for end-pairing uncertainty; LRGASP-type benchmark papers were not examined; I did not search patents or conference-only abstracts; PMID for Zheng 2022 not retrieved (DOI and PMC id are confirmed).
