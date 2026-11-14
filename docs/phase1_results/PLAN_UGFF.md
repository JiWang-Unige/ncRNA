# PLAN_UGFF — selective unified GFF3 (CONVERGENCE_20261112 item 5)

> Written 2026-11-14, before any merge or scoring. Engineering deliverable, not a pre-registered claim. The two
> rule parameters (`t_int`, the CpG cut-off) are set on chr18 (development) by the rules below and then applied
> unchanged to chr20 and chr9. No model is re-read on chr20 or chr9: their lnc layers are the predictions already
> read and reported (chr20 v1.9 `P15/chr20`, chr9 v2.0 `V20/read/chr9`).

## 1. What the file contains

One GFF3 per chromosome, `outputs/PHASE1/UGFF/<chrom>/unified.gff3`, three layers:

| Layer | Source | Feature types |
|---|---|---|
| Coding | Tiberius `mammalia_nosofttmasking_v2`, run as in PREREG_BASELINES (chr20 existing; chr18, chr9 new, `ugff/tiberius_queue.sh`) | `gene` → `mRNA` → `exon`, `CDS` (CDS only: Tiberius predicts no UTR) |
| lncRNA | our shipped model (seed 42 run: chr20 `g1b-42-ke`, v1.9, raw confidence ≥ 0.6; chr18 / chr9 `v20-g1b-42-ke`, v2.0, calibrated probability ≥ 0.36) | `lncRNA_gene`, and only when the structure is resolved (§3) `lnc_RNA` → `exon` |
| Small ncRNA | Phase 0 0.6 Rfam / tRNAscan with the P1.6 integration rules (chr20 only; no small-RNA run exists for chr18 / chr9) | `ncRNA_gene` → `ncRNA` → `exon` |

Our own mRNA predictions are not written: coding structure comes from Tiberius (our coding L2 is far weaker).

## 2. Coding guardrail

An lnc gene is dropped when any of its exons overlaps a Tiberius CDS on the same strand. Fixed a priori, not tuned.
Reported on all three chromosomes: number dropped, and among them how many hit a v50 lnc gene, only a v50 mRNA
gene, or nothing; lnc L1 Sn and strict Pr before and after.

## 3. Selective structure

- Intron confidence = mean over the intron's bases of the summed lnc intron-state posteriors (the intron term of
  `chrom.gene_confidence`), from the shipped run's `post.npy`.
- A transcript (`lnc_RNA` + exons) is written only for a multi-exon lnc gene whose every intron has confidence ≥
  `t_int`. All other lnc genes (single-exon, or any intron below `t_int`) are written as `lncRNA_gene` only, with
  `structure=unresolved` and the decoded exons kept as an attribute (`decoded_exons=`), so no transcript model is
  asserted that the evidence does not support.
- **Rule for `t_int` (chr18, three v2.0 runs):** grid 0.50–0.95 step 0.05; intron precision = share of written
  lnc introns that exactly match an intron of a v50 lncRNA transcript on the same strand. `t_int` = the smallest
  grid value whose three-run mean intron precision is ≥ 0.50 with ≥ 20 written introns in every run. If no value
  qualifies, `t_int = none`: every lnc gene is written at locus level.
- Error cost reported on chr18 for forced (all decoded exons written as one transcript, as `pred.gff3` does) versus
  selective output: written introns, intron precision (lnc introns; any v50 intron), and exact intron chains.

## 4. Start annotation

- `promoter_cpg_oe` = CpG observed/expected of ±500 bp around the predicted 5′ end (as in END / DEDEND).
- `start_class=CpG_rich` when it is ≥ the upper quartile of the same statistic over annotated v50 lncRNA gene
  starts on chr18 (computed once, stored); otherwise `other`. Every lnc gene also carries `end_precision=low`.
- Development check on chr18 (descriptive): among lnc predictions that hit a v50 lnc gene, share whose 5′ end is
  within 1 kb of a 5′ end of a hit gene's transcripts, CpG_rich versus other.

## 5. Other attributes and checks

- `confidence` (the score the threshold was applied to), `support=k/3` (runs of the same model family with a
  kept, guardrail-passing lnc gene overlapping on the same strand).
- Small RNA overlapping an lnc gene's decoded exons on the same strand: `Note=overlaps_lnc:<ID>` (P1.6 rule).
- Validity: unique IDs, every `Parent` defined, children inside parents, sorted by start; the file is re-read with
  `ruler.read_pred_gff3` and the lnc layer re-scored to confirm the round trip.
- `SUMMARY.json` per chromosome: counts per layer, guardrail effect, resolved transcripts, Tiberius mRNA L1 / L2.
