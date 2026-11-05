#!/bin/bash
# S2-AG analysis on ont, CPU only (PLAN_S2_AG.md sections 2-3): waits until ag_colab_run.sh has fetched all four
# chromosomes (sha256 checked), then offsets on chr20 -> score files for chr18 / chr9 / chr7 -> s1_analyze.py per
# chromosome (model "alphagenome") -> RNA-seq segments (threshold on chr18; reads chr18 / chr20 / chr9).
cd ~/ncrna-annot || exit 1
D=outputs/PHASE1/S2AG
R=$D/auto/raw
PY=~/envs/generanno/bin/python
mkdir -p $D/logs
log() { echo "$(date -u +%FT%TZ) $*" >> $D/post_queue.log; }
until [ -e $R/chr20.fetched ] && [ -e $R/chr18.fetched ] && [ -e $R/chr9.fetched ] && [ -e $R/chr7.fetched ]; do
  grep -q "END wall_s" $D/auto/run.log 2>/dev/null && ! { [ -e $R/chr20.fetched ] && [ -e $R/chr18.fetched ] && [ -e $R/chr9.fetched ] && [ -e $R/chr7.fetched ]; } && { log "Colab run ended without all four chromosomes"; exit 1; }
  sleep 120
done
A="nice -n 19 taskset -c 2-31 env OMP_NUM_THREADS=8 $PY"
log START
$A scripts/phase1/s2ag/ag_post.py --chrom chr20 --pick-offsets --out $D/ag_chr20.npy > $D/logs/post20.log 2>&1 || { log "post20 FAIL"; exit 1; }
for c in chr18 chr9 chr7; do
  $A scripts/phase1/s2ag/ag_post.py --chrom $c --offsets $D/ag_chr20.npy.json --out $D/ag_$c.npy > $D/logs/post_$c.log 2>&1 || { log "post $c FAIL"; exit 1; }
done
log "scores done"
S1="$A scripts/phase1/s1/s1_analyze.py --judge alphagenome"
$S1 --chrom chr20 --scores alphagenome=$D/ag_chr20.npy --out $D/S1AG_chr20.json > $D/logs/s1_chr20.log 2>&1; log "s1 chr20 rc=$?"
$S1 --chrom chr18 --codes data/p1/chr18_codes.npy --snaptron "data/batch2/snaptron/chr18_*.tsv" \
  --scores alphagenome=$D/ag_chr18.npy --out $D/S1AG_chr18.json > $D/logs/s1_chr18.log 2>&1; log "s1 chr18 rc=$?"
$S1 --chrom chr7 --codes data/splice_own/codes_chr7.npy --snaptron "data/batch2/snaptron/chr7_*.tsv" \
  --scores alphagenome=$D/ag_chr7.npy --out $D/S1AG_chr7.json > $D/logs/s1_chr7.log 2>&1; log "s1 chr7 rc=$?"
$S1 --chrom chr9 --codes data/p1/chr9_codes.npy --snaptron "data/v20_chr9/snaptron/chr9_*.tsv" \
  --scores alphagenome=$D/ag_chr9.npy --out $D/S1AG_chr9.json > $D/logs/s1_chr9.log 2>&1; log "s1 chr9 rc=$?"
$A scripts/phase1/s2ag/ag_segments.py --out $D/SEGMENTS.json > $D/logs/segments.log 2>&1; log "segments rc=$?"
log QUEUE_END
