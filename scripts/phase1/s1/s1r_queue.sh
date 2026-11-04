#!/bin/bash
# S1R queue on ont (PLAN_BATCH2.md A; card batch2-downloads 2026-11-04).  Waits for the mouse queue (one ncRNA GPU
# queue at a time) and for the Snaptron downloads, then: SpliceAI / own S10k on chr7, Pangolin on chr20 (offsets
# picked there) / chr18 / chr7, one GPU job at a time behind gpu_gate.sh; then s1_analyze.py on chr18 (judge own)
# and chr7 (judge spliceai), chr20 again with Pangolin added (new output file), and the R4 second half on chr18.
cd ~/ncrna-annot || exit 1
. scripts/phase1/gpu_gate.sh
D=outputs/PHASE1/S1
PY=~/envs/generanno/bin/python
mkdir -p $D/logs
log() { echo "$(date -u +%FT%TZ) $*" >> $D/s1r_queue.log; }
until grep -q FETCH_END data/batch2/fetch.log 2>/dev/null && grep -qE "QUEUE_END|FAIL" outputs/PHASE1/P17_MOUSE/queue.log 2>/dev/null; do sleep 300; done
step() {  # name done-file command...
  local name=$1 done=$2 rc
  shift 2
  if [ -s "$done" ]; then log "$name skip"; return 0; fi
  while true; do
    gpu_wait >> $D/s1r_queue.log
    log "$name start"
    nice -n 19 taskset -c 2-31 env OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 NCRNA_GPU_FRAC=0.25 "$@" > $D/logs/$name.log 2>&1
    rc=$?
    [ $rc -eq 75 ] && { log "$name gate busy, retry"; sleep 120; continue; }
    log "$name rc=$rc"
    return $rc
  done
}
G50=data/human/gencode.v50.annotation.gtf.gz
log START
step sai7 $D/sai_chr7.npy $PY scripts/phase1/splice_diag/sai_score.py --codes data/splice_own/codes_chr7.npy --truth $G50 --chrom chr7 --out $D/sai_chr7.npy || exit 1
step own7 $D/own_chr7.npy $PY scripts/phase1/splice_diag/own_splice.py score --codes data/splice_own/codes_chr7.npy --truth $G50 --chrom chr7 --out $D/own_chr7.npy || exit 1
step pang20 $D/pang_chr20.npy $PY scripts/phase1/s1/pangolin_score.py --codes data/p1/chr20_codes.npy --chrom chr20 --pick-offsets --out $D/pang_chr20.npy || exit 1
step pang18 $D/pang_chr18.npy $PY scripts/phase1/s1/pangolin_score.py --codes data/p1/chr18_codes.npy --chrom chr18 --offsets $D/pang_chr20.npy.json --out $D/pang_chr18.npy || exit 1
step pang7 $D/pang_chr7.npy $PY scripts/phase1/s1/pangolin_score.py --codes data/splice_own/codes_chr7.npy --chrom chr7 --offsets $D/pang_chr20.npy.json --out $D/pang_chr7.npy || exit 1
log GPU_DONE
A="nice -n 19 taskset -c 2-31 env OMP_NUM_THREADS=8 $PY scripts/phase1/s1/s1_analyze.py"
$A --chrom chr18 --codes data/p1/chr18_codes.npy --snaptron "data/batch2/snaptron/chr18_*.tsv" --judge own \
  --scores spliceai=outputs/PHASE1/SPLICE_DIAG/sai_chr18.npy own=outputs/PHASE1/SPLICE_OWN/own_chr18.npy \
  vipsania=outputs/PHASE1/VIPSANIA/chr18/post6.npy:2 pangolin=$D/pang_chr18.npy --out $D/S1R_chr18.json > $D/logs/analyze18.log 2>&1
log "analyze18 rc=$?"
$A --chrom chr7 --codes data/splice_own/codes_chr7.npy --snaptron "data/batch2/snaptron/chr7_*.tsv" --judge spliceai \
  --scores spliceai=$D/sai_chr7.npy own=$D/own_chr7.npy pangolin=$D/pang_chr7.npy --out $D/S1R_chr7.json > $D/logs/analyze7.log 2>&1
log "analyze7 rc=$?"
$A --chrom chr20 --scores pangolin=$D/pang_chr20.npy --judge pangolin --out $D/S1R_chr20_pangolin.json > $D/logs/analyze20p.log 2>&1
log "analyze20p rc=$?"
nice -n 19 $PY scripts/phase1/s1/s1_r4b.py --out $D/R4B_chr18.json > $D/logs/r4b.log 2>&1
log "r4b rc=$? QUEUE_END"
