#!/bin/bash
# P1.6 zero-shot panel queue on ont (user decision 2026-11-03, card p15-next): chr18 reproduction of the three
# main runs (bit for bit vs the P1.5 batch-2 read) -> the four panel chromosomes -> scoring.  Runs of one target in
# parallel; CPU nice 19 on cores 2-31, 2 threads per job; GPU gate before every job (exit 75 = busy -> retry).
cd ~/ncrna-annot || exit 1
. scripts/phase1/gpu_gate.sh
PY=~/envs/generanno/bin/python
C=configs/p16_v19.json; O=outputs/PHASE1/P16; mkdir -p $O
RUNS="g1b-42-ke g1b-123-ke g1b-456-ke"
log() { echo "$(date -u +%FT%TZ) $*" >> $O/queue.log; }
run_t() {  # <target> [<run>]
  local rc
  while true; do
    gpu_wait
    nice -n 19 taskset -c 2-31 env OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMBA_NUM_THREADS=2 NCRNA_GPU_FRAC=0.15 \
      $PY scripts/phase1/p16_read.py --config $C --target $1 ${2:+--runs $2} >> $O/read_$1${2:+_$2}.log 2>&1
    rc=$?
    [ $rc -eq 75 ] && { sleep 120; continue; }
    log "$1 ${2:-all} rc=$rc"; return $rc
  done
}
log START
for r in $RUNS; do run_t chr18 $r & sleep 60; done; wait
run_t chr18   # every run skipped (read_meta.json exists); rewrites REPRO.json / REPRO_OK from all three
[ -e $O/chr18/REPRO_OK ] || { log REPRO_FAIL; exit 1; }
log REPRO_OK
for t in worm_I fly_2L zebrafish_22 chicken_5; do
  for r in $RUNS; do run_t $t $r & sleep 60; done; wait
done
nice -n 19 taskset -c 2-31 env PYTHONHASHSEED=0 $PY scripts/phase1/p16_score.py --config $C > $O/score.log 2>&1
log "score rc=$?"
log QUEUE_END
