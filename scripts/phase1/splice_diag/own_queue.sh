#!/bin/bash
# Own splice net + chr18 re-diagnosis on ont (user decision 2026-11-03, card splice-next option 1;
# docs/phase1_results/PLAN_SPLICE_OWN.md): train s2k then s10k -> select on chr22 -> score chr18 ->
# splice_diag.py w in {0.5, 1, 2} x 3 main runs (up to 4 at a time) -> summarize.py -> verdict.
# CPU: nice 19, cores 2-31, 2 threads per job; GPU: gpu_gate.sh before every job (exit 75 = gate busy -> retry).
cd ~/ncrna-annot || exit 1
. scripts/phase1/gpu_gate.sh
D=outputs/PHASE1/SPLICE_OWN
PY=~/envs/generanno/bin/python
OS=scripts/phase1/splice_diag/own_splice.py
mkdir -p $D/diag
log() { echo "$(date -u +%FT%TZ) $*" >> $D/queue.log; }
run1() {  # <name> <script> <args...>
  local n=$1 rc; shift
  while true; do
    gpu_wait
    nice -n 19 taskset -c 2-31 env NCRNA_THREADS=2 OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=2 MKL_NUM_THREADS=2 \
      $PY "$@" >> $D/$n.log 2>&1
    rc=$?
    [ $rc -eq 75 ] && { sleep 120; continue; }
    log "$n rc=$rc"
    return $rc
  done
}
log START
[ -f data/splice_own/prep.json ] || { log "no prep"; exit 1; }
run1 train_s2k $OS train --variant s2k || exit 1
run1 train_s10k $OS train --variant s10k || exit 1
$PY $OS select > $D/select.log 2>&1 || { log SELECT_FAIL; exit 1; }
log "selected $(grep -o '"selected": "[a-z0-9]*"' $D/SELECT.json)"
run1 score_chr18 $OS score --codes data/p1/chr18_codes.npy --truth outputs/PHASE1/DIAG_L2/chr18_v50.gtf --chrom chr18 \
  --out $D/own_chr18.npy || exit 1
for w in 1 2 0.5; do for r in g1b-42-ke g1b-123-ke g1b-456-ke; do
  while [ "$(jobs -rp | wc -l)" -ge 4 ]; do wait -n; done
  run1 diag_${r}_w$w scripts/phase1/splice_diag/splice_diag.py --run $r --w $w --sai $D/own_chr18.npy --out $D/diag/${r}_w$w &
  sleep 90  # let the job take its GPU memory before the next gate check
done; done
wait
$PY scripts/phase1/splice_diag/summarize.py --dir $D/diag > $D/summary.md 2> $D/summary.err || { log SUMMARY_FAIL; exit 1; }
$PY $OS verdict --dir $D/diag > $D/verdict.log 2>&1 || { log VERDICT_FAIL; exit 1; }
log "QUEUE_END effective=$(grep -o '"effective": [a-z]*' $D/diag/VERDICT.json | tail -1)"
