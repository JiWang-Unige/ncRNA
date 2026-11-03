#!/bin/bash
# chr18 splice diagnostic queue on ont (user decision 2026-11-03): w = 0 reproduction first (post sha256 must equal
# the P1.5 chr18 read), then w in {0.5, 1, 2} x 3 main runs and one truth-site oracle (seed 42), up to 4 at a time.
# CPU: nice 19, cores 2-31, 2 threads per job; GPU: gpu_gate.sh before every job (exit 75 = gate busy -> retry).
cd ~/ncrna-annot || exit 1
. scripts/phase1/gpu_gate.sh
D=outputs/PHASE1/SPLICE_DIAG
PY=~/envs/generanno/bin/python
run1() {  # <name> <splice_diag.py args...>
  local n=$1 rc; shift
  while true; do
    gpu_wait
    nice -n 19 taskset -c 2-31 env NCRNA_THREADS=2 OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=2 MKL_NUM_THREADS=2 \
      $PY scripts/phase1/splice_diag/splice_diag.py --out $D/$n "$@" >> $D/$n.log 2>&1
    rc=$?
    [ $rc -eq 75 ] && { sleep 120; continue; }
    echo "$(date -u +%FT%TZ) $n rc=$rc" >> $D/queue.log
    return $rc
  done
}
echo "$(date -u +%FT%TZ) START" >> $D/queue.log
run1 g1b-42-ke_w0 --run g1b-42-ke --w 0 || exit 1
grep -q '"post_matches_p15": true' $D/g1b-42-ke_w0.log || { echo "$(date -u +%FT%TZ) REPRO_FAIL" >> $D/queue.log; exit 1; }
echo "$(date -u +%FT%TZ) REPRO_OK" >> $D/queue.log
JOBS=()
for w in 1 2 0.5; do for r in g1b-42-ke g1b-123-ke g1b-456-ke; do JOBS+=("${r}_w$w --run $r --w $w"); done; done
JOBS+=("g1b-42-ke_oracle --run g1b-42-ke --w 1 --oracle")
for j in "${JOBS[@]}"; do
  while [ "$(jobs -rp | wc -l)" -ge 4 ]; do wait -n; done
  run1 $j &
  sleep 90  # let the job take its GPU memory before the next gate check
done
wait
echo "$(date -u +%FT%TZ) QUEUE_END" >> $D/queue.log
