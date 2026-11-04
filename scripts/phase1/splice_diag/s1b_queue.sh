#!/bin/bash
# S1b queue on ont (PLAN_S1B.md, card subjournal-portfolio 2026-11-04): structure-failure decomposition on chr18.
# A0 = v1.9 reproduction (post sha256 must equal the P1.5 chr18 read), A1 = lnc gene spans known (--span-oracle),
# A2 = truth splice sites (--oracle --w 1, as SPLICE_DIAG), A3 = both; 3 main runs each; at most 2 jobs at a time
# (GPU fraction 0.15 each, <= 30% together); gpu_gate.sh before every job (exit 75 = gate busy -> retry).
cd ~/ncrna-annot || exit 1
. scripts/phase1/gpu_gate.sh
D=outputs/PHASE1/S1B
PY=~/envs/generanno/bin/python
mkdir -p $D
run1() {  # <name> <splice_diag.py args...>
  local n=$1 rc; shift
  while true; do
    gpu_wait
    nice -n 19 taskset -c 2-31 env NCRNA_THREADS=2 OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=2 MKL_NUM_THREADS=2 NCRNA_GPU_FRAC=0.15 \
      $PY scripts/phase1/splice_diag/splice_diag.py --out $D/$n "$@" >> $D/$n.log 2>&1
    rc=$?
    [ $rc -eq 75 ] && { sleep 120; continue; }
    echo "$(date -u +%FT%TZ) $n rc=$rc" >> $D/queue.log
    return $rc
  done
}
echo "$(date -u +%FT%TZ) START" >> $D/queue.log
run1 g1b-42-ke_A0 --run g1b-42-ke --w 0 || exit 1
grep -q '"post_matches_p15": true' $D/g1b-42-ke_A0.log || { echo "$(date -u +%FT%TZ) REPRO_FAIL" >> $D/queue.log; exit 1; }
echo "$(date -u +%FT%TZ) REPRO_OK" >> $D/queue.log
JOBS=()
for r in g1b-42-ke g1b-123-ke g1b-456-ke; do
  JOBS+=("${r}_A1 --run $r --w 0 --span-oracle" "${r}_A2 --run $r --w 1 --oracle" "${r}_A3 --run $r --w 1 --oracle --span-oracle")
done
for j in "${JOBS[@]}"; do
  while [ "$(jobs -rp | wc -l)" -ge 2 ]; do wait -n; done
  run1 $j &
  sleep 90  # let the job take its GPU memory before the next gate check
done
wait
cmp -s <(python3 -c "import json;print(json.dumps(json.load(open('$D/g1b-42-ke_A2/genes.json'))))") \
       <(python3 -c "import json;print(json.dumps(json.load(open('outputs/PHASE1/SPLICE_DIAG/g1b-42-ke_oracle/genes.json'))))") \
  && echo "$(date -u +%FT%TZ) A2_REPRO_OK" >> $D/queue.log || echo "$(date -u +%FT%TZ) A2_REPRO_FAIL" >> $D/queue.log
echo "$(date -u +%FT%TZ) QUEUE_END" >> $D/queue.log
