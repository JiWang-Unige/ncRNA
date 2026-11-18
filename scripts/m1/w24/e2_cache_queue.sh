#!/bin/bash
# Run e2_cache.py for chr18 through the GPU gate (exit 75 = busy -> wait and retry).  CPU nice 19, cores 2-31.
cd ~/ncrna-annot || exit 1
. scripts/phase1/gpu_gate.sh
O=outputs/M1/W24/cache; mkdir -p $O
while true; do
  gpu_wait
  nice -n 19 taskset -c 2-31 env NCRNA_GPU_FRAC=0.25 OMP_NUM_THREADS=4 ~/envs/generanno/bin/python scripts/m1/w24/e2_cache.py --target chr18 >> $O/cache.log 2>&1
  rc=$?
  [ $rc -eq 75 ] && { sleep 120; continue; }
  echo "$(date -u +%FT%TZ) CACHE_END rc=$rc" >> $O/cache.log
  break
done
