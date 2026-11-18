#!/bin/bash
# chr18 development grid for E2: jobs from a file (one "run cond frac seed beta tol" per line), 4 at a time, CPU only.
#   e2_grid.sh <abs jobs.txt> <abs log>      (stderr of every job kept in <log>.err)
J=$1; LOG=$2
cd $HOME/ncrna-annot/scripts/m1/w24 || exit 1
PY=$HOME/envs/generanno/bin/python
n=0
while read -r run cond frac seed beta tol; do
  [ -z "$run" ] && continue
  ( nice -n 19 taskset -c 2-31 $PY e2_decode.py --target chr18 --run $run --cond $cond --frac $frac --seed $seed \
      --beta $beta --tol $tol --procs 6 > $LOG.$n.out 2>> $LOG.err; echo "$run $cond $frac $seed $beta $tol rc=$?" >> $LOG ) &
  n=$((n + 1))
  while [ "$(jobs -rp | wc -l)" -ge 4 ]; do sleep 5; done
done < $J
wait
echo "GRID_END $(date -u +%FT%TZ)" >> $LOG
