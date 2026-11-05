#!/bin/bash
# Resume of v20_chr9_queue.sh after an ont reboot during the section 4.2 replications: only the steps after the
# chr9 read and scoring (which are not repeated): remaining S1b decodes (a decode with result.json is skipped),
# S3, verdicts, chr9 feature cache deletion.  Same commands and settings as the queue.
cd ~/ncrna-annot || exit 1
. scripts/phase1/gpu_gate.sh
O=outputs/PHASE1/V20; REP=$O/rep; PY=~/envs/generanno/bin/python; export PYTHONPATH=src
C=configs/v20.json; N="nice -n 19 taskset -c 2-31"
log() { echo "$(date -u +%FT%TZ) $*" >> $O/chr9_queue.log; }
gstep() {
  local name=$1 done=$2 rc; shift 2
  [ -s "$done" ] && { log "$name skip"; return 0; }
  while true; do
    gpu_wait >> $O/chr9_queue.log
    log "$name start"
    $N env OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 "$@" > $REP/logs/$name.log 2>&1; rc=$?
    [ $rc -eq 75 ] && { sleep 120; continue; }
    log "$name rc=$rc"; return $rc
  done
}
log "RESUME after ont reboot (~2026-11-05T04:17Z); chr9 read, scoring, R-S1 already done and not repeated"
for r in v20-g1b-42-ke v20-g1b-123-ke v20-g1b-456-ke; do
  for arm in "A1 --w 0 --span-oracle" "A2 --w 1 --oracle" "A3 --w 1 --oracle --span-oracle"; do
    set -- $arm; a=$1; shift
    while [ "$(jobs -rp | wc -l)" -ge 2 ]; do wait -n; done
    gstep s1b_${r}_$a $REP/S1B/${r}_$a/result.json env NCRNA_THREADS=2 NUMBA_NUM_THREADS=2 NCRNA_GPU_FRAC=0.15 \
      $PY scripts/phase1/splice_diag/splice_diag.py --config $C --target chr9 --run $r "$@" --out $REP/S1B/${r}_$a &
    sleep 90
  done
done
wait
for r in v20-g1b-42-ke v20-g1b-123-ke v20-g1b-456-ke; do for a in A1 A2 A3; do [ -s $REP/S1B/${r}_$a/result.json ] || { log "S1B $r $a missing - stop"; exit 1; }; done; done
log "R-S1b decoded"
$N env OMP_NUM_THREADS=4 $PY scripts/phase1/s3/s3_features.py --config $C --chrom chr9 --v20 --out $REP/S3_chr9.json > $REP/logs/s3.log 2>&1 || { log S3_FAIL; exit 1; }
$PY scripts/phase1/v20/v20_replicate.py > $REP/logs/replicate.log 2>&1 || { log REPLICATE_FAIL; exit 1; }
log "replications done"
rm -f data/p1_feat/chr9.npy && log "chr9 feature cache deleted"
log QUEUE_END
