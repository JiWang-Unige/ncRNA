#!/bin/bash
# Ops wrapper (not protocol code): finish the preview queue from rat_12 on. v30_test_queue.sh is not resumable past a
# deleted test cache (it would re-extract mouse_chr19, already read). Same steps as its group-C loop, verbatim
# commands; extraction pinned to cpus 16-31 (random CPU faults on cpu 0).
cd ~/ncrna-annot || exit 1
. scripts/phase1/gpu_gate.sh
PY=~/envs/generanno/bin/python
C=configs/v30_preview.json; SFX=_preview; STAGE=preview; O=outputs/PHASE1/V30
log() { echo "$(date -u +%FT%TZ) $*" >> $O/test_queue$SFX.log; }
RUNS=$(python3 -c "import json; print(\" \".join(json.load(open(\"$C\"))[\"targets\"][\"chr18\"][\"runs\"]))")
run_t() {
  local rc
  while true; do
    gpu_wait
    nice -n 19 taskset -c 2-31 env OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMBA_NUM_THREADS=2 NCRNA_GPU_FRAC=0.15 \
      $PY scripts/phase1/p16_read.py --config $C --target $1 --runs $2 >> $O/read_$1_$2.log 2>&1
    rc=$?
    [ $rc -eq 75 ] && { sleep 120; continue; }
    log "$1 $2 rc=$rc"; return $rc
  done
}
read_all() {
  for r in $RUNS; do
    while [ "$(jobs -rp | wc -l)" -ge 2 ]; do wait -n; done
    run_t $1 $r & sleep 90
  done
  wait
  for r in $RUNS; do [ -s $O/read/$1/$r/read_meta.json ] || { log "READ_FAIL $1 $r"; exit 1; }; done
}
log "RESUME (logs/v30_preview_tail.sh) from rat_12; mouse_chr19 already read"
for t in rat_12 cow_18 xtrop_10; do
  if [ ! -s data/p1_feat/$t.npy ]; then
    free=$(df --output=avail -BG / | tail -1 | tr -dc 0-9)
    [ "$free" -ge 150 ] || { log "STOP: only ${free} GB free before extracting $t"; exit 1; }
    rm -f data/p1_feat/$t.cover.json data/p1_shards/ont/${t}_t*.json
    NAMES=$t CHROM_ONLY=1 NCRNA_GPU_FRAC=0.28 taskset -c 16-31 nice -n 19 bash scripts/phase1/p1_feat_ont.sh >> $O/feat_$t$SFX.log 2>&1
    grep -q "ALL DONE fails=0" $O/feat_$t$SFX.log || { log "FEAT_FAIL $t"; exit 1; }
  fi
  [ -s data/p1_feat/$t.npy ] && echo "$(date -u +%FT%TZ) $STAGE $t $(sha256sum data/p1_feat/$t.npy | cut -c1-64)" >> $O/feat_sha.log
  read_all $t
  rm -f data/p1_feat/$t.npy data/p1_feat/$t.cover.json
  log "$t read, feature cache deleted"
done
nice -n 19 taskset -c 2-31 env PYTHONHASHSEED=0 OMP_NUM_THREADS=8 $PY scripts/phase1/v30/v30_score.py --config $C > $O/score$SFX.log 2>&1
log "score rc=$?"
log QUEUE_END
