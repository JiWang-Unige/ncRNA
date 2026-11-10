#!/bin/bash
# v3.0 test reads on ont (protocol section 18, execution details 9-10).  Started by hand once configs/v30.json and
# v30_score.py are committed (the commit sha goes into outputs/PHASE1/V30/CONFIG_COMMIT).  chr18 reproduction of all
# nine runs (REPRO_OK) -> group B (feature caches kept from P1.6) -> mouse_chr19 and group C (features extracted,
# read by the nine runs, deleted) -> v30_score.py.  p16_read.py unchanged; GPU through gpu_gate.sh, <= 2 reads at a time
# (0.15 each), feature extraction 0.28; CPU nice 19 on cores 2-31.  Stops before an extraction when < 150 GB are free.
# Section 18.1: v30_test_queue.sh preview -> configs/v30_preview.json (7 runs; CONFIG_COMMIT_preview, test_queue_preview.log);
# v30_test_queue.sh [final] -> configs/v30.json (runs already read under the preview config are skipped by p16_read.py).
# The sha256 of every extracted test feature cache is logged (feat_sha.log) before it is deleted.
cd ~/ncrna-annot || exit 1
. scripts/phase1/gpu_gate.sh
PY=~/envs/generanno/bin/python
STAGE=${1:-final}; O=outputs/PHASE1/V30
if [ "$STAGE" = preview ]; then C=configs/v30_preview.json; SFX=_preview; else C=configs/v30.json; SFX=; fi
log() { echo "$(date -u +%FT%TZ) $*" >> $O/test_queue$SFX.log; }
if [ "$STAGE" = final ] && [ -e $O/test_queue_preview.log ]; then  # chr18/REPRO_OK is per config: never overlap
  until grep -qE "QUEUE_END|FAIL|STOP" $O/test_queue_preview.log; do sleep 600; done
fi
[ -s $O/CONFIG_COMMIT$SFX ] || { log "no CONFIG_COMMIT$SFX - $C not committed"; exit 1; }
RUNS=$(python3 -c "import json; print(' '.join(json.load(open('$C'))['targets']['chr18']['runs']))")
run_t() {  # <target> <run>
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
read_all() {  # <target>: the nine runs, two at a time
  for r in $RUNS; do
    while [ "$(jobs -rp | wc -l)" -ge 2 ]; do wait -n; done
    run_t $1 $r & sleep 90
  done
  wait
  for r in $RUNS; do [ -s $O/read/$1/$r/read_meta.json ] || { log "READ_FAIL $1 $r"; exit 1; }; done
}
log "START $STAGE config $(sha256sum $C | cut -c1-12) commit $(cat $O/CONFIG_COMMIT$SFX)"
read_all chr18
nice -n 19 $PY scripts/phase1/p16_read.py --config $C --target chr18 >> $O/read_chr18_all$SFX.log 2>&1
[ -e $O/read/chr18/REPRO_OK ] || { log REPRO_FAIL; exit 1; }
log REPRO_OK
for t in chicken_5 zebrafish_22 fly_2L worm_I; do read_all $t; log "$t read"; done
for t in mouse_chr19 rat_12 cow_18 xtrop_10; do
  if [ ! -s data/p1_feat/$t.npy ]; then
    free=$(df --output=avail -BG / | tail -1 | tr -dc 0-9)
    [ "$free" -ge 150 ] || { log "STOP: only ${free} GB free before extracting $t"; exit 1; }
    # a cache deleted after an earlier read (mouse_chr19, S5) leaves its cover file and shard markers behind
    rm -f data/p1_feat/$t.cover.json data/p1_shards/ont/${t}_t*.json
    NAMES=$t CHROM_ONLY=1 NCRNA_GPU_FRAC=0.28 bash scripts/phase1/p1_feat_ont.sh >> $O/feat_$t$SFX.log 2>&1
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
