#!/bin/bash
# M1 weeks 2-4 (docs/m1/PLAN_M1_W24.md): DNA-only v2.0 reads of the test parts under configs/m1_w24.json, plus Tiberius on
# chr12 / chr13.  chr18 reproduction first -> Tiberius chr12, chr13 -> per part (w24_chr12a/b, w24_chr13a/b, w24_chr20):
# extract features (from scratch, cores 2-31, one retry), read the three runs, delete the feature cache.
# GPU through gpu_gate.sh (TE first); CPU nice 19 on cores 2-31.  No WTC11 data and no truth is read here.
cd ~/ncrna-annot || exit 1
. scripts/phase1/gpu_gate.sh
PY=~/envs/generanno/bin/python
C=configs/m1_w24.json; O=outputs/M1/W24; mkdir -p $O/tiberius $O/read
log() { echo "$(date -u +%FT%TZ) $*" >> $O/queue.log; }
RUNS="v20-g1b-123-ke v20-g1b-42-ke v20-g1b-456-ke"
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
read_all() {
  for r in $RUNS; do
    while [ "$(jobs -rp | wc -l)" -ge 2 ]; do wait -n; done
    run_t $1 $r & sleep 90
  done
  wait
  for r in $RUNS; do [ -s $O/read/$1/$r/read_meta.json ] || { log "READ_FAIL $1 $r"; exit 1; }; done
}
log "START config $(sha256sum $C | cut -c1-12)"
if [ ! -e $O/read/chr18/REPRO_OK ]; then
  read_all chr18
  nice -n 19 $PY scripts/phase1/p16_read.py --config $C --target chr18 >> $O/read_chr18_all.log 2>&1
  [ -e $O/read/chr18/REPRO_OK ] || { log REPRO_FAIL; exit 1; }
  rm -f $O/read/chr18/*/post.npy   # identical to outputs/PHASE1/V20/read/chr18 (sha256 checked)
  log REPRO_OK
fi
for c in chr12 chr13; do
  if [ ! -s $O/tiberius/$c.gtf ]; then
    gpu_wait
    (cd ext/tiberius && nice -n 19 taskset -c 2-31 env TF_FORCE_GPU_ALLOW_GROWTH=true OMP_NUM_THREADS=4 \
      ~/envs/tiberius/bin/python tiberius.py --genome ~/ncrna-annot/data/human/$c.fa \
      --model_cfg mammalia_nosofttmasking_v2 --batch_size 2 --out ~/ncrna-annot/$O/tiberius/$c.gtf) > $O/tiberius/$c.run.log 2>&1
    log "tiberius $c rc=$?"
    [ -s $O/tiberius/$c.gtf ] || { log "TIBERIUS_FAIL $c"; exit 1; }
  fi
done
for t in w24_chr20 w24_chr12a w24_chr12b w24_chr13a w24_chr13b; do
  done_all=1; for r in $RUNS; do [ -s $O/read/$t/$r/read_meta.json ] || done_all=0; done
  [ $done_all = 1 ] && { log "$t already read"; continue; }
  if [ ! -e $O/feat_$t.complete ]; then   # a partial cache is never reused
    ok=0
    for attempt in 1 2; do
      rm -f data/p1_feat/$t.npy data/p1_feat/$t.cover.json data/p1_shards/ont/${t}_t*.json
      free=$(df --output=avail -BG / | tail -1 | tr -dc 0-9)
      [ "$free" -ge 120 ] || { log "STOP: only ${free} GB free before extracting $t"; exit 1; }
      taskset -c 2-31 env NAMES=$t CHROM_ONLY=1 NCRNA_GPU_FRAC=0.28 bash scripts/phase1/p1_feat_ont.sh >> $O/feat_$t.log 2>&1
      grep "ALL DONE" $O/feat_$t.log | tail -1 | grep -q "fails=0" && { ok=1; break; }
      log "feature extraction of $t failed (attempt $attempt)"
    done
    [ $ok = 1 ] || { log "FEAT_FAIL $t"; exit 1; }
    touch $O/feat_$t.complete
  fi
  echo "$(date -u +%FT%TZ) $t $(sha256sum data/p1_feat/$t.npy | cut -c1-64)" >> $O/feat_sha.log
  read_all $t
  rm -f data/p1_feat/$t.npy data/p1_feat/$t.cover.json
  log "$t read, feature cache deleted"
done
log QUEUE_END
