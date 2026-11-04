#!/bin/bash
# S5 mouse chr19 zero-shot on ont (PLAN_BATCH2.md B; card batch2-downloads 2026-11-04): codes -> features ->
# config -> chr18 reproduction of the three main runs -> mouse_chr19 reads -> scoring -> delete the feature cache.
# Same pipeline as P1.6 (p16_read.py unchanged).  GPU: gpu_gate.sh before every job, <= 2 jobs at a time
# (0.15 each), feature extraction capped at 0.28; CPU nice 19 on cores 2-31.
cd ~/ncrna-annot || exit 1
. scripts/phase1/gpu_gate.sh
PY=~/envs/generanno/bin/python
C=configs/p17_mouse.json; O=outputs/PHASE1/P17_MOUSE; mkdir -p $O data/mouse
RUNS="g1b-42-ke g1b-123-ke g1b-456-ke"
log() { echo "$(date -u +%FT%TZ) $*" >> $O/queue.log; }
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
log START
if [ ! -s data/mouse/chr19.fa ]; then
  zcat data/batch2/mouse/mm39_chr19.fa.gz > data/mouse/chr19.fa && (cd data/mouse && md5sum chr19.fa > MD5SUMS)
fi
[ -s data/p1/mouse_chr19_codes.npy ] || nice -n 19 $PY scripts/phase1/prep_panel.py --only mouse_chr19 >> $O/prep.log 2>&1 || { log PREP_FAIL; exit 1; }
log "prep ok"
NAMES=mouse_chr19 CHROM_ONLY=1 NCRNA_GPU_FRAC=0.28 bash scripts/phase1/p1_feat_ont.sh >> $O/feat.log 2>&1
grep -q "ALL DONE fails=0" $O/feat.log || { log FEAT_FAIL; exit 1; }
log "features ok"
$PY scripts/phase1/p16_config.py --commit "$(cat $O/COMMIT 2>/dev/null || echo unknown)" --previous configs/p15_v19b.json \
  --out $C --only mouse_chr19 --out-dir $O >> $O/config.log 2>&1 || { log CONFIG_FAIL; exit 1; }
for t in chr18 mouse_chr19; do
  for r in $RUNS; do
    while [ "$(jobs -rp | wc -l)" -ge 2 ]; do wait -n; done
    run_t $t $r & sleep 90
  done
  wait
  if [ $t = chr18 ]; then
    nice -n 19 $PY scripts/phase1/p16_read.py --config $C --target chr18 >> $O/read_chr18_all.log 2>&1
    [ -e $O/chr18/REPRO_OK ] || { log REPRO_FAIL; exit 1; }
    log REPRO_OK
  fi
done
nice -n 19 taskset -c 2-31 env PYTHONHASHSEED=0 $PY scripts/phase1/p16_score.py --config $C > $O/score.log 2>&1
log "score rc=$?"
[ -s $O/mouse_chr19/SCORE.json ] && rm -f data/p1_feat/mouse_chr19.npy && log "feature cache deleted"
log QUEUE_END
