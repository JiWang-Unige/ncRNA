#!/bin/bash
# v2.0 chr18 development on ont (protocol section 17, detail 5), after the retraining queue: dev config -> chr18 reads
# of the three retrained runs (p15_read.py through run_formal_ont.sh: GPU gate, lock, yield retries) -> calibration
# and tau (CPU).  The hexamer table is built beforehand (v20_calib.py hexamer, CPU, training chromosomes only).
cd ~/ncrna-annot || exit 1
O=outputs/PHASE1/V20; PY=~/envs/generanno/bin/python; export PYTHONPATH=src
log() { echo "$(date -u +%FT%TZ) $*" >> $O/dev_queue.log; }
until grep -qE "TRAIN_END" $O/train_queue.log 2>/dev/null; do sleep 600; done
for s in 42 123 456; do [ -s outputs/PHASE1/FORMAL/v20-g1b-$s-ke/ep12.pt ] || { log "v20-g1b-$s-ke not trained - stop"; exit 1; }; done
log START
$PY scripts/phase1/v20/v20_config.py dev --commit "$(cat $O/COMMIT)" >> $O/dev_config.log 2>&1 || { log CONFIG_FAIL; exit 1; }
NCRNA_CPUS=2-31 NCRNA_NICE=19 NCRNA_GPU_FRAC=0.30 OMP_NUM_THREADS=4 bash scripts/phase1/run_formal_ont.sh v20dev-chr18 --script p15_read.py \
  --config configs/v20_dev.json --target chr18 >> logs/formal_v20dev-chr18.log 2>&1 || { log READ_FAIL; exit 1; }
log "chr18 reads ok"
[ -s data/p1v20/hexamer_table.json ] || nice -n 19 $PY scripts/phase1/v20/v20_calib.py hexamer >> $O/hexamer.log 2>&1 || { log HEXAMER_FAIL; exit 1; }
nice -n 19 taskset -c 2-31 env OMP_NUM_THREADS=8 $PY scripts/phase1/v20/v20_calib.py calibrate > $O/calibrate.log 2>&1 || { log CALIB_FAIL; exit 1; }
log DEV_END
