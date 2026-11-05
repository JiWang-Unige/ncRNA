#!/bin/bash
# v3.0 chr18 development on ont (protocol section 18, execution detail 8), after the training queue: dev config ->
# human chr18 reads of the six new runs (p15_read.py through run_formal_ont.sh: GPU gate, lock, yield retries) ->
# per-arm calibration and tau (CPU, v2.0 procedure and section 17.1 rule).  Nothing outside human chr18 is read.
cd ~/ncrna-annot || exit 1
O=outputs/PHASE1/V30; PY=~/envs/generanno/bin/python; export PYTHONPATH=src
log() { echo "$(date -u +%FT%TZ) $*" >> $O/dev_queue.log; }
until grep -qE "TRAIN_END|prep failed" $O/train_queue.log 2>/dev/null; do sleep 600; done
for a in h6 hm; do for s in 42 123 456; do
  [ -s outputs/PHASE1/FORMAL/v30$a-g1b-$s-ke/ep12.pt ] || { log "v30$a-g1b-$s-ke not trained - stop"; exit 1; }
done; done
log START
$PY scripts/phase1/v30/v30_dev.py config --commit "$(cat $O/COMMIT)" >> $O/dev_config.log 2>&1 || { log CONFIG_FAIL; exit 1; }
NCRNA_CPUS=2-31 NCRNA_NICE=19 NCRNA_GPU_FRAC=0.30 OMP_NUM_THREADS=4 bash scripts/phase1/run_formal_ont.sh v30dev-chr18 --script p15_read.py \
  --config configs/v30_dev.json --target chr18 >> logs/formal_v30dev-chr18.log 2>&1 || { log READ_FAIL; exit 1; }
log "chr18 reads ok"
for a in h6 hm; do
  nice -n 19 taskset -c 2-31 env OMP_NUM_THREADS=8 $PY scripts/phase1/v30/v30_dev.py arm $a > $O/dev_$a.log 2>&1 || { log "CALIB_FAIL $a"; exit 1; }
  log "arm $a: $(tail -1 $O/dev_$a.log)"
done
log DEV_END
