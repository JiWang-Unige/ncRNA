#!/bin/bash
# v3.0 section 18.1 item 2 on ont: the failed seed-123 runs of H6 and HM are replaced by seed 789 (recipe of
# v30_train_queue.sh: run_formal_ont.sh, GPU gate incl. TE containers, yield / crash retries, fraction 0.30) ->
# failed-run check -> chr18 reads of the replacements (configs/v30_dev2.json, out_dir V30/dev) -> per-arm calibration
# and tau on the FAILCHECK.json run sets (V30/dev/<arm>).  If a 789 run fails the check, its arm keeps seeds 42 / 456
# (no further retraining).  Nothing outside human chr18 is read.  configs/v30.json is built and committed by hand after DEV2_END.
cd ~/ncrna-annot || exit 1
O=outputs/PHASE1/V30; PY=~/envs/generanno/bin/python; export PYTHONPATH=src
log() { echo "$(date -u +%FT%TZ) $*" >> $O/retrain_queue.log; }
log START
for arm in h6 hm; do
  n=v30$arm-g1b-789-ke
  if [ -s outputs/PHASE1/FORMAL/$n/ep12.pt ]; then log "$n skip (trained)"; continue; fi
  NCRNA_CPUS=2-31 NCRNA_NICE=19 NCRNA_GPU_FRAC=0.30 bash scripts/phase1/run_formal_ont.sh $n --seed 789 --data v30$arm --accum 2 >> logs/formal_$n.log 2>&1
  log "$n rc=$?"
  [ -s outputs/PHASE1/FORMAL/$n/ep12.pt ] || { log "$n not trained - stop"; exit 1; }
done
log TRAIN_END
$PY scripts/phase1/v30/v30_failcheck.py > $O/failcheck.log 2>&1 || { log FAILCHECK_FAIL; exit 1; }
log "failcheck: $(python3 -c "import json; d=json.load(open('$O/FAILCHECK.json')); print({a: v['runs'] for a, v in d['arms'].items()})")"
R2=$(python3 -c "import json; d=json.load(open('$O/FAILCHECK.json')); print(','.join(n for a in d['arms'].values() for n in a['runs'] if '-789-' in n))")
if [ -n "$R2" ]; then
  $PY scripts/phase1/v30/v30_dev.py config --commit "$(cat $O/COMMIT2)" --runs "$R2" --out configs/v30_dev2.json >> $O/dev_config.log 2>&1 || { log CONFIG2_FAIL; exit 1; }
  NCRNA_CPUS=2-31 NCRNA_NICE=19 NCRNA_GPU_FRAC=0.30 OMP_NUM_THREADS=4 bash scripts/phase1/run_formal_ont.sh v30dev2-chr18 --script p15_read.py \
    --config configs/v30_dev2.json --target chr18 >> logs/formal_v30dev2-chr18.log 2>&1 || { log READ2_FAIL; exit 1; }
  log "chr18 reads of $R2 ok"
fi
for a in h6 hm; do
  nice -n 19 taskset -c 2-31 env OMP_NUM_THREADS=8 $PY scripts/phase1/v30/v30_dev.py arm $a --stage final > $O/dev_final_$a.log 2>&1 || { log "CALIB_FAIL $a"; exit 1; }
  log "arm $a: $(tail -1 $O/dev_final_$a.log)"
done
log DEV2_END
