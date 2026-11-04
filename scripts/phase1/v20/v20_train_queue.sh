#!/bin/bash
# v2.0 retraining on ont (protocol section 17, detail 4): after the data prep, the three main G1B runs on the v2.0
# window set, one at a time through run_formal_ont.sh (GPU gate, lock, yield / crash retries), same settings as
# g1b-*-ke (v50uev labels, 12 epochs, 2 x 4 accumulation, every epoch kept, GPU fraction 0.30).
cd ~/ncrna-annot || exit 1
O=outputs/PHASE1/V20; mkdir -p $O logs
log() { echo "$(date -u +%FT%TZ) $*" >> $O/train_queue.log; }
until grep -qE "PREP_END|FAIL" $O/prep_queue.log 2>/dev/null; do sleep 300; done
grep -q PREP_END $O/prep_queue.log || { log "prep failed - not training"; exit 1; }
log START
for s in 42 123 456; do
  n=v20-g1b-$s-ke
  if [ -s outputs/PHASE1/FORMAL/$n/ep12.pt ]; then log "$n skip (trained)"; continue; fi
  NCRNA_CPUS=2-31 NCRNA_NICE=19 NCRNA_GPU_FRAC=0.30 bash scripts/phase1/run_formal_ont.sh $n --seed $s --data v20 --accum 2 >> logs/formal_$n.log 2>&1
  log "$n rc=$?"
done
log TRAIN_END
