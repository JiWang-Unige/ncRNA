#!/bin/bash
# Unified GFF3 (docs/phase1_results/PLAN_UGFF.md): Tiberius coding layer for chr18 (development) and chr9 (held-out;
# Tiberius is an external tool run with its published settings, nothing is tuned on chr9).  Same call as the chr20
# baseline (baselines_queue.sh): mammalia_nosofttmasking_v2, batch 2.  GPU via gpu_gate.sh; CPU nice 19, cores 2-31.
cd ~/ncrna-annot || exit 1
. scripts/phase1/gpu_gate.sh
D=outputs/PHASE1/BASELINES/tiberius
log() { echo "$(date -u +%FT%TZ) $*" >> outputs/PHASE1/UGFF/queue.log; }
mkdir -p outputs/PHASE1/UGFF
log START
for c in chr18 chr9; do
  [ -s $D/$c.gtf ] && { log "$c exists"; continue; }
  gpu_wait
  (cd ext/tiberius && nice -n 19 taskset -c 2-31 env TF_FORCE_GPU_ALLOW_GROWTH=true OMP_NUM_THREADS=4 \
    ~/envs/tiberius/bin/python tiberius.py --genome ~/ncrna-annot/data/human/$c.fa \
    --model_cfg mammalia_nosofttmasking_v2 --batch_size 2 --out ~/ncrna-annot/$D/$c.gtf) > $D/run_$c.log 2>&1
  log "tiberius $c rc=$?"
done
log DONE
