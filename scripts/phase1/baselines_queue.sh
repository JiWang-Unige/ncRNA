#!/bin/bash
# chr20 baselines on ont (card baselines-chr20; docs/phase1_results/PREREG_BASELINES.md): after the own splice net has
# finished training -> Tiberius (mammalia_nosofttmasking_v2, batch 2) -> SegmentNT inference -> decoding + scoring.
# One read of chr20 each.  GPU via gpu_gate.sh (exit 75 = busy -> retry); CPU nice 19, cores 2-31.
cd ~/ncrna-annot || exit 1
. scripts/phase1/gpu_gate.sh
D=outputs/PHASE1/BASELINES
mkdir -p $D/tiberius $D/segnt
log() { echo "$(date -u +%FT%TZ) $*" >> $D/queue.log; }
log START
until [ -f outputs/PHASE1/SPLICE_OWN/s2k/DONE ] && [ -f outputs/PHASE1/SPLICE_OWN/s10k/DONE ]; do sleep 300; done
log "splice training done"
if [ ! -s $D/tiberius/chr20.gtf ]; then
  gpu_wait
  (cd ext/tiberius && nice -n 19 taskset -c 2-31 env TF_FORCE_GPU_ALLOW_GROWTH=true OMP_NUM_THREADS=4 \
    ~/envs/tiberius/bin/python tiberius.py --genome ~/ncrna-annot/data/human/chr20.fa \
    --model_cfg mammalia_nosofttmasking_v2 --batch_size 2 --out ~/ncrna-annot/$D/tiberius/chr20.gtf) \
    > $D/tiberius/run.log 2>&1
  log "tiberius rc=$?"
fi
nice -n 19 ~/envs/generanno/bin/python scripts/phase1/tiberius_score.py > $D/tiberius/score.log 2>&1
log "tiberius score rc=$?"
while [ ! -f $D/segnt/infer.json ]; do
  gpu_wait
  nice -n 19 taskset -c 2-31 env OMP_NUM_THREADS=2 HF_HUB_OFFLINE=1 ~/envs/generanno/bin/python scripts/phase1/segnt_chr20.py infer \
    >> $D/segnt/infer.log 2>&1
  rc=$?
  [ $rc -eq 75 ] && { sleep 120; continue; }
  log "segnt infer rc=$rc"; [ $rc -eq 0 ] || exit 1
done
nice -n 19 taskset -c 2-31 ~/envs/generanno/bin/python scripts/phase1/segnt_chr20.py decode > $D/segnt/decode.log 2>&1
log "segnt decode rc=$? QUEUE_END"
