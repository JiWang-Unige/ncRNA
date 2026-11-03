#!/bin/bash
# chr18 Vipsania / Tiberius coding-evidence diagnostic on ont (card vipsania-ref option 1; REVIEW_VIPSANIA.md section 5):
# Vipsania chr18 posterior -> Vipsania chr18 / chr20 annotation -> Tiberius chr18 (batch 1) -> CPU analysis (vip_diag.py).
# Every GPU step waits on gpu_gate.sh and is capped at 7000 MiB; per-process GPU memory is sampled every 30 s into
# logs/gpu_mem_<step>.log.  CPU: nice 19, cores 2-31.
cd ~/ncrna-annot || exit 1
. scripts/phase1/gpu_gate.sh
D=outputs/PHASE1/VIPSANIA
VP=~/envs/vipsania_ncrna/bin/python
mkdir -p $D/logs $D/chr18 $D/chr20
log() { echo "$(date -u +%FT%TZ) $*" >> $D/queue.log; }
step() {  # name done-file command...
  local name=$1 done=$2 rc mon
  shift 2
  if [ -s "$done" ]; then log "$name skip"; return 0; fi
  gpu_wait >> $D/queue.log
  log "$name start"
  (while sleep 30; do echo "$(date -u +%T) $(nvidia-smi --query-compute-apps=pid,used_memory --format=csv,noheader | tr '\n' ' ')"; done) \
    >> $D/logs/gpu_mem_$name.log &
  mon=$!
  nice -n 19 taskset -c 2-31 env OMP_NUM_THREADS=4 TF_NUM_INTRAOP_THREADS=4 TF_NUM_INTEROP_THREADS=2 "$@" > $D/logs/$name.log 2>&1
  rc=$?
  kill $mon
  log "$name rc=$rc"
  return $rc
}
log START
step post18 $D/chr18/post.json $VP scripts/phase1/vipsania/vip_run.py post --chrom chr18 || exit 1
step ann18 $D/chr18/annotate.json $VP scripts/phase1/vipsania/vip_run.py annotate --chrom chr18 || exit 1
step ann20 $D/chr20/annotate.json $VP scripts/phase1/vipsania/vip_run.py annotate --chrom chr20 || exit 1
step tib18 $D/chr18/tiberius.ok bash -c "cd ext/tiberius && ~/envs/tiberius/bin/python ~/ncrna-annot/scripts/phase1/vipsania/tf_cap.py 7000 \
  tiberius.py --genome ~/ncrna-annot/data/human/chr18.fa --model_cfg mammalia_nosofttmasking_v2 --batch_size 1 \
  --out ~/ncrna-annot/$D/chr18/tiberius.gtf && touch ~/ncrna-annot/$D/chr18/tiberius.ok" || exit 1
log GPU_DONE
nice -n 19 taskset -c 2-31 env OMP_NUM_THREADS=8 ~/envs/generanno/bin/python scripts/phase1/vipsania/vip_diag.py > $D/logs/diag.log 2>&1
log "diag rc=$? QUEUE_END"
