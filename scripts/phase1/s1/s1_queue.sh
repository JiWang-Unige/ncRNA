#!/bin/bash
# S1 queue on ont (PLAN_S1.md, card subjournal-portfolio 2026-11-04): after the S1b queue ends, chr20 splice scores of
# SpliceAI, the own S10k net and Vipsania (one at a time, gpu_gate.sh before each, <= 30% GPU), then the CPU analysis.
cd ~/ncrna-annot || exit 1
. scripts/phase1/gpu_gate.sh
D=outputs/PHASE1/S1
PY=~/envs/generanno/bin/python
VP=~/envs/vipsania_ncrna/bin/python
mkdir -p $D/logs
log() { echo "$(date -u +%FT%TZ) $*" >> $D/queue.log; }
until grep -q QUEUE_END outputs/PHASE1/S1B/queue.log 2>/dev/null; do sleep 300; done
step() {  # name done-file command...
  local name=$1 done=$2 rc
  shift 2
  if [ -s "$done" ]; then log "$name skip"; return 0; fi
  while true; do
    gpu_wait >> $D/queue.log
    log "$name start"
    nice -n 19 taskset -c 2-31 env OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 TF_NUM_INTRAOP_THREADS=4 TF_NUM_INTEROP_THREADS=2 "$@" > $D/logs/$name.log 2>&1
    rc=$?
    [ $rc -eq 75 ] && { log "$name gate busy, retry"; sleep 120; continue; }
    log "$name rc=$rc"
    return $rc
  done
}
log START
step sai20 $D/sai_chr20.npy env NCRNA_GPU_FRAC=0.25 $PY scripts/phase1/splice_diag/sai_score.py --codes data/p1/chr20_codes.npy \
  --truth data/human/gencode.v50.annotation.gtf.gz --chrom chr20 --out $D/sai_chr20.npy || exit 1
step own20 $D/own_chr20.npy env NCRNA_GPU_FRAC=0.25 $PY scripts/phase1/splice_diag/own_splice.py score --codes data/p1/chr20_codes.npy \
  --truth data/human/gencode.v50.annotation.gtf.gz --chrom chr20 --out $D/own_chr20.npy || exit 1
step vip20 outputs/PHASE1/VIPSANIA/chr20/post.json env VIP_GPU_MB=6800 $VP scripts/phase1/vipsania/vip_run.py post --chrom chr20 || exit 1
log GPU_DONE
nice -n 19 taskset -c 2-31 env OMP_NUM_THREADS=8 $PY scripts/phase1/s1/s1_analyze.py \
  --scores spliceai=$D/sai_chr20.npy own=$D/own_chr20.npy vipsania=outputs/PHASE1/VIPSANIA/chr20/post6.npy:2 \
  --out $D/S1.json > $D/logs/analyze.log 2>&1
log "analyze rc=$? QUEUE_END"
