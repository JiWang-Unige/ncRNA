#!/bin/bash
# PLAN_END stages 0 and 1 on ont, one after another through run_formal_ont.sh (GPU gate, lock, yield / resume).
#   nohup bash scripts/phase1/ends/end_queue.sh > outputs/PHASE1/ENDS/queue.log 2>&1 &
cd ${NCRNA_ROOT:-~/ncrna-annot}
export NCRNA_CPUS=${NCRNA_CPUS:-2-31} NCRNA_NICE=${NCRNA_NICE:-10}
R=scripts/phase1/run_formal_ont.sh
bash $R stage0 --script ends/dump_bnd.py || exit 1
bash $R end-ann-42 --script ends/train_end_probe.py --arm ann --seed 42 || exit 1
bash $R end-ev-42 --script ends/train_end_probe.py --arm ev --seed 42 || exit 1
echo "QUEUE_END $(date -u +%FT%TZ)"
