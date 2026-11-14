#!/bin/bash
# PLAN_DEDEND (frozen 8deaaff) on ont: Puffin on chr8/9/18/20, PolyaID on chr18/20, then metrics (CPU).
#   nohup bash scripts/phase1/dedend/dedend_queue.sh > outputs/PHASE1/DEDEND/queue.log 2>&1 &
set -u
cd ${NCRNA_ROOT:-~/ncrna-annot}
export NCRNA_CPUS=${NCRNA_CPUS:-2-31} NCRNA_NICE=${NCRNA_NICE:-10}
R=scripts/phase1/run_formal_ont.sh
log() { echo "$(date -u +%Y-%m-%dT%H:%MZ) $*"; }
NCRNA_PY=~/envs/puffin/bin/python bash $R puffin --script dedend/puffin_run.py --chroms chr8,chr9,chr18,chr20 || exit 1
NCRNA_PY=~/envs/polya/bin/python bash $R polya --script dedend/polya_run.py --chroms chr18,chr20 || exit 1
PYTHONPATH=src nice -n 10 ~/envs/generanno/bin/python scripts/phase1/dedend/dedend_metrics.py || exit 1
log "QUEUE_END"
