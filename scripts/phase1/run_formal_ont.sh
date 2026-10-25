#!/bin/bash
# Wrapper for train_formal_ont.py on ont: wait for the GPU gate, run, and on exit code 75 (yielded to a TE / other
# user's job) wait and restart from last.pt.  One ncRNA GPU run at a time (lock).
#   nohup bash scripts/phase1/run_formal_ont.sh <name> [train_formal_ont.py args...] >> logs/formal_<name>.log 2>&1 &
set -u
cd ~/ncrna-annot
NAME=$1; shift
. scripts/phase1/gpu_gate.sh
LOCK=logs/formal_gpu.lock
log() { echo "$(date -u +%Y-%m-%dT%H:%MZ) $*"; }
while [ -e $LOCK ] && kill -0 "$(cat $LOCK)" 2>/dev/null; do log "another ncRNA GPU run (pid $(cat $LOCK)) active - waiting"; sleep 300; done
echo $$ > $LOCK; trap 'rm -f $LOCK' EXIT
while true; do
  gpu_wait
  log "start $NAME $*"
  PYTHONPATH=src ~/envs/generanno/bin/python -u scripts/phase1/train_formal_ont.py --name "$NAME" "$@"
  rc=$?
  if [ $rc -eq 0 ]; then log "DONE $NAME"; exit 0; fi
  if [ $rc -eq 75 ]; then log "yielded GPU ($NAME); will resume from last.pt"; sleep 120; continue; fi
  log "FAILED $NAME rc=$rc"; exit $rc
done
