#!/bin/bash
# Wrapper for train_formal_ont.py on ont: wait for the GPU gate, run, and on exit code 75 (yielded to a TE / other
# user's job) wait and restart from last.pt.  One ncRNA GPU run at a time (lock).  PYTHONFAULTHANDLER: a crash
# in native code (e.g. the segfault of 2026-10-29) still leaves a Python traceback in the log.
# Random crashes on ont (9 unrelated programs segfaulted 10-23..29, nearly all on CPU 0/1 = core 0 of the i9-13900K):
# a run killed by SIGSEGV / SIGABRT / SIGBUS (rc 139 / 134 / 135) is retried up to NCRNA_CRASH_RETRIES (default 2) times
# (training resumes from last.pt as after a yield); NCRNA_CPUS (e.g. 2-31) pins the run to those CPUs with taskset.
#   nohup bash scripts/phase1/run_formal_ont.sh <name> [train_formal_ont.py args...] >> logs/formal_<name>.log 2>&1 &
set -u
cd ${NCRNA_ROOT:-~/ncrna-annot}
NAME=$1; shift
. scripts/phase1/gpu_gate.sh
LOCK=logs/formal_gpu.lock
log() { echo "$(date -u +%Y-%m-%dT%H:%MZ) $*"; }
while [ -e $LOCK ] && kill -0 "$(cat $LOCK)" 2>/dev/null; do log "another ncRNA GPU run (pid $(cat $LOCK)) active - waiting"; sleep 300; done
echo $$ > $LOCK; trap 'rm -f $LOCK' EXIT
PIN=""; [ -n "${NCRNA_CPUS:-}" ] && PIN="taskset -c $NCRNA_CPUS"
crashes=0
while true; do
  gpu_wait
  log "start $NAME $*"
  PYTHONFAULTHANDLER=1 PYTHONPATH=src $PIN nice -n ${NCRNA_NICE:-0} ${NCRNA_PY:-$HOME/envs/generanno/bin/python} -u scripts/phase1/train_formal_ont.py --name "$NAME" "$@"
  rc=$?
  if [ $rc -eq 0 ]; then log "DONE $NAME"; exit 0; fi
  if [ $rc -eq 75 ]; then log "yielded GPU ($NAME); will resume from last.pt"; sleep 120; continue; fi
  if { [ $rc -eq 139 ] || [ $rc -eq 134 ] || [ $rc -eq 135 ]; } && [ $crashes -lt ${NCRNA_CRASH_RETRIES:-2} ]; then
    crashes=$((crashes + 1)); log "CRASHED $NAME rc=$rc - retry $crashes/${NCRNA_CRASH_RETRIES:-2} (resumes from last.pt)"; sleep 60; continue
  fi
  log "FAILED $NAME rc=$rc"; exit $rc
done
