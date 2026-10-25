#!/bin/bash
# ncRNA GPU gate on ont (user rule 2026-10-26: "只要不阻塞TE项目，随意推进任务").
# Source it:  . scripts/phase1/gpu_gate.sh ;  gpu_ok && echo go ;  gpu_wait   (blocks, logs every 10 min)
# ncRNA may run when every GPU process belongs to the user (ji) and is NOT a TE-project process, and enough memory
# is free.  It waits for: any TE process (cwd/args under ab-initio-TE or TE-GLM, env te_benchmark/earlgrey/
# repeatmasker), any other user's process (shared lab server), or free memory < NCRNA_NEED_MB (default 10000).
: "${NCRNA_NEED_MB:=10000}"
TE_PAT='ab-initio-TE|TE-GLM|te_benchmark|earlgrey|repeatmasker|RepeatMasker'
gpu_why_busy() {  # prints the reason and returns 0 when busy; returns 1 when ncRNA may run
  local me pid user where free
  me=$(id -un)
  for pid in $(nvidia-smi --query-compute-apps=pid --format=csv,noheader 2>/dev/null); do
    user=$(ps -o user= -p "$pid" 2>/dev/null | tr -d ' ')
    [ -z "$user" ] && continue
    if [ "$user" != "$me" ]; then echo "other user $user (pid $pid)"; return 0; fi
    where="$(readlink /proc/$pid/cwd 2>/dev/null) $(ps -o args= -p "$pid" 2>/dev/null)"
    if echo "$where" | grep -qE "$TE_PAT"; then echo "TE process pid $pid: $(echo "$where" | cut -c1-120)"; return 0; fi
  done
  free=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)
  if [ "$free" -lt "$NCRNA_NEED_MB" ]; then echo "only ${free} MiB free (< $NCRNA_NEED_MB)"; return 0; fi
  return 1
}
gpu_ok() { ! gpu_why_busy > /dev/null; }
gpu_wait() {
  local n=0 why
  while why=$(gpu_why_busy); do
    [ $((n % 10)) -eq 0 ] && echo "$(date -u +%Y-%m-%dT%H:%MZ) GPU gate: $why - waiting"
    n=$((n + 1)); sleep 60
  done
}
