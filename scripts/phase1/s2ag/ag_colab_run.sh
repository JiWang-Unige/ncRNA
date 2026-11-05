#!/bin/bash
# S2-AG (PLAN_S2_AG.md): run ag_colab_job.py in a Colab VM and fetch its outputs to ont as each part completes.
#   ag_colab_run.sh smoke|full|auto|cputest <deadline_s> <gpu[,gpu...]|CPU> <cu_cap> <chrom>...
# auto = smoke on the first chromosome, then, in the same session and only if the smoke check passes (fast path equal
# to predict_sequence), the full run on the longest prefix of the chromosome list (given in priority order) whose
# projected CU fits in what is left of cu_cap; the first chromosome must fit, otherwise stop.
# cputest = the smoke steps on a CPU runtime with 16,384 bp windows (install, model load, track mapping, fetch).
# Assignment: tries each listed accelerator in turn; with ASSIGN_HOURS=h (default 0: one round) repeats every 20 min
# for up to h hours (premium GPUs are often out of capacity; a failed attempt costs nothing).
# Transport: everything after assignment goes over ssh / scp through the Colab CLI's ssh bridge (~/.ssh/colab.conf,
# host colab-<session>, key ~/.ssh/colab_ed25519).  Not `colab exec` (a poll hung for 42 min on 11-05 while the VM
# kept billing) and not `colab upload / download` (Jupyter contents API, refuses files of ~100 MB and more).
# Every remote call has a hard timeout; 3 failed polls in a row end the run; a watchdog kills the run at
# deadline + 15 min.  Caps, counted from the moment a session is READY: wall clock (deadline_s) and CU (cu_cap).
# On exit, whatever happened (also on TERM / INT / HUP): stop the session, unassign orphans, log the usage.
set -u
MODE=$1; DEADLINE_S=$2; GPU=$3; CAP=$4; shift 4; CHROMS="$*"
cd ~/ncrna-annot || exit 1
S=agrun
C="timeout 900 $HOME/bin/colab_long --auth oauth2"
P=$HOME/.local/share/uv/tools/google-colab-cli/bin/python
D=outputs/PHASE1/S2AG/$MODE
mkdir -p $D/raw
TMP=$(mktemp -d)
SSH="ssh -F $HOME/.ssh/colab.conf colab-$S"
SCP="scp -q -F $HOME/.ssh/colab.conf"
T0=$SECONDS
WD=
log() { echo "$(date -u +%FT%TZ) $*" | tee -a $D/run.log; }
usage_line() { timeout 120 $HOME/bin/colab_long --auth oauth2 usage 2>&1; }
bal() { usage_line | sed -n 's/.*Current balance: \([0-9.]*\).*/\1/p' | head -1; }
rate() { usage_line | sed -n 's/.*Usage rate: \([0-9.]*\).*/\1/p' | head -1; }
spent() { python3 -c "print(round($B0 - ${1:-$B0}, 2))"; }
B0=$(bal)
[ -n "$B0" ] || { echo "ABORT: no balance"; exit 1; }
cleanup() {
  trap - EXIT TERM INT HUP
  [ -n "$WD" ] && kill $WD 2>/dev/null
  for f in job.log job.out job.1.log job.1.out; do timeout 60 $SCP colab-$S:/content/$f $D/$f >/dev/null 2>&1; done
  timeout 120 $HOME/bin/colab_long --auth oauth2 stop -s $S 2>&1 | tail -1
  timeout 300 $P $HOME/bin/colab_cleanup.py --orphans
  log "after: $(usage_line | grep -E 'balance|rate|assign' | tr '\n' ' ')"
  log "END wall_s $((SECONDS - T0)) spent_cu $(spent "$(bal)")"
  rm -rf $TMP
}
trap cleanup EXIT
trap 'log "signal: exiting"; exit 1' TERM INT HUP
log "START mode=$MODE gpu=$GPU cap=$CAP chroms=$CHROMS balance=$B0"

W=data/alphagenome/hf_all_folds.tar
[ -s $W ] || tar -cf $W -C data/alphagenome/hf_all_folds _CHECKPOINT_METADATA _METADATA manifest.ocdbt d ocdbt.process_0
WSHA=$(sha256sum $W | cut -d' ' -f1)
job_args() {  # mode chrom...
  python3 - "$WSHA" "$@" <<'EOF'
import json, sys
c = json.load(open("configs/s2ag.json"))
a = {"mode": sys.argv[2], "chroms": sys.argv[3:], "ag_commit": c["ag_commit"], "ckpt_tar_sha256": sys.argv[1],
     "codes_sha256": {k: c["codes_sha256"][k] for k in sys.argv[3:]}}
if a["mode"] == "cputest":
    a.update(cpu=True, window=16384, smoke_windows=1)
print(json.dumps(a))
EOF
}

ok=0
A_END=$((SECONDS + ${ASSIGN_HOURS:-0} * 3600))
while [ $ok = 0 ]; do
  for g in ${GPU//,/ }; do
    if [ "$g" = CPU ]; then acc=; else acc="--gpu $g"; fi
    if $C new -s $S $acc 2>&1 | tail -1 | grep -q READY; then ok=1; GPU=$g; break; fi
    log "new $g failed: $(usage_line | grep -i assign | tr -d '\n')"
    timeout 300 $P $HOME/bin/colab_cleanup.py --orphans >/dev/null
  done
  [ $ok = 1 ] || [ $SECONDS -ge $A_END ] && break
  sleep 1200
done
[ $ok = 1 ] || { log "ABORT: no session on $GPU"; exit 1; }
T0=$SECONDS
( sleep $((DEADLINE_S + 900)); kill -TERM $$ ) & WD=$!
log "session READY on $GPU; rate: $(rate)/h"
timeout 120 $SSH true || { log "ABORT: ssh to the session failed"; exit 1; }
t=$SECONDS; timeout 1800 $SCP $W colab-$S:/content/ag_ckpt.tar; rc=$?
log "weights scp rc=$rc s $((SECONDS - t)) bytes $(stat -c %s $W)"
[ $rc = 0 ] || { log "ABORT: weights upload failed"; exit 1; }
timeout 120 $SCP scripts/phase1/s2ag/ag_colab_job.py colab-$S:/content/job.py || { log "ABORT: job upload failed"; exit 1; }

fetch() {  # tag: download <tag>.done, then every file it lists, checking sha256
  local tag=$1 f want got t
  [ -e $D/raw/$tag.fetched ] && return 0
  timeout 120 $SCP colab-$S:/content/out/$tag.done $D/raw/$tag.done.json >/dev/null 2>&1 || return 1
  for f in $(python3 -c "import json,sys; print(' '.join(json.load(open(sys.argv[1]))['files']))" $D/raw/$tag.done.json); do
    want=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['files'][sys.argv[2]])" $D/raw/$tag.done.json $f)
    t=$SECONDS
    timeout 1800 $SCP colab-$S:/content/out/$f $D/raw/$f >/dev/null 2>&1
    got=$(sha256sum $D/raw/$f 2>/dev/null | cut -d' ' -f1)
    log "fetch $f bytes $(stat -c %s $D/raw/$f 2>/dev/null) s $((SECONDS - t)) sha_ok $([ "$got" = "$want" ] && echo 1 || echo 0)"
    [ "$got" = "$want" ] || return 1
  done
  touch $D/raw/$tag.fetched
}

POLL='tail -n 3 /content/job.log 2>/dev/null || echo nolog; pgrep -f "python3 [j]ob.py" >/dev/null && echo JOB_ALIVE; echo DONE_TAGS $(cd /content/out 2>/dev/null && ls *.done 2>/dev/null | sed "s/\.done$//")'
run_job() {  # mode chrom...: launch, poll every 60 s, fetch finished parts; returns 0 only on the job's own end line
  local fails=0 OUT B rc=1
  job_args "$@" > $TMP/job_args.json
  timeout 120 $SCP $TMP/job_args.json colab-$S:/content/job_args.json || { log "ABORT: args upload failed"; return 1; }
  timeout 120 $SSH 'cd /content; for f in job.log job.out; do [ -f $f ] && mv $f ${f%.*}.1.${f##*.}; done; nohup python3 job.py > job.out 2>&1 < /dev/null & echo launched' \
    || { log "ABORT: launch failed"; return 1; }
  log "launched $*"
  while [ $((SECONDS - T0)) -lt "$DEADLINE_S" ]; do
    sleep 60
    OUT=$(timeout 90 $SSH "$POLL" 2>&1)
    if ! echo "$OUT" | grep -q DONE_TAGS; then
      fails=$((fails + 1)); log "poll failed ($fails): $(echo "$OUT" | tail -1 | cut -c1-200)"
      [ $fails -ge 3 ] && { log "ABORT: session unreachable"; return 1; }
      continue
    fi
    fails=0
    B=$(bal)
    log "poll $((SECONDS - T0)) s spent_cu $(spent "$B") | $(echo "$OUT" | grep -v 'DONE_TAGS\|JOB_ALIVE' | tail -1 | cut -c1-300)"
    for tag in $(echo "$OUT" | sed -n 's/^DONE_TAGS //p'); do fetch $tag || log "fetch $tag failed (retry next poll)"; done
    if echo "$OUT" | grep -q '"total_s"'; then
      log "JOB_END $(echo "$OUT" | grep total_s | tail -1 | cut -c1-300)"
      echo "$OUT" | grep total_s | tail -1 | grep -q FAILED || rc=0
      break
    fi
    if ! echo "$OUT" | grep -q JOB_ALIVE; then log "ABORT: job process gone without an end line"; break; fi
    if python3 -c "import sys; sys.exit(0 if $B0 - ${B:-$B0} >= $CAP else 1)"; then log "ABORT: CU cap $CAP reached"; return 1; fi
  done
  [ $((SECONDS - T0)) -ge "$DEADLINE_S" ] && log "DEADLINE hit"
  OUT=$(timeout 90 $SSH "$POLL" 2>&1)
  for tag in $(echo "$OUT" | sed -n 's/^DONE_TAGS //p'); do fetch $tag || log "fetch $tag failed (final)"; done
  return $rc
}

if [ "$MODE" != auto ]; then
  run_job "$MODE" $CHROMS
  exit
fi
run_job smoke ${CHROMS%% *} || exit 1
[ -e $D/raw/smoke.fetched ] || { log "STOP: smoke outputs not fetched"; exit 1; }
# go / stop for the full run, written before any smoke output was seen
R=$(rate)
python3 - $D/raw/smoke.done.json "${R:-0}" "$(spent "$(bal)")" "$CAP" $CHROMS > $TMP/check.txt <<'EOF'
import json, math, sys
s, rate, spent, cap = json.load(open(sys.argv[1])), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
L = {"chr20": 64444167, "chr18": 80373285, "chr9": 138394717, "chr7": 159345973}
d = s["fast_vs_public"]
same = "public_error" not in d and d["splice_max_abs"] <= 1e-3 and d["rna_max_abs"] <= 1e-3 * max(1.0, d["rna_max_value"])
per = sum(s["window_s"]) / len(s["window_s"])
pick, hours = [], 0.0
for c in sys.argv[5:]:
    h = sum(math.ceil(L[x] / 524288) for x in pick + [c]) * per * 1.25 / 3600 + 0.25   # +25%, 15 min for I/O
    if spent + rate * h > cap:
        break
    pick, hours = pick + [c], h
go = same and rate > 0 and bool(pick)
print("GO" if go else "STOP", json.dumps({"fast_equals_public": same, "rate_per_h": rate, "window_s": round(per, 2),
                                          "chroms": pick, "projected_h": round(hours, 2),
                                          "projected_cu": round(rate * hours, 2), "spent_cu": spent, "cap": cap}))
print("CHROMS", " ".join(pick))
EOF
log "CHECK $(head -1 $TMP/check.txt)"
grep -q '^GO' $TMP/check.txt || exit 1
run_job full $(sed -n 's/^CHROMS //p' $TMP/check.txt)
