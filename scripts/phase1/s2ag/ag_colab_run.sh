#!/bin/bash
# S2-AG (PLAN_S2_AG.md): run ag_colab_job.py in a Colab GPU VM and fetch its outputs to ont as each part completes.
#   ag_colab_run.sh smoke|full|auto <deadline_s> <gpu[,gpu...]> <cu_cap> <chrom>...
# auto = smoke on the first chromosome, then, in the same session and only if the smoke check passes (fast path equal
# to predict_sequence; projected CU of the full run within what is left of cu_cap), the full run on all chromosomes.
# Assignment: tries each listed GPU in turn; with ASSIGN_HOURS=h (default 0: one round) repeats every 20 min for up to
# h hours (premium GPUs are often out of capacity; a failed attempt costs nothing).  Uses ~/bin/colab_long (the stock CLI
# with a 600 s HTTP read timeout instead of 120 s).
# Caps, counted from the moment a session is READY: wall clock (deadline_s) and compute units spent (cu_cap).
# On exit, whatever happened: stop the session, unassign orphans (colab_cleanup.py --orphans), log the usage.
set -u
MODE=$1; DEADLINE_S=$2; GPU=$3; CAP=$4; shift 4; CHROMS="$*"
cd ~/ncrna-annot || exit 1
S=agrun
C="$HOME/bin/colab_long --auth oauth2"
P=$HOME/.local/share/uv/tools/google-colab-cli/bin/python
D=outputs/PHASE1/S2AG/$MODE
mkdir -p $D/raw
TMP=$(mktemp -d)
T0=$SECONDS
log() { echo "$(date -u +%FT%TZ) $*" | tee -a $D/run.log; }
bal() { $C usage 2>&1 | sed -n 's/.*Current balance: \([0-9.]*\).*/\1/p' | head -1; }
rate() { $C usage 2>&1 | sed -n 's/.*Usage rate: \([0-9.]*\).*/\1/p' | head -1; }
spent() { python3 -c "print(round($B0 - ${1:-$B0}, 2))"; }
B0=$(bal)
[ -n "$B0" ] || { echo "ABORT: no balance"; exit 1; }
cleanup() {
  for f in job.log job.out job.1.log job.1.out; do $C download -s $S /content/$f $D/$f >/dev/null 2>&1; done
  $C stop -s $S 2>&1 | tail -1
  $P $HOME/bin/colab_cleanup.py --orphans
  log "after: $($C usage 2>&1 | grep -E 'balance|rate|assign' | tr '\n' ' ')"
  log "END wall_s $((SECONDS - T0)) spent_cu $(spent "$(bal)")"
  rm -rf $TMP
}
trap cleanup EXIT
log "START mode=$MODE gpu=$GPU cap=$CAP chroms=$CHROMS balance=$B0"

W=data/alphagenome/hf_all_folds.tar
[ -s $W ] || tar -cf $W -C data/alphagenome/hf_all_folds _CHECKPOINT_METADATA _METADATA manifest.ocdbt d ocdbt.process_0
job_args() {  # mode chrom...
  python3 - "$@" <<'EOF'
import json, sys
c = json.load(open("configs/s2ag.json"))
print(json.dumps({"mode": sys.argv[1], "chroms": sys.argv[2:], "ag_commit": c["ag_commit"],
                  "codes_sha256": {k: c["codes_sha256"][k] for k in sys.argv[2:]}}))
EOF
}
cat > $TMP/launch.py <<'EOF'
import os, subprocess
for f in ("job.log", "job.out"):  # a previous job's log moves aside, so its "total_s" line cannot end the next poll
    if os.path.exists(f"/content/{f}"):
        os.replace(f"/content/{f}", f"/content/{f[:-4]}.1.{f[-3:]}")
subprocess.Popen("nohup python3 /content/job.py > /content/job.out 2>&1 &", shell=True); print("launched")
EOF
cat > $TMP/poll.py <<'EOF'
import glob, os
p = "/content/job.log"
print("\n".join(open(p).read().splitlines()[-3:]) if os.path.exists(p) else "nolog")
print("DONE_TAGS", " ".join(sorted(os.path.basename(f)[:-5] for f in glob.glob("/content/out/*.done"))))
EOF

ok=0
A_END=$((SECONDS + ${ASSIGN_HOURS:-0} * 3600))
while [ $ok = 0 ]; do
  for g in ${GPU//,/ }; do
    if $C new -s $S --gpu "$g" 2>&1 | tail -1 | grep -q READY; then ok=1; GPU=$g; break; fi
    log "new $g failed: $($C usage 2>&1 | grep -i assign | tr -d '\n')"
    $P $HOME/bin/colab_cleanup.py --orphans >/dev/null
  done
  [ $ok = 1 ] || [ $SECONDS -ge $A_END ] && break
  sleep 1200
done
[ $ok = 1 ] || { log "ABORT: no session on $GPU"; exit 1; }
T0=$SECONDS
log "session READY on $GPU; rate: $(rate)/h"
t=$SECONDS; $C upload -s $S $W /content/ag_ckpt.tar 2>&1 | tail -1; log "weights upload_s $((SECONDS - t)) bytes $(stat -c %s $W)"
$C upload -s $S scripts/phase1/s2ag/ag_colab_job.py /content/job.py 2>&1 | tail -1

fetch() {  # tag: download <tag>.done, then every file it lists, checking sha256
  local tag=$1 f want got t
  [ -e $D/raw/$tag.fetched ] && return 0
  $C download -s $S /content/out/$tag.done $D/raw/$tag.done.json >/dev/null 2>&1 || return 1
  for f in $(python3 -c "import json,sys; print(' '.join(json.load(open(sys.argv[1]))['files']))" $D/raw/$tag.done.json); do
    want=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['files'][sys.argv[2]])" $D/raw/$tag.done.json $f)
    t=$SECONDS
    $C download -s $S /content/out/$f $D/raw/$f >/dev/null 2>&1
    got=$(sha256sum $D/raw/$f 2>/dev/null | cut -d' ' -f1)
    log "fetch $f bytes $(stat -c %s $D/raw/$f 2>/dev/null) s $((SECONDS - t)) sha_ok $([ "$got" = "$want" ] && echo 1 || echo 0)"
    [ "$got" = "$want" ] || return 1
  done
  touch $D/raw/$tag.fetched
}

run_job() {  # mode chrom...: launch, poll every 60 s, fetch finished parts; returns 0 only on the job's own end line
  local fails=0 OUT B rc=1
  job_args "$@" > $TMP/job_args.json
  $C upload -s $S $TMP/job_args.json /content/job_args.json 2>&1 | tail -1
  $C exec -s $S --timeout 60 -f $TMP/launch.py 2>&1 | tail -1
  log "launched $*"
  while [ $((SECONDS - T0)) -lt "$DEADLINE_S" ]; do
    sleep 60
    OUT=$($C exec -s $S --timeout 60 -f $TMP/poll.py 2>&1)
    if ! echo "$OUT" | grep -q DONE_TAGS; then
      fails=$((fails + 1)); log "poll failed ($fails): $(echo "$OUT" | tail -1 | cut -c1-200)"
      [ $fails -ge 5 ] && { log "ABORT: session unreachable"; return 1; }
      continue
    fi
    fails=0
    B=$(bal)
    log "poll $((SECONDS - T0)) s spent_cu $(spent "$B") | $(echo "$OUT" | grep -v DONE_TAGS | tail -1 | cut -c1-300)"
    for tag in $(echo "$OUT" | sed -n 's/^DONE_TAGS //p'); do fetch $tag || log "fetch $tag failed (retry next poll)"; done
    if echo "$OUT" | grep -q '"total_s"'; then
      log "JOB_END $(echo "$OUT" | grep total_s | tail -1 | cut -c1-300)"
      echo "$OUT" | grep total_s | tail -1 | grep -q FAILED || rc=0
      break
    fi
    if python3 -c "import sys; sys.exit(0 if $B0 - ${B:-$B0} >= $CAP else 1)"; then log "ABORT: CU cap $CAP reached"; return 1; fi
  done
  [ $((SECONDS - T0)) -ge "$DEADLINE_S" ] && log "DEADLINE hit"
  OUT=$($C exec -s $S --timeout 60 -f $TMP/poll.py 2>&1)
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
python3 - $D/raw/smoke.done.json "${R:-0}" "$(spent "$(bal)")" "$CAP" > $TMP/check.txt <<'EOF'
import json, sys
s, rate, spent, cap = json.load(open(sys.argv[1])), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
d = s["fast_vs_public"]
same = "public_error" not in d and d["splice_max_abs"] <= 1e-3 and d["rna_max_abs"] <= 1e-3 * max(1.0, d["rna_max_value"])
hours = s["est_predict_min"] * 1.25 / 60 + 0.25          # +25% and 15 min for downloads, saving and fetching
need = rate * hours
go = same and rate > 0 and spent + need <= cap
print("GO" if go else "STOP", json.dumps({"fast_equals_public": same, "rate_per_h": rate, "projected_h": round(hours, 2),
                                          "projected_cu": round(need, 2), "spent_cu": spent, "cap": cap}))
EOF
log "CHECK $(cat $TMP/check.txt)"
grep -q '^GO' $TMP/check.txt || exit 1
run_job full $CHROMS
