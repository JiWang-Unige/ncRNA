#!/bin/bash
# S2-AG (PLAN_S2_AG.md): run ag_colab_job.py in a Colab GPU VM and fetch its outputs to ont as each part completes.
#   ag_colab_run.sh smoke|full <deadline_s> <gpu> <cu_cap> <chrom>...
# Caps: wall clock (deadline_s) and compute units spent by this run (cu_cap, measured from the balance at start).
# On exit, whatever happened: stop the session, unassign orphans (colab_cleanup.py --orphans), log the usage.
set -u
MODE=$1; DEADLINE_S=$2; GPU=$3; CAP=$4; shift 4; CHROMS="$*"
cd ~/ncrna-annot || exit 1
S=agrun
C="$HOME/.local/bin/colab --auth oauth2"
P=$HOME/.local/share/uv/tools/google-colab-cli/bin/python
D=outputs/PHASE1/S2AG/$MODE
mkdir -p $D/raw
TMP=$(mktemp -d)
T0=$SECONDS
log() { echo "$(date -u +%FT%TZ) $*" | tee -a $D/run.log; }
bal() { $C usage 2>&1 | sed -n 's/.*Current balance: \([0-9.]*\).*/\1/p' | head -1; }
spent() { python3 -c "print(round($B0 - ${1:-$B0}, 2))"; }
B0=$(bal)
[ -n "$B0" ] || { echo "ABORT: no balance"; exit 1; }
cleanup() {
  $C download -s $S /content/job.log $D/job.log >/dev/null 2>&1
  $C download -s $S /content/job.out $D/job.out >/dev/null 2>&1
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
python3 - "$MODE" $CHROMS > $TMP/job_args.json <<'EOF'
import json, sys
c = json.load(open("configs/s2ag.json"))
print(json.dumps({"mode": sys.argv[1], "chroms": sys.argv[2:], "ag_commit": c["ag_commit"],
                  "codes_sha256": {k: c["codes_sha256"][k] for k in sys.argv[2:]}}))
EOF
cat > $TMP/launch.py <<'EOF'
import subprocess; subprocess.Popen("nohup python3 /content/job.py > /content/job.out 2>&1 &", shell=True); print("launched")
EOF
cat > $TMP/poll.py <<'EOF'
import glob, os
p = "/content/job.log"
print("\n".join(open(p).read().splitlines()[-3:]) if os.path.exists(p) else "nolog")
print("DONE_TAGS", " ".join(sorted(os.path.basename(f)[:-5] for f in glob.glob("/content/out/*.done"))))
EOF

ok=0
for try in 1 2; do
  if $C new -s $S --gpu "$GPU" 2>&1 | tail -1 | grep -q READY; then ok=1; break; fi
  log "new failed (try $try)"
  $P $HOME/bin/colab_cleanup.py --all >/dev/null
done
[ $ok = 1 ] || { log "ABORT: no $GPU session"; exit 1; }
log "session READY; rate: $($C usage 2>&1 | grep -i 'rate')"
t=$SECONDS; $C upload -s $S $W /content/ag_ckpt.tar 2>&1 | tail -1; log "weights upload_s $((SECONDS - t)) bytes $(stat -c %s $W)"
$C upload -s $S scripts/phase1/s2ag/ag_colab_job.py /content/job.py 2>&1 | tail -1
$C upload -s $S $TMP/job_args.json /content/job_args.json 2>&1 | tail -1
$C exec -s $S --timeout 60 -f $TMP/launch.py 2>&1 | tail -1

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

fails=0
while [ $((SECONDS - T0)) -lt "$DEADLINE_S" ]; do
  sleep 60
  OUT=$($C exec -s $S --timeout 60 -f $TMP/poll.py 2>&1)
  if ! echo "$OUT" | grep -q DONE_TAGS; then
    fails=$((fails + 1)); log "poll failed ($fails): $(echo "$OUT" | tail -1 | cut -c1-200)"
    [ $fails -ge 5 ] && { log "ABORT: session unreachable"; break; }
    continue
  fi
  fails=0
  B=$(bal)
  log "poll $((SECONDS - T0)) s spent_cu $(spent "$B") | $(echo "$OUT" | grep -v DONE_TAGS | tail -1 | cut -c1-300)"
  for tag in $(echo "$OUT" | sed -n 's/^DONE_TAGS //p'); do fetch $tag || log "fetch $tag failed (retry next poll)"; done
  if echo "$OUT" | grep -q '"total_s"'; then log "JOB_END $(echo "$OUT" | grep total_s | tail -1 | cut -c1-300)"; break; fi
  if python3 -c "import sys; sys.exit(0 if $B0 - ${B:-$B0} >= $CAP else 1)"; then log "ABORT: CU cap $CAP reached"; break; fi
done
[ $((SECONDS - T0)) -ge "$DEADLINE_S" ] && log "DEADLINE hit"
OUT=$($C exec -s $S --timeout 60 -f $TMP/poll.py 2>&1)
for tag in $(echo "$OUT" | sed -n 's/^DONE_TAGS //p'); do fetch $tag || log "fetch $tag failed (final)"; done
