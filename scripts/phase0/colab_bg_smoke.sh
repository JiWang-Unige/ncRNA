#!/bin/bash
# Colab route test: start a script in the VM in the background, poll its log with short execs.
# Hard wall-clock cap (DEADLINE_S); always stops the session and unassigns orphans on exit.
#   colab_bg_smoke.sh <local_script.py> [deadline_s] [gpu]
set -u
SCRIPT=$1; DEADLINE_S=${2:-1080}; GPU=${3:-T4}; S=bgsmoke
C="$HOME/.local/bin/colab --auth oauth2"
P=$HOME/.local/share/uv/tools/google-colab-cli/bin/python
T0=$SECONDS
cleanup() {
  $C stop -s $S 2>&1 | tail -1
  $P $HOME/bin/colab_cleanup.py --orphans
  $C usage 2>&1 | tail -3
  echo "wall_s $((SECONDS - T0))"
}
trap cleanup EXIT
$C usage 2>&1 | head -1
ok=0
for try in 1 2; do
  if $C new -s $S --gpu "$GPU" 2>&1 | tail -1 | grep -q READY; then ok=1; break; fi
  echo "new failed (try $try); releasing any half-assigned VM"
  $P $HOME/bin/colab_cleanup.py --all >/dev/null
done
[ $ok = 1 ] || { echo "ABORT: could not get a $GPU session"; exit 1; }

$C upload -s $S "$SCRIPT" /content/job.py 2>&1 | tail -1
echo 'import subprocess; subprocess.Popen("nohup python3 /content/job.py > /content/job.log 2>&1 &", shell=True); print("launched")' \
  > /tmp/colab_launch_$$.py
$C exec -s $S --timeout 60 -f /tmp/colab_launch_$$.py 2>&1 | tail -1
echo 'import os; p="/content/job.log"; print(open(p).read()[-2500:] if os.path.exists(p) else "nolog")' \
  > /tmp/colab_poll_$$.py
while [ $((SECONDS - T0)) -lt "$DEADLINE_S" ]; do
  sleep 60
  LOG=$($C exec -s $S --timeout 60 -f /tmp/colab_poll_$$.py 2>&1)
  echo "--- poll at $((SECONDS - T0)) s: $(echo "$LOG" | tail -1 | cut -c1-120)"
  if echo "$LOG" | grep -q "total_s"; then echo "$LOG"; echo "DONE"; break; fi
done
[ $((SECONDS - T0)) -ge "$DEADLINE_S" ] && { echo "DEADLINE hit; last log:"; echo "$LOG" | tail -20; }
rm -f /tmp/colab_launch_$$.py /tmp/colab_poll_$$.py
