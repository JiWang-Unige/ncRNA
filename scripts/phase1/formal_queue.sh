#!/bin/bash
# Protocol v1.8 (section 15) queue on ont.  usage: nohup bash scripts/phase1/formal_queue.sh queues/formal_v18.txt >> logs/formal_queue.log 2>&1 &
# Lines (in order):  <name> <train_formal_ont.py args...>   |   check-frag <name>   |   await <file>   |   # comment   |   END
# await <file>: block until that file (relative to the root) exists, e.g. a marker that the pilot check on ont passed;
# a file <file>.fail stops the queue.
# Runs one entry at a time through run_formal_ont.sh (GPU gate: waits for TE / other users, resumes after yielding).
# State: <queue>.state ("<name> done|failed").  A failed run or a triggered fragmentation check stops the queue.
# The queue file is re-read after every entry, so entries can be appended while it runs; without END it waits for more.
set -u
cd ${NCRNA_ROOT:-~/ncrna-annot}
Q=$1; S=$Q.state; touch $S
log() { echo "$(date -u +%Y-%m-%dT%H:%MZ) $*"; }
state() { awk -v n="$1" '$1==n{print $2}' $S | tail -1; }
while true; do
  next=""
  while read -r name rest; do
    [ -z "$name" ] && continue
    case $name in \#*) continue;; END) next=END; break;; esac
    key=$name; [ "$name" = check-frag ] && key="check-frag:$rest"; [ "$name" = await ] && key="await:$rest"
    [ -n "$(state $key)" ] && continue
    next="$name $rest"; break
  done < $Q
  if [ "$next" = END ]; then log "QUEUE END"; exit 0; fi
  if [ -z "$next" ]; then sleep 600; continue; fi
  set -- $next
  if [ "$1" = await ]; then
    if [ -e "$2.fail" ]; then echo "await:$2 failed" >> $S; log "await $2: .fail marker present - queue stopped"; exit 2; fi
    if [ -e "$2" ]; then echo "await:$2 done" >> $S; log "await $2 satisfied"; continue; fi
    sleep 600; continue
  fi
  if [ "$1" = check-frag ]; then
    r=outputs/PHASE1/FORMAL/$2/result.json
    n=$(python3 -c "import json;print(json.load(open('$r'))['S3']['n_pred_lnc'])")
    if [ "$n" -gt $((5 * 929)) ]; then
      echo "check-frag:$2 failed" >> $S; log "FRAGMENTATION_RULE fired: $2 n_pred_lnc=$n > 5x929 - queue stopped for discussion (protocol section 15 P1.3)"; exit 2
    fi
    echo "check-frag:$2 done" >> $S; log "check-frag $2 ok (n_pred_lnc=$n <= 4645)"; continue
  fi
  name=$1; shift
  log "run $name $*"
  if bash scripts/phase1/run_formal_ont.sh "$name" "$@" >> logs/formal_$name.log 2>&1; then
    echo "$name done" >> $S; log "done $name"
  else
    echo "$name failed" >> $S; log "FAILED $name (see logs/formal_$name.log) - queue stopped"; exit 1
  fi
done
