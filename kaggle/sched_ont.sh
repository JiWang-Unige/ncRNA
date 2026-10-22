#!/bin/bash
# Autonomous castwj scheduler on ont (protocol v1.6/v1.7 runs).
# usage: nohup bash kaggle/sched_ont.sh kaggle/queue_<name>.txt > logs/kaggle_sched_<name>.log 2>&1 &
# queue lines: "<slug> [dep1,dep2,...]" in priority order.  Max 2 concurrent GPU sessions (T4).  On the weekly-quota
# error it sleeps 30 min and retries.  A kernel ending in ERROR/CANCEL is re-pushed once (infrastructure-error rule);
# a second failure marks it failed and its dependents are skipped.  Outputs of completed kernels are downloaded to
# outputs/PHASE1/P1G/<slug>.  State is kept in <queue>.state so the script can be restarted.
cd ~/ncrna-annot
Q=$1; S=$Q.state; touch $S
declare -A STATE
while read -r k v; do [ -n "$k" ] && STATE[$k]=$v; done < $S
save() { : > $S; for k in "${!STATE[@]}"; do echo "$k ${STATE[$k]}" >> $S; done; }
st() { ~/bin/kaggle-as castwj kernels status castwj/$1 2>&1 | tail -1 | grep -oE "RUNNING|QUEUED|COMPLETE|ERROR|CANCEL" | head -1; }
push() { ~/bin/kaggle-as castwj kernels push -p kaggle/$1 --accelerator NvidiaTeslaT4 2>&1 | tail -1; }
log() { echo "$(date -u +%FT%H:%MZ) $*"; }
while true; do
  active=0; alldone=1
  while read -r slug deps; do
    [ -z "$slug" ] && continue
    s=${STATE[$slug]:-pending}
    case $s in
      pushed|retried)
        r=$(st $slug)
        if [ "$r" = COMPLETE ]; then
          STATE[$slug]=done; save; log "done $slug"
          mkdir -p outputs/PHASE1/P1G/$slug && ~/bin/kaggle-as castwj kernels output castwj/$slug -p outputs/PHASE1/P1G/$slug > /dev/null 2>&1
        elif [ "$r" = ERROR ] || [ "$r" = CANCEL ]; then
          if [ $s = retried ]; then STATE[$slug]=failed; save; log "FAILED twice $slug"
          else out=$(push $slug); log "retry $slug: $out"; echo "$out" | grep -q successfully && { STATE[$slug]=retried; save; active=$((active+1)); }; fi
        else active=$((active+1)); fi ;;
    esac
    [[ "${STATE[$slug]:-pending}" =~ ^(done|failed)$ ]] || alldone=0
  done < $Q
  [ $alldone = 1 ] && { log "ALL_DONE"; exit 0; }
  quota=0
  while read -r slug deps; do
    [ -z "$slug" ] && continue
    [ "${STATE[$slug]:-pending}" = pending ] || continue
    [ $active -ge 2 ] && break
    ok=1; bad=0
    for d in ${deps//,/ }; do
      [ "${STATE[$d]:-pending}" = done ] || ok=0
      [ "${STATE[$d]:-pending}" = failed ] && bad=1
    done
    if [ $bad = 1 ]; then STATE[$slug]=failed; save; log "skip $slug (dependency failed)"; continue; fi
    [ $ok = 1 ] || continue
    out=$(push $slug); log "push $slug: $out"
    if echo "$out" | grep -q successfully; then STATE[$slug]=pushed; save; active=$((active+1)); sleep 30
    elif echo "$out" | grep -qi quota; then quota=1; break
    else break; fi
  done < $Q
  if [ $quota = 1 ]; then sleep 1800; else sleep 300; fi
done
