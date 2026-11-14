#!/bin/bash
# PLAN_UGFF: after Tiberius chr18 -> dev (RULES.json) -> build chr18, chr20; after Tiberius chr9 -> build chr9.  CPU only.
cd ~/ncrna-annot || exit 1
Q=outputs/PHASE1/UGFF/queue.log
log() { echo "$(date -u +%FT%TZ) $*" >> $Q; }
PY="nice -n 19 taskset -c 2-31 ~/envs/generanno/bin/python scripts/phase1/ugff/ugff.py"
until grep -q "tiberius chr18 rc=" $Q; do sleep 60; done
$PY dev > outputs/PHASE1/UGFF/dev.log 2>&1; log "dev rc=$?"
for c in chr18 chr20; do $PY build --chrom $c > outputs/PHASE1/UGFF/build_$c.log 2>&1; log "build $c rc=$?"; done
until grep -q "tiberius chr9 rc=" $Q; do sleep 60; done
$PY build --chrom chr9 > outputs/PHASE1/UGFF/build_chr9.log 2>&1; log "build chr9 rc=$?"
log UGFF_DONE
