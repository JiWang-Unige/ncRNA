#!/bin/bash
# v2.0 data prep on ont (protocol section 17, details 1-3): replacement windows -> v1.2 labels and events for the
# new window list -> bit-exact check of the unchanged rows -> features of the 154 new windows (gpu_gate.sh before
# each shard, GPU fraction 0.28) -> data/p1_feat_v20/train_new.npy.  CPU nice 19 on cores 2-31.
cd ~/ncrna-annot || exit 1
. scripts/phase1/gpu_gate.sh
PY=~/envs/generanno/bin/python; export PYTHONPATH=src
O=outputs/PHASE1/V20; mkdir -p $O data/p1v20 data/p1_feat_v20
log() { echo "$(date -u +%FT%TZ) $*" >> $O/prep_queue.log; }
N="nice -n 19 taskset -c 2-31"
log START
[ -s data/p1v20/replaced.json ] || $N $PY scripts/phase1/v20/prep_v20_windows.py > $O/prep_windows.log 2>&1 || { log WINDOWS_FAIL; exit 1; }
log "windows ok"
[ -s data/p1v20/train_bits_v50uevlnc.npy ] || $N $PY scripts/phase1/prep_screen.py --windows-train data/p1v20/train_windows.tsv --splits train --out data/p1v20 > $O/prep_labels.log 2>&1 || { log LABELS_FAIL; exit 1; }
log "labels ok"
[ -s data/p1v20/train_events.npy ] || $N $PY scripts/phase1/prep_events.py --windows-train data/p1v20/train_windows.tsv --splits train --out data/p1v20 > $O/prep_events.log 2>&1 || { log EVENTS_FAIL; exit 1; }
log "events ok"
$N $PY scripts/phase1/v20/check_v20_data.py > $O/check.log 2>&1 || { log CHECK_FAIL; exit 1; }
log "check ok (unchanged rows bit-identical to data/p1v12)"
STG=data/p1_shards/v20stage; mkdir -p $STG
for r in 0:100 100:154; do
  a=${r%%:*}; b=${r##*:}
  [ -e data/p1_feat_v20/done_$a.json ] && continue
  gpu_wait >> $O/prep_queue.log
  rm -rf $STG/*
  NCRNA_GPU_FRAC=0.28 $N $PY scripts/phase1/extract_p1_shard.py --out $STG --kind train --name train_new --codes data/p1v20/train_new_codes.npy --range $a:$b > $O/feat_$a.log 2>&1 || { log FEAT_FAIL $a; exit 1; }
  $N $PY scripts/phase1/assemble_p1_shards.py --shards $STG --feat-root data/p1_feat_v20 --codes-dir data/p1v20 >> $O/feat_$a.log 2>&1 || { log ASSEMBLE_FAIL $a; exit 1; }
  cp $STG/*.json data/p1_feat_v20/done_$a.json; rm -rf $STG/*
  log "features $a:$b ok"
done
log PREP_END
