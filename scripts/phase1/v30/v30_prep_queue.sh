#!/bin/bash
# v3.0 data prep on ont (protocol section 18, execution details 3-6): waits for the stage-1 downloads, then
# regression check of the extended label scripts on the v2.0 windows -> mm39 chr19 identity -> mouse U-ev (srav1m
# lifted to mm39) -> human / mouse windows -> v1.2 labels and events -> summary -> features of both window sets
# (gpu_gate.sh before every shard, GPU fraction 0.28) -> data/p1_feat_v30/{human_new,mouse}.npy.  CPU nice 19, cores 2-31.
cd ~/ncrna-annot || exit 1
. scripts/phase1/gpu_gate.sh
PY=~/envs/generanno/bin/python; export PYTHONPATH=src
O=outputs/PHASE1/V30; mkdir -p $O data/p1v30 data/p1_feat_v30
log() { echo "$(date -u +%FT%TZ) $*" >> $O/prep_queue.log; }
N="nice -n 19 taskset -c 2-31"
V=scripts/phase1/v30
VM39=data/batch2/mouse/gencode.vM39.annotation.gtf.gz
until grep -qE "STAGE1 done|FAIL" data/v30/fetch.log 2>/dev/null; do sleep 120; done
grep -q "STAGE1 done" data/v30/fetch.log || { log "downloads failed"; exit 1; }
log START
if [ ! -s data/p1v30/REGRESS.json ]; then
  $N $PY scripts/phase1/prep_screen.py --windows-train data/p1v20/train_windows.tsv --splits train --out data/p1v30/regress > $O/regress_labels.log 2>&1 &&
  $N $PY scripts/phase1/prep_events.py --windows-train data/p1v20/train_windows.tsv --splits train --out data/p1v30/regress > $O/regress_events.log 2>&1 &&
  $N $PY $V/check_v30_data.py regress > $O/regress_check.log 2>&1 || { log REGRESS_FAIL; exit 1; }
  rm -rf data/p1v30/regress
fi
log "regression ok (default options reproduce data/p1v20 bit for bit)"
$N $PY $V/check_v30_data.py mm39 > $O/mm39_check.log 2>&1 || { log MM39_FAIL; exit 1; }
log "mm39 chr19 = S5 chr19"
[ -s data/v30/junctions_mm39/LIFT_STATS.json ] || $N $PY $V/mouse_uev.py > $O/mouse_uev.log 2>&1 || { log UEV_FAIL; exit 1; }
log "mouse U-ev junctions lifted"
[ -s data/p1v30/human_new_windows.tsv ] || $N $PY $V/prep_v30_windows.py human > $O/windows_human.log 2>&1 || { log WINDOWS_FAIL human; exit 1; }
[ -s data/p1v30/mouse_windows.tsv ] || $N $PY $V/prep_v30_windows.py mouse > $O/windows_mouse.log 2>&1 || { log WINDOWS_FAIL mouse; exit 1; }
log "windows ok"
H="--windows-train data/p1v30/human_new_windows.tsv --splits train --out data/p1v30/human_new"
M="--windows-train data/p1v30/mouse_windows.tsv --splits train --out data/p1v30/mouse --genome data/v30/mouse/mm39.fa.gz --gtf $VM39 --junction-cache data/v30/junctions_mm39"
[ -s data/p1v30/human_new/train_bits_v50uev.npy ] || $N $PY scripts/phase1/prep_screen.py $H > $O/labels_human.log 2>&1 || { log LABELS_FAIL human; exit 1; }
[ -s data/p1v30/human_new/train_events.npy ] || $N $PY scripts/phase1/prep_events.py $H > $O/events_human.log 2>&1 || { log EVENTS_FAIL human; exit 1; }
[ -s data/p1v30/mouse/train_bits_v50uev.npy ] || $N $PY scripts/phase1/prep_screen.py $M > $O/labels_mouse.log 2>&1 || { log LABELS_FAIL mouse; exit 1; }
[ -s data/p1v30/mouse/train_events.npy ] || $N $PY scripts/phase1/prep_events.py $M > $O/events_mouse.log 2>&1 || { log EVENTS_FAIL mouse; exit 1; }
$N $PY $V/check_v30_data.py summary > $O/summary.log 2>&1 || { log SUMMARY_FAIL; exit 1; }
log "labels and events ok"
STG=data/p1_shards/v30stage; mkdir -p $STG
for part in human_new mouse; do
  for a in $(seq 0 250 2750); do
    b=$((a + 250))
    [ -e data/p1_feat_v30/done_${part}_$a.json ] && continue
    gpu_wait >> $O/prep_queue.log
    rm -rf $STG/*
    NCRNA_GPU_FRAC=0.28 $N $PY scripts/phase1/extract_p1_shard.py --out $STG --kind train --name $part --codes data/p1v30/$part/train_codes.npy --range $a:$b > $O/feat_${part}_$a.log 2>&1 || { log FEAT_FAIL $part $a; exit 1; }
    $N $PY scripts/phase1/assemble_p1_shards.py --shards $STG --feat-root data/p1_feat_v30 --codes-dir data/p1v30/$part >> $O/feat_${part}_$a.log 2>&1 || { log ASSEMBLE_FAIL $part $a; exit 1; }
    cp $STG/*.json data/p1_feat_v30/done_${part}_$a.json; rm -rf $STG/*
    log "features $part $a:$b ok"
  done
done
log PREP_END
