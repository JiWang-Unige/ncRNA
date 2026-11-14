#!/bin/bash
# P1.1 formal feature cache on ont (decision formal-features-ont, 2026-10-23): GENERanno L17+L23, fp16, bf16 compute.
# Shard by shard (train 100 windows, chromosomes 200 tiles): before EVERY shard pass the GPU gate (gpu_gate.sh:
# wait for TE processes and other users; the user's own coding jobs may share); each shard is < 1 min of GPU.
# Resumable: a shard whose .json sits in data/p1_shards/ont is skipped.  Shards are assembled into data/p1_feat
# (same layout as extract_p1.py) right away and their .npy deleted, so peak disk ~ cache + one shard.
#   nohup bash scripts/phase1/p1_feat_ont.sh > logs/p1_feat_ont.log 2>&1 &
set -u
cd ~/ncrna-annot
PY=~/envs/generanno/bin/python; export PYTHONPATH=src
DONE=${P1_DONE:-data/p1_shards/ont}; STG=data/p1_shards/stage; FEAT=data/p1_feat; mkdir -p $DONE $STG $FEAT logs
LOCK=logs/p1_feat_ont.lock
if [ -e $LOCK ] && kill -0 "$(cat $LOCK)" 2>/dev/null; then echo "another runner (pid $(cat $LOCK)) is active"; exit 1; fi
echo $$ > $LOCK; trap 'rm -f $LOCK' EXIT
TRAIN_STEP=${TRAIN_STEP:-100}; TILE_STEP=${TILE_STEP:-200}
log() { echo "$(date -u +%Y-%m-%dT%H:%MZ) $*"; }
. scripts/phase1/gpu_gate.sh   # TE / other users block; own coding jobs may share (rule 2026-10-26)
export NCRNA_GPU_FRAC=${NCRNA_GPU_FRAC:-0.4}   # cap this process at ~10 GB of the 24 GB card
wait_idle() { gpu_wait; }
run_shard() {  # stem, extract args...
  local stem=$1; shift
  [ -e $DONE/$stem.json ] && { log "skip $stem (done)"; return 0; }
  wait_idle
  rm -rf $STG/*
  log "extract $stem"
  $PY scripts/phase1/extract_p1_shard.py --out $STG "$@" > $STG/extract.log 2>&1 || { log "EXTRACT FAILED $stem"; tail -5 $STG/extract.log; return 1; }
  [ -e $STG/$stem.json ] || { log "stem mismatch: expected $stem, got $(ls $STG/*.json)"; return 1; }
  $PY scripts/phase1/assemble_p1_shards.py --shards $STG --feat-root $FEAT --codes-dir data/p1 > $STG/assemble.log 2>&1 || { log "ASSEMBLE FAILED $stem"; tail -5 $STG/assemble.log; return 1; }
  mv $STG/$stem.json $DONE/; rm -f $STG/$stem.npy
  log "done $stem ($(grep -o '"window_s": [0-9.]*\|"tile_s": [0-9.]*' $DONE/$stem.json))"
}
T0=$SECONDS; fails=0
CHROM_ONLY=${CHROM_ONLY:-0}   # 1: chromosomes only (P1.6 panel: NAMES="chicken_5 zebrafish_22 fly_2L worm_I")
NAMES=${NAMES:-"chr18 chr20 t2t_chr20"}
if [ "$CHROM_ONLY" != 1 ]; then
# 1) training windows
N=$($PY -c "import numpy as np;print(len(np.load('data/p1/train_codes.npy',mmap_mode='r')))")
for ((a=0; a<N; a+=TRAIN_STEP)); do b=$((a+TRAIN_STEP)); [ $b -gt $N ] && b=$N
  run_shard $(printf "train_%05d_%05d" $a $b) --kind train --codes data/p1/train_codes.npy --range $a:$b || fails=$((fails+1)); done
# 1b) checkpoint-selection windows (300 chr18 windows of P1-S; starts are not on the 6-mer grid -> own features)
V=data/p1v12/val_codes.npy
if [ -e $V ]; then NV=$($PY -c "import numpy as np;print(len(np.load('$V',mmap_mode='r')))")
  for ((a=0; a<NV; a+=TRAIN_STEP)); do b=$((a+TRAIN_STEP)); [ $b -gt $NV ] && b=$NV
    run_shard $(printf "val_%05d_%05d" $a $b) --kind train --name val --codes $V --range $a:$b || fails=$((fails+1)); done
else log "WARN no $V - validation windows skipped"; fi
fi
# 2) chromosomes (default: chr18 validation, chr20 test, T2T chr20 if its codes exist)
for name in $NAMES; do
  C=data/p1/${name}_codes.npy; [ -e $C ] || { log "WARN no codes for $name ($C) - skipped"; continue; }
  NT=$($PY -c "import numpy as np,sys;sys.path.insert(0,'src');from ncrna_annot.features import generanno as G;print(len(G.chromosome_tiles(len(np.load('$C',mmap_mode='r')))))")
  for ((k=0; k<NT; k+=TILE_STEP)); do k1=$((k+TILE_STEP)); [ $k1 -gt $NT ] && k1=$NT
    run_shard $(printf "%s_t%05d_%05d" $name $k $k1) --kind chrom --name $name --codes $C --tiles $k:$k1 || fails=$((fails+1)); done
done
[ "$CHROM_ONLY" = 1 ] && { du -sh $FEAT | log "cache size: $(cat)"; log "ALL DONE fails=$fails wall_s=$((SECONDS-T0))"; exit 0; }
$PY - <<PY
import json, glob, os, time
js = [json.load(open(f)) for f in glob.glob("$DONE/*.json")]
meta = dict(layers=[17, 23], dim=4096, dtype="float16", compute=sorted({j.get("compute") for j in js}),
            gpu=sorted({j.get("gpu") for j in js}), shards=len(js), train_windows=sum(j["shape"][0] for j in js if j["kind"] == "train"),
            chrom_tiles={n: sum(len(j["tiles"]) for j in js if j["kind"] == "chrom" and j["name"] == n) for n in ("chr18", "chr20", "t2t_chr20")},
            finished=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
json.dump(meta, open("$FEAT/meta.json", "w"), indent=1); print(json.dumps(meta))
PY
du -sh $FEAT | log "cache size: $(cat)"
log "ALL DONE fails=$fails wall_s=$((SECONDS-T0))"
