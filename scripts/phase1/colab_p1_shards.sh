#!/bin/bash
# Colab route for P1.1 features (decision colab-topup, 2026-10-22): extract one shard in a Colab VM, download it
# to ont, release the VM.  Runs ON ont (Colab CLI + token live there).  No secrets enter the VM.
#   colab_p1_shards.sh train A:B [GPU] [DEADLINE_S]          e.g. colab_p1_shards.sh train 0:20 T4
#   colab_p1_shards.sh chrom K0:K1 [GPU] [DEADLINE_S] [NAME]  e.g. colab_p1_shards.sh chrom 0:10 T4 1800 chr18
# Output: data/p1_shards/colab/<stem>.{npy,json}; then on ont:
#   python3 scripts/phase1/assemble_p1_shards.py --shards data/p1_shards/colab [--ref-dir data/p1_shards/ont]
set -u
KIND=$1; RANGE=$2; GPU=${3:-T4}; DEADLINE_S=${4:-1800}; NAME=${5:-chr18}
ROOT=$HOME/ncrna-annot; S=p1shard
C="$HOME/.local/bin/colab --auth oauth2"
P=$HOME/.local/share/uv/tools/google-colab-cli/bin/python
OUT=$ROOT/data/p1_shards/colab; mkdir -p "$OUT"
STG=$(mktemp -d /tmp/p1shard.XXXX); T0=$SECONDS
cleanup() {
  $C stop -s $S 2>&1 | tail -1
  $P $HOME/bin/colab_cleanup.py --orphans
  $C usage 2>&1 | tail -3
  rm -rf "$STG"; echo "wall_s $((SECONDS - T0))"
}
trap cleanup EXIT
cd $ROOT
tar czf $STG/src.tgz --exclude=__pycache__ src
# stage only the codes the shard needs (train: the A:B windows; chrom: the bp span of tiles K0:K1)
ARGS=$(python3 - "$KIND" "$RANGE" "$NAME" "$STG" <<'PY'
import sys, numpy as np
sys.path.insert(0, "src"); from ncrna_annot.features import generanno as G
kind, rng, name, stg = sys.argv[1:5]; a, b = map(int, rng.split(":"))
if kind == "train":
    c = np.load("data/p1/train_codes.npy", mmap_mode="r"); np.save(f"{stg}/codes.npy", np.asarray(c[a:b]))
    print(f"--kind train --range 0:{b-a} --index-offset {a}")
else:
    c = np.load(f"data/p1/{name}_codes.npy", mmap_mode="r"); st = G.chromosome_tiles(len(c))
    lo, hi = st[a], st[b-1] + G.WIN; np.save(f"{stg}/codes.npy", np.asarray(c[lo:hi]))
    print(f"--kind chrom --name {name} --tiles {a}:{b} --chrom-len {len(c)} --bp-offset {lo}")
PY
)
echo "args: $ARGS  ($(du -h $STG/codes.npy | cut -f1) codes)"
$C usage 2>&1 | head -1
ok=0
for try in 1 2; do
  if $C new -s $S --gpu "$GPU" 2>&1 | tail -1 | grep -q READY; then ok=1; break; fi
  echo "new failed (try $try); releasing any half-assigned VM"; $P $HOME/bin/colab_cleanup.py --all >/dev/null
done
[ $ok = 1 ] || { echo "ABORT: could not get a $GPU session"; exit 1; }
for f in $STG/src.tgz $STG/codes.npy scripts/phase1/extract_p1_shard.py; do $C upload -s $S "$f" /content/$(basename $f) 2>&1 | tail -1; done
cat > $STG/launch.py <<PY
import subprocess, tarfile
tarfile.open("/content/src.tgz").extractall("/content/code")
subprocess.Popen("nohup bash -c 'pip install -q transformers==4.57.6 && python3 /content/extract_p1_shard.py --src /content/code/src --codes /content/codes.npy --out /content/shards $ARGS' > /content/job.log 2>&1 &", shell=True)
print("launched")
PY
$C exec -s $S --timeout 60 -f $STG/launch.py 2>&1 | tail -1
echo 'import os; p="/content/job.log"; print(open(p).read()[-2500:] if os.path.exists(p) else "nolog")' > $STG/poll.py
while [ $((SECONDS - T0)) -lt "$DEADLINE_S" ]; do
  sleep 60
  LOG=$($C exec -s $S --timeout 60 -f $STG/poll.py 2>&1)
  echo "--- poll at $((SECONDS - T0)) s: $(echo "$LOG" | tail -1 | cut -c1-120)"
  if echo "$LOG" | grep -q "total_s"; then echo "$LOG" | tail -3; break; fi
  if echo "$LOG" | grep -q -E "Traceback|Error"; then echo "$LOG" | tail -20; echo "JOB FAILED"; exit 2; fi
done
[ $((SECONDS - T0)) -ge "$DEADLINE_S" ] && { echo "DEADLINE hit; last log:"; echo "$LOG" | tail -20; exit 3; }
STEM=$(echo "$LOG" | grep -o '"stem": "[^"]*"' | tail -1 | cut -d'"' -f4)
echo "stem $STEM"
TD=$SECONDS
for ext in json npy; do $C download -s $S /content/shards/$STEM.$ext $OUT/$STEM.$ext 2>&1 | tail -1; done
DT=$((SECONDS - TD)); MB=$(( $(stat -c %s $OUT/$STEM.npy) / 1048576 ))
echo "download ${MB} MB in ${DT} s = $(( MB / (DT>0?DT:1) )) MB/s"
python3 - $OUT/$STEM.npy $OUT/$STEM.json <<'PY'
import sys, json, hashlib
h = hashlib.sha256(open(sys.argv[1], "rb").read()).hexdigest(); m = json.load(open(sys.argv[2]))
print("sha256", "OK" if h == m["sha256"] else "MISMATCH", "| compute", m["compute"], "| gpu", m["gpu"], "| per-unit s", m.get("window_s", m.get("tile_s")))
PY
