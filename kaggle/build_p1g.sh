#!/bin/bash
# v1.6 (protocol section 13) G1B arm on Kaggle castwj.  Runs ON ont.  usage: bash kaggle/build_p1g.sh [new|version]
# dataset ncrna-p1g-data (src + P1-G0 event targets).  Training and chr18 evaluation are separate kernels because
# p1s-data + p1s2-data + the 13.7 GB chr18 feature output no longer start together (empty-log ERROR, 2026-10-22).
#   ncrna-p1g-g1b-<tag>       train  (p1s-data + p1g-data)                         -> best.pt
#   ncrna-p1g-g1b-<tag>-eval  eval   (p1s2-data + p1g-data + chr18feat + train kernel output) -> result.json, GFF3
set -euo pipefail
ACC=castwj
cd ~/ncrna-annot
DS=kaggle/ds_p1g_$ACC; rm -rf $DS && mkdir -p $DS
cp data/p1g/train_events.npy data/p1g/val_events.npy data/p1g/events_stats.json $DS/
tar czf $DS/src.tgz src
printf '{"title": "ncrna-p1g-data", "id": "%s/ncrna-p1g-data", "licenses": [{"name": "CC0-1.0"}]}' $ACC > $DS/dataset-metadata.json
HDR='import subprocess, sys
subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "numba"])
import glob, json, os, tarfile
def find(name, hint=""):
    c = [x for x in glob.glob(f"/kaggle/input/**/{name}", recursive=True) if hint in x]
    return sorted(c)[0]
tarfile.open(find("src.tgz", "p1g-data")).extractall("/tmp/code"); sys.path.insert(0, "/tmp/code/src")
import numpy as np, torch
from ncrna_annot.train import screen_g1b'
TRAIN='tr = dict(feats=np.load(find("feat_train.npy"), mmap_mode="r"), codes=np.load(find("train_codes.npy")),
          bits=np.load(find("train_bits_v50uev.npy")), events=np.load(find("train_events.npy"), mmap_mode="r"))
va = dict(feats=np.load(find("feat_val.npy"), mmap_mode="r"), codes=np.load(find("val_codes.npy")),
          bits=np.load(find("val_bits_v50uev.npy")))'
CH='ch = dict(feats=np.load(find("chr18_feat.npy"), mmap_mode="r"), codes=np.load(find("chr18_codes.npy")),
          truth_gtf=find("chr18_v50_full.gtf"), chrom="chr18")'
meta() {  # slug datasets kernel_sources
  printf '{"id": "%s/%s", "title": "%s", "code_file": "kernel.py", "language": "python", "kernel_type": "script", "is_private": "true", "enable_gpu": "true", "enable_tpu": "false", "enable_internet": "true", "machine_shape": "", "dataset_sources": %s, "competition_sources": [], "kernel_sources": %s, "model_sources": []}' $ACC $1 $1 "$2" "$3"
}
TRDS="[\"$ACC/ncrna-p1s-data\", \"$ACC/ncrna-p1g-data\"]"
EVDS="[\"$ACC/ncrna-p1s2-data\", \"$ACC/ncrna-p1g-data\"]"
for tag in smoke-run 42 123 456; do
  K=kaggle/ncrna-p1g-g1b-$tag; mkdir -p $K
  if [ $tag = smoke-run ]; then
    cat > $K/kernel.py <<PY
$HDR
$TRAIN
tr = {k: np.asarray(v[:64]) for k, v in tr.items()}
va = {k: np.asarray(v[:8]) for k, v in va.items()}
screen_g1b.run(0, tr, va, None, "/kaggle/working/p1g_g1b_smoke", epochs=1)
print("TRAIN_DONE", flush=True)
PY
  else
    cat > $K/kernel.py <<PY
$HDR
$TRAIN
screen_g1b.run($tag, tr, va, None, "/kaggle/working/p1g_g1b_$tag")
print("TRAIN_DONE", flush=True)
PY
  fi
  meta ncrna-p1g-g1b-$tag "$TRDS" "[]" > $K/kernel-metadata.json
  E=kaggle/ncrna-p1g-g1b-$tag-eval; mkdir -p $E
  if [ $tag = smoke-run ]; then
    cat > $E/kernel.py <<PY
$HDR
$CH
ch.update(feats=ch["feats"][:500000], codes=ch["codes"][:3000000], gff_out="/kaggle/working/smoke_chr18.gff3")
screen_g1b.evaluate_checkpoint(find("best.pt", "ncrna-p1g-g1b-smoke-run/"), 512, ch, "/kaggle/working/p1g_g1b_smoke")
print("EVAL_DONE", flush=True)
PY
  else
    cat > $E/kernel.py <<PY
$HDR
$CH
ch["gff_out"] = "/kaggle/working/g1b_${tag}_chr18.gff3"
screen_g1b.evaluate_checkpoint(find("best.pt", "ncrna-p1g-g1b-$tag/"), 512, ch, "/kaggle/working/p1g_g1b_$tag")
print("EVAL_DONE", flush=True)
PY
  fi
  meta ncrna-p1g-g1b-$tag-eval "$EVDS" "[\"$ACC/ncrna-p1s2-chr18feat\", \"$ACC/ncrna-p1g-g1b-$tag\"]" > $E/kernel-metadata.json
done
for k in kaggle/ncrna-p1g-g1b-*; do python3 -m py_compile $k/kernel.py; python3 -c "import json,sys; json.load(open(sys.argv[1]))" $k/kernel-metadata.json; done
if [ "${1:-}" = version ]; then ~/bin/kaggle-as $ACC datasets version -p $DS -m "src update" -q; fi
ls -d kaggle/ncrna-p1g-g1b-*; du -sh $DS
