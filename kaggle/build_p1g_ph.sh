#!/bin/bash
# v1.7 (protocol section 14) G2-PH arms on Kaggle castwj.  Runs ON ont.  usage: bash kaggle/build_p1g_ph.sh [G2PH|G2PH1]
#   ncrna-p1g-<arm>-<tag>       train (p1s-data + p1g-data; code inline)          -> best.pt
#   ncrna-p1g-<arm>-<tag>-eval  chr18 evaluation (p1s2-data + chr18 features + train output; code inline)
# tags: smoke-run (64 windows, 1 epoch; eval on the first 1.2 Mb), 42, 123, 456.
set -euo pipefail
ACC=castwj
ARM=${1:-G2PH}
SLUG=$(echo $ARM | tr 'A-Z' 'a-z')
cd ~/ncrna-annot
tar czf /tmp/p1g_src_ph.tgz --exclude=__pycache__ src
SRC64=$(base64 -w0 /tmp/p1g_src_ph.tgz)
HDR="import subprocess, sys
subprocess.check_call([sys.executable, '-m', 'pip', 'install', '-q', 'numba'])
import base64, glob, io, json, os, tarfile
def find(name, hint=''):
    c = [x for x in glob.glob(f'/kaggle/input/**/{name}', recursive=True) if hint in x]
    return sorted(c)[0]
tarfile.open(fileobj=io.BytesIO(base64.b64decode('$SRC64'))).extractall('/tmp/code'); sys.path.insert(0, '/tmp/code/src')
import numpy as np, torch
from ncrna_annot.train import screen_g2ph"
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
EVDS="[\"$ACC/ncrna-p1s2-data\"]"
for tag in smoke-run 42 123 456; do
  K=kaggle/ncrna-p1g-$SLUG-$tag; mkdir -p $K
  if [ $tag = smoke-run ]; then
    SEED=0; SUB='tr = {k: np.asarray(v[:64]) for k, v in tr.items()}
va = {k: np.asarray(v[:8]) for k, v in va.items()}'; EP=1
  else
    SEED=$tag; SUB=''; EP=12
  fi
  cat > $K/kernel.py <<PY
$HDR
$TRAIN
$SUB
screen_g2ph.run('$ARM', $SEED, tr, va, '/kaggle/working/p1g_${SLUG}_$tag', epochs=$EP)
print('TRAIN_DONE', flush=True)
PY
  meta ncrna-p1g-$SLUG-$tag "$TRDS" "[]" > $K/kernel-metadata.json
  E=kaggle/ncrna-p1g-$SLUG-$tag-eval; mkdir -p $E
  if [ $tag = smoke-run ]; then
    CUT='ch.update(feats=ch["feats"][:200000], codes=ch["codes"][:1200000])'
  else
    CUT=''
  fi
  cat > $E/kernel.py <<PY
$HDR
$CH
$CUT
ch["gff_out"] = "/kaggle/working/${SLUG}_${tag}_chr18.gff3"
screen_g2ph.evaluate_checkpoint(find("best.pt", "ncrna-p1g-$SLUG-$tag/"), "$ARM", 512, ch, "/kaggle/working/p1g_${SLUG}_$tag")
print("EVAL_DONE", flush=True)
PY
  meta ncrna-p1g-$SLUG-$tag-eval "$EVDS" "[\"$ACC/ncrna-p1s2-chr18feat\", \"$ACC/ncrna-p1g-$SLUG-$tag\"]" > $E/kernel-metadata.json
done
for k in kaggle/ncrna-p1g-$SLUG-*; do python3 -m py_compile $k/kernel.py; python3 -c "import json,sys; json.load(open(sys.argv[1]))" $k/kernel-metadata.json; done
ls -d kaggle/ncrna-p1g-$SLUG-*
