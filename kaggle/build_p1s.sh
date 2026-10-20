#!/bin/bash
# P1-S proxy screening on Kaggle (protocol v1.1 section 8). Runs ON ont. usage: bash kaggle/build_p1s.sh <account>
set -euo pipefail
ACC=${1:-castwj}
cd ~/ncrna-annot
DS=kaggle/ds_p1s_$ACC
rm -rf $DS && mkdir -p $DS
ln data/p05_feat/feat_train.npy $DS/feat_train.npy; ln data/p05_feat/feat_val.npy $DS/feat_val.npy
for f in data/p1s_v12/*.npy data/p1s_v12/val_windows.tsv data/p1s_v12/val_truth_v50.gtf data/p1s_v12/repair_stats.json; do cp $f $DS/; done
tar czf $DS/src.tgz src
printf '{"title": "ncrna-p1s-data", "id": "%s/ncrna-p1s-data", "licenses": [{"name": "CC0-1.0"}]}' $ACC > $DS/dataset-metadata.json
for spec in r:v50uev:0:argmax s2:v50uevlnc:1: s4:v50f:0: s5:v41uev:0:; do
  IFS=: read arm bits lnc modes <<< "$spec"
  K=kaggle/ncrna-p1s-$arm-$ACC; mkdir -p $K
  cat > $K/kernel.py <<PY
import glob, json, os, subprocess, sys, tarfile
subprocess.check_call([sys.executable, '-m', 'pip', 'install', '-q', 'numba'])
IN = os.path.dirname(glob.glob('/kaggle/input/**/feat_train.npy', recursive=True)[0])
tarfile.open(f'{IN}/src.tgz').extractall('/tmp/code'); sys.path.insert(0, '/tmp/code/src')
import numpy as np
from ncrna_annot.train import screen
tr = dict(feats=np.load(f'{IN}/feat_train.npy', mmap_mode='r'), codes=np.load(f'{IN}/train_codes.npy'),
          bits=np.load(f'{IN}/train_bits_$bits.npy'))
va = dict(feats=np.load(f'{IN}/feat_val.npy', mmap_mode='r'), codes=np.load(f'{IN}/val_codes.npy'),
          starts=[int(l.split()[1]) for l in open(f'{IN}/val_windows.tsv')], truth_gtf=f'{IN}/val_truth_v50.gtf', chrom='chr18')
modes = tuple(m for m in '$modes'.split(',') if m)
res = screen.run('${arm^^}', tr, va, '/kaggle/working/p1s_$arm', lnc_only=bool($lnc), extra_modes=modes)
print('RESULT', json.dumps(res), flush=True)
PY
  printf '{"id": "%s/ncrna-p1s-%s", "title": "ncrna-p1s-%s", "code_file": "kernel.py", "language": "python", "kernel_type": "script", "is_private": "true", "enable_gpu": "true", "enable_tpu": "false", "enable_internet": "true", "machine_shape": "", "dataset_sources": ["%s/ncrna-p1s-data"], "competition_sources": [], "kernel_sources": [], "model_sources": []}' $ACC $arm $arm $ACC > $K/kernel-metadata.json
done
du -shL $DS; ls kaggle | grep p1s
