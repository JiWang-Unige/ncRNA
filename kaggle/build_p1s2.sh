#!/bin/bash
# P1-S2 (protocol v1.4 section 11) on Kaggle castwj. Runs ON ont. usage: bash kaggle/build_p1s2.sh
set -euo pipefail
ACC=castwj
cd ~/ncrna-annot
DS=kaggle/ds_p1s2_$ACC; rm -rf $DS && mkdir -p $DS
cp data/p05_feat/pca.npz data/p1/chr18_codes.npy $DS/
zcat data/human/gencode.v50.annotation.gtf.gz | awk '$1=="chr18"' > $DS/chr18_v50_full.gtf
tar czf $DS/src.tgz src
printf '{"title": "ncrna-p1s2-data", "id": "%s/ncrna-p1s2-data", "licenses": [{"name": "CC0-1.0"}]}' $ACC > $DS/dataset-metadata.json
HDR='import glob, json, os, subprocess, sys, tarfile
def find(name, hint=""):
    c = [x for x in glob.glob(f"/kaggle/input/**/{name}", recursive=True) if hint in x]
    return sorted(c)[0]
tarfile.open(find("src.tgz", "p1s2-data")).extractall("/tmp/code"); sys.path.insert(0, "/tmp/code/src")
import numpy as np, torch'
# 1) chr18 proxy-feature extraction
K=kaggle/ncrna-p1s2-chr18feat; mkdir -p $K
cat > $K/kernel.py <<PY
import subprocess, sys
subprocess.check_call([sys.executable, '-m', 'pip', 'install', '-q', 'transformers==4.57.6'])
$HDR
from ncrna_annot.features import generanno as G
pca = np.load(find('pca.npz')); mu = torch.tensor(pca['mean'], device='cuda'); comp = torch.tensor(pca['comp'], device='cuda')
tok, model = G.load()
def proj_window(codes):
    ids = torch.tensor([tok(G.codes_to_str(codes))['input_ids']], device='cuda')
    with torch.no_grad(), torch.autocast('cuda', dtype=G.half_dtype()):
        h = model(input_ids=ids, output_hidden_states=True).hidden_states[23][0, 1:].float()
    return ((h - mu) @ comp.T).half().cpu().numpy()
# consistency check vs Phase 0 cache (A10G bf16)
tc = np.load(find('train_codes.npy'), mmap_mode='r'); tf = np.load(find('feat_train.npy'), mmap_mode='r')
sims = []
for i in range(20):
    a = proj_window(np.asarray(tc[i])).astype(np.float32); b = np.asarray(tf[i], np.float32)
    sims.append(float(np.median((a * b).sum(1) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1) + 1e-9))))
json.dump({'median_token_cosine_per_window': sims, 'median': float(np.median(sims))}, open('/kaggle/working/consistency.json', 'w'), indent=1)
print('CONSISTENCY', np.median(sims), flush=True)
c = np.load(find('chr18_codes.npy')); starts = G.chromosome_tiles(len(c))
out = np.lib.format.open_memmap('/kaggle/working/chr18_feat.npy', 'w+', np.float16, (len(c) // 6, 512))
for k, s in enumerate(starts):
    lo, hi = G.kept_tokens(s, starts); out[s // 6 + lo:s // 6 + hi] = proj_window(c[s:s + G.WIN])[lo:hi]
    if k % 500 == 0: print('tile', k, len(starts), flush=True)
out.flush(); print('DONE', flush=True)
PY
printf '{"id": "%s/ncrna-p1s2-chr18feat", "title": "ncrna-p1s2-chr18feat", "code_file": "kernel.py", "language": "python", "kernel_type": "script", "is_private": "true", "enable_gpu": "true", "enable_tpu": "false", "enable_internet": "true", "machine_shape": "", "dataset_sources": ["%s/ncrna-p1s-data", "%s/ncrna-p1s2-data"], "competition_sources": [], "kernel_sources": [], "model_sources": []}' $ACC $ACC $ACC > $K/kernel-metadata.json
# 2) training/evaluation kernels: arm x seed
for spec in r:v50uev:0 s2:v50uevlnc:1 s4:v50f:0; do
  IFS=: read arm bits lnc <<< "$spec"
  for seed in 42 123 456; do
    K=kaggle/ncrna-p1s2-$arm-$seed; mkdir -p $K
    cat > $K/kernel.py <<PY
import subprocess, sys
subprocess.check_call([sys.executable, '-m', 'pip', 'install', '-q', 'numba'])
$HDR
from ncrna_annot.train import screen2
tr = dict(feats=np.load(find('feat_train.npy'), mmap_mode='r'), codes=np.load(find('train_codes.npy')), bits=np.load(find('train_bits_$bits.npy')))
va = dict(feats=np.load(find('feat_val.npy'), mmap_mode='r'), codes=np.load(find('val_codes.npy')), bits=np.load(find('val_bits_$bits.npy')))
ch = dict(feats=np.load(find('chr18_feat.npy'), mmap_mode='r'), codes=np.load(find('chr18_codes.npy')), truth_gtf=find('chr18_v50_full.gtf'), chrom='chr18')
screen2.run('${arm^^}', $seed, tr, va, ch, '/kaggle/working/p1s2_${arm}_${seed}', lnc_only=bool($lnc))
PY
    printf '{"id": "%s/ncrna-p1s2-%s-%s", "title": "ncrna-p1s2-%s-%s", "code_file": "kernel.py", "language": "python", "kernel_type": "script", "is_private": "true", "enable_gpu": "true", "enable_tpu": "false", "enable_internet": "true", "machine_shape": "", "dataset_sources": ["%s/ncrna-p1s-data", "%s/ncrna-p1s2-data"], "competition_sources": [], "kernel_sources": ["%s/ncrna-p1s2-chr18feat"], "model_sources": []}' $ACC $arm $seed $arm $seed $ACC $ACC $ACC > $K/kernel-metadata.json
  done
done
ls -d kaggle/ncrna-p1s2-*; du -sh $DS
