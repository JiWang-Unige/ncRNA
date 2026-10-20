#!/bin/bash
# P1.2 smoke on Kaggle: dataset (src + 20 windows + chr18 first 1 Mb + chr18 v50 GTF) and a script kernel.
# Runs ON ont (data/p1 lives there). usage: bash kaggle/build_p1_smoke.sh
set -euo pipefail
cd ~/ncrna-annot
DS=kaggle/ds_p1_smoke; K=kaggle/ncrna-p1-smoke
rm -rf $DS && mkdir -p $DS $K
cp -r src $DS/src
~/envs/genatator/bin/python - <<'PY'
import numpy as np
d = 'data/p1'; o = 'kaggle/ds_p1_smoke'
np.save(f'{o}/train_codes.npy', np.load(f'{d}/train_codes.npy', mmap_mode='r')[:20])
np.save(f'{o}/train_bits.npy', np.load(f'{d}/train_bits_v50uev.npy', mmap_mode='r')[:20])
np.save(f'{o}/chr18_codes_1Mb.npy', np.load(f'{d}/chr18_codes.npy', mmap_mode='r')[:1_000_000])
PY
zcat data/human/gencode.v50.annotation.gtf.gz | awk '$1=="chr18"' | gzip > $DS/chr18_v50.gtf.gz
tar czf $DS/src.tgz -C $DS src && rm -rf $DS/src
printf '{"title": "ncrna-p1-smoke", "id": "jiunige/ncrna-p1-smoke", "licenses": [{"name": "CC0-1.0"}]}' > $DS/dataset-metadata.json
cat > $K/kernel.py <<'PY'
import glob, json, os, subprocess, sys, tarfile
subprocess.check_call([sys.executable, '-m', 'pip', 'install', '-q', 'transformers==4.57.6', 'numba'])
IN = os.path.dirname(glob.glob('/kaggle/input/**/train_codes.npy', recursive=True)[0])
tarfile.open(f'{IN}/src.tgz').extractall('/tmp/code'); sys.path.insert(0, '/tmp/code/src')
import numpy as np, torch
from ncrna_annot.features import generanno as G
from ncrna_annot.train import loop
tok, model = G.load()
codes = np.load(f'{IN}/train_codes.npy'); bits = np.load(f'{IN}/train_bits.npy')
feats = np.concatenate([G.features(tok, model, list(codes[i:i + 2])).cpu().numpy() for i in range(0, len(codes), 2)])
vc = np.load(f'{IN}/chr18_codes_1Mb.npy'); starts = G.chromosome_tiles(len(vc)); vf = np.zeros(((len(vc) // 6), G.DIM), np.float16)
for s in starts:
    lo, hi = G.kept_tokens(s, starts); vf[s // 6 + lo:s // 6 + hi] = G.features(tok, model, [vc[s:s + G.WIN]])[0, lo:hi].cpu().numpy()
del model; torch.cuda.empty_cache()
best = loop.train(feats, codes, bits, {'feats': vf, 'codes': vc[:len(vf) * 6], 'truth_gtf': f'{IN}/chr18_v50.gtf.gz', 'chrom': 'chr18'},
                  '/kaggle/working/smoke', epochs=1, batch=4)
print('BEST', best, flush=True)
PY
printf '{"id": "jiunige/ncrna-p1-smoke-run", "title": "ncrna-p1-smoke-run", "code_file": "kernel.py", "language": "python", "kernel_type": "script", "is_private": "true", "enable_gpu": "true", "enable_tpu": "false", "enable_internet": "true", "machine_shape": "", "dataset_sources": ["jiunige/ncrna-p1-smoke"], "competition_sources": [], "kernel_sources": [], "model_sources": []}' > $K/kernel-metadata.json
du -sh $DS
