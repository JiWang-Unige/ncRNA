#!/bin/bash
# v1.6 G2-0 (protocol section 13): candidate scores from the frozen G1B checkpoint of the same seed, for the
# 2,000 training and 300 validation windows (GPU).  Runs ON ont.  usage: bash kaggle/build_p1g_g2c.sh
# Output per seed: {train,val}_cand_scores.npy float32 [n, 24576, 8] (decode/g2cand.py channel order).
set -euo pipefail
ACC=castwj
cd ~/ncrna-annot
tar czf /tmp/p1g_src_g2c.tgz --exclude=__pycache__ src
SRC64=$(base64 -w0 /tmp/p1g_src_g2c.tgz)
for seed in 42 123 456; do
  K=kaggle/ncrna-p1g-g2c-$seed; mkdir -p $K
  cat > $K/kernel.py <<PY
import base64, glob, io, sys, tarfile, time
def find(name, hint=''):
    c = [x for x in glob.glob(f'/kaggle/input/**/{name}', recursive=True) if hint in x]
    return sorted(c)[0]
tarfile.open(fileobj=io.BytesIO(base64.b64decode('$SRC64'))).extractall('/tmp/code'); sys.path.insert(0, '/tmp/code/src')
import numpy as np, torch
from ncrna_annot.decode import g2cand
from ncrna_annot.grammar import g1
from ncrna_annot.heads.c0 import C0BHead
ck = torch.load(find('best.pt', 'ncrna-p1g-g1b-$seed/'), map_location='cuda')
head, trans, coef = C0BHead(512).cuda(), g1.G1Transitions().cuda(), g1.BoundaryCoef().cuda()
head.load_state_dict(ck['head']); trans.load_state_dict(ck['trans']); coef.load_state_dict(ck['coef']); head.eval()
T = trans().detach().double()
import copy; coef64 = copy.deepcopy(coef).double()
for sp in ('val', 'train'):
    feats = np.load(find(f'feat_{sp}.npy'), mmap_mode='r'); codes = np.load(find(f'{sp}_codes.npy'))
    out = np.lib.format.open_memmap(f'/kaggle/working/{sp}_cand_scores.npy', 'w+', np.float32, (len(codes), codes.shape[1], 8))
    t0 = time.time()
    for i in range(0, len(codes), 4):
        f = torch.as_tensor(np.asarray(feats[i:i + 4], np.float32), device='cuda')
        c = torch.as_tensor(np.asarray(codes[i:i + 4]).astype(np.int64), device='cuda')
        with torch.no_grad():
            p, b = head(f, c)
            edges = g1.boundary_edges(b.double(), coef64)
        gates = torch.as_tensor(np.stack([g1.gate_mask(x) for x in np.asarray(codes[i:i + 4])]), device='cuda').double()
        m = g1.marginals(p.double(), T, gates, edges=edges)   # float64: true-site posteriors go far below 1e-7
        for j in range(len(f)):
            out[i + j] = g2cand.scores(m[j].cpu().numpy(), b[j].float().cpu().numpy())
        if i % 400 == 0: print(sp, i, round(time.time() - t0), flush=True)
    out.flush()
print('G2C_DONE', flush=True)
PY
  printf '{"id": "%s/ncrna-p1g-g2c-%s", "title": "ncrna-p1g-g2c-%s", "code_file": "kernel.py", "language": "python", "kernel_type": "script", "is_private": "true", "enable_gpu": "true", "enable_tpu": "false", "enable_internet": "true", "machine_shape": "", "dataset_sources": ["%s/ncrna-p1s-data"], "competition_sources": [], "kernel_sources": ["%s/ncrna-p1g-g1b-%s"], "model_sources": []}' $ACC $seed $seed $ACC $ACC $seed > $K/kernel-metadata.json
  python3 -m py_compile $K/kernel.py
done
ls -d kaggle/ncrna-p1g-g2c-*
