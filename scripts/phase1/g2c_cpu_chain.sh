#!/bin/bash
# G2-0 candidate scores on ont CPU for one G1B seed: 4 shards per split, then merge.  usage: g2c_cpu_chain.sh <seed>
cd ~/ncrna-annot
S=$1; CK=outputs/PHASE1/P1G/g1b-$S/p1g_g1b_$S/best.pt; O=data/p1g/g2c_$S
for sp in val train; do
  for i in 0 1 2 3; do PYTHONPATH=src nice -n 10 ~/envs/genatator/bin/python -u scripts/phase1/g2c_cpu.py $CK $O $sp $i/4 2 > logs/g2c_${S}_${sp}_$i.log 2>&1 & done
  wait
  PYTHONPATH=src ~/envs/genatator/bin/python - $O $sp <<PY
import sys, numpy as np
o, sp = sys.argv[1:3]
parts = [np.load(f"{o}/{sp}_cand_scores.part{i}.npy") for i in range(4)]
n = sum(len(p) for p in parts)
out = np.zeros((n,) + parts[0].shape[1:], np.float32)
for i, p in enumerate(parts):
    out[i::4] = p
np.save(f"{o}/{sp}_cand_scores.npy", out)
print("merged", sp, out.shape, flush=True)
PY
done
echo G2C_CPU_DONE
