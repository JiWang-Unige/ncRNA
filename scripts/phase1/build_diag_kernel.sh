#!/bin/bash
# v1.5 L2 diagnostic on Kaggle castwj (GPU potentials + numba S3 decode + diag_l2). Runs ON ont.
# usage: bash scripts/phase1/build_diag_kernel.sh <arm> <seed>   (e.g. r 123)
set -euo pipefail
ACC=castwj; ARM=$1; SEED=$2
cd ~/ncrna-annot
K=kaggle/ncrna-diag-l2-$ARM-$SEED; mkdir -p $K
python3 - "$ARM" "$SEED" > $K/kernel.py <<'PY'
import sys
arm, seed = sys.argv[1:3]
diag = open("scripts/phase1/diag_l2.py").read()
print(f'''import subprocess, sys, glob, json, os, tarfile, time
subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "numba"])
def find(name, hint=""):
    c = [x for x in glob.glob(f"/kaggle/input/**/{{name}}", recursive=True) if hint in x]
    return sorted(c)[0]
tarfile.open(find("src.tgz", "p1s2-data")).extractall("/tmp/code"); sys.path.insert(0, "/tmp/code/src")
open("/tmp/code/diag_l2.py", "w").write({diag!r})
import numpy as np, torch
from ncrna_annot.decode import chrom as C
from ncrna_annot.decode import genes as dec
from ncrna_annot.grammar import g1
from ncrna_annot.heads.c0 import C0Head
t0 = time.time()
feats = np.load(find("chr18_feat.npy"), mmap_mode="r"); codes = np.load(find("chr18_codes.npy"))
head, trans = C0Head(feats.shape[-1]).cuda(), g1.G1Transitions()
st = torch.load(find("best.pt", "p1s2-{arm}-{seed}"), map_location="cpu")
head.load_state_dict(st["head"]); trans.load_state_dict(st["trans"])
L = len(feats) * 6
pot = dec.chromosome_potentials(head, feats, codes[:L], device="cuda")
print("potentials s", round(time.time() - t0), flush=True)
genes = C.decode_blocks(pot, codes[:L], trans().detach().double().numpy(), 20, 50)
genes = [g for g in genes if g["type"] != "l" or sum(b - a + 1 for a, b in g["exons"]) >= 200]
out = "/kaggle/working/{arm}{seed}_chr18.gff3"
dec.write_gff3(genes, "chr18", out)
print("genes", len(genes), "lnc", sum(g["type"] == "l" for g in genes), "s", round(time.time() - t0), flush=True)
r = subprocess.run([sys.executable, "/tmp/code/diag_l2.py", find("chr18_v50_full.gtf"), "chr18", out,
                    "/kaggle/working/diag_{arm}{seed}.json"], env=dict(os.environ, PYTHONPATH="/tmp/code/src"),
                   capture_output=True, text=True)
open("/kaggle/working/diag_{arm}{seed}.out", "w").write(r.stdout + r.stderr)
print(r.stdout[-3000:], r.stderr[-3000:], flush=True)
print("DIAG_DONE" if r.returncode == 0 else "DIAG_FAILED", flush=True)
''')
PY
printf '{"id": "%s/ncrna-diag-l2-%s-%s", "title": "ncrna-diag-l2-%s-%s", "code_file": "kernel.py", "language": "python", "kernel_type": "script", "is_private": "true", "enable_gpu": "true", "enable_tpu": "false", "enable_internet": "true", "machine_shape": "", "dataset_sources": ["%s/ncrna-p1s2-data"], "competition_sources": [], "kernel_sources": ["%s/ncrna-p1s2-chr18feat", "%s/ncrna-p1s2-%s-%s"], "model_sources": []}' $ACC $ARM $SEED $ARM $SEED $ACC $ACC $ACC $ARM $SEED > $K/kernel-metadata.json
python3 -m py_compile $K/kernel.py && echo BUILT $K
