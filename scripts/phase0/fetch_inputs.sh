#!/bin/bash
# Phase 0 inputs on ont (CPU + network only):
#  - GENATATOR pipeline repo + 4 stage models (HF snapshot, pinned revision recorded)
#  - hg38 chr20/chr21 FASTA (UCSC) and GENCODE v41 / v50 basic+comprehensive GTF (EBI)
set -euo pipefail
ROOT=${ROOT:-$HOME/ncrna-annot}
PY=${PY:-$HOME/envs/genatator/bin/python}
D=$ROOT/data/human
mkdir -p "$D" "$ROOT/models" "$ROOT/logs"
exec > >(tee -a "$ROOT/logs/fetch_inputs.log") 2>&1
echo "== $(date -Is) fetch inputs"

"$PY" - "$ROOT/models" <<'EOF'
import json, sys
from huggingface_hub import snapshot_download, HfApi
out = sys.argv[1]
repos = ["AIRI-Institute/genatator-pipeline",
         "AIRI-Institute/genatator-moderngena-base-multispecies-edge-model",
         "AIRI-Institute/genatator-moderngena-base-multispecies-region-model",
         "AIRI-Institute/genatator-caduceus-ps-multispecies-transcript-type",
         "AIRI-Institute/genatator-caduceus-ps-multispecies-segmentation"]
api, rec = HfApi(), {}
for r in repos:
    sha = api.model_info(r).sha
    path = snapshot_download(r, revision=sha, local_dir=f"{out}/{r.split('/')[1]}")
    rec[r] = {"revision": sha, "path": path}
    print(r, sha)
json.dump(rec, open(f"{out}/MODEL_REVISIONS.json", "w"), indent=1)
EOF

cd "$D"
for c in chr20 chr21; do
  [ -s $c.fa ] || { curl -sSfL https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/$c.fa.gz | gunzip > $c.fa; }
done
G=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human
for v in 41 50; do
  f=gencode.v$v.annotation.gtf.gz
  [ -s $f ] || curl -sSfLO $G/release_$v/$f
done
md5sum chr20.fa chr21.fa gencode.v*.gtf.gz > MD5SUMS
echo "== $(date -Is) done"
