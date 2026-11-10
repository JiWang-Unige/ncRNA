#!/bin/bash
# resume after FEAT_FAIL rat_12 (random SystemError/segfaults on CPU 0 -> pin everything to E-cores 16-31, nice 19)
cd ~/ncrna-annot; O=outputs/PHASE1/V30; PY=~/envs/generanno/bin/python; export PYTHONPATH=src
log() { echo "$(date -u +%FT%TZ) $*" >> $O/test_queue_preview.log; }
log "RETRY rat_12 failed shards pinned to cpus 16-31 (random SystemError / segfaults on cpu 0; i9-13900K)"
NAMES=rat_12 CHROM_ONLY=1 NCRNA_GPU_FRAC=0.28 bash scripts/phase1/p1_feat_ont.sh >> $O/feat_rat_12_preview_retry.log 2>&1
grep -q "ALL DONE fails=0" $O/feat_rat_12_preview_retry.log || { log "FEAT_FAIL rat_12 retry"; exit 1; }
# determinism / silent-corruption check: re-extract one shard done before (t00200_00400) and compare sha256
. scripts/phase1/gpu_gate.sh; gpu_wait
CK=data/p1_shards/check_rat; rm -rf $CK; mkdir -p $CK
NCRNA_GPU_FRAC=0.28 $PY scripts/phase1/extract_p1_shard.py --out $CK --kind chrom --name rat_12 --codes data/p1/rat_12_codes.npy --tiles 200:400 > $CK/extract.log 2>&1
a=$(python3 -c "import json;print(json.load(open(\"data/p1_shards/ont/rat_12_t00200_00400.json\"))[\"sha256\"])")
b=$(sha256sum $CK/rat_12_t00200_00400.npy | cut -c1-64)
if [ "$a" = "$b" ]; then log "SHARD_RECHECK_OK rat_12_t00200_00400 $b"; else log "SHARD_RECHECK_MISMATCH rat_12_t00200_00400 old $a new $b"; rm -rf $CK; exit 1; fi
rm -rf $CK
exec bash scripts/phase1/v30/v30_test_queue.sh preview
