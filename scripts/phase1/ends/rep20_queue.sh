#!/bin/bash
# PLAN_REP20 (frozen 2026-11-14): end analyses on chr20 with existing models, on ont.
#   1) re-extract chr20 GENERanno features (p1_feat_ont.sh, own done-dir; GPU gate per shard)
#   2) stage 0: v2.0 boundary logits on chr20;  3) the four END probes (ep04.pt), forward only
#   4) end_metrics on chr20 (CPU);  5) delete the chr20 feature cache again (disk: ont is ~98% full)
#   nohup bash scripts/phase1/ends/rep20_queue.sh > outputs/PHASE1/ENDS/rep20_queue.log 2>&1 &
set -u
cd ${NCRNA_ROOT:-~/ncrna-annot}
export NCRNA_CPUS=${NCRNA_CPUS:-2-31} NCRNA_NICE=${NCRNA_NICE:-10}
log() { echo "$(date -u +%Y-%m-%dT%H:%MZ) $*"; }
R=scripts/phase1/run_formal_ont.sh
F=data/p1_feat/chr20.npy
free_gb=$(df -BG --output=avail data/p1_feat | tail -1 | tr -dc 0-9)
if [ ! -e $F ] && [ "$free_gb" -lt 120 ]; then log "ABORT: only ${free_gb} GB free (< 120 GB, PLAN_REP20 section 3)"; exit 1; fi
log "free ${free_gb} GB"
P1_DONE=data/p1_shards/ont_rep20 CHROM_ONLY=1 NAMES=chr20 bash scripts/phase1/p1_feat_ont.sh || { log "feature extraction failed"; exit 1; }
# identity check against the shards recorded at the P1.5 read
~/envs/generanno/bin/python - <<'PY'
import json, glob, os
old = {os.path.basename(f): json.load(open(f)).get("sha256") for f in glob.glob("data/p1_shards/ont/chr20_t*.json")}
new = {os.path.basename(f): json.load(open(f)).get("sha256") for f in glob.glob("data/p1_shards/ont_rep20/chr20_t*.json")}
same = sum(old.get(k) == v for k, v in new.items())
rep = {"old": len(old), "new": len(new), "sha256_identical": same}
print("SHA_CHECK", json.dumps(rep))
json.dump(rep, open("outputs/PHASE1/ENDS/rep20_sha_check.json", "w"))
PY
bash $R stage0_chr20 --script ends/dump_bnd.py --chrom chr20 || exit 1
for p in ann-42 ev-42 ann-123 ev-123; do
  arm=${p%-*}; seed=${p#*-}
  bash $R end-$p --script ends/train_end_probe.py --arm $arm --seed $seed --forward-only --chrom chr20 || exit 1
done
PYTHONPATH=src nice -n 10 ~/envs/generanno/bin/python scripts/phase1/ends/end_metrics.py --chrom chr20 || exit 1
rm -f $F && log "deleted $F"
log "QUEUE_END"
