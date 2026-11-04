#!/bin/bash
# v2.0 chr9 read on ont (protocol section 17 + amendment 1; configs/v20.json committed before this runs):
# chr18 reproduction (reads + REPRO_OK) -> chr9 features -> chr9 single read -> v2.0 scoring -> section 4.2
# replications (Snaptron chr9; SpliceAI / own / Pangolin chr9 scores; s1_analyze; S1b A1-A3 x 3 runs; S3) -> verdicts
# -> chr9 feature cache deleted.  GPU: gpu_gate.sh before every job (TE containers included), <= 30% in total.
cd ~/ncrna-annot || exit 1
. scripts/phase1/gpu_gate.sh
O=outputs/PHASE1/V20; REP=$O/rep; PY=~/envs/generanno/bin/python; export PYTHONPATH=src
C=configs/v20.json; N="nice -n 19 taskset -c 2-31"
mkdir -p $REP/logs $REP/S1B data/v20_chr9/snaptron
log() { echo "$(date -u +%FT%TZ) $*" >> $O/chr9_queue.log; }
gstep() {  # name done-file command... (one GPU job, gate before, retry on exit 75)
  local name=$1 done=$2 rc; shift 2
  [ -s "$done" ] && { log "$name skip"; return 0; }
  while true; do
    gpu_wait >> $O/chr9_queue.log
    log "$name start"
    $N env OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 "$@" > $REP/logs/$name.log 2>&1; rc=$?
    [ $rc -eq 75 ] && { sleep 120; continue; }
    log "$name rc=$rc"; return $rc
  done
}
log START
# 1) chr18 reproduction
NCRNA_CPUS=2-31 NCRNA_NICE=19 NCRNA_GPU_FRAC=0.30 bash scripts/phase1/run_formal_ont.sh v20read-chr18 --script p15_read.py --config $C --target chr18 >> logs/formal_v20read-chr18.log 2>&1 || { log READ18_FAIL; exit 1; }
$N $PY scripts/phase1/p15_score.py --config $C --target chr18 > $O/repro_chr18.log 2>&1
[ -e $O/read/chr18/REPRO_OK ] || { log REPRO_FAIL; exit 1; }
log REPRO_OK
# 2) chr9 features (cache only with >= 250 GB free, protocol section 17 detail 6)
free=$(df -BG --output=avail ~ | tail -1 | tr -dc 0-9)
[ "$free" -ge 250 ] || { log "only ${free} GB free - stop (stream mode not implemented)"; exit 1; }
NAMES=chr9 CHROM_ONLY=1 NCRNA_GPU_FRAC=0.28 bash scripts/phase1/p1_feat_ont.sh >> $O/feat_chr9.log 2>&1
grep -q "ALL DONE fails=0" $O/feat_chr9.log || { log FEAT9_FAIL; exit 1; }
log "chr9 features ok"
# 3) the single chr9 read
NCRNA_CPUS=2-31 NCRNA_NICE=19 NCRNA_GPU_FRAC=0.30 bash scripts/phase1/run_formal_ont.sh v20read-chr9 --script p15_read.py --config $C --target chr9 >> logs/formal_v20read-chr9.log 2>&1 || { log READ9_FAIL; exit 1; }
log "chr9 read ok"
# 4) scoring (section 4.1)
$N env OMP_NUM_THREADS=8 $PY scripts/phase1/v20/v20_score.py --config $C --target chr9 > $O/score_chr9.log 2>&1 || { log SCORE_FAIL; exit 1; }
log "chr9 scored"
# 5) section 4.2 replications
F="snaptron_id,chromosome,start,end,length,strand,annotated,left_motif,right_motif,samples_count,coverage_sum"
for s in $(seq 0 1000000 138394717); do
  e=$((s + 1000000)); [ $s -eq 0 ] && s=1
  f=data/v20_chr9/snaptron/chr9_$s.tsv; [ -s $f ] && continue
  for i in 1 2 3; do nice -n 19 curl -sSfL --retry 2 --max-time 1800 -o $f.part "https://snaptron.cs.jhu.edu/gtexv2/snaptron?regions=chr9:$s-$e&rfilter=samples_count>:5&fields=$F" && mv $f.part $f && break; sleep 30; done
  sleep 2
done
$PY - <<'PY'
import glob, hashlib, json, os, time
fs = sorted(glob.glob("data/v20_chr9/snaptron/*.tsv"))
json.dump({"fetched": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "query": "gtexv2, samples_count >= 5, 1 Mb windows",
           "files": {os.path.basename(f): {"bytes": os.path.getsize(f), "sha256": hashlib.sha256(open(f, "rb").read()).hexdigest()} for f in fs}},
          open("data/v20_chr9/MANIFEST.json", "w"), indent=1)
print(len(fs), "files")
PY
log "snaptron chr9 ok"
G50=data/human/gencode.v50.annotation.gtf.gz
gstep sai9 $REP/sai_chr9.npy env NCRNA_GPU_FRAC=0.25 $PY scripts/phase1/splice_diag/sai_score.py --codes data/p1/chr9_codes.npy --truth $G50 --chrom chr9 --out $REP/sai_chr9.npy || exit 1
gstep own9 $REP/own_chr9.npy env NCRNA_GPU_FRAC=0.25 $PY scripts/phase1/splice_diag/own_splice.py score --codes data/p1/chr9_codes.npy --truth $G50 --chrom chr9 --out $REP/own_chr9.npy || exit 1
gstep pang9 $REP/pang_chr9.npy env NCRNA_GPU_FRAC=0.25 $PY scripts/phase1/s1/pangolin_score.py --codes data/p1/chr9_codes.npy --chrom chr9 --offsets outputs/PHASE1/S1/pang_chr20.npy.json --out $REP/pang_chr9.npy || exit 1
$N env OMP_NUM_THREADS=8 $PY scripts/phase1/s1/s1_analyze.py --chrom chr9 --codes data/p1/chr9_codes.npy --snaptron "data/v20_chr9/snaptron/chr9_*.tsv" --judge own \
  --scores spliceai=$REP/sai_chr9.npy own=$REP/own_chr9.npy pangolin=$REP/pang_chr9.npy --out $REP/S1_chr9.json > $REP/logs/s1.log 2>&1 || { log S1_FAIL; exit 1; }
log "R-S1 analysed"
for r in v20-g1b-42-ke v20-g1b-123-ke v20-g1b-456-ke; do
  for arm in "A1 --w 0 --span-oracle" "A2 --w 1 --oracle" "A3 --w 1 --oracle --span-oracle"; do
    set -- $arm; a=$1; shift
    while [ "$(jobs -rp | wc -l)" -ge 2 ]; do wait -n; done
    gstep s1b_${r}_$a $REP/S1B/${r}_$a/result.json env NCRNA_THREADS=2 NUMBA_NUM_THREADS=2 NCRNA_GPU_FRAC=0.15 \
      $PY scripts/phase1/splice_diag/splice_diag.py --config $C --target chr9 --run $r "$@" --out $REP/S1B/${r}_$a &
    sleep 90
  done
done
wait
log "R-S1b decoded"
$N env OMP_NUM_THREADS=4 $PY scripts/phase1/s3/s3_features.py --config $C --chrom chr9 --v20 --out $REP/S3_chr9.json > $REP/logs/s3.log 2>&1 || { log S3_FAIL; exit 1; }
$PY scripts/phase1/v20/v20_replicate.py > $REP/logs/replicate.log 2>&1 || { log REPLICATE_FAIL; exit 1; }
log "replications done"
rm -f data/p1_feat/chr9.npy && log "chr9 feature cache deleted"
log QUEUE_END
