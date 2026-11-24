#!/bin/bash
# AUTH post-freeze test queue (docs/auth/PLAN_AUTH.md section 7 step 5).  Run ONLY after the author has approved the
# freeze and outputs/AUTH/FROZEN has been written; every Python step re-checks FROZEN and refuses if RESULT.json exists.
#   1. FROZEN check + soft-masked test chromosome sequences (IsoQuant reference)
#   2. rep1 (ENCFF373TKM) StringTie -L / IsoQuant de novo per test chromosome, as scripts/m1/w24/dev_longread.sh
#   3. auth_score.py --mode test (the one-time read) -> outputs/AUTH/RESULT.json
# CPU only, nice 19, cores 2-31, THREADS (default 8) per tool; yields to TE jobs.
#   setsid nohup bash /home/ji/ncrna-annot/scripts/auth/auth_queue.sh > /home/ji/ncrna-annot/outputs/AUTH/queue.log 2>&1 < /dev/null &
set -euo pipefail
R=/home/ji/ncrna-annot
PY=/home/ji/envs/generanno/bin/python
E=/home/ji/envs/lrtools/bin
T=${THREADS:-8}
cd $R
echo "[$(date -u +%FT%TZ)] queue start"
nice -n 19 taskset -c 2-31 $PY $R/scripts/auth/auth_score.py --mode prep
B=$R/data/wtc11/encode/ENCFF373TKM.bam
[ -s $B.bai ]
for C in chr2 chr6 chr10 chr16 chr22; do
  O=$R/outputs/AUTH/longread/$C
  if [ -s $O/isoquant/rep1/rep1.transcript_models.gtf ] && [ -s $O/stringtie.gtf ]; then
    echo "[$(date -u +%FT%TZ)] $C long reads already done"; continue
  fi
  mkdir -p $O
  nice -n 19 taskset -c 2-31 $E/samtools view -b -@ 4 $B $C > $O/rep1.bam
  nice -n 19 taskset -c 2-31 $E/samtools index $O/rep1.bam
  echo "[$(date -u +%FT%TZ)] $C reads $($E/samtools view -c $O/rep1.bam)"
  nice -n 19 taskset -c 2-31 $E/stringtie -L -p $T -o $O/stringtie.gtf $O/rep1.bam
  echo "[$(date -u +%FT%TZ)] $C stringtie transcripts $(awk '$3=="transcript"' $O/stringtie.gtf | wc -l)"
  rm -rf $O/isoquant
  nice -n 19 taskset -c 2-31 $E/python $E/isoquant.py --reference $R/outputs/AUTH/seq/$C.fa --bam $O/rep1.bam \
    --data_type pacbio_ccs --threads $T --output $O/isoquant --prefix rep1 > $O/isoquant.log 2>&1
  [ -s $O/isoquant/rep1/rep1.transcript_models.gtf ]
  echo "[$(date -u +%FT%TZ)] $C isoquant done"
done
nice -n 19 taskset -c 2-31 $PY $R/scripts/auth/auth_score.py --mode test
echo "[$(date -u +%FT%TZ)] queue done"
