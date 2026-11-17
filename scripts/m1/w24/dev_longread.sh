#!/bin/bash
# M1 weeks 2-4: full-transcriptome baselines on one chromosome from PacBio rep1 (ENCFF373TKM), reference-free
# (no annotation given, as our method): StringTie -L and IsoQuant de novo.  CPU only, nice 19, cores 2-31.
#   dev_longread.sh chr18 -> outputs/M1/W24/longread/<chrom>/{stringtie.gtf, isoquant/...}
set -eu
C=$1
cd ~/ncrna-annot
E=~/envs/lrtools/bin
O=outputs/M1/W24/longread/$C; mkdir -p $O
B=data/wtc11/encode/ENCFF373TKM.bam
[ -s $B.bai ] || nice -n 19 $E/samtools index -@ 4 $B
nice -n 19 $E/samtools view -b -@ 4 $B $C > $O/rep1.bam
nice -n 19 $E/samtools index $O/rep1.bam
echo "reads $($E/samtools view -c $O/rep1.bam)"
nice -n 19 taskset -c 2-31 $E/stringtie -L -p 8 -o $O/stringtie.gtf $O/rep1.bam
echo "stringtie transcripts $(awk '$3=="transcript"' $O/stringtie.gtf | wc -l)"
rm -rf $O/isoquant
nice -n 19 taskset -c 2-31 $E/python $E/isoquant.py --reference data/human/$C.fa --bam $O/rep1.bam --data_type pacbio_ccs \
  --threads 8 --output $O/isoquant --prefix rep1 > $O/isoquant.log 2>&1
ls $O/isoquant/rep1/ | head -20
