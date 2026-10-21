#!/bin/bash
cd ~/ncrna-annot
scripts/phase1/kaggle_output_resume.py castwj castwj/ncrna-p1s2-chr18feat chr18_feat.npy data/p1s2_chr18/chr18_feat.npy && sleep 0 || { echo DOWNLOAD_FAILED; exit 1; }
F=data/p1s2_chr18/chr18_feat.npy
ls -la $F
[ $(stat -c %s $F) -ge 10998470000 ] || { echo DOWNLOAD_INCOMPLETE; exit 1; }
O=outputs/PHASE1/DIAG_L2; mkdir -p $O
PYTHONPATH=src nice -n 10 ~/envs/genatator/bin/python -u scripts/phase1/diag_decode.py outputs/PHASE1/P1S2/r-123/p1s2_r_123/best.pt $F data/p1/chr18_codes.npy chr18 $O/r123_chr18.gff3 0 && zcat data/human/gencode.v50.annotation.gtf.gz | awk "\$1==\"chr18\"" > $O/chr18_v50.gtf && PYTHONPATH=src nice -n 10 ~/envs/genatator/bin/python -u scripts/phase1/diag_l2.py $O/chr18_v50.gtf chr18 $O/r123_chr18.gff3 $O/diag_r123.json > $O/diag_r123.out && echo DIAG_DONE
