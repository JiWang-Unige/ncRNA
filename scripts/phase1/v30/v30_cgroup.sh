#!/bin/bash
# v3.0 C group (protocol 18 detail 2): after the stage-1 downloads, choose the test chromosomes by rule
# (select_c_chroms.py -> data/v30/C_CHROMS.json), then download their FASTA (fetch_v30.sh stage2).
cd ~/ncrna-annot
until grep -qE "STAGE1 done|FAIL" data/v30/fetch.log; do sleep 60; done
grep -q "STAGE1 done" data/v30/fetch.log || exit 1
nice -n 19 ~/envs/generanno/bin/python scripts/phase1/v30/select_c_chroms.py > outputs/PHASE1/V30/select_c.log 2>&1 || { echo "$(date -u +%FT%TZ) SELECT_FAIL" >> data/v30/fetch.log; exit 1; }
nice -n 19 scripts/phase1/v30/fetch_v30.sh stage2 > data/v30_fetch_stage2.out 2>&1
