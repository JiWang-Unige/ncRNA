#!/bin/bash
# v3.0 test sequences (protocol section 18, execution detail 9), CPU only, no model involved: the C-group FASTA
# (Ensembl 116, data/v30/ensembl) -> data/<species>/<chrom>.fa with the header ">chrom" and MD5SUMS, as
# fetch_species.sh did for the P1.6 panel -> prep_panel.py codes (data/p1/<name>_codes.npy).
cd ~/ncrna-annot || exit 1
PY=~/envs/generanno/bin/python
names=()
for short in rat cow xtrop; do
  c=$(python3 -c "import json; print(json.load(open('data/v30/C_CHROMS.json'))['$short']['chrom'])")
  f=$(ls data/v30/ensembl/${short}_*.dna.*.$c.fa.gz)
  mkdir -p data/$short
  if [ ! -s data/$short/$c.fa ]; then
    zcat "$f" | sed '1s/^>\([^ ]*\).*/>\1/' > data/$short/$c.fa
    (cd data/$short && md5sum $c.fa > MD5SUMS)
  fi
  names+=("${short}_$c")
done
nice -n 19 $PY scripts/phase1/prep_panel.py --only "$(IFS=,; echo "${names[*]}")"
