#!/bin/bash
# Phase 0.6 small-ncRNA tools on ont: Infernal + tRNAscan-SE (bioconda) and current Rfam CMs.
set -euo pipefail
ROOT=${ROOT:-$HOME/ncrna-annot}
ENV=${ENV:-$HOME/envs/smallrna}
MM=$(command -v micromamba || echo "$HOME/micromamba/bin/micromamba")
[ -x "$MM" ] || MM=$(ls "$HOME"/micromamba/micromamba "$HOME"/bin/micromamba 2>/dev/null | head -1)
exec > >(tee -a "$ROOT/logs/install_smallrna_env.log") 2>&1
echo "== $(date -Is) micromamba=$MM env=$ENV"
[ -x "$ENV/bin/cmscan" ] || "$MM" create -y -p "$ENV" -c conda-forge -c bioconda infernal trnascan-se
"$ENV/bin/cmscan" -h | sed -n 2p
"$ENV/bin/tRNAscan-SE" -h 2>&1 | sed -n 2p || true

R=$ROOT/data/rfam; mkdir -p "$R"; cd "$R"
F=https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT
[ -s Rfam.cm ] || { curl -sSfL $F/Rfam.cm.gz | gunzip > Rfam.cm; }
[ -s Rfam.clanin ] || curl -sSfLO $F/Rfam.clanin
[ -s Rfam.cm.i1m ] || "$ENV/bin/cmpress" Rfam.cm
curl -sSfL $F/README | grep -iE 'release|version' | head -3 > RFAM_VERSION.txt || true
md5sum Rfam.cm Rfam.clanin > MD5SUMS
echo "== $(date -Is) done; Rfam $(cat RFAM_VERSION.txt | head -1)"
