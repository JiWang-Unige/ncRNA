#!/bin/bash
# Phase 0.1b inputs: Ensembl (release 116) chromosome FASTA + full GTF for the animal panel.
# Chromosomes fixed in PROTOCOL v1.1 (chosen by size before any result).
set -euo pipefail
ROOT=${ROOT:-$HOME/ncrna-annot}
E=https://ftp.ensembl.org/pub
exec > >(tee -a "$ROOT/logs/fetch_species.log") 2>&1
REL=${REL:-116}   # pinned; current_gtf/ is not served any more
echo "== $(date -Is) Ensembl release $REL"

fetch() {  # species_dir  out_name  chrom...
  local sp=$1 name=$2; shift 2
  local D=$ROOT/data/$name; mkdir -p "$D"; cd "$D"
  local gtf; gtf=$(curl -s $E/release-$REL/gtf/$sp/ | grep -oE 'href="[^"]+\.[0-9]+\.gtf\.gz"' | grep -v -E 'abinitio|chr\.gtf|chr_patch' | head -1 | cut -d'"' -f2)
  [ -s "$gtf" ] || curl -sSfLO "$E/release-$REL/gtf/$sp/$gtf"
  local listing; listing=$(curl -s $E/release-$REL/fasta/$sp/dna/)
  for c in "$@"; do
    local f; f=$(echo "$listing" | grep -oE "href=\"[^\"]+dna\.(primary_assembly|chromosome)\.$c\.fa\.gz\"" | head -1 | cut -d'"' -f2)
    [ -s "$c.fa" ] || { curl -sSfL "$E/release-$REL/fasta/$sp/dna/$f" | gunzip | sed '1s/^>\([^ ]*\).*/>\1/' > "$c.fa"; }
    echo "$name $c $f $(grep -v '>' $c.fa | tr -d '\n' | wc -c) bp"
  done
  echo "$name gtf $gtf"
  md5sum ./*.fa "$gtf" > MD5SUMS
}
fetch gallus_gallus chicken 5
fetch danio_rerio zebrafish 22
fetch drosophila_melanogaster fly 2L
fetch caenorhabditis_elegans worm I
echo "== $(date -Is) done"
