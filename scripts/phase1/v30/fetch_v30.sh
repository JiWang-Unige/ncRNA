#!/bin/bash
# v3.0 downloads (protocol section 18, execution detail 1), ont, CPU / network only.
#   fetch_v30.sh stage1   mm39 genome, mm10 chrom sizes, mm10 -> mm39 chain, Snaptron srav1m junctions (mm10, filter as
#                         evidence.SNAP) for the mouse training chromosomes, Ensembl 116 GTFs of rat / cow / X. tropicalis
#   fetch_v30.sh stage2   Ensembl 116 FASTA of the C-group chromosomes in data/v30/C_CHROMS.json (chosen by rule first)
# Every file's sha256 goes into data/v30/MANIFEST.json (appended per stage).
set -uo pipefail
cd ~/ncrna-annot || exit 1
B=data/v30
mkdir -p $B/mouse $B/srav1m $B/chains $B/ensembl
log() { echo "$(date -u +%FT%TZ) $*" | tee -a $B/fetch.log; }
get() {  # relative_out url
  [ -s "$B/$1" ] && { log "have $1"; return 0; }
  curl -sSfL --retry 3 --max-time 3600 -o "$B/$1.part" "$2" && mv "$B/$1.part" "$B/$1" && log "got $1 $(stat -c %s $B/$1)" \
    || { log "FAIL $1 $2"; return 1; }
}
manifest() {  # stage file...
  python3 - "$@" <<'EOF'
import hashlib, json, os, sys, time
p = "data/v30/MANIFEST.json"
m = json.load(open(p)) if os.path.exists(p) else {"card": "v30-freeze", "files": {}}
for f in sys.argv[2:]:
    h = hashlib.sha256(open(f, "rb").read()).hexdigest()
    m["files"][f] = {"sha256": h, "bytes": os.path.getsize(f), "stage": sys.argv[1]}
m.setdefault("fetched", {})[sys.argv[1]] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
json.dump(m, open(p, "w"), indent=1)
EOF
}
E=https://ftp.ensembl.org/pub/release-116
SPECIES="rattus_norvegicus:rat bos_taurus:cow xenopus_tropicalis:xtrop"
TRAIN_CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chrX"
F="chromosome,start,end,strand,left_motif,right_motif,samples_count,coverage_sum"

case "${1:-}" in
stage1)
  log "STAGE1 start"
  get mouse/mm39.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz || exit 1
  get mouse/mm10.chrom.sizes https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes || exit 1
  get chains/mm10ToMm39.over.chain.gz https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm39.over.chain.gz || exit 1
  for c in $TRAIN_CHROMS; do
    n=$(awk -v c=$c '$1 == c {print $2}' $B/mouse/mm10.chrom.sizes)
    get srav1m/$c.tsv "https://snaptron.cs.jhu.edu/srav1m/snaptron?regions=$c:1-$n&rfilter=samples_count>:99&rfilter=coverage_sum>:499&fields=$F" || exit 1
  done
  for sp in $SPECIES; do
    d=${sp%%:*}; n=${sp##*:}
    g=$(curl -s $E/gtf/$d/ | grep -oE 'href="[^"]+\.116\.gtf\.gz"' | grep -v -E 'abinitio|chr\.gtf|chr_patch' | head -1 | cut -d'"' -f2)
    [ -n "$g" ] || { log "FAIL no GTF listed for $d"; exit 1; }
    get ensembl/${n}_$g "$E/gtf/$d/$g" || exit 1
  done
  manifest stage1 $B/mouse/mm39.fa.gz $B/mouse/mm10.chrom.sizes $B/chains/mm10ToMm39.over.chain.gz $B/srav1m/*.tsv $B/ensembl/*.gtf.gz
  log "STAGE1 done"
  ;;
stage2)
  [ -s $B/C_CHROMS.json ] || { log "FAIL stage2 needs data/v30/C_CHROMS.json"; exit 1; }
  log "STAGE2 start"
  files=()
  for sp in $SPECIES; do
    d=${sp%%:*}; n=${sp##*:}
    c=$(python3 -c "import json,sys; print(json.load(open('$B/C_CHROMS.json'))['$n']['chrom'])")
    f=$(curl -s $E/fasta/$d/dna/ | grep -oE "href=\"[^\"]+dna\.(primary_assembly|chromosome)\.$c\.fa\.gz\"" | head -1 | cut -d'"' -f2)
    [ -n "$f" ] || { log "FAIL no FASTA listed for $d $c"; exit 1; }
    get ensembl/${n}_$f "$E/fasta/$d/dna/$f" || exit 1
    files+=("$B/ensembl/${n}_$f")
  done
  manifest stage2 "${files[@]}"
  log "STAGE2 done"
  ;;
*) echo "usage: $0 stage1|stage2"; exit 2 ;;
esac
