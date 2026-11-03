#!/bin/bash
# Bio-validation evidence download on ont (card bioval-data, user 2026-11-03; docs/phase1_results/PREREG_BIOVAL.md,
# frozen in 7c7bc88 before this runs).  Downloads only, then size + sha256 into data/bioval/MANIFEST.json; nothing
# here reads the file contents.  nice 19.  PolyASite's TLS certificate has expired (noted in the review), so that one
# file is fetched with -k; its sha256 is recorded like every other file.
cd ~/ncrna-annot || exit 1
R=data/bioval/raw
mkdir -p $R/snaptron
get() {  # <name> <url> [curl opts]
  local n=$1 u=$2 i; shift 2
  for i in 1 2 3; do
    nice -n 19 curl -sSfL --retry 2 --max-time 1800 "$@" -o "$R/$n.part" "$u" && { mv "$R/$n.part" "$R/$n"; echo "ok $n"; return 0; }
    sleep 30
  done
  rm -f "$R/$n.part"; echo "FAIL $n"; return 1
}
get fantom5_cage_peaks.bed.gz "https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz"
get polyasite2_GRCh38.96.tsv.gz "https://polyasite.unibas.ch/download/atlas/2.0/GRCh38.96/atlas.clusters.2.0.GRCh38.96.tsv.gz" -k
get gtex_v9_flair_filter_transcripts.gtf.gz "https://storage.googleapis.com/adult-gtex/long-read-data/v9/long-read-RNA-seq/flair_filter_transcripts.gtf.gz"
get chess3.1.3.GRCh38.gtf.gz "https://github.com/chess-genome/chess/releases/download/v.3.1.3/chess3.1.3.GRCh38.gtf.gz"
get lncbook2.1_GRCh38.gtf.gz "https://ngdc.cncb.ac.cn/lncbook/files/lncRNA_LncBookv2.1_GRCh38.gtf.gz" --http1.1
get rnacentral_homo_sapiens.GRCh38.bed.gz "https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/genome_coordinates/bed/homo_sapiens.GRCh38.bed.gz"
get rnacentral_release_notes.txt "https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/release_notes.txt"
F="snaptron_id,chromosome,start,end,length,strand,annotated,left_motif,right_motif,samples_count,coverage_sum"
for s in $(seq 0 1000000 64000000); do
  e=$((s + 1000000)); [ $s -eq 0 ] && s=1
  get snaptron/chr20_$s.tsv "https://snaptron.cs.jhu.edu/gtexv2/snaptron?regions=chr20:$s-$e&rfilter=samples_count>:5&fields=$F"
  sleep 2
done
python3 - <<'EOF'
import glob, hashlib, json, os, time
R = "data/bioval/raw"
files = sorted(f for f in glob.glob(f"{R}/**/*", recursive=True) if os.path.isfile(f))
man = {"fetched": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "prereg_commit": "7c7bc88",
       "files": {os.path.relpath(f, R): {"bytes": os.path.getsize(f), "sha256": hashlib.sha256(open(f, "rb").read()).hexdigest()}
                 for f in files}}
json.dump(man, open("data/bioval/MANIFEST.json", "w"), indent=1)
print(len(files), "files", sum(v["bytes"] for v in man["files"].values()), "bytes")
EOF
