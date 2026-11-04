#!/bin/bash
# Batch-2 downloads on ont (card batch2-downloads, user 2026-11-04 "1，全部批准"; SCREEN_SUBJOURNAL.md).  Downloads only,
# each source in its own directory under data/batch2/, then size + sha256 into data/batch2/MANIFEST.json.  Nothing here
# runs or unpacks the downloaded files.  nice 19.  SpliceTransformer is NOT fetched: its weights are only on Google
# Drive, which the card did not cover (asked separately).
cd ~/ncrna-annot || exit 1
B=data/batch2
mkdir -p $B/snaptron $B/pangolin $B/mouse $B/chains
get() {  # <path under B> <url> [curl opts]
  local n=$1 u=$2 i; shift 2
  [ -s "$B/$n" ] && { echo "have $n"; return 0; }
  for i in 1 2 3; do
    nice -n 19 curl -sSfL --retry 2 --max-time 1800 "$@" -o "$B/$n.part" "$u" && { mv "$B/$n.part" "$B/$n"; echo "ok $n"; return 0; }
    sleep 30
  done
  rm -f "$B/$n.part"; echo "FAIL $n"; return 1
}
get pangolin/Pangolin-5cf94b8db938c658391b4305cd7ce33297d44ff7.tar.gz "https://github.com/tkzeng/Pangolin/archive/5cf94b8db938c658391b4305cd7ce33297d44ff7.tar.gz"
get mouse/mm39_chr19.fa.gz "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/chromosomes/chr19.fa.gz"
get mouse/gencode.vM39.annotation.gtf.gz "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M39/gencode.vM39.annotation.gtf.gz"
get chains/GCF_016699485.2ToHg38.over.chain.gz "https://hgdownload.soe.ucsc.edu/hubs/GCF/016/699/485/GCF_016699485.2/liftOver/GCF_016699485.2ToHg38.over.chain.gz"
get chains/GCF_016699485.2.chromAlias.txt "https://hgdownload.soe.ucsc.edu/hubs/GCF/016/699/485/GCF_016699485.2/GCF_016699485.2.chromAlias.txt"
get chains/danRer11ToHg38.over.chain.gz "https://hgdownload.soe.ucsc.edu/goldenPath/danRer11/liftOver/danRer11ToHg38.over.chain.gz"
# same Snaptron query and filter as bioval_fetch.sh (gtexv2, samples_count >= 5), 1 Mb windows
F="snaptron_id,chromosome,start,end,length,strand,annotated,left_motif,right_motif,samples_count,coverage_sum"
for spec in chr18:80373285 chr7:159345973; do
  c=${spec%%:*}; L=${spec##*:}
  for s in $(seq 0 1000000 $L); do
    e=$((s + 1000000)); [ $s -eq 0 ] && s=1
    get snaptron/${c}_$s.tsv "https://snaptron.cs.jhu.edu/gtexv2/snaptron?regions=$c:$s-$e&rfilter=samples_count>:5&fields=$F"
    sleep 2
  done
done
python3 - <<'PY'
import glob, hashlib, json, os, time
B = "data/batch2"
files = sorted(f for f in glob.glob(f"{B}/**/*", recursive=True) if os.path.isfile(f) and not f.endswith((".part", "MANIFEST.json")))
man = {"fetched": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "card": "batch2-downloads",
       "files": {os.path.relpath(f, B): {"bytes": os.path.getsize(f), "sha256": hashlib.sha256(open(f, "rb").read()).hexdigest()}
                 for f in files}}
json.dump(man, open(f"{B}/MANIFEST.json", "w"), indent=1)
print(len(files), "files", sum(v["bytes"] for v in man["files"].values()), "bytes")
PY
echo FETCH_END
