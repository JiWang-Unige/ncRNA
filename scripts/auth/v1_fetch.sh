#!/usr/bin/env bash
# AUTH v1 step 2 (docs/auth/v1/FEASIBILITY_V1.md section 5; author approved 2026-11-24): small public files for the
# GM12878 / H1 descriptive check.  No BAM.  ENCODE md5 from the portal metadata, FANTOM5 md5 from hg38_latest/basic/00MD5SUM.txt.
#   bash scripts/auth/v1_fetch.sh            -> data/v1/{gm12878,h1}/encode, data/v1/fantom5, data/v1/FETCH.tsv
set -euo pipefail
ROOT=${NCRNA_ROOT:-$HOME/ncrna-annot}
D=$ROOT/data/v1
F5=https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest
F5MD5=$ROOT/data/fantom5_stageB/00MD5SUM.txt
mkdir -p $D/gm12878/encode $D/h1/encode $D/fantom5
LOG=$D/FETCH.tsv
[ -f $LOG ] || printf "file\tsource\tmd5_expected\tmd5_got\tbytes\tutc\n" > $LOG

get() {   # dest url md5
  local dest=$1 url=$2 md5=$3
  if [ -s "$dest" ] && [ "$(md5sum < "$dest" | cut -d' ' -f1)" = "$md5" ]; then echo "have $dest"; return; fi
  curl -sSL --retry 3 --retry-delay 5 -o "$dest.part" "$url"
  local got; got=$(md5sum < "$dest.part" | cut -d' ' -f1)
  if [ -n "$md5" ] && [ "$got" != "$md5" ]; then echo "MD5 MISMATCH $dest $got != $md5" >&2; exit 1; fi
  mv "$dest.part" "$dest"
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "${dest#$ROOT/}" "$url" "${md5:-none}" "$got" "$(stat -c %s "$dest")" "$(date -u +%FT%TZ)" >> $LOG
  echo "ok $dest"
}

enc() {   # cell accession ext md5
  get $D/$1/encode/$2.$3 "https://www.encodeproject.org/files/$2/@@download/$2.$3" $4
}

# GM12878: PacBio Sequel II (ENCSR838WFC) reps 1-2 = truth; QuantSeq (ENCSR275KDL, 1 rep, 2 files); RAMPAGE TSS (ENCSR000AEI)
enc gm12878 ENCFF871ODI gtf.gz b8853e2217d92824b944c614e7be6622   # transcriptome annotations rep [1]
enc gm12878 ENCFF993PKH gtf.gz 50819f3646970a62b39777e9af319077   # transcriptome annotations rep [2]
enc gm12878 ENCFF079JHG tsv bb1f5ba6bc220e816524919428dca814   # transcript quantifications rep [1]
enc gm12878 ENCFF386OOG tsv ccd15fd0b6b01eeb508e08daf0e4dc42   # transcript quantifications rep [2]
enc gm12878 ENCFF259EZO bed.gz 37220abb8b39b91d90c317d41073bff2   # polyA sites rep [1]
enc gm12878 ENCFF448CNL bed.gz 329ebb22cd1127c16088d740afdefa77   # polyA sites rep [1]
enc gm12878 ENCFF871RLN bed.gz 1edb700cfd3d14a4dce71244020862ba   # transcription start sites rep [1, 2]
enc gm12878 ENCFF591KXZ bed.gz a11182fa0ffd739c5f35f675d9fcc355   # transcription start sites rep [1]
enc gm12878 ENCFF579JKU bed.gz 004855da1ad6259834b82786beaad368   # transcription start sites rep [2]
# H1: LRGASP PacBio (ENCSR271KEJ) reps 1-3; LRGASP QuantSeq (ENCSR813AOB) reps 1-3
enc h1 ENCFF721HVH gtf.gz fd870491c4c564ecaa707f11e7902df2   # transcriptome annotations rep [1]
enc h1 ENCFF814SOM gtf.gz 9d80aeac07aa65cd1e94926e36e86a94   # transcriptome annotations rep [2]
enc h1 ENCFF236CGS gtf.gz 19ca5c67730ed24889bd1967bad446b3   # transcriptome annotations rep [3]
enc h1 ENCFF580WQM tsv db8923f23707ffc442f6702dac2c66dd   # transcript quantifications rep [1]
enc h1 ENCFF563JZI tsv c1c7dc762fcf2a35829567d760ad49a0   # transcript quantifications rep [2]
enc h1 ENCFF014CCY tsv 901c5c335d6914c37a42aebfec074fdc   # transcript quantifications rep [3]
enc h1 ENCFF970YYF bed.gz 1f2bb3eacc9e7fbfc3f3b8f9d1e580b9   # polyA sites rep [1]
enc h1 ENCFF138ZWV bed.gz e9045c35100a635e7d550bbfe8df3359   # polyA sites rep [2]
enc h1 ENCFF215NMJ bed.gz 0a9de47743884a6bb8f4e9c4851f5849   # polyA sites rep [3]

# FANTOM5 hg38 CTSS (3 biological replicates each); md5 from 00MD5SUM.txt (paths there are single-encoded)
f5() {  # cell dir encoded_name short
  local rel="$2/$(python3 -c "import sys,urllib.parse;print(urllib.parse.unquote(sys.argv[1]))" "$3")"
  local md5; md5=$(awk -v p="$rel" '$2==p{print $1}' $F5MD5)
  [ -n "$md5" ] || { echo "no md5 for $rel" >&2; exit 1; }
  get $D/fantom5/$1_$4.ctss.bed.gz "$F5/basic/$2/$3" $md5
}
f5 gm12878 human.cell_line.hCAGE "B%2520lymphoblastoid%2520cell%2520line%253a%2520GM12878%2520ENCODE%252c%2520biol_rep1.CNhs12331.10821-111C2.hg38.nobarcode.ctss.bed.gz" rep1
f5 gm12878 human.cell_line.hCAGE "B%2520lymphoblastoid%2520cell%2520line%253a%2520GM12878%2520ENCODE%252c%2520biol_rep2.CNhs12332.10822-111C3.hg38.nobarcode.ctss.bed.gz" rep2
f5 gm12878 human.cell_line.hCAGE "B%2520lymphoblastoid%2520cell%2520line%253a%2520GM12878%2520ENCODE%252c%2520biol_rep3.CNhs12333.10823-111C4.hg38.nobarcode.ctss.bed.gz" rep3
f5 h1 human.timecourse.hCAGE "H1%2520embryonic%2520stem%2520cells%2520differentiation%2520to%2520CD34%252b%2520HSC%252c%2520day00%252c%2520biol_rep1.CNhs14067.13523-145F4.hg38.nobarcode.ctss.bed.gz" rep1
f5 h1 human.timecourse.hCAGE "H1%2520embryonic%2520stem%2520cells%2520differentiation%2520to%2520CD34%252b%2520HSC%252c%2520day00%252c%2520biol_rep2.CNhs14068.13524-145F5.hg38.nobarcode.ctss.bed.gz" rep2
f5 h1 human.timecourse.hCAGE "H1%2520embryonic%2520stem%2520cells%2520differentiation%2520to%2520CD34%252b%2520HSC%252c%2520day00%252c%2520biol_rep3.CNhs13964.13525-145F6.hg38.nobarcode.ctss.bed.gz" rep3

# FANTOM5 hg38 DPI peak set (no published md5: md5 recorded in FETCH.tsv)
get $D/fantom5/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz "$F5/extra/CAGE_peaks/hg38_fair%2Bnew_CAGE_peaks_phase1and2.bed.gz" ""
echo done
