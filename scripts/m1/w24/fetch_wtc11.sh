#!/bin/bash
# M1 weeks 2-4: download WTC11 files listed in manifest.tsv (url, name, md5 or -, size) into data/wtc11, verify.
D=~/ncrna-annot/data/wtc11; mkdir -p $D/encode $D/geo; cd $D
log() { echo "$(date -u +%FT%TZ) $*" >> fetch.log; }
log START
while IFS=$'\t' read -r url name md5 size; do
  sub=encode; [[ $url == *ncbi* ]] && sub=geo
  f=$sub/$name
  if [ ! -s $f ] || [ "$(stat -c %s $f)" != "$size" ]; then
    curl -fsSL --retry 5 --retry-delay 20 --max-time 7200 -o $f.part "$url" && mv $f.part $f
  fi
  got_size=$(stat -c %s $f 2>/dev/null); got_md5=$(md5sum $f 2>/dev/null | cut -d' ' -f1)
  ok=1; [ "$got_size" = "$size" ] || ok=0; [ "$md5" = - ] || [ "$md5" = "$got_md5" ] || ok=0
  log "$name size $got_size/$size md5 $got_md5 ok=$ok"
done < manifest.tsv
log "END fails=$(grep -c 'ok=0' fetch.log)"
