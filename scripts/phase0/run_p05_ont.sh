#!/bin/bash
# Run the 0.5 heads on ont once (a) the feature pull from Modal finished and (b) the GPU has no
# other compute process (lab rule: one GPU job at a time). Polls every 5 min.
cd ~/ncrna-annot
F=data/p05_feat
until [ -s $F/feat_chr20.npy ] && ! pgrep -f "modal volume get" >/dev/null && \
      [ -z "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader)" ]; do
  echo "$(date -Is) waiting (feat_chr20: $(ls -la $F/feat_chr20.npy 2>/dev/null | awk '{print $5}'), gpu apps: $(nvidia-smi --query-compute-apps=pid --format=csv,noheader | wc -l))"
  sleep 300
done
echo "$(date -Is) start"
~/envs/genatator/bin/python -u scripts/phase0/train_05_local.py
echo "$(date -Is) done rc=$?"
