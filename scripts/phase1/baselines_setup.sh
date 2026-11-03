#!/bin/bash
# chr20 baselines setup on ont (card baselines-chr20; docs/phase1_results/PREREG_BASELINES.md, 9608a3a).
# Tiberius: clone + uv env (Python 3.12, TensorFlow per its requirements).  SegmentNT: HF snapshot at a fixed revision.
# Network + CPU only, nice 19; no GPU.
cd ~/ncrna-annot || exit 1
mkdir -p ext
UV=~/.local/bin/uv
echo "$(date -u +%FT%TZ) tiberius clone"
[ -d ext/tiberius ] || nice -n 19 git clone -q https://github.com/Gaius-Augustus/Tiberius.git ext/tiberius || exit 1
git -C ext/tiberius rev-parse HEAD > ext/tiberius.commit
[ -d ~/envs/tiberius ] || $UV venv -q ~/envs/tiberius --python 3.12 || exit 1
echo "$(date -u +%FT%TZ) tiberius install"
(cd ext/tiberius && VIRTUAL_ENV=~/envs/tiberius nice -n 19 $UV pip install -q "tensorflow[and-cuda]>=2.17,<2.18" && \
  VIRTUAL_ENV=~/envs/tiberius nice -n 19 $UV pip install -q ".[from_source]") || echo "TIBERIUS_INSTALL_FAIL"
~/envs/tiberius/bin/python -c "import tensorflow as tf; print('tf', tf.__version__)" 2>&1 | tail -1
echo "$(date -u +%FT%TZ) segment_nt download"
nice -n 19 ~/envs/generanno/bin/python - <<'PY'
from huggingface_hub import snapshot_download
p = snapshot_download("InstaDeepAI/segment_nt", revision="1048ad869036bb8e50a15d82b323ebeb52e8ab92",
                      local_dir="ext/segment_nt",
                      allow_patterns=["config.json", "*.py", "pytorch_model.bin", "tokenizer_config.json", "vocab.txt",
                                      "special_tokens_map.json"])
print("segment_nt at", p)
PY
sha256sum ext/segment_nt/pytorch_model.bin > ext/segment_nt.sha256
echo "$(date -u +%FT%TZ) SETUP_END"
