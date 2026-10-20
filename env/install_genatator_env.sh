#!/bin/bash
# Build the GENATATOR-pipeline inference env on ont (uv venv, no GPU needed to install).
# Pins follow AIRI-Institute/genatator-pipeline environment.yml (torch 2.2.2+cu121,
# mamba-ssm 2.2.2, causal-conv1d 1.4.0, flash-attn 2.6.3, transformers 4.57.6).
# CUDA extensions come from the upstream GitHub release wheels (no local compile).
set -euo pipefail
ENV=${ENV:-$HOME/envs/genatator}
UV=$HOME/.local/bin/uv
LOG=${LOG:-$HOME/ncrna-annot/logs/install_genatator_env.log}
mkdir -p "$(dirname "$LOG")"
exec > >(tee -a "$LOG") 2>&1
echo "== $(date -Is) install into $ENV"

[ -x "$ENV/bin/python" ] || $UV venv --python 3.11 "$ENV"
PIP="$UV pip install --python $ENV/bin/python"

$PIP torch==2.2.2+cu121 torchvision==0.17.2+cu121 torchaudio==2.2.2+cu121 \
  --index-url https://download.pytorch.org/whl/cu121
$PIP packaging==26.0 ninja==1.13.0 psutil==7.2.2 einops==0.8.2 numpy==1.26.4 triton==2.2.0

GH=https://github.com
$PIP "$GH/Dao-AILab/causal-conv1d/releases/download/v1.4.0/causal_conv1d-1.4.0+cu122torch2.2cxx11abiFALSE-cp311-cp311-linux_x86_64.whl"
$PIP "$GH/state-spaces/mamba/releases/download/v2.2.2/mamba_ssm-2.2.2+cu122torch2.2cxx11abiFALSE-cp311-cp311-linux_x86_64.whl"
$PIP "$GH/Dao-AILab/flash-attention/releases/download/v2.6.3/flash_attn-2.6.3+cu123torch2.2cxx11abiFALSE-cp311-cp311-linux_x86_64.whl"

# docker/requirements.txt of the pipeline repo, minus flask (web service only)
$PIP biopython==1.86 huggingface-hub==0.36.2 hf-xet==1.4.2 safetensors==0.7.0 \
  tokenizers==0.22.2 transformers==4.57.6 scipy==1.17.1 pyfaidx~=0.7.1 pyyaml==6.0.3 \
  regex==2026.2.28 tqdm==4.67.3 rich==14.3.3 h5py pandas pytest

"$ENV/bin/python" - <<'EOF'
import torch, transformers, mamba_ssm, causal_conv1d, flash_attn
print("torch", torch.__version__, "cuda", torch.version.cuda, "available", torch.cuda.is_available())
print("transformers", transformers.__version__, "mamba_ssm", mamba_ssm.__version__,
      "causal_conv1d", causal_conv1d.__version__, "flash_attn", flash_attn.__version__)
EOF
echo "== $(date -Is) env OK"
