# Colab VM smoke: the runtime ships Python 3.13 / torch 2.11, so build the pinned
# GENATATOR stack (py3.11, torch 2.2.2+cu121, mamba-ssm 2.2.2) in a separate uv venv.
import subprocess, sys, time
t0 = time.time()
def sh(c):
    r = subprocess.run(c, shell=True, capture_output=True, text=True)
    print(f"$ {c[:160]}\n{(r.stdout + r.stderr)[-1200:]}", flush=True)
V = "/content/py311"
GH = "https://github.com"
sh("nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv")
sh(f"{sys.executable} -m pip install -q uv && uv venv -q --python 3.11 {V}")
UV = f"uv pip install -q --python {V}/bin/python"
sh(f"{UV} torch==2.2.2 --index-url https://download.pytorch.org/whl/cu121")
sh(f"{UV} numpy==1.26.4 einops==0.8.2 triton==2.2.0 packaging ninja setuptools transformers==4.57.6 "
   f"{GH}/Dao-AILab/causal-conv1d/releases/download/v1.4.0/causal_conv1d-1.4.0+cu122torch2.2cxx11abiFALSE-cp311-cp311-linux_x86_64.whl "
   f"{GH}/state-spaces/mamba/releases/download/v2.2.2/mamba_ssm-2.2.2+cu122torch2.2cxx11abiFALSE-cp311-cp311-linux_x86_64.whl")
print("install_s", round(time.time() - t0, 1), flush=True)
sh(f"""{V}/bin/python -c "
import torch, mamba_ssm
from mamba_ssm import Mamba
m = Mamba(d_model=64).cuda()
y = m(torch.randn(2, 4096, 64, device='cuda'))
print('torch', torch.__version__, 'cuda', torch.version.cuda, 'mamba', mamba_ssm.__version__, 'out', tuple(y.shape), 'finite', bool(torch.isfinite(y).all()))
" """)
print("total_s", round(time.time() - t0, 1), flush=True)
