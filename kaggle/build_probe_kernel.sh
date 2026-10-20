#!/bin/bash
# Build a Kaggle script kernel that runs probe_backbones.py for the given backbones.
# usage: build_probe_kernel.sh <slug> <n_train> <n_eval> <epochs> <backbone>...
set -euo pipefail
slug=$1 ntr=$2 nev=$3 ep=$4; shift 4
here=$(cd "$(dirname "$0")" && pwd); d=$here/$slug; mkdir -p $d
head -n -2 $here/../scripts/phase0/probe_backbones.py > $d/kernel.py
cat >> $d/kernel.py <<PY
if __name__ == "__main__":
    import subprocess, sys, traceback
    # Kaggle ships transformers 5.x; NTv2/GENERanno remote code needs the 4.x API (pin = ont env)
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "transformers==4.57.6"])
    for b in $(printf '"%s",' "$@"):
        sys.argv = ["probe", "--backbone", b, "--data-dir", "/tmp/probe_data",
                    "--out", f"/kaggle/working/probe_{b}.json", "--n-train", "$ntr",
                    "--n-eval", "$nev", "--epochs", "$ep"]
        try:
            main()
        except Exception:
            traceback.print_exc()
PY
cat > $d/kernel-metadata.json <<JS
{"id": "jiunige/$slug", "title": "$slug", "code_file": "kernel.py", "language": "python",
 "kernel_type": "script", "is_private": "true", "enable_gpu": "true", "enable_tpu": "false",
 "enable_internet": "true", "machine_shape": "", "dataset_sources": [], "competition_sources": [],
 "kernel_sources": [], "model_sources": []}
JS
echo $d
