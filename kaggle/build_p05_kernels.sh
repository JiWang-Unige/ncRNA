#!/bin/bash
# Build one Kaggle script kernel per 0.5 arm (3 seeds each) from scripts/phase0/train_05_core.py.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
for arm in F Uev Upu; do
  slug=ncrna-p05-$(echo $arm | tr A-Z a-z); d=$here/$slug; mkdir -p $d
  cat $here/../scripts/phase0/train_05_core.py > $d/kernel.py
  cat >> $d/kernel.py <<PY


if __name__ == "__main__":
    import glob
    IN = os.path.dirname(glob.glob("/kaggle/input/**/feat_train.npy", recursive=True)[0])
    for s in (42, 123, 456):
        t = time.time()
        r = run("$arm", s, IN, IN, "/kaggle/working/out")
        print(json.dumps({"arm": "$arm", "seed": s, "best": r[2], "s": round(time.time() - t)}), flush=True)
PY
  cat > $d/kernel-metadata.json <<JS
{"id": "jiunige/$slug", "title": "$slug", "code_file": "kernel.py", "language": "python",
 "kernel_type": "script", "is_private": "true", "enable_gpu": "true", "enable_tpu": "false",
 "enable_internet": "false", "machine_shape": "", "dataset_sources": ["jiunige/ncrna-p05-features"],
 "competition_sources": [], "kernel_sources": [], "model_sources": []}
JS
  echo $d
done
