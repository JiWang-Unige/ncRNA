"""Phase 0.4 full probes on Modal (one A10G container per backbone, run in parallel).

Reuses the GENATATOR image (torch 2.2.2, mamba-ssm 2.2.2, transformers 4.57.6) so Caduceus loads.
  MODAL_PROFILE=ji-wang-unige ~/envs/modal/bin/modal run -d probe_modal.py \
      --backbones ntv2,generanno_base,generanno_cds,caduceus_ps
"""
import os
import subprocess

import modal

from genatator_modal import image as base_image

image = (base_image.add_local_file(os.path.join(os.path.dirname(__file__), "probe_backbones.py"),
                                   "/root/probe_backbones.py")
         .add_local_python_source("genatator_modal"))  # the container re-imports this module
app = modal.App("ncrna-probe", image=image)
vol = modal.Volume.from_name("ncrna-annot", create_if_missing=True)
hf_cache = modal.Volume.from_name("ncrna-hf-cache", create_if_missing=True)


@app.function(gpu="A10G", cpu=4, memory=32768, timeout=6 * 3600,
              volumes={"/vol": vol, "/root/.cache/huggingface": hf_cache})
def probe(backbone: str, n_train: int = 300, n_eval: int = 300, epochs: int = 2):
    out = "/vol/outputs/PHASE0/0.4_probe/full"
    os.makedirs(out, exist_ok=True)
    with open(f"{out}/probe_{backbone}.log", "w") as log:
        rc = subprocess.call(["python", "-u", "/root/probe_backbones.py", "--backbone", backbone,
                              "--data-dir", "/vol/data/probe", "--out", f"{out}/probe_{backbone}.json",
                              "--n-train", str(n_train), "--n-eval", str(n_eval), "--epochs", str(epochs)],
                             stdout=log, stderr=subprocess.STDOUT)
    vol.commit()
    hf_cache.commit()
    return backbone, rc


@app.local_entrypoint()
def main(backbones: str = "ntv2,generanno_base,generanno_cds,caduceus_ps"):
    for b, rc in probe.map(backbones.split(",")):
        print(b, "rc", rc)
