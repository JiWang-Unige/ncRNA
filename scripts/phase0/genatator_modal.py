"""Run the official GENATATOR pipeline on Modal (used when the ont GPU is busy).

Same pins as env/install_genatator_env.sh. Inputs live on the Modal volume `ncrna-annot`:
  /vol/models/<repo-name>/ (snapshots pinned in MODEL_REVISIONS.json), /vol/data/human/<chrom>.fa
Run from ont:  MODAL_PROFILE=ji-wang-unige ~/envs/modal/bin/modal run -d genatator_modal.py \
                 --fasta data/human/chr21.fa --tag PHASE0/0.1a_genatator/chr21
"""
import json
import os
import shutil
import subprocess
import threading
import time

import modal

GH = "https://github.com"
image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("procps")
    .pip_install("torch==2.2.2", "torchvision==0.17.2", "torchaudio==2.2.2",
                 index_url="https://download.pytorch.org/whl/cu121")
    .pip_install("packaging==26.0", "ninja==1.13.0", "psutil==7.2.2", "einops==0.8.2",
                 "numpy==1.26.4", "triton==2.2.0")
    .pip_install(
        f"{GH}/Dao-AILab/causal-conv1d/releases/download/v1.4.0/causal_conv1d-1.4.0+cu122torch2.2cxx11abiFALSE-cp311-cp311-linux_x86_64.whl",
        f"{GH}/state-spaces/mamba/releases/download/v2.2.2/mamba_ssm-2.2.2+cu122torch2.2cxx11abiFALSE-cp311-cp311-linux_x86_64.whl",
        f"{GH}/Dao-AILab/flash-attention/releases/download/v2.6.3/flash_attn-2.6.3+cu123torch2.2cxx11abiFALSE-cp311-cp311-linux_x86_64.whl")
    .pip_install("biopython==1.86", "huggingface-hub==0.36.2", "safetensors==0.7.0",
                 "tokenizers==0.22.2", "transformers==4.57.6", "scipy==1.17.1", "pyfaidx~=0.7.1",
                 "pyyaml==6.0.3", "regex==2026.2.28", "tqdm==4.67.3", "rich==14.3.3", "h5py")
    .add_local_file(os.path.join(os.path.dirname(__file__), "genatator_cli.py"), "/root/genatator_cli.py")
)
app = modal.App("ncrna-genatator", image=image)
vol = modal.Volume.from_name("ncrna-annot", create_if_missing=True)
STAGES = {"edge_model_path": "genatator-moderngena-base-multispecies-edge-model",
          "region_model_path": "genatator-moderngena-base-multispecies-region-model",
          "transcript_type_model_path": "genatator-caduceus-ps-multispecies-transcript-type",
          "segmentation_model_path": "genatator-caduceus-ps-multispecies-segmentation"}


@app.function(gpu="A10G", cpu=8, memory=32768, timeout=8 * 3600, volumes={"/vol": vol})
def run(fasta: str, tag: str, ram_gib: int = 24):
    out = f"/vol/outputs/{tag}"
    work = "/tmp/genatator"  # checkpoints on local disk: the pipeline needs POSIX locks
    os.makedirs(out, exist_ok=True)
    os.makedirs(work, exist_ok=True)
    repo = "/vol/models/genatator-pipeline"
    cfg = json.load(open(f"{repo}/config.json"))
    for k, d in STAGES.items():
        cfg[k] = f"/vol/models/{d}"
    json.dump(cfg, open(f"{out}/config.json", "w"), indent=1)
    cmd = ["python", "-u", "/root/genatator_cli.py", "--fasta", f"/vol/{fasta}",
           "--output-gff", f"{work}/pred.gff3", "--config", f"{out}/config.json",
           "--devices", "0", "--checkpoint-dir", f"{work}/ckpt",
           "--cpu-ram-limit-gib-per-node", str(ram_gib), "--local-files-only"]
    stop = threading.Event()

    def heartbeat():  # make the log visible on the volume while the run is going
        while not stop.wait(300):
            shutil.copy(f"{work}/run.log", f"{out}/run.log")
            subprocess.run(["nvidia-smi", "--query-gpu=memory.used,utilization.gpu",
                            "--format=csv,noheader"], stdout=open(f"{out}/gpu.log", "a"))
            vol.commit()

    t0 = time.time()
    threading.Thread(target=heartbeat, daemon=True).start()
    with open(f"{work}/run.log", "w") as log:
        rc = subprocess.call(cmd, cwd=repo, stdout=log, stderr=subprocess.STDOUT)
    stop.set()
    for f in ("run.log", "pred.gff3", "pred.gff3.summary.json", "pred.gff3.summary.txt"):
        if os.path.exists(f"{work}/{f}"):
            shutil.copy(f"{work}/{f}", f"{out}/{f}")
    json.dump({"rc": rc, "elapsed_s": round(time.time() - t0, 1), "cmd": cmd, "gpu": "A10G"},
              open(f"{out}/modal_run.json", "w"), indent=1)
    vol.commit()
    return rc


@app.local_entrypoint()
def main(fasta: str, tag: str, ram_gib: int = 24):
    print("rc", run.remote(fasta, tag, ram_gib))
