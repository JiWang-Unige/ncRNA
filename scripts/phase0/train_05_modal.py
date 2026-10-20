"""0.5 heads (PROTOCOL v2) on Modal; logic in train_05_core.py (shared with train_05_local.py).

  MODAL_PROFILE=ji-wang-unige modal run -d train_05_modal.py [--smoke] [--arms F,Uev,Upu] [--seeds 42,123,456]
"""
import modal

from genatator_modal import image as base_image

image = base_image.add_local_python_source("genatator_modal", "train_05_core")
app = modal.App("ncrna-p05-train", image=image)
vol = modal.Volume.from_name("ncrna-annot")


@app.function(gpu="A10G", cpu=4, memory=49152, timeout=6 * 3600, volumes={"/vol": vol})
def run(arm: str, seed: int, smoke: bool = False):
    import train_05_core
    feat = "/vol/data/p05" + ("/smoke" if smoke else "")
    out = "/vol/outputs/PHASE0/0.5_unknown" + ("_smoke" if smoke else "")
    r = train_05_core.run(arm, seed, "/vol/data/p05", feat, out, smoke)
    vol.commit()
    return r


@app.local_entrypoint()
def main(smoke: bool = False, arms: str = "F,Uev,Upu", seeds: str = "42,123,456"):
    jobs = [(a, int(s)) for a in arms.split(",") for s in seeds.split(",")]
    for r in run.starmap([(a, s, smoke) for a, s in jobs]):
        print(r)
