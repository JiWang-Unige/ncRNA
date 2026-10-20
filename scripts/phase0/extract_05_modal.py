"""0.5 feature cache (PROTOCOL v2) on Modal A10G.

Frozen GENERanno-1.2b-cds, hidden_states[23] (6-mer tokens), PCA -> 512 (fit on 200 train windows),
fp16 to /vol/data/p05/feat_{train,val}.npy [n, 4096, 512] and feat_chr20.npy [len//6, 512].
chr20 tiling: 24,576-bp windows, step 20,472 (multiple of 6), keep the central 342 + 3412 tokens.
  MODAL_PROFILE=ji-wang-unige modal run -d extract_05_modal.py [--smoke]
"""
import json
import os
import time

import modal

from genatator_modal import image as base_image

image = base_image.add_local_python_source("genatator_modal")
app = modal.App("ncrna-p05-extract", image=image)
vol = modal.Volume.from_name("ncrna-annot")
hf_cache = modal.Volume.from_name("ncrna-hf-cache", create_if_missing=True)
REPO, LAYER, K, WIN, TOK = "GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview", 23, 512, 24576, 6
NT, STEP, MARGIN = 4096, 20472, 342  # tokens per window, bp step, kept-margin tokens


@app.function(gpu="A10G", cpu=8, memory=65536, timeout=8 * 3600,
              volumes={"/vol": vol, "/root/.cache/huggingface": hf_cache})
def extract(smoke: bool = False):
    import numpy as np
    import torch
    from transformers import AutoModel, AutoTokenizer
    D = "/vol/data/p05" + ("/smoke" if smoke else "")
    src = "/vol/data/p05"
    os.makedirs(D, exist_ok=True)
    t0 = time.time()
    tok = AutoTokenizer.from_pretrained(REPO, trust_remote_code=True)
    model = AutoModel.from_pretrained(REPO, trust_remote_code=True).cuda().eval()

    def feats(seqs):
        ids = torch.tensor([tok(s)["input_ids"] for s in seqs], device="cuda")
        assert ids.shape[1] == NT + 1, ids.shape  # BOS + 4096 six-mers
        with torch.no_grad(), torch.autocast("cuda", dtype=torch.bfloat16):
            h = model(input_ids=ids, output_hidden_states=True).hidden_states[LAYER]
        return h[:, 1:].float()

    train = [l.strip() for l in open(f"{src}/seqs_train.txt")]
    val = [l.strip() for l in open(f"{src}/seqs_val.txt")]
    if smoke:
        train, val = train[:20], val[:8]
    # PCA on 200 training windows (subsample 256 tokens each)
    rng = np.random.default_rng(0)
    X = torch.cat([feats([s])[0][torch.as_tensor(rng.choice(NT, 256, replace=False))] for s in train[:200]])
    mu = X.mean(0)
    U, S, V = torch.linalg.svd(X - mu, full_matrices=False)
    comp = V[:K]  # [K, 2048]
    evr = float((S[:K] ** 2).sum() / (S ** 2).sum())
    np.savez(f"{D}/pca.npz", mean=mu.cpu().numpy(), comp=comp.cpu().numpy(), evr=evr)
    proj = lambda h: ((h - mu) @ comp.T).half().cpu().numpy()
    for split, seqs in (("train", train), ("val", val)):
        out = np.lib.format.open_memmap(f"{D}/feat_{split}.npy", "w+", np.float16, (len(seqs), NT, K))
        for i in range(0, len(seqs), 4):
            out[i:i + 4] = proj(feats(seqs[i:i + 4]))
        out.flush()
        print(split, "done", round(time.time() - t0), flush=True)
    seq = "".join(l.strip() for l in open("/vol/data/human/chr20.fa") if not l.startswith(">")).upper()
    seq = "".join(c if c in "ACGT" else "A" for c in seq)
    L = len(seq) if not smoke else 1_000_000
    ntok = L // TOK
    out = np.lib.format.open_memmap(f"{D}/feat_chr20.npy", "w+", np.float16, (ntok, K))
    starts = list(range(0, L - WIN + 1, STEP))
    last = ((L - WIN) // TOK) * TOK
    if starts[-1] != last:
        starts.append(last)
    for j in range(0, len(starts), 4):
        batch = starts[j:j + 4]
        P = proj(feats([seq[s:s + WIN] for s in batch]))
        for s, p in zip(batch, P):
            t0w = s // TOK
            lo = 0 if s == 0 else MARGIN
            hi = NT if s == starts[-1] else NT - MARGIN
            out[t0w + lo:t0w + hi] = p[lo:hi]
    out.flush()
    json.dump({"layer": LAYER, "k": K, "pca_evr": evr, "n_train": len(train), "n_val": len(val),
               "chr20_tokens": ntok, "chr20_windows": len(starts), "elapsed_s": round(time.time() - t0, 1),
               "smoke": smoke}, open(f"{D}/extract_meta.json", "w"), indent=1)
    vol.commit()
    hf_cache.commit()
    return evr


@app.local_entrypoint()
def main(smoke: bool = False):
    print("pca evr", extract.remote(smoke))
