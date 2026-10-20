"""Frozen GENERanno-1.2b-cds features (Phase 1 v1 section 2): hidden_states 17 and 23 concatenated.

One 4,096-d vector per 6-mer token (no BOS). Windows of 24,576 bp = 4,096 tokens; non-ACGT -> 'A'
so the 6-mer grid stays aligned (as in Phase 0). Whole chromosomes are tiled with step 20,472 bp
and the central tokens kept (margin 342 tokens), as in Phase 0 0.5.
"""
import numpy as np
import torch

REPO = "GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview"
LAYERS = (17, 23)
WIN, TOK, NT, STEP, MARGIN = 24576, 6, 4096, 20472, 342
DIM = 2048 * len(LAYERS)
ACGT = np.array(list("ACGTA"))  # code 4 (N) -> A


def half_dtype():
    """bf16 on Ampere+ (e.g. RTX 4090, A10G), fp16 on older GPUs (e.g. Kaggle T4, no native bf16)."""
    return torch.bfloat16 if torch.cuda.get_device_capability()[0] >= 8 else torch.float16


def load(device="cuda"):
    from transformers import AutoModel, AutoTokenizer
    tok = AutoTokenizer.from_pretrained(REPO, trust_remote_code=True)
    model = AutoModel.from_pretrained(REPO, trust_remote_code=True, torch_dtype=half_dtype()).to(device).eval()
    return tok, model


def codes_to_str(codes):
    return "".join(ACGT[np.asarray(codes, np.int64)])


def features(tok, model, codes_batch, dtype=None):
    """codes_batch: list of int8 arrays of length 24,576 -> float16 tensor [B, 4096, DIM] (on device)."""
    dev = next(model.parameters()).device
    if dtype is None:
        dtype = half_dtype()
    ids = torch.tensor([tok(codes_to_str(c))["input_ids"] for c in codes_batch], device=dev)
    assert ids.shape[1] == NT + 1, ids.shape
    with torch.no_grad(), torch.autocast("cuda", dtype=dtype):
        hs = model(input_ids=ids, output_hidden_states=True).hidden_states
    return torch.cat([hs[k][:, 1:] for k in LAYERS], -1).half()


def chromosome_tiles(length):
    """window starts (bp) covering a chromosome; last start aligned to the 6-mer grid."""
    starts = list(range(0, length - WIN + 1, STEP))
    last = ((length - WIN) // TOK) * TOK
    if not starts or starts[-1] != last:
        starts.append(last)
    return starts


def kept_tokens(start, starts):
    """[lo, hi) token range of window `start` that is written to the chromosome array."""
    lo = 0 if start == starts[0] else MARGIN
    hi = NT if start == starts[-1] else NT - MARGIN
    return lo, hi
