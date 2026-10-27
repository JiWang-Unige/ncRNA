"""Chromosome-wide work arrays (potentials, posteriors) in RAM, or - when NCRNA_SCRATCH is set (win3090, where the
CPU side must stay under ~12 GB next to TE's EarlGrey) - as unlinked memory-mapped files on disk.  Same dtype,
same values: only where the bytes live changes."""
import os
import uuid

import numpy as np


def zeros(shape, dtype):
    d = os.environ.get("NCRNA_SCRATCH")
    if not d:
        return np.zeros(shape, dtype)
    os.makedirs(d, exist_ok=True)
    p = os.path.join(d, f"arr_{os.getpid()}_{uuid.uuid4().hex}.npy")
    a = np.lib.format.open_memmap(p, "w+", np.dtype(dtype), tuple(shape))  # zero-filled sparse file
    os.unlink(p)  # space is released when the mapping goes away
    return a


def astype(x, dtype, block=1 << 22):
    if not os.environ.get("NCRNA_SCRATCH"):
        return x.astype(dtype)
    out = zeros(x.shape, dtype)
    for i in range(0, len(x), block):
        out[i:i + block] = x[i:i + block]
    return out
