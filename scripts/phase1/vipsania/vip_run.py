#!/usr/bin/env python3
"""Vipsania on ont (card vipsania-ref option 1, 2026-11-03; design in docs/phase1_results/REVIEW_VIPSANIA.md section 5).

  vip_run.py annotate --chrom chr18|chr20  -> outputs/PHASE1/VIPSANIA/<chrom>/vipsania.gff3 + annotate.json
  vip_run.py post --chrom chr18            -> outputs/PHASE1/VIPSANIA/<chrom>/post6.npy (float16 [L, 6]) + post.json
      columns: cod+, cod-, acc+, don+, acc-, don-   (positions 0-based, both strands in forward coordinates)
      cod = E0 + E1 + E2 + START + EI0-2 + IE0-2 + STOP (every state but IR and the six intron states)
      don = EI0 + EI1 + EI2 (exon -> intron), acc = IE0 + IE1 + IE2 (intron -> exon)
      state order from bricks2marble.tf.hmm.tools.state_names(2):
      IR, I00 I10 I20 I01 I11 I21, E0 E1 E2, START, EI0 EI1 EI2, IE0 IE1 IE2, STOP

Model etb1go6q (Vertebrata, Vipsania 1.1.0, pretrained, no finetuning) from ext/vipsania_models (sha256 in
outputs/PHASE1/VIPSANIA/model_sha256.txt); input data/human/<chrom>.fa extracted from data/human/hg38.fa.gz as is
(soft-masked; Vipsania reads lowercase as repeats).  annotate = the package's own annotate_model (Viterbi, T = 200,000).
post = the HMM posterior (HMMMode.POSTERIOR) on windows of T = 200,000 placed every 150,000 bp with 25,000 bp context on
each side (N padding at the chromosome ends); only the central 150,000 bp of each window are kept.

GPU memory is capped at VIP_GPU_MB (default 7000 MiB, about 30% of the 4090) by a logical-device limit set before
anything touches the GPU; batch B = 2.  The caller (vip_queue.sh) waits on gpu_gate.sh before every GPU step.
Run with ~/envs/vipsania_ncrna/bin/python.
"""
import argparse
import json
import math
import os
import sys
import time

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
MODEL, MODEL_DIR = "etb1go6q", f"{ROOT}/ext/vipsania_models"
T, STEP, B = 200_000, 150_000, 2
MARGIN = (T - STEP) // 2
GPU_MB = int(os.environ.get("VIP_GPU_MB", "7000"))

import numpy as np  # noqa: E402
import tensorflow as tf  # noqa: E402

_gpus = tf.config.list_physical_devices("GPU")
if not _gpus:
    sys.exit("no GPU visible to TensorFlow")
tf.config.set_logical_device_configuration(_gpus[0], [tf.config.LogicalDeviceConfiguration(memory_limit=GPU_MB)])


def peak_gib():
    try:
        return round(tf.config.experimental.get_memory_info("GPU:0")["peak"] / 1024 ** 3, 2)
    except Exception:  # noqa: BLE001  (not every allocator keeps statistics)
        return None


def read_fasta(path, chrom):
    seq, keep = [], False
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                keep = line[1:].split()[0] == chrom
                continue
            if keep:
                seq.append(line.strip())
    return "".join(seq)


def encode(seq):
    """bricks2marble nucleotide codes: A C G T = 0-3, N (and anything else) = 4, a c g t = 5-8 (repeat-masked)."""
    lut = np.full(256, 4, np.int8)
    for i, c in enumerate("ACGT"):
        lut[ord(c)] = i
    for i, c in enumerate("acgt"):
        lut[ord(c)] = 5 + i
    return lut[np.frombuffer(seq.encode(), np.uint8)]


def annotate(chrom, out):
    from vipsania.cli.annotate import annotate_model
    t0 = time.time()
    annotate_model(MODEL, f"{ROOT}/data/human/{chrom}.fa", output=f"{out}/vipsania.gff3", model_dir=MODEL_DIR,
                   T=T, B=B, include=[chrom])
    json.dump({"chrom": chrom, "model": MODEL, "T": T, "B": B, "gpu_mb_cap": GPU_MB, "peak_gib": peak_gib(),
               "seconds": round(time.time() - t0)}, open(f"{out}/annotate.json", "w"), indent=1)


def post(chrom, out):
    import vipsania
    from bricks2marble.struct.fasta import one_hot
    from bricks2marble.tf.hmm.tools import state_names
    from hidten import HMMMode
    from vipsania.data.mask import pseudo_mask
    names = state_names(2)
    cod = [i for i, n in enumerate(names) if n != "IR" and not n.startswith("I") or n.startswith("IE")]
    don = [names.index(f"EI{k}") for k in range(3)]
    acc = [names.index(f"IE{k}") for k in range(3)]
    assert len(names) == 18 and len(cod) == 11, (names, cod)
    t0 = time.time()
    seq = encode(read_fasta(f"{ROOT}/data/human/{chrom}.fa", chrom))
    L = len(seq)
    V = vipsania.create_model(MODEL, build=True, compile=False, load=True, id_parent_folder=MODEL_DIR,
                              override_config={"model": {"hmm": {"repeats_emitter": None, "repeats_at_borders": None}}})
    par = math.isqrt(T)
    while T % par:
        par -= 1
    V.set_options(parallel=par, tree_depth=(T - 1).bit_length() - 1)
    V.toggle_inference(HMMMode.POSTERIOR, hmm=-1, enable=True)
    call = tf.function(lambda x: V(x), jit_compile=True)
    res = np.lib.format.open_memmap(f"{out}/post6.npy.part", mode="w+", dtype=np.float16, shape=(L, 6))
    starts = list(range(0, L, STEP))
    pad = np.full(T, 4, np.int8)
    for i in range(0, len(starts), B):
        group = starts[i:i + B]
        batch = []
        for s in group:
            w = pad.copy()
            a, b = max(s - MARGIN, 0), min(s - MARGIN + T, L)
            w[a - (s - MARGIN):b - (s - MARGIN)] = seq[a:b]
            batch.append(w)
        while len(batch) < B:
            batch.append(pad)
        x = pseudo_mask(tf.convert_to_tensor(one_hot(np.stack(batch), repeats="track", N="track"), dtype=tf.float32))
        p = call(x).numpy()  # (B, T, 2, 18): forward strand, reverse strand (already in forward coordinates)
        assert p.shape == (B, T, 2, 18), p.shape
        for j, s in enumerate(group):
            e = min(s + STEP, L)
            q = p[j, MARGIN:MARGIN + (e - s)]
            res[s:e] = np.stack([q[:, 0, cod].sum(-1), q[:, 1, cod].sum(-1), q[:, 0, acc].sum(-1), q[:, 0, don].sum(-1),
                                 q[:, 1, acc].sum(-1), q[:, 1, don].sum(-1)], -1).astype(np.float16)
        if i % (20 * B) == 0:
            print(json.dumps({"window": i, "of": len(starts), "s": round(time.time() - t0), "peak_gib": peak_gib()}),
                  flush=True)
    res.flush()
    del res
    os.replace(f"{out}/post6.npy.part", f"{out}/post6.npy")
    json.dump({"chrom": chrom, "L": L, "model": MODEL, "T": T, "step": STEP, "B": B, "parallel": par,
               "cols": ["cod+", "cod-", "acc+", "don+", "acc-", "don-"], "cod_states": [names[k] for k in cod],
               "gpu_mb_cap": GPU_MB, "peak_gib": peak_gib(), "seconds": round(time.time() - t0)},
              open(f"{out}/post.json", "w"), indent=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=("annotate", "post"))
    ap.add_argument("--chrom", required=True, choices=("chr18", "chr20"))
    a = ap.parse_args()
    out = f"{ROOT}/outputs/PHASE1/VIPSANIA/{a.chrom}"
    os.makedirs(out, exist_ok=True)
    (annotate if a.what == "annotate" else post)(a.chrom, out)
    print("DONE", a.what, a.chrom, flush=True)


if __name__ == "__main__":
    main()
