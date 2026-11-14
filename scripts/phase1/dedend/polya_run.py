#!/usr/bin/env python3
"""PLAN_DEDEND: PolyaID / PolyaStrength (Stroup & Ji 2023; github.com/zhejilab/PolyaModelsHuman, commit d70bcd0) at every
position within +-20 kb of the canonical TES of every lncRNA / protein-coding gene, on the gene's strand.

Window for position p (0-based) = 240 nt with p at index 120, as human_pipeline.get_window: + strand genome[p-120:p+120],
- strand reverse complement of genome[p-119:p+121].  One-hot A,C,G,T (isolearn OneHotEncoder order), N = all zero.
Architectures copied from human_pipeline.make_polyaid_model / make_polyastrength_model, weights loaded with load_weights.
PolyaStrength is computed only where PolyaID >= 0.75 (PLAN_DEDEND section 3).

  polya_run.py --name polya --chroms chr18,chr20 [--smoke]
  -> outputs/PHASE1/DEDEND/<name>/<chrom>_<strand>.npz : pos (int64), polyaid (float16), strength (float32, NaN if not computed)
"""
import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
PM = os.path.expanduser(os.environ.get("POLYA_DIR", "~/tools/dedend/PolyaModelsHuman"))
sys.path.insert(0, os.path.join(ROOT, "src"))
HALF, BATCH, CUT = 20_000, 4096, 0.75


def gate_busy():
    r = subprocess.run(["bash", "-c", f". {ROOT}/scripts/phase1/gpu_gate.sh; gpu_why_busy"], capture_output=True, text=True)
    why = r.stdout.strip() if r.returncode == 0 else ""
    return "" if why.startswith("only ") else why


def models():
    import tensorflow as tf
    gpus = tf.config.list_physical_devices("GPU")
    if gpus:
        tf.config.set_logical_device_configuration(gpus[0], [tf.config.LogicalDeviceConfiguration(memory_limit=7000)])
    from keras import Input
    from keras.models import Model
    from keras.layers import Dense, Dropout, Flatten, Bidirectional, Conv1D, MaxPooling1D, LSTM, ReLU
    inp = Input(shape=(240, 4))
    x = ReLU()(Conv1D(512, 8, padding="valid", strides=1, name="bin_conv")(inp))
    x = MaxPooling1D(pool_size=3, strides=3, name="bin_pool")(x)
    x = Dropout(0.7)(x)
    x = Bidirectional(LSTM(units=128, return_sequences=True, name="bin_lstm"))(x)
    x = Dropout(0.7)(x)
    x = Flatten(name="bin_flatten")(x)
    b = ReLU()(Dense(256, name="bin_dense1")(x)); b = Dropout(0.7)(b)
    b = ReLU()(Dense(256, name="bin_dense2")(b)); b = Dropout(0.7)(b)
    b = ReLU()(Dense(128, name="bin_dense3")(b)); b = Dropout(0.7)(b)
    b = ReLU()(Dense(64, name="bin_dense4")(b))
    b = Dense(1, activation="sigmoid", name="bin_predictions")(b)
    p = ReLU()(Dense(256, name="prob_dense1")(x)); p = Dropout(0.1)(p)
    p = ReLU()(Dense(256, name="prob_dense2")(p)); p = Dropout(0.1)(p)
    p = ReLU()(Dense(128, name="prob_dense3")(p)); p = Dropout(0.1)(p)
    p = ReLU()(Dense(64, name="prob_dense4")(p))
    p = Dense(50, activation="softmax", name="prob_predictions")(p)
    pid = Model(inputs=[inp], outputs=[b, p], name="model")
    pid.load_weights(f"{PM}/resources/published_models/PolyaID.h5")
    inp2 = Input(shape=(240, 4))
    y = ReLU()(Conv1D(64, 8, padding="valid", strides=1, name="seq_conv")(inp2))
    y = MaxPooling1D(pool_size=3, strides=3, name="bin_pool")(y)
    y = Dropout(0.1)(y)
    y = Bidirectional(LSTM(units=16, return_sequences=True, name="bin_lstm"))(y)
    y = Dropout(0.1)(y)
    y = Flatten(name="bin_flatten")(y)
    s = ReLU()(Dense(128, name="bin_dense1")(y)); s = Dropout(0.1)(s)
    s = ReLU()(Dense(64, name="bin_dense2")(s))
    s = Dense(1, activation="linear", name="bin_predictions")(s)
    pst = Model(inputs=inp2, outputs=s, name="model")
    pst.load_weights(f"{PM}/resources/published_models/PolyaStrength.h5")
    return pid, pst


def codes_of(path):
    with open(path) as fh:
        s = "".join(l.strip() for l in fh if not l.startswith(">")).upper()
    lut = np.full(256, 4, np.uint8)
    for i, ch in enumerate("ACGT"):
        lut[ord(ch)] = i
    return lut[np.frombuffer(s.encode(), np.uint8)]


EYE = np.vstack([np.eye(4, dtype=np.float32), np.zeros((1, 4), np.float32)])


def windows(codes, pos, strand):
    if strand == "+":
        idx = pos[:, None] + np.arange(-120, 120)[None]
        return EYE[codes[idx]]
    idx = pos[:, None] + np.arange(120, -120, -1)[None]       # p+120 ... p-119, reversed order
    c = codes[idx]
    c = np.where(c < 4, 3 - c, 4)                              # complement A<->T, C<->G
    return EYE[c]


def tes_rows(chrom):
    sys.path.insert(0, os.path.join(ROOT, "scripts/phase1/ends"))
    sys.path.insert(0, os.path.join(ROOT, "scripts/phase1"))
    sys.path.insert(0, os.path.join(ROOT, "scripts/phase1/s3"))
    from ncrna_annot.data import labels
    import prep_events as PE
    tx = labels.read_transcripts(f"{ROOT}/data/human/gencode.v50.annotation.gtf.gz", {chrom})[chrom]
    return [(s, PE.tx_events(s, ex)[f"TES{s}"][0]) for k, s, ex, canon in tx if canon]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--chroms", default="chr18,chr20")
    ap.add_argument("--smoke", action="store_true", help="first 30 TES of chr18 only, out dir <name>_smoke")
    a = ap.parse_args()
    out = f"{ROOT}/outputs/PHASE1/DEDEND/{a.name}{'_smoke' if a.smoke else ''}"
    os.makedirs(out, exist_ok=True)
    pid, pst = models()
    for c in (["chr18"] if a.smoke else a.chroms.split(",")):
        codes = codes_of(f"{ROOT}/data/human/{c}.fa")
        L = len(codes)
        rows = tes_rows(c)
        if a.smoke:
            rows = rows[:30]
        for strand in "+-":
            dst = f"{out}/{c}_{strand}.npz"
            if os.path.exists(dst):
                continue
            m = np.zeros(L, bool)
            for s, x in rows:
                if s == strand:
                    m[max(120, x - HALF):min(L - 121, x + HALF + 1)] = True
            pos = np.flatnonzero(m)
            part = f"{out}/{c}_{strand}.part.npy"
            prob = np.load(part) if os.path.exists(part) else np.full(len(pos), np.nan, np.float32)
            t0 = time.time()
            todo = np.flatnonzero(np.isnan(prob))
            for i in range(0, len(todo), BATCH * 64):
                why = gate_busy()
                if why:
                    np.save(part, prob)
                    print(json.dumps({"yield_gpu": why, "chrom": c, "strand": strand, "done": int((~np.isnan(prob)).sum())}), flush=True)
                    sys.exit(75)
                blk = todo[i:i + BATCH * 64]
                for j in range(0, len(blk), BATCH):
                    sel = blk[j:j + BATCH]
                    prob[sel] = pid.predict_on_batch(windows(codes, pos[sel], strand))[0][:, 0]
                np.save(part, prob)
                print(json.dumps({"chrom": c, "strand": strand, "done": int(i + len(blk)), "of": len(todo),
                                  "s": round(time.time() - t0)}), flush=True)
            strength = np.full(len(pos), np.nan, np.float32)
            hi = np.flatnonzero(prob >= CUT)
            for j in range(0, len(hi), BATCH):
                sel = hi[j:j + BATCH]
                strength[sel] = pst.predict_on_batch(windows(codes, pos[sel], strand))[:, 0]
            np.savez(dst, pos=pos, polyaid=prob.astype(np.float16), strength=strength)
            os.remove(part)
            print(json.dumps({"chrom": c, "strand": strand, "n_pos": len(pos), "n_polyaid_ge_0.75": len(hi),
                              "s": round(time.time() - t0)}), flush=True)


if __name__ == "__main__":
    main()
