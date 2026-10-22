"""Shard bookkeeping (extract_p1_shard.py + assemble_p1_shards.py) with stub features: every kept token of a
chromosome lands at its own row, train windows keep their global index, sha256 and cover files are written."""
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SCR = os.path.join(HERE, "..", "scripts", "phase1")
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from ncrna_annot.features import generanno as G  # noqa: E402


def run(*args):
    subprocess.check_call([sys.executable] + list(args))


def test_chrom_and_train_shards(tmp_path):
    L = 2 * G.STEP + G.WIN + 1234  # tiles at 0, STEP, 2*STEP plus the grid-aligned last one
    rng = np.random.default_rng(1)
    chrom = rng.integers(0, 4, L, dtype=np.int8)
    train = rng.integers(0, 4, (5, G.WIN), dtype=np.int8)
    d = tmp_path / "p1"; d.mkdir()
    np.save(d / "chrX_codes.npy", chrom); np.save(d / "train_codes.npy", train)
    starts = G.chromosome_tiles(L)
    assert len(starts) == 4
    sh = tmp_path / "sh"
    # chrom shard from a bp slice (tiles 1:4), as colab_p1_shards.sh stages it
    lo = starts[1]
    np.save(tmp_path / "slice.npy", chrom[lo:starts[3] + G.WIN])
    run(f"{SCR}/extract_p1_shard.py", "--cpu-stub", "--kind", "chrom", "--name", "chrX", "--codes", str(tmp_path / "slice.npy"),
        "--tiles", "1:4", "--chrom-len", str(L), "--bp-offset", str(lo), "--out", str(sh))
    run(f"{SCR}/extract_p1_shard.py", "--cpu-stub", "--kind", "train", "--codes", str(d / "train_codes.npy"),
        "--range", "0:2", "--index-offset", "3", "--out", str(sh))
    m = json.load(open(sh / "chrX_t00001_00004.json"))
    rows = sum(hi - lo for _, _, lo, hi, _, _ in m["tiles"])
    assert m["shape"] == [rows, G.DIM] and m["tiles"][-1][3] == G.NT  # last tile keeps its right margin
    feat = tmp_path / "feat"
    run(f"{SCR}/assemble_p1_shards.py", "--shards", str(sh), "--feat-root", str(feat), "--codes-dir", str(d))
    arr = np.load(feat / "chrX.npy", mmap_mode="r")
    assert arr.shape == (L // G.TOK, G.DIM)
    shard = np.load(sh / "chrX_t00001_00004.npy")
    for k, s, lo_, hi_, r0, r1 in m["tiles"]:
        assert np.array_equal(arr[s // G.TOK + lo_:s // G.TOK + hi_], shard[r0:r1])
    assert np.abs(np.asarray(arr[:starts[1] // G.TOK + G.MARGIN], np.float32)).sum() == 0  # tile 0 untouched
    tr = np.load(feat / "train.npy", mmap_mode="r")
    assert tr.shape == (5, G.NT, G.DIM) and np.array_equal(tr[3:5], np.load(sh / "train_00003_00005.npy"))
    assert json.load(open(feat / "train.cover.json")) == [[3, 5]]
    # a corrupted shard is refused
    with open(sh / "train_00003_00005.npy", "r+b") as fh:
        fh.seek(200); fh.write(b"\x00\x01")
    r = subprocess.run([sys.executable, f"{SCR}/assemble_p1_shards.py", "--shards", str(sh), "--feat-root", str(feat),
                        "--codes-dir", str(d), "--check-only"], capture_output=True, text=True)
    assert r.returncode == 1 and "SHA256 MISMATCH" in r.stdout
