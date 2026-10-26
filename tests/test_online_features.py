"""Online features equal the on-disk cache built by extract_p1_shard.py + assemble_p1_shards.py (stand-in backbone)."""
import numpy as np

from ncrna_annot.features import generanno as G
from ncrna_annot.features import online


def fake_fn(windows):
    # deterministic, position-dependent stand-in: token t of a window gets a value from its 6-mer and position
    out = np.empty((len(windows), G.NT, 8), np.float16)
    for i, w in enumerate(windows):
        k = np.asarray(w[:G.NT * G.TOK], np.int64).reshape(G.NT, G.TOK) @ (4 ** np.arange(G.TOK))
        out[i] = ((k[:, None] * np.arange(1, 9)) % 997 / 997 + np.arange(G.NT)[:, None] / G.NT).astype(np.float16)
    return out


def test_online_chromosome_matches_stitched_cache(monkeypatch):
    monkeypatch.setattr(G, "DIM", 8)
    rng = np.random.default_rng(0)
    L = 3 * G.STEP + G.WIN + 777
    codes = rng.integers(0, 4, L).astype(np.int8)
    starts = G.chromosome_tiles(L)
    ref = np.zeros((L // G.TOK, 8), np.float16)  # extract_p1.py semantics: sequential writes of kept tokens
    for s in starts:
        lo, hi = G.kept_tokens(s, starts)
        ref[s // G.TOK + lo:s // G.TOK + hi] = fake_fn([codes[s:s + G.WIN]])[0, lo:hi]
    ch = online.OnlineChromosome(codes, fake_fn, max_tiles=2)
    assert ch.shape == ref.shape and len(ch) == len(ref)
    assert np.array_equal(ch[0:len(ch)], ref)
    for a, b in [(0, 5000), (3000, 3500), (len(ref) - 4500, len(ref)), (7000, 7001)]:
        assert np.array_equal(ch[a:b], ref[a:b])


def test_online_windows_with_and_without_cache(monkeypatch):
    monkeypatch.setattr(G, "DIM", 8)
    rng = np.random.default_rng(1)
    codes = rng.integers(0, 4, (5, G.WIN)).astype(np.int8)
    ref = fake_fn(list(codes))
    for cache_all in (False, True):
        w = online.OnlineWindows(codes, fake_fn, cache_all=cache_all)
        assert w.shape == (5, G.NT, 8) and len(w) == 5
        assert np.array_equal(w[[3, 1]], ref[[3, 1]]) and np.array_equal(w[1:4], ref[1:4])
        assert np.array_equal(w[[1, 3, 4]], ref[[1, 3, 4]]) and np.array_equal(w[2], ref[2])
