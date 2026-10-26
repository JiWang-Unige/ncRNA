"""Backbone features computed on the fly (win3090: no 100 GB cache on the family PC; protocol v1.8 section 15).

OnlineWindows(codes)[idx]   -> float16 [len(idx), 4096, DIM], identical to data/p1_feat/<set>.npy[idx]
OnlineChromosome(codes)[a:b] -> float16 [b - a, DIM], identical to the stitched data/p1_feat/<chrom>.npy[a:b]
(same tiles, kept-token ranges and overwrite order as extract_p1.py / assemble_p1_shards.py).  `fn(list of int8
windows) -> float16 array [n, 4096, DIM]` is the backbone call; tests pass a deterministic stand-in.
"""
from collections import OrderedDict

import numpy as np

from . import generanno as G


def backbone_fn(tok, model, batch=4):
    def fn(windows):
        out = []
        for i in range(0, len(windows), batch):
            out.append(G.features(tok, model, windows[i:i + batch]).cpu().numpy())
        return np.concatenate(out) if out else np.zeros((0, G.NT, G.DIM), np.float16)
    return fn


class OnlineWindows:
    def __init__(self, codes, fn, cache_all=False):
        self.codes, self.fn, self.cache = codes, fn, ({} if cache_all else None)
        self.shape = (len(codes), G.NT, G.DIM)
        self.dtype = np.float16

    def __len__(self):
        return self.shape[0]

    def __getitem__(self, idx):
        if isinstance(idx, slice):
            idx = np.arange(len(self))[idx]
        elif np.isscalar(idx):
            return self[np.array([int(idx)])][0]
        idx = np.asarray(idx, np.int64)
        out = np.empty((len(idx), G.NT, G.DIM), np.float16)
        todo = [j for j, k in enumerate(idx) if self.cache is None or int(k) not in self.cache]
        if todo:
            f = self.fn([np.asarray(self.codes[int(idx[j])]) for j in todo])
            for n, j in enumerate(todo):
                out[j] = f[n]
                if self.cache is not None:
                    self.cache[int(idx[j])] = f[n]
        if self.cache is not None:
            for j, k in enumerate(idx):
                if j not in todo:
                    out[j] = self.cache[int(k)]
        return out


class OnlineChromosome:
    def __init__(self, codes, fn, max_tiles=64):
        self.codes, self.fn, self.max_tiles = codes, fn, max_tiles
        L = len(codes)
        self.starts = G.chromosome_tiles(L)
        self.ranges = []  # per tile: (row_lo, row_hi, tok_lo, tok_hi)
        for s in self.starts:
            lo, hi = G.kept_tokens(s, self.starts)
            self.ranges.append((s // G.TOK + lo, s // G.TOK + hi, lo, hi))
        self.shape = (L // G.TOK, G.DIM)
        self.dtype = np.float16
        self._lru = OrderedDict()
        self.computed = 0

    def __len__(self):
        return self.shape[0]

    def _tile(self, k):
        if k in self._lru:
            self._lru.move_to_end(k)
            return self._lru[k]
        s = self.starts[k]
        _, _, lo, hi = self.ranges[k]
        t = self.fn([np.asarray(self.codes[s:s + G.WIN])])[0, lo:hi]
        self.computed += 1
        self._lru[k] = t
        if len(self._lru) > self.max_tiles:
            self._lru.popitem(last=False)
        return t

    def __getitem__(self, sl):
        if not isinstance(sl, slice) or sl.step not in (None, 1):
            raise TypeError("OnlineChromosome supports contiguous slices only")
        a, b, _ = sl.indices(len(self))
        out = np.zeros((max(0, b - a), G.DIM), np.float16)
        for k, (r0, r1, _, _) in enumerate(self.ranges):  # increasing tile order = overwrite order of the cache
            if r1 <= a or r0 >= b:
                continue
            t = self._tile(k)
            x0, x1 = max(a, r0), min(b, r1)
            out[x0 - a:x1 - a] = t[x0 - r0:x1 - r0]
        return out
