"""v2.0 (protocol section 17, execution detail 3): read-only view of the training features in which the replaced
rows come from a small new cache and all other rows from the v1.9 cache (no 100 GB copy)."""
import numpy as np


class Overlay:
    def __init__(self, base, new, rows):
        """base [N, ...] memmap; new [len(rows), ...]; rows: indices of base replaced by new[0], new[1], ..."""
        self.base, self.new = base, new
        self.map = {int(r): i for i, r in enumerate(rows)}
        assert len(self.map) == len(new) and base.shape[1:] == new.shape[1:]
        self.shape, self.dtype = base.shape, base.dtype

    def __len__(self):
        return self.shape[0]

    def _row(self, i):
        i = int(i)
        return self.new[self.map[i]] if i in self.map else self.base[i]

    def __getitem__(self, idx):
        if isinstance(idx, (int, np.integer)):
            return np.asarray(self._row(idx))
        if isinstance(idx, slice):
            idx = range(*idx.indices(self.shape[0]))
        return np.stack([np.asarray(self._row(i)) for i in idx])
