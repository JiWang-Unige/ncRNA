"""v3.0 (protocol section 18, execution detail 7): read-only row-wise concatenation of training arrays (v2.0 windows
followed by the 3,000 new human or mouse windows) without copying the ~100 GB feature caches."""
import numpy as np


class Concat:
    def __init__(self, parts):
        """parts: arrays / memmaps / Overlay views with equal trailing shapes and dtypes."""
        self.parts = list(parts)
        assert all(p.shape[1:] == self.parts[0].shape[1:] and p.dtype == self.parts[0].dtype for p in self.parts)
        self.starts = np.cumsum([0] + [len(p) for p in self.parts])
        self.shape = (int(self.starts[-1]),) + tuple(self.parts[0].shape[1:])
        self.dtype = self.parts[0].dtype

    def __len__(self):
        return self.shape[0]

    def _row(self, i):
        i = int(i)
        if i < 0:
            i += self.shape[0]
        if not 0 <= i < self.shape[0]:
            raise IndexError(i)
        k = int(np.searchsorted(self.starts, i, "right") - 1)
        return self.parts[k][i - self.starts[k]]

    def __getitem__(self, idx):
        if isinstance(idx, (int, np.integer)):
            return np.asarray(self._row(idx))
        if isinstance(idx, slice):
            idx = range(*idx.indices(self.shape[0]))
        return np.stack([np.asarray(self._row(i)) for i in idx])
