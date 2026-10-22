"""C0 head: aligned 6-mer token features + DNA one-hot -> per-base G1 state potentials.

Adapted from coding-rna pipeline/m28/src/m28/core.py::SharedContext (commit c20e6dc), output layer
changed to the 17 G1 states; widths scaled for 4,096-d input. Not imported from the frozen repo.
"""
import torch
from torch import nn

from ..grammar import g1


class C0Head(nn.Module):
    def __init__(self, feature_dim, hidden=256, dna_hidden=32):
        super().__init__()
        self.project = nn.Linear(feature_dim, hidden)
        self.context = nn.GRU(hidden, hidden // 2, batch_first=True, bidirectional=True)
        self.dna = nn.Sequential(nn.Conv1d(5, dna_hidden, 9, padding=4), nn.GELU(),
                                 nn.Conv1d(dna_hidden, dna_hidden, 5, padding=2), nn.GELU())
        self.fuse = nn.Sequential(nn.Linear(hidden + dna_hidden, hidden), nn.GELU(), nn.LayerNorm(hidden))
        self.out = nn.Linear(hidden, g1.S)

    def hidden(self, token_features, codes):
        """fused per-base representation [B, L, hidden]."""
        B, L = codes.shape
        if token_features.shape[1] * 6 != L:
            raise ValueError("expect one feature vector per six bases")
        ctx, _ = self.context(self.project(token_features.float()))
        ctx = ctx.repeat_interleave(6, dim=1)
        dna = nn.functional.one_hot(codes.clamp(0, 4), 5).float()
        local = self.dna(dna.transpose(1, 2)).transpose(1, 2)
        return self.fuse(torch.cat([ctx, local], -1))

    def forward(self, token_features, codes):
        """token_features [B, T, D] (T = L/6), codes int64 [B, L] -> potentials [B, L, S] (fp32)."""
        return self.out(self.hidden(token_features, codes))


class C0BHead(C0Head):
    """G1B (protocol v1.6 section 13): C0 plus a boundary head -> (potentials [B, L, S], logits [B, L, 4])
    with logits in g1.BOUNDARY order (TSS+, TES+, TSS-, TES-)."""

    def __init__(self, feature_dim, hidden=256, dna_hidden=32):
        super().__init__(feature_dim, hidden, dna_hidden)
        self.bnd = nn.Linear(hidden, len(g1.BOUNDARY))

    def forward(self, token_features, codes):
        h = self.hidden(token_features, codes)
        return self.out(h), self.bnd(h)
