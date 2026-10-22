"""Whole-chromosome decoding and gene confidence (protocol v1.4 section 11).

decode_blocks: duration-constrained Viterbi (S3) in 1 Mb blocks with 100 kb padding on both sides;
  a gene is kept from the block whose core [s, e) contains its first base (bounded memory).
chromosome_marginals: forward-backward posteriors on 24,576-bp tiles (step 20,472, central part kept).
gene_confidence: mean posterior of the gene's decoded states over its whole span
  (exon bases: P(kE s); intron bases: P(kL s) + P(kI s) + P(kR s)).
"""
import copy

import numpy as np
import torch

from ..features import generanno as G
from ..grammar import g1


def decode_blocks(pot, codes, T, min_exon=20, min_intron=50, block=1_000_000, pad=100_000, mode="s3", bnd=None, coef=None):
    """pot float32 [L, S] (already masked for lnc-only arms), codes int8 [L] -> genes (1-based exons).
    G1B: bnd float [L, 4] boundary logits + coef (g1.BoundaryCoef) add gene entry / exit edge bonuses (S3 only)."""
    L = len(pot)
    out, coef_cpu = [], None
    for s in range(0, L, block):
        e = min(s + block, L)
        a, b = max(s - pad, 0), min(e + pad, L)
        gates = g1.gate_mask(np.asarray(codes[a:b]))
        if mode == "s3":
            edges = None
            if bnd is not None and coef_cpu is None:
                coef_cpu = copy.deepcopy(coef).double().cpu()
            if bnd is not None:
                with torch.no_grad():
                    edges = g1.boundary_edges(torch.as_tensor(np.asarray(bnd[a:b]), dtype=torch.float64), coef_cpu)
            path, _ = g1.viterbi_duration(pot[a:b], T, gates, None, min_exon, min_intron, edges=edges)
        else:  # "argmax" (S1): per-base argmax of posteriors passed in as `pot`
            path = np.asarray(pot[a:b]).argmax(-1).astype(np.int8)
        for gene in g1.path_to_genes(path, offset=a):
            first = gene["exons"][0][0] - 1  # 0-based
            if s <= first < e:
                out.append(gene)
    return out


def chromosome_marginals(emit_fn, L, T, codes, device="cuda", coef=None):
    """emit_fn(start, end) -> torch [end-start, S] potentials on `device` (G1B: a tuple (potentials,
    boundary logits [end-start, 4]) and coef given).  Returns float16 [L, S]."""
    post = np.zeros((L, g1.S), np.float16)
    starts = G.chromosome_tiles(L)
    for k, s in enumerate(starts):
        e = emit_fn(s, s + G.WIN)
        edges = None
        if coef is not None:
            e, b = e
            with torch.no_grad():
                edges = g1.boundary_edges(b[None].to(e.dtype), coef)
        e = e[None]
        gates = torch.as_tensor(g1.gate_mask(np.asarray(codes[s:s + G.WIN])), device=device)[None]
        m = g1.marginals(e, T, gates, edges=edges)[0].float().cpu().numpy()
        lo = 0 if k == 0 else G.MARGIN * G.TOK
        hi = G.WIN if k == len(starts) - 1 else G.WIN - G.MARGIN * G.TOK
        post[s + lo:s + hi] = m[lo:hi]
    return post


def gene_confidence(gene, post):
    k, s = gene["type"], gene["strand"]
    E = g1.IDX[f"{k}E{s}"]
    intron = [g1.IDX[f"{k}{p}{s}"] for p in ("L", "I", "R")]
    ex = gene["exons"]
    vals = []
    for a, b in ex:
        vals.append(post[a - 1:b, E].astype(np.float64))
    for (_, b1), (a2, _) in zip(ex, ex[1:]):
        vals.append(post[b1:a2 - 1][:, intron].astype(np.float64).sum(1))
    v = np.concatenate(vals)
    return float(v.mean()) if len(v) else 0.0
