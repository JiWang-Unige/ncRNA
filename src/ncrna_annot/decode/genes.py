"""Whole-chromosome decoding: head potentials over tiled token features -> Viterbi -> genes -> GFF3."""
import numpy as np
import torch

from ..grammar import g1

CHUNK, MARGIN = 4096, 256  # tokens per head call, context margin on both sides


def chromosome_potentials(head, feats, codes, device="cuda"):
    """feats: array-like [ntok, D] (memmap ok), codes int8 [L] -> float32 [ntok*6, S] potentials."""
    ntok = len(feats)
    out = np.zeros((ntok * 6, g1.S), np.float32)
    head.eval()
    with torch.no_grad():
        for s in range(0, ntok, CHUNK - 2 * MARGIN):
            lo, hi = max(s - MARGIN, 0), min(s + CHUNK - MARGIN, ntok)
            f = torch.as_tensor(np.asarray(feats[lo:hi], np.float32), device=device)[None]
            c = torch.as_tensor(np.asarray(codes[lo * 6:hi * 6], np.int64), device=device)[None]
            p = head(f, c)[0].cpu().numpy()
            a, b = s, min(s + CHUNK - 2 * MARGIN, ntok)
            out[a * 6:b * 6] = p[(a - lo) * 6:(b - lo) * 6]
    return out


def decode(potentials, codes, T):
    """-> list of genes (dict type/strand/exons 1-based closed) for the covered prefix of the chromosome."""
    n = len(potentials)
    gates = g1.gate_mask(np.asarray(codes[:n]))
    path, _ = g1.viterbi(potentials, T, gates)
    return g1.path_to_genes(path)


def write_gff3(genes, chrom, path, source="ncrna_annot", min_len=0):
    kind = {"l": "lnc_RNA", "m": "mRNA"}
    with open(path, "w") as fh:
        fh.write("##gff-version 3\n")
        for i, g in enumerate(genes):
            ex = g["exons"]
            if sum(b - a + 1 for a, b in ex) < min_len:
                continue
            gid = f"{source}_{chrom}_{i}"
            a, b = ex[0][0], ex[-1][1]
            fh.write(f"{chrom}\t{source}\tgene\t{a}\t{b}\t.\t{g['strand']}\t.\tID={gid}\n")
            fh.write(f"{chrom}\t{source}\t{kind[g['type']]}\t{a}\t{b}\t.\t{g['strand']}\t.\tID={gid}.t1;Parent={gid}\n")
            for j, (x, y) in enumerate(ex):
                fh.write(f"{chrom}\t{source}\texon\t{x}\t{y}\t.\t{g['strand']}\t.\tID={gid}.t1.e{j};Parent={gid}.t1\n")
