"""Phase 1 labels: per-base allowed-state bitmasks for the G1 grammar (protocol v1 section 1-2).

Each base gets a uint32 bitmask over g1.STATES (bit i = state i allowed).
  annotated on strand s:   exon -> {kE s};  intron -> {kL s, kI s, kR s}
      precedence inside a strand: exon > intron, lncRNA > mRNA (same as Phase 0 0.5)
  annotated on both strands: union of both strands' sets
  unannotated:  {IG}, or in a U-ev evidence interval {IG} + lnc states of both strands
Gene sets come from a GTF (GENCODE gene_type / Ensembl gene_biotype); only lncRNA and
protein_coding transcripts are used (other biotypes are background, like the protocol).
"""
import gzip

import numpy as np

from ..grammar import g1

LNC, PC = "lncRNA", "protein_coding"


def bit(name):
    return np.uint32(1) << np.uint32(g1.IDX[name])


IG_BIT = bit("IG")
LNC_ALL = np.uint32(0)
for _s in g1.STRANDS:
    for _p in g1.PARTS:
        LNC_ALL |= bit(f"l{_p}{_s}")


def read_transcripts(gtf, chroms):
    """{chrom: [(type 'm'|'l', strand, sorted exons 0-based half-open)]}"""
    tx = {}
    op = gzip.open if gtf.endswith(".gz") else open
    with op(gtf, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 9 or f[0] not in chroms or f[2] != "exon":
                continue
            a = f[8]
            gt = (a.split('gene_type "')[1] if 'gene_type "' in a else a.split('gene_biotype "')[1]).split('"')[0]
            if gt not in (LNC, PC):
                continue
            tid = a.split('transcript_id "')[1].split('"')[0]
            key = (f[0], tid)
            tx.setdefault(key, ["l" if gt == LNC else "m", f[6], []])[2].append((int(f[3]) - 1, int(f[4])))
    out = {}
    for (c, _), (k, s, ex) in tx.items():
        out.setdefault(c, []).append((k, s, sorted(ex)))
    return out


def gene_spans(gtf, chroms):
    """all genes (any biotype) as merged [start, end) spans per chrom, for 'outside annotation' tests."""
    spans = {}
    op = gzip.open if gtf.endswith(".gz") else open
    with op(gtf, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 9 or f[0] not in chroms or f[2] != "gene":
                continue
            spans.setdefault(f[0], []).append((int(f[3]) - 1, int(f[4])))
    return spans


def strand_classes(n, txs, strand):
    """uint8 per base for one strand: 0 none, 1 mRNA intron, 2 mRNA exon, 3 lnc intron, 4 lnc exon."""
    c = np.zeros(n, np.uint8)
    mine = [(k, ex) for k, s, ex in txs if s == strand]
    for k, ex in mine:                                   # 1. mRNA introns
        if k == "m":
            for (_, e1), (s2, _) in zip(ex, ex[1:]):
                seg = c[e1:s2]
                seg[seg == 0] = 1
    for k, ex in mine:                                   # 2. mRNA exons
        if k == "m":
            for a, b in ex:
                c[a:b] = 2
    for k, ex in mine:                                   # 3. lnc introns (over background / mRNA intron)
        if k == "l":
            for (_, e1), (s2, _) in zip(ex, ex[1:]):
                seg = c[e1:s2]
                seg[(seg == 0) | (seg == 1)] = 3
    for k, ex in mine:                                   # 4. lnc exons (over everything)
        if k == "l":
            for a, b in ex:
                c[a:b] = 4
    return c


def allowed_bits(n, txs, unknown=None):
    """uint32 [n] allowed-state bitmask."""
    out = np.zeros(n, np.uint32)
    for s in g1.STRANDS:
        c = strand_classes(n, txs, s)
        for code, k, parts in ((1, "m", "LIR"), (2, "m", "E"), (3, "l", "LIR"), (4, "l", "E")):
            m = np.uint32(0)
            for p in parts:
                m |= bit(f"{k}{p}{s}")
            out[c == code] |= m
    un = out == 0
    out[un] = IG_BIT
    if unknown is not None:
        out[un & unknown] |= LNC_ALL
    return out


def bits_to_bool(b):
    """uint32 [..] -> bool [.., S]"""
    return ((b[..., None] >> np.arange(g1.S, dtype=np.uint32)) & 1).astype(bool)
