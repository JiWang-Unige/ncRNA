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
    """{chrom: [(type 'm'|'l', strand, sorted exons 0-based half-open, is_canonical)]}
    is_canonical = the transcript carries tag "Ensembl_canonical" (one per gene in GENCODE/Ensembl)."""
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
            rec = tx.setdefault(key, ["l" if gt == LNC else "m", f[6], [], False])
            rec[2].append((int(f[3]) - 1, int(f[4])))
            rec[3] = rec[3] or 'tag "Ensembl_canonical"' in a
    out = {}
    for (c, _), (k, s, ex, canon) in tx.items():
        out.setdefault(c, []).append((k, s, sorted(ex), canon))
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
    mine = [(t[0], t[2]) for t in txs if t[1] == strand]
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


ALL_BITS = np.uint32((1 << g1.S) - 1)


def _type_strand_bits(k, s):
    m = IG_BIT
    for p in g1.PARTS:
        m |= bit(f"{k}{p}{s}")
    return m


def allowed_bits(n, txs, unknown=None, canonical_only=True):
    """uint32 [n] allowed-state bitmask (protocol v1.2).

    canonical_only: deterministic labels come from each gene's Ensembl_canonical transcript only;
    exons of the other isoforms that the canonical labels do not call exon (same type and strand)
    become unknown for that gene: {IG} + all states of its type/strand.
    """
    main = [t for t in txs if len(t) < 4 or t[3]] if canonical_only else txs
    out = np.zeros(n, np.uint32)
    cls = {}
    for s in g1.STRANDS:
        c = strand_classes(n, main, s)
        cls[s] = c
        for code, k, parts in ((1, "m", "LIR"), (2, "m", "E"), (3, "l", "LIR"), (4, "l", "E")):
            m = np.uint32(0)
            for p in parts:
                m |= bit(f"{k}{p}{s}")
            out[c == code] |= m
    un = out == 0
    out[un] = IG_BIT
    if unknown is not None:
        out[un & unknown] |= LNC_ALL
    if canonical_only:
        for t in txs:
            if len(t) < 4 or t[3]:
                continue
            k, s, ex = t[0], t[1], t[2]
            exon_code = 4 if k == "l" else 2
            for a, b in ex:
                seg = cls[s][a:b]
                sel = np.flatnonzero(seg != exon_code) + a
                out[sel] |= _type_strand_bits(k, s)
    return out


def repair(bits, codes, flank=50, max_iter=20):
    """Make the observation feasible under G1: wherever the constrained Viterbi (zero potentials) needs a
    forbidden step, allow all states within +-flank bp; iterate.  Returns (bits, relaxed_bp, n_violations)."""
    T = g1.G1Transitions()().detach().double().numpy()  # structure only (all free parameters 0)
    bits = np.array(bits, np.uint32, copy=True)
    gates = g1.gate_mask(np.asarray(codes))
    first_violations = None
    for _ in range(max_iter):
        al = bits_to_bool(bits)
        x = np.where(al, gates, g1.NEG)
        path, sc = g1.viterbi(np.zeros((len(bits), g1.S)), T, gates, al)
        if sc > g1.NEG / 2:
            break
        p = path.astype(np.int64)
        bad_e = x[np.arange(len(p)), p] < g1.NEG / 2
        bad_t = np.concatenate([[False], T[p[:-1], p[1:]] < g1.NEG / 2])
        bad = np.flatnonzero(bad_e | bad_t)
        if first_violations is None:
            first_violations = len(bad)
        for t in bad:
            bits[max(t - flank, 0):t + flank] |= ALL_BITS
    else:
        raise RuntimeError("label repair did not converge")
    relaxed = int((bits == ALL_BITS).sum())
    return bits, relaxed, first_violations or 0


def bits_to_bool(b):
    """uint32 [..] -> bool [.., S]"""
    return ((b[..., None] >> np.arange(g1.S, dtype=np.uint32)) & 1).astype(bool)
