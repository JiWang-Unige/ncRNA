"""U-ev unknown intervals (Phase 0 protocol v2 definition, reused unchanged by Phase 1 v1).

Snaptron srav3h (recount3, 2021) junctions with samples_count >= 100, coverage_sum >= 500,
canonical motif, length <= 500 kb, no overlap with any annotated gene span of the reference
annotation; unknown = [start - 500, end + 500] merged, minus all annotated gene spans.
"""
import urllib.request

import numpy as np

SNAP = ("https://snaptron.cs.jhu.edu/srav3h/snaptron?regions={c}:1-{n}&rfilter=samples_count>:99"
        "&rfilter=coverage_sum>:499&fields=chromosome,start,end,strand,left_motif,right_motif,samples_count,coverage_sum")
CANON = {"GT-AG", "CT-AC", "GC-AG", "CT-GC", "AT-AC", "GT-AT"}


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def subtract(iv, cut):
    cut, out = merge(cut), []
    for s, e in iv:
        cur = s
        for cs, ce in cut:
            if ce <= cur or cs >= e:
                continue
            if cs > cur:
                out.append([cur, cs])
            cur = max(cur, ce)
        if cur < e:
            out.append([cur, e])
    return out


def fetch_junctions(chrom, n, cache=None):
    if cache:
        try:
            return open(cache).read().splitlines()
        except FileNotFoundError:
            pass
    rows = urllib.request.urlopen(SNAP.format(c=chrom, n=n), timeout=900).read().decode().splitlines()[1:]
    if cache:
        open(cache, "w").write("\n".join(rows))
    return rows


def unknown_mask(chrom, n, spans, cache=None):
    """bool [n] U-ev mask given annotated gene spans (list of [s, e) 0-based)."""
    gs = merge(spans)
    starts = np.array([s for s, _ in gs]) if gs else np.zeros(0, int)
    keep = []
    for r in fetch_junctions(chrom, n, cache):
        f = r.split("\t")
        s, e = int(f[2]) - 1, int(f[3])
        if f"{f[5].upper()}-{f[6].upper()}" not in CANON or e - s > 500_000:
            continue
        i = np.searchsorted(starts, e) - 1
        if i >= 0 and gs[i][1] > s:
            continue
        keep.append([max(s - 500, 0), min(e + 500, n)])
    m = np.zeros(n, bool)
    for s, e in subtract(merge(keep), gs):
        m[s:e] = True
    return m
