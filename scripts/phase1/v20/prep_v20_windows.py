#!/usr/bin/env python3
"""v2.0 (protocol section 17, execution detail 1): replace the 154 v1.9 training windows on chr9 in place.

data/p1v20/train_windows.tsv   the 3,000 P1.1 windows in the same order; rows on chr9 replaced by a new window of the
                               same kind (lnc / pc / rand) drawn with prep_p1.py's rule from TRAIN without chr9
                               (chromosome by length; lnc / pc centred on a random v50 exon midpoint of that class,
                               rand uniform; > 10% N rejected), numpy default_rng(1), rows replaced in index order
data/p1v20/replaced.json       replaced row indices, old and new windows
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/..")
import prep_p1 as P  # noqa: E402

OLD = f"{P.ROOT}/data/p1/train_windows.tsv"
OUT = f"{P.ROOT}/data/p1v20"


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = [l.rstrip("\n").split("\t") for l in open(OLD)]
    rows = [(c, int(s), k) for c, s, k in rows]
    assert len(rows) == P.NTRAIN
    repl = [i for i, (c, _, _) in enumerate(rows) if c == "chr9"]
    assert len(repl) == 154, len(repl)
    train = [c for c in P.TRAIN if c != "chr9"]
    seqs = P.read_genome(set(train))
    tx = P.labels.read_transcripts(P.GTF["v50"], set(train))
    pos = {"lnc": {}, "pc": {}}
    for c in train:
        for kind, k in (("lnc", "l"), ("pc", "m")):
            # labels.read_transcripts now also returns is_canonical (v1.2); the transcript set and order are the same
            # as when prep_p1.py sampled the v1.9 windows, so the exon-midpoint pool is unchanged
            ex = [(a, b) for kk, _, exs, *_ in tx.get(c, []) if kk == k for a, b in exs]
            pos[kind][c] = np.array([(a + b) // 2 for a, b in ex], np.int64)
    codes = {c: P.g1.encode(seqs[c]) for c in train}
    lens = np.array([len(seqs[c]) for c in train], float)
    rng = np.random.default_rng(1)
    new = []
    for i in repl:
        kind = rows[i][2]
        while True:
            c = train[rng.choice(len(train), p=lens / lens.sum())]
            if kind == "rand":
                st = int(rng.integers(0, len(seqs[c]) - P.WIN))
            else:
                if not len(pos[kind][c]):
                    continue
                st = int(rng.choice(pos[kind][c])) - P.WIN // 2
            st = min(max(st, 0), len(seqs[c]) - P.WIN)
            if (codes[c][st:st + P.WIN] == 4).sum() > P.WIN // 10:
                continue
            break
        new.append((i, rows[i], (c, st, kind)))
        rows[i] = (c, st, kind)
    with open(f"{OUT}/train_windows.tsv", "w") as fh:
        fh.writelines(f"{c}\t{st}\t{k}\n" for c, st, k in rows)
    json.dump({"rule": "protocol section 17 detail 1; numpy default_rng(1)", "n_replaced": len(new),
               "replaced": [{"row": i, "old": list(o), "new": list(n)} for i, o, n in new],
               "kinds": {k: sum(1 for _, o, _ in new if o[2] == k) for k in ("lnc", "pc", "rand")},
               "new_chroms": {c: sum(1 for _, _, n in new if n[0] == c) for c in train}}, open(f"{OUT}/replaced.json", "w"), indent=1)
    print(json.dumps({"replaced": len(new), "kinds": {k: sum(1 for _, o, _ in new if o[2] == k) for k in ("lnc", "pc", "rand")}}), flush=True)


if __name__ == "__main__":
    main()
