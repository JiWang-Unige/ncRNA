#!/usr/bin/env python3
"""P1-S S3 (protocol v1.3 section 10): decode-time minimum durations on R's best checkpoint (CPU on ont).

usage: decode_s3.py <best.pt> <out_dir> [min_exon min_intron] [lnc_only 0|1]
"""
import json
import os
import sys

import numpy as np
import torch

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ncrna_annot.decode import genes as dec  # noqa: E402
from ncrna_annot.grammar import g1  # noqa: E402
from ncrna_annot.heads.c0 import C0Head  # noqa: E402
from ncrna_annot.train import screen  # noqa: E402

ROOT = os.path.expanduser("~/ncrna-annot")


def main():
    ck, out = sys.argv[1], sys.argv[2]
    me, mi = (int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 else (20, 50)
    lnc_only = len(sys.argv) > 5 and sys.argv[5] == "1"
    os.makedirs(out, exist_ok=True)
    torch.set_num_threads(4)
    feats = np.load(f"{ROOT}/data/p05_feat/feat_val.npy", mmap_mode="r")
    codes = np.load(f"{ROOT}/data/p1s_v12/val_codes.npy")
    starts = [int(l.split()[1]) for l in open(f"{ROOT}/data/p1s_v12/val_windows.tsv")]
    head = C0Head(feats.shape[-1])
    trans = g1.G1Transitions()
    st = torch.load(ck, map_location="cpu")
    head.load_state_dict(st["head"])
    trans.load_state_dict(st["trans"])
    head.eval()
    T = trans().detach().double().numpy()
    genes = []
    with torch.no_grad():
        for i in range(len(codes)):
            e = head(torch.as_tensor(np.asarray(feats[i], np.float32))[None],
                     torch.as_tensor(codes[i].astype(np.int64))[None])[0].numpy()
            if lnc_only:
                e[:, screen.M_STATES] = g1.NEG
            path, _ = g1.viterbi_duration(e, T, g1.gate_mask(codes[i]), None, me, mi)
            genes += g1.path_to_genes(path, offset=starts[i])
    seen, uniq = set(), []
    for g in genes:
        k = (g["type"], g["strand"], tuple(map(tuple, g["exons"])))
        if k not in seen:
            seen.add(k)
            uniq.append(g)
    keep = [g for g in uniq if g["type"] != "l" or sum(b - a + 1 for a, b in g["exons"]) >= 200]
    dec.write_gff3(uniq, "chr18", f"{out}/val_s3_nofilter.gff3")
    dec.write_gff3(keep, "chr18", f"{out}/val_s3.gff3")
    truth = f"{ROOT}/data/p1s_v12/val_truth_v50.gtf"
    res = {"min_exon": me, "min_intron": mi, "checkpoint": ck, "lnc_only": lnc_only,
           "S3": screen.score(f"{out}/val_s3.gff3", truth, "chr18"),
           "S3_without_length_filter": screen.score(f"{out}/val_s3_nofilter.gff3", truth, "chr18")}
    json.dump(res, open(f"{out}/s3.json", "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
