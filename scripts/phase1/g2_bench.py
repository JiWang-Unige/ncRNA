#!/usr/bin/env python3
"""G2-0 throughput check (protocol v1.6 section 13): lattice sizes and DP cost on real candidates.

usage: g2_bench.py <cand_dir> <gate.json> <out.json> [device] [n_windows] [viterbi_windows]
Uses the thresholds from g2_gate.py on validation windows with random potentials (cost does not depend on values):
node counts, time of segment scores + log Z + observed log Z + backward surrogate per window (level 3, batch 8),
and a Viterbi on concatenated windows (default ~1 Mb).  Extrapolates T4 hours per seed: 12 epochs x 2,000 windows + val.
"""
import json
import os
import sys
import time

import numpy as np
import torch

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ncrna_annot.decode import g2cand  # noqa: E402
from ncrna_annot.grammar import g1, g2  # noqa: E402

ROOT = os.path.expanduser("~/ncrna-annot")


def main():
    cdir, gate, out = sys.argv[1:4]
    dev = sys.argv[4] if len(sys.argv) > 4 else "cpu"
    nw = int(sys.argv[5]) if len(sys.argv) > 5 else 8
    theta = np.array([json.load(open(gate))["theta"][c] for c in g2cand.CHANNELS])
    sc = np.load(f"{cdir}/val_cand_scores.npy", mmap_mode="r")
    codes = np.load(f"{ROOT}/data/p1s_v12/val_codes.npy", mmap_mode="r")
    bits = np.load(f"{ROOT}/data/p1s_v12/val_bits_v50uev.npy", mmap_mode="r")
    L = sc.shape[1]
    rng = np.random.default_rng(0)
    idx = rng.choice(len(sc), nw, replace=False)
    lats = [g2.build_lattice(L, g2cand.candidates(np.asarray(sc[i], np.float32), theta)) for i in idx]
    nodes = [len(l["pos"]) for l in lats]
    from ncrna_annot.data import labels
    prm = g2.G2Params(level=3).to(dev)
    lat = g2.collate(lats, dev)
    x = torch.randn(nw, L, g1.S, device=dev, requires_grad=True)
    bnd = torch.randn(nw, L, 4, device=dev, requires_grad=True)
    h = torch.randn(nw, L, 256, device=dev, requires_grad=True)
    gates = torch.as_tensor(np.stack([g1.gate_mask(np.asarray(codes[i])) for i in idx]), device=dev)
    al = torch.as_tensor(labels.bits_to_bool(np.asarray(bits[idx])), device=dev)
    if dev == "cuda":
        torch.cuda.synchronize()
    t0 = time.time()
    SEG, IG = g2.segment_scores(x, gates, bnd, lat, prm, h=h)
    SEGo, IGo = g2.segment_scores(x, gates, bnd, lat, prm, h=h, allowed=al)
    t1 = time.time()
    loss, feas = g2.nll(SEG, IG, SEGo, IGo, lat)
    t2 = time.time()
    loss.sum().backward()
    if dev == "cuda":
        torch.cuda.synchronize()
    t3 = time.time()
    per_win = (t3 - t0) / nw
    res = {"device": dev, "n_windows": nw, "nodes_per_window": nodes, "nodes_median": float(np.median(nodes)),
           "K": int(lats[0]["SI"].shape[1]), "feasible_obs": feas.tolist(),
           "s_scores": t1 - t0, "s_dp": t2 - t1, "s_backward": t3 - t2, "s_per_window_train_step": per_win,
           "est_hours_per_seed_12ep": per_win * (12 * 2000 + 12 * 300) / 3600}
    # 1 Mb Viterbi: concatenate ~40 windows' candidates on a 1 Mb axis (cost proxy)
    nb = int(sys.argv[6]) if len(sys.argv) > 6 else 1_000_000 // L
    cat = [[] for _ in range(8)]
    for w in range(min(nb, len(sc))):
        cs = g2cand.candidates(np.asarray(sc[w], np.float32), theta)
        for c in range(8):
            cat[c].extend((cs[c] + w * L).tolist())
    LL = min(nb, len(sc)) * L
    t0 = time.time()
    latb = g2.collate([g2.build_lattice(LL, cat)], "cpu")
    with torch.no_grad():
        xb = torch.randn(1, LL, g1.S)
        SEGb, IGb = g2.segment_scores(xb, torch.zeros(1, LL, g1.S), torch.randn(1, LL, 4), latb,
                                      g2.G2Params(level=2), min_len=(20, 50))
        g2.viterbi(SEGb, IGb, latb)
    res.update({"viterbi_Mb": LL / 1e6, "viterbi_nodes": int(latb["pos"].shape[1]), "viterbi_s": time.time() - t0})
    json.dump(res, open(out, "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
