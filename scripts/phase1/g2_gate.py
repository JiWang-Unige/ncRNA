#!/usr/bin/env python3
"""G2-0 feasibility gate (protocol v1.6 section 13), CPU on ont.

usage: g2_gate.py <cand_dir with {train,val}_cand_scores.npy> <out.json> [splice_target]
Thresholds on the training windows (splice recall >= 99% exact, ends recall >= 95% within +-50 bp), applied to the
300 validation windows.  Gate: recall of lnc canonical introns (L and R both candidates) >= 90%.
If it fails, the protocol allows ONE re-run with splice_target = 0.995.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ncrna_annot.decode import g2cand  # noqa: E402

ROOT = os.path.expanduser("~/ncrna-annot")
K_LOOKBACK = 256


def event_values(sc, ex):
    vals = {c: [] for c in range(8)}
    for w in np.unique(ex["ev_w"]):
        sel = ex["ev_w"] == w
        ch = np.array([g2cand.EV_TO_CAND[int(c)] for c in ex["ev_ch"][sel]])
        v = g2cand.true_event_values(np.asarray(sc[w], np.float32), ex["ev_pos"][sel], ch)
        for c, x in zip(ch, v):
            vals[c].append(x)
    return vals


def main():
    cdir, out = sys.argv[1:3]
    if len(sys.argv) > 3:
        g2cand.TARGET["splice"] = float(sys.argv[3])
    sc = {sp: np.load(f"{cdir}/{sp}_cand_scores.npy", mmap_mode="r") for sp in ("train", "val")}
    ex = {sp: dict(np.load(f"{ROOT}/data/p1g/{sp}_exact.npz")) for sp in ("train", "val")}
    theta = g2cand.choose_theta(event_values(sc["train"], ex["train"]))
    res = {"targets": dict(g2cand.TARGET), "theta": dict(zip(g2cand.CHANNELS, theta.tolist()))}
    # validation: per-channel recall, density, lnc intron recall, intron-chain ceiling
    V = ex["val"]
    vv = event_values(sc["val"], V)
    res["val_event_recall"] = {g2cand.CHANNELS[c]: float(np.mean(np.asarray(vv[c]) >= theta[c])) if vv[c] else None for c in range(8)}
    n_cand = np.zeros(8)
    intron_ok = np.zeros(len(V["in_w"]), bool)
    reach_ok = np.zeros(len(V["in_w"]), bool)
    nodes = []
    for w in range(len(sc["val"])):
        cands = g2cand.candidates(np.asarray(sc["val"][w], np.float32), theta)
        n_cand += [len(x) for x in cands]
        nodes.append(len(set(np.concatenate([cands[c] + (1 if c in (1, 3, 5, 7) else 0) for c in range(8)]).tolist())) + 2)
        sel = np.flatnonzero(V["in_w"] == w)
        for i in sel:
            s = V["in_s"][i]
            Lc, Rc = (cands[0], cands[1]) if s == 0 else (cands[2], cands[3])
            intron_ok[i] = V["in_L"][i] in set(Lc.tolist()) and V["in_R"][i] in set(Rc.tolist())
            # look-back reach: the intron's start node must be among the last K intron-start nodes before its end node
            n_between = int(np.sum((Lc >= V["in_L"][i]) & (Lc <= V["in_R"][i])))
            reach_ok[i] = intron_ok[i] and n_between <= K_LOOKBACK
    res["val_candidates_per_kb"] = dict(zip(g2cand.CHANNELS, (n_cand / (len(sc["val"]) * sc["val"].shape[1] / 1000)).tolist()))
    res["val_nodes_per_window"] = {"median": float(np.median(nodes)), "q90": float(np.quantile(nodes, 0.9))}
    for k, name in ((1, "lnc"), (0, "mrna")):
        m = V["in_k"] == k
        res[f"val_{name}_intron_recall"] = {"n": int(m.sum()), "recall": float(intron_ok[m].mean()) if m.any() else None,
                                            "recall_and_reachable_K256": float(reach_ok[m].mean()) if m.any() else None}
        # chain ceiling: complete multi-exon canonical transcripts inside a window with all introns recalled
        txs = [t for t in range(len(V["tx_w"])) if V["tx_k"][t] == k and V["tx_inside"][t] and V["tx_n_intron"][t] > 0]
        ok = [bool(np.all(intron_ok[V["in_tx"] == t])) for t in txs]
        res[f"val_{name}_chain_ceiling"] = {"n": len(txs), "frac": float(np.mean(ok)) if ok else None}
    res["gate_lnc_intron_recall_ge_0.90"] = bool(res["val_lnc_intron_recall"]["recall"] >= 0.90)
    json.dump(res, open(out, "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
