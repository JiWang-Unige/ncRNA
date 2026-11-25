#!/usr/bin/env python
"""B0-lite small 1D CNN (PROTOCOL_B0.md section 5): input = +-25 nt one-hot + reference reactivity track
(the reference track is only available at observed positions; unobserved are 0 with a mask channel).
Same GroupKFold(5) by gene as b0_lite.py. Run with ~/envs/te_benchmark/bin/python on the ont GPU."""
import os, sys, json
import numpy as np, pandas as pd, torch, torch.nn as nn
from sklearn.model_selection import GroupKFold

ROOT = os.path.expanduser("~/ncrna-annot")
OUT = f"{ROOT}/outputs/context_pilot/B0"
W = 25
torch.manual_seed(0); np.random.seed(0)


def windows(df, tx_seq, ref_tracks):
    X = np.zeros((len(df), 6, 2 * W + 1), np.float32)
    lut = {"A": 0, "C": 1, "G": 2, "T": 3}
    for k, (t, i) in enumerate(zip(df.tid.values, df.pos.values)):
        s = tx_seq[t]; r = ref_tracks[t]
        for j in range(-W, W + 1):
            p = i + j
            if 0 <= p < len(s):
                c = lut.get(s[p])
                if c is not None:
                    X[k, c, j + W] = 1
                if p < len(r) and not np.isnan(r[p]):
                    X[k, 4, j + W] = r[p]; X[k, 5, j + W] = 1
    return X


class Net(nn.Module):
    def __init__(self, cin):
        super().__init__()
        self.f = nn.Sequential(nn.Conv1d(cin, 64, 7, padding=3), nn.ReLU(), nn.Conv1d(64, 64, 5, padding=2), nn.ReLU(),
                               nn.Conv1d(64, 32, 3, padding=1), nn.ReLU(), nn.AdaptiveAvgPool1d(1), nn.Flatten(), nn.Linear(32, 1))

    def forward(self, x):
        return self.f(x).squeeze(-1)


def train_pred(Xtr, ytr, Xte, dev, epochs=8):
    net = Net(Xtr.shape[1]).to(dev); opt = torch.optim.Adam(net.parameters(), 1e-3)
    Xt, yt = torch.tensor(Xtr), torch.tensor(ytr, dtype=torch.float32)
    for _ in range(epochs):
        perm = torch.randperm(len(Xt))
        for b in range(0, len(Xt), 2048):
            ix = perm[b:b + 2048]
            opt.zero_grad(); loss = ((net(Xt[ix].to(dev)) - yt[ix].to(dev)) ** 2).mean(); loss.backward(); opt.step()
    net.eval(); out = []
    with torch.no_grad():
        for b in range(0, len(Xte), 8192):
            out.append(net(torch.tensor(Xte[b:b + 8192]).to(dev)).cpu().numpy())
    return np.concatenate(out)


if __name__ == "__main__":
    task = sys.argv[1]
    sys.path.insert(0, os.path.dirname(__file__))
    # inputs exported by b0_lite.py (rnafate env) as npz/pickle-free files to avoid pandas version mismatch
    z = np.load(f"{OUT}/{task}_cnn_in.npz", allow_pickle=True)
    df = pd.DataFrame({"tid": z["tid"], "pos": z["pos"], "gid": z["gid"], "delta": z["delta"]})
    tx_seq = z["tx_seq"].item(); ref_tracks = z["ref_tracks"].item()
    X = windows(df, tx_seq, ref_tracks)
    y = df.delta.values.astype(np.float32)
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    pE, pS = np.zeros(len(df)), np.zeros(len(df))
    for tr, te in GroupKFold(5).split(df, groups=df.gid.values):
        pE[te] = train_pred(X[tr], y[tr], X[te], dev)
        pS[te] = train_pred(X[tr][:, :4], y[tr], X[te][:, :4], dev)
    np.save(f"{OUT}/{task}_cnn_E.npy", pE); np.save(f"{OUT}/{task}_cnn_S.npy", pS)
    print("saved CNN predictions; run b0_lite.py cnnsum", task)
