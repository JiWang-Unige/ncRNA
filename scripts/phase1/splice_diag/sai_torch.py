"""SpliceAI-10k (Jaganathan et al. 2019, package spliceai 1.3.1) inference in PyTorch, for the chr18 splice
diagnostic only (not part of any frozen model).

The Keras graph and weights are exported once with a Keras 2.2.4 environment (sai_export.py: <i>.json = model.to_json(),
<i>.npz = layer weights, ref<i>.npy = Keras output on a fixed random input); this module rebuilds the same graph
(InputLayer, Conv1D same/valid with dilation, BatchNormalization eps 1e-3 inference mode, ReLU, Add, Cropping1D,
softmax) and is checked against the Keras reference (check()).  Output channels: 0 neither, 1 acceptor, 2 donor.
The network is fully convolutional with a receptive field of exactly +-5000 nt, and the last 5000 nt on each side are
cropped, so scoring an input of length n + 10000 gives the same n outputs as SpliceAI's own 5000-nt blocks.
"""
import json

import numpy as np
import torch
import torch.nn.functional as F

CTX = 5000


class KerasGraph(torch.nn.Module):
    def __init__(self, cfg_json, npz):
        super().__init__()
        cfg = json.load(open(cfg_json))["config"]
        w = np.load(npz)
        self.layers = []
        self.params = torch.nn.ParameterDict()
        for l in cfg["layers"]:
            c, name = l["config"], l["name"]
            ins = [n[0] for n in l["inbound_nodes"][0]] if l["inbound_nodes"] else []
            op = {"cls": l["class_name"], "name": name, "ins": ins}
            if op["cls"] == "Conv1D":
                k, d = c["kernel_size"][0], c["dilation_rate"][0]
                op.update(pad=(d * (k - 1) // 2) if c["padding"] == "same" else 0, dil=d, act=c["activation"])
                assert c["padding"] in ("same", "valid") and k % 2 == 1 and c["strides"] == [1]
                self.params[name + "_w"] = torch.nn.Parameter(torch.as_tensor(w[f"{name}/0"]).permute(2, 1, 0).contiguous(), False)
                self.params[name + "_b"] = torch.nn.Parameter(torch.as_tensor(w[f"{name}/1"]), False)
            elif op["cls"] == "BatchNormalization":
                g, b, m, v = (torch.as_tensor(w[f"{name}/{j}"]).double() for j in range(4))
                scale = g / torch.sqrt(v + c["epsilon"])
                self.params[name + "_s"] = torch.nn.Parameter(scale.float(), False)
                self.params[name + "_t"] = torch.nn.Parameter((b - m * scale).float(), False)
            elif op["cls"] == "Cropping1D":
                op["crop"] = tuple(c["cropping"])
            elif op["cls"] == "Activation":
                op["act"] = c["activation"]
            else:
                assert op["cls"] in ("InputLayer", "Add"), op["cls"]
            self.layers.append(op)
        self.out = cfg["output_layers"][0][0]

    def forward(self, x):
        """x [B, 4, L] one-hot (A, C, G, T; N = 0) -> [B, 3, L - 10000] class probabilities."""
        t = {}
        for op in self.layers:
            cls, a = op["cls"], [t[i] for i in op["ins"]]
            if cls == "InputLayer":
                y = x
            elif cls == "Conv1D":
                y = F.conv1d(a[0], self.params[op["name"] + "_w"], self.params[op["name"] + "_b"], padding=op["pad"],
                             dilation=op["dil"])
                y = torch.softmax(y, 1) if op["act"] == "softmax" else y
                assert op["act"] in ("linear", "softmax")
            elif cls == "BatchNormalization":
                y = a[0] * self.params[op["name"] + "_s"][:, None] + self.params[op["name"] + "_t"][:, None]
            elif cls == "Activation":
                assert op["act"] == "relu"
                y = torch.relu(a[0])
            elif cls == "Add":
                y = a[0] + a[1]
            elif cls == "Cropping1D":
                lo, hi = op["crop"]
                y = a[0][..., lo:a[0].shape[-1] - hi]
            t[op["name"]] = y
        return t[self.out]


def load_models(export_dir, device="cuda"):
    return [KerasGraph(f"{export_dir}/spliceai{i}.json", f"{export_dir}/spliceai{i}.npz").to(device).eval()
            for i in range(1, 6)]


def one_hot(codes):
    """int8 codes (A0 C1 G2 T3, other 4) -> float32 [4, L]."""
    c = np.asarray(codes)
    x = np.zeros((4, len(c)), np.float32)
    ok = c < 4
    x[c[ok], np.flatnonzero(ok)] = 1.0
    return x


@torch.no_grad()
def score(models, codes, strand, block=200_000, device="cuda"):
    """Ensemble-mean (acceptor, donor) probabilities for every base of `codes` on one strand -> float32 [L, 2].
    strand "-": the reverse complement is scored and mapped back, so out[p] refers to genomic position p
    (acceptor = first exon base, donor = last exon base, in transcript orientation, as SpliceAI labels them)."""
    L = len(codes)
    x = one_hot(codes)
    if strand == "-":
        x = x[::-1, ::-1].copy()  # complement (A<->T, C<->G = reversed channel order) and reverse
    xp = np.zeros((4, L + 2 * CTX), np.float32)
    xp[:, CTX:CTX + L] = x
    out = np.zeros((L, 2), np.float32)
    for s in range(0, L, block):
        e = min(s + block, L)
        xt = torch.as_tensor(xp[None, :, s:e + 2 * CTX], device=device)
        y = sum(m(xt) for m in models) / len(models)
        out[s:e] = y[0, 1:].T.float().cpu().numpy()
    return out[::-1].copy() if strand == "-" else out


def check(export_dir, device="cuda"):
    """max |torch - keras| over the 5 models on the exported reference input."""
    x = np.load(f"{export_dir}/ref_input.npy")[0].T  # [4, L]
    xt = torch.as_tensor(x[None], device=device)
    worst = 0.0
    for i, m in enumerate(load_models(export_dir, device), 1):
        with torch.no_grad():
            y = m(xt)[0].T.cpu().numpy()
        worst = max(worst, float(np.abs(y - np.load(f"{export_dir}/ref{i}.npy")[0]).max()))
    return worst
