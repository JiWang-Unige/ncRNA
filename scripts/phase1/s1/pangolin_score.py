#!/usr/bin/env python3
"""S1R (PLAN_BATCH2.md A): Pangolin (tkzeng/Pangolin 5cf94b8, unpacked to ext/pangolin) splice scores for every base
of a chromosome, both strands, in the s1_analyze.py column order (acc+, don+, acc-, don-).
Pangolin's 12 standard models (final.{1,2,3}.{0,2,4,6}.3.v2: heart, liver, brain, testis x 3 replicates) give a
P(splice site) per position that does not tell donors from acceptors; per tissue the 3 replicates are averaged and
the site score is the maximum over the 4 tissues; the strand's score fills both of its columns.  The position
convention is fixed by --offsets: chosen once on chr20 (offset -2..+2 per column with the best AP on v50 mRNA sites
among motif-allowed positions, as for Vipsania on chr18) and reused on the other chromosomes.

  pangolin_score.py --codes data/p1/chr20_codes.npy --chrom chr20 --pick-offsets --out outputs/PHASE1/S1/pang_chr20.npy
  pangolin_score.py --codes data/p1/chr18_codes.npy --chrom chr18 --offsets outputs/PHASE1/S1/pang_chr20.npy.json --out ...
"""
import argparse
import importlib.util
import json
import os
import sys
import time

import numpy as np
import torch

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1/splice_diag"))
CTX = 5000
CHAN = (1, 4, 7, 10)          # P(splice) output channel of tissue model 0 / 2 / 4 / 6
TISSUES = (0, 2, 4, 6)


def load_models(device):
    spec = importlib.util.spec_from_file_location("pangolin_model", f"{ROOT}/ext/pangolin/pangolin/model.py")
    pm = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pm)
    models = {}
    for t in TISSUES:
        for j in (1, 2, 3):
            m = pm.Pangolin(pm.L, pm.W, pm.AR)
            m.load_state_dict(torch.load(f"{ROOT}/ext/pangolin/pangolin/models/final.{j}.{t}.3.v2", map_location=device))
            models.setdefault(t, []).append(m.to(device).eval())
    return models


@torch.no_grad()
def score_strand(models, codes, strand, device, block=100_000):
    import sai_torch
    L = len(codes)
    x = sai_torch.one_hot(codes)
    if strand == "-":
        x = x[::-1, ::-1].copy()
    xp = np.zeros((4, L + 2 * CTX), np.float32)
    xp[:, CTX:CTX + L] = x
    out = np.zeros(L, np.float32)
    for s in range(0, L, block):
        e = min(s + block, L)
        xt = torch.as_tensor(xp[None, :, s:e + 2 * CTX], device=device)
        per = [torch.stack([m(xt)[0, CHAN[k]] for m in models[t]]).mean(0) for k, t in enumerate(TISSUES)]
        out[s:e] = torch.stack(per).max(0).values.float().cpu().numpy()
    return out[::-1].copy() if strand == "-" else out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--codes", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--v50", default="data/human/gencode.v50.annotation.gtf.gz")
    ap.add_argument("--pick-offsets", action="store_true")
    ap.add_argument("--offsets")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    torch.cuda.set_per_process_memory_fraction(float(os.environ.get("NCRNA_GPU_FRAC", "0.25")), 0)
    import own_splice
    import sai_score
    import splice_diag
    from ncrna_annot.eval import ruler
    t0 = time.time()
    codes = np.asarray(np.load(f"{ROOT}/{a.codes}", mmap_mode="r"))
    models = load_models("cuda")
    raw = {s: score_strand(models, codes, s, "cuda") for s in ("+", "-")}
    meta = {"chrom": a.chrom, "pangolin_commit": "5cf94b8db938c658391b4305cd7ce33297d44ff7", "score_s": round(time.time() - t0)}
    base = np.stack([raw["+"], raw["+"], raw["-"], raw["-"]], 1)          # acc+, don+, acc-, don-
    truth = ruler.read_gtf(f"{ROOT}/{a.v50}", a.chrom)
    ts = {cls: sai_score.truth_sites(truth, cls) for cls in ("mrna", "lnc")}
    if a.pick_offsets:
        ok = splice_diag.site_allowed(codes)
        offs = []
        for col in range(4):
            cand = np.flatnonzero(ok[:, col])
            pos = ts["mrna"][col]
            keep = cand[~np.isin(cand, np.setdiff1d(ts["lnc"][col], pos))]
            lab = np.isin(keep, pos)
            aps = {o: own_splice.average_precision(np.roll(base[:, col], o)[keep], lab) for o in range(-2, 3)}
            offs.append(max(aps, key=lambda o: aps[o]))
            meta.setdefault("offset_ap", {})[own_splice.COLS[col]] = {str(o): v for o, v in aps.items()}
    else:
        offs = json.load(open(f"{ROOT}/{a.offsets}"))["offsets"]
    meta["offsets"] = [int(o) for o in offs]
    sc = np.stack([np.roll(base[:, c], offs[c]) for c in range(4)], 1).astype(np.float16)
    out = f"{ROOT}/{a.out}"
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.save(out, sc)
    meta["ap"] = own_splice.site_ap(np.asarray(sc, np.float32), codes, lambda cls, col: ts[cls][col])
    meta["total_s"] = round(time.time() - t0)
    json.dump(meta, open(out + ".json", "w"), indent=1)
    print(json.dumps({k: meta[k] for k in ("chrom", "offsets", "score_s", "total_s")}), flush=True)


if __name__ == "__main__":
    main()
