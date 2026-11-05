#!/usr/bin/env python3
"""S2-AG analysis A input (PLAN_S2_AG.md section 2): AlphaGenome SPLICE_SITES per base -> the s1_analyze.py score file.

Input: <raw>/ag_<chrom>_splice.npy [L, 4] float16 (acc+, don+, acc-, don-, no offsets), written by ag_colab_job.py and
fetched with its sha256 checked.  Position offsets as for Pangolin (pangolin_score.py): chosen once on chr20 (offset
-2..+2 per column with the best AP on v50 mRNA sites among motif-allowed positions, lnc-only sites left out) and reused
on the other chromosomes.

  ag_post.py --chrom chr20 --pick-offsets --out outputs/PHASE1/S2AG/ag_chr20.npy
  ag_post.py --chrom chr18 --offsets outputs/PHASE1/S2AG/ag_chr20.npy.json --out outputs/PHASE1/S2AG/ag_chr18.npy
"""
import argparse
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1/splice_diag"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--raw", default="outputs/PHASE1/S2AG/auto/raw")
    ap.add_argument("--v50", default="data/human/gencode.v50.annotation.gtf.gz")
    ap.add_argument("--pick-offsets", action="store_true")
    ap.add_argument("--offsets")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    import own_splice
    import sai_score
    import splice_diag
    from ncrna_annot.eval import ruler
    t0 = time.time()
    cfg = json.load(open(f"{ROOT}/configs/s2ag.json"))
    codes = np.asarray(np.load(f"{ROOT}/{cfg['codes'][a.chrom]}", mmap_mode="r"))
    base = np.load(f"{ROOT}/{a.raw}/ag_{a.chrom}_splice.npy").astype(np.float32)
    assert base.shape == (len(codes), 4), (base.shape, len(codes))
    done = json.load(open(f"{ROOT}/{a.raw}/{a.chrom}.done.json"))
    meta = {"chrom": a.chrom, "ag_commit": cfg["ag_commit"], "weights_revision": cfg["weights"]["hf_revision"],
            "raw_sha256": done["files"][f"ag_{a.chrom}_splice.npy"], "splice_cols": done.get("splice_cols")}
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
    print(json.dumps({k: meta[k] for k in ("chrom", "offsets", "total_s")}), flush=True)


if __name__ == "__main__":
    main()
