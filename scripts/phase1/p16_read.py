#!/usr/bin/env python3
"""P1.6 zero-shot panel read (descriptive): the three frozen v1.9 main runs -> posteriors, MEAX genes and GFF3.

  p16_read.py --config configs/p16_v19.json --target chr18 [--runs g1b-42-ke]     reproduction of the P1.5 chr18 read
  p16_read.py --config configs/p16_v19.json --target chicken_5 | zebrafish_22 | fly_2L | worm_I

Same code path as p15_read.py (its posteriors(), the config's decoding, lnc < 200 nt dropped, gene confidence,
GFF3 of the mRNA genes + lnc genes with confidence >= tau*); outputs <out_dir>/<target>/<run>/ with post.npy,
genes.json, pred.gff3, read_meta.json (written last; a run that has it is skipped).  chr18: post.npy and genes.json
must have the sha256 of the P1.5 batch-2 chr18 read; when all main runs match, <out_dir>/chr18/REPRO_OK is written
(config sha256), which the panel targets require.  No truth is read here.  Exit 75 when the GPU gate is busy.
"""
import argparse
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p15_read as R  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--target", required=True)
    ap.add_argument("--runs", default="")
    a = ap.parse_args()
    cfg_path = a.config if os.path.isabs(a.config) else f"{ROOT}/{a.config}"
    cfg, cfg_sha = json.load(open(cfg_path)), R.sha256(cfg_path)
    out_dir = f"{ROOT}/{cfg['out_dir']}"
    assert a.target in cfg["targets"], a.target
    if a.target != "chr18":
        ok = f"{out_dir}/chr18/REPRO_OK"
        if not os.path.exists(ok) or json.load(open(ok)).get("config_sha256") != cfg_sha:
            sys.exit(f"{a.target} refused: no chr18 reproduction for this config ({ok})")
    R.check_code(cfg)
    import torch
    torch.cuda.set_per_process_memory_fraction(float(os.environ.get("NCRNA_GPU_FRAC", "0.30")), 0)
    from ncrna_annot.decode import chrom
    from ncrna_annot.grammar import g1
    from ncrna_annot.train import screen2 as s2
    feats, codes, tg = R.load_target(cfg, a.target)
    d = cfg["decode"]
    Tg = g1.G1Transitions()().detach().double().numpy()
    for name in [r for r in a.runs.split(",") if r] or tg["runs"]:
        run, out = cfg["runs"][name], f"{out_dir}/{a.target}/{name}"
        if os.path.exists(f"{out}/read_meta.json"):
            print(json.dumps({"skip": name, "done": True}), flush=True)
            continue
        if R.gate_busy():
            sys.exit(75)
        os.makedirs(out, exist_ok=True)
        bad = [p for p, h in zip(run["ckpts"], run["ckpt_sha256"]) if R.sha256(f"{ROOT}/{p}") != h]
        if bad:
            sys.exit(f"{name}: checkpoint hash mismatch {bad}")
        ck = s2.load_ckpt([f"{ROOT}/{p}" for p in run["ckpts"]], "cuda")
        t0 = time.time()
        meta = {"name": name, "target": a.target, "config": a.config, "config_sha256": cfg_sha, "run": run,
                "host": os.uname().nodename, "torch_threads": torch.get_num_threads(),
                "started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
        post, L = R.posteriors(run, ck, feats, codes, "cuda")
        meta["post_s"] = round(time.time() - t0)
        np.save(f"{out}/post.npy", post)
        t1 = time.time()
        genes = chrom.decode_blocks(post, codes[:L], Tg, d["min_exon"], d["min_intron"], mode="mea", gamma=d["gamma"],
                                    level=d["level"], log_gain=False)
        genes = [g for g in genes if g["type"] != "l" or sum(b - a_ + 1 for a_, b in g["exons"]) >= d["lnc_min_len"]]
        for g in genes:
            g["conf"] = chrom.gene_confidence(g, post)
        meta["decode_s"] = round(time.time() - t1)
        json.dump([{"type": g["type"], "strand": g["strand"], "exons": [list(x) for x in g["exons"]], "conf": g["conf"]}
                   for g in genes], open(f"{out}/genes.json", "w"))
        kept = [g for g in genes if g["type"] != "l" or g["conf"] >= d["tau"]]
        R.write_gff3(kept, tg["seqid"], f"{out}/pred.gff3", "ncrna_annot_v19")
        meta.update(L=L, n_genes=len(genes), n_pred_lnc_all=sum(g["type"] == "l" for g in genes),
                    n_pred_lnc_kept=sum(g["type"] == "l" for g in kept), n_pred_mrna=sum(g["type"] == "m" for g in kept),
                    post_sha256=R.sha256(f"{out}/post.npy"), genes_sha256=R.sha256(f"{out}/genes.json"),
                    read_s=round(time.time() - t0), mb_per_h=round(L / 1e6 / max(time.time() - t0, 1) * 3600, 1))
        json.dump(meta, open(f"{out}/read_meta.json", "w"), indent=1)
        print(json.dumps({"done": name, **{k: meta[k] for k in ("post_s", "decode_s", "n_pred_lnc_kept", "n_pred_mrna")}}),
              flush=True)
        del post
    if a.target == "chr18":
        rep = {}
        for name, ref in cfg["chr18_reference"].items():
            m = f"{out_dir}/chr18/{name}/read_meta.json"
            m = json.load(open(m)) if os.path.exists(m) else {}
            rep[name] = {"post_equal": m.get("post_sha256") == ref["post_sha256"], "genes_equal": m.get("genes_sha256") == ref["genes_sha256"]}
        ok = all(v["post_equal"] and v["genes_equal"] for v in rep.values())
        res = {"config_sha256": cfg_sha, "ok": ok, "runs": rep, "time": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
        json.dump(res, open(f"{out_dir}/chr18/REPRO.json", "w"), indent=1)
        if ok:
            json.dump(res, open(f"{out_dir}/chr18/REPRO_OK", "w"), indent=1)
        print(json.dumps(res), flush=True)
    print("RUN_DONE", flush=True)


if __name__ == "__main__":
    main()
