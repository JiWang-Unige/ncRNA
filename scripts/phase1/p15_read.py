#!/usr/bin/env python3
"""Protocol v1.9 section 16, P1.5 read: frozen runs -> posteriors, MEAX genes and GFF3 on one target chromosome.

  p15_read.py --name p15-chr18 --config configs/p15_v19.json --target chr18 [--runs g1b-42-ke,g1b-123-ke]
  p15_read.py --name p15-chr20 --config configs/p15_v19.json --target chr20
  p15_read.py --name p15-t2t --config configs/p15_v19.json --target t2t_chr20
(queue: `<name> --script p15_read.py --config ... --target ...`, see run_formal_ont.sh)

Per run (default: every run the config lists for the target): the sha256 of its epoch weights and of the code files
must match the config; the weights are averaged (screen2.load_ckpt); posteriors by the float64 CPU forward-backward
on the cached features, i.e. the code path of train_formal_ont.py --eval-only --exact-post (same emissions, same
tiling); exon-level MEA at gamma* with the grammar's allowed transitions (as dev_decode.py); lnc genes < 200 nt
dropped; gene confidence = chrom.gene_confidence.  Writes <out>/<target>/<run>/, <out> = the config's "out_dir"
(default outputs/PHASE1/P15; the second read batch of 2026-11-01 uses its own directory):
  post.npy        float16 [L, S] posteriors
  genes.json      every decoded gene with its confidence (before tau*)
  pred.gff3       mRNA genes + lnc genes with confidence >= tau* (seqid of the target, `confidence` attribute)
  read_meta.json  hashes, settings, timings (written last: a run with read_meta.json is skipped on a restart)
AB6 (arm G2PH) is decoded and scored by screen_g2ph's own evaluation (result.json + pred.gff3, no tau*).
No truth is read here except inside that G2PH evaluation; scoring is p15_score.py (CPU).
chr20 / t2t_chr20 are refused until p15_score.py --target chr18 has written <out>/chr18/REPRO_OK for this config
(its sha256).  Exit 75 when the GPU gate is busy before a run (run_formal_ont.sh restarts the job).
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def gate_busy():
    r = subprocess.run(["bash", "-c", f". {ROOT}/scripts/phase1/gpu_gate.sh; gpu_why_busy"], capture_output=True, text=True)
    return r.stdout.strip() if r.returncode == 0 else ""


def check_code(cfg):
    bad = [p for p, h in cfg["code_files"].items() if sha256(f"{ROOT}/{p}") != h]
    if bad:
        sys.exit(f"code on disk differs from the config: {bad}")


def load_target(cfg, target):
    t = cfg["targets"][target]
    cover = sorted(json.load(open(f"{ROOT}/{t['cover']}")))
    assert cover[0][0] == 0 and all(a == b for (_, b), (a, _) in zip(cover, cover[1:])) and cover[-1][1] == t["tiles"], \
        f"{target}: feature cache incomplete {cover}"
    feats = np.load(f"{ROOT}/{t['feats']}", mmap_mode="r")
    codes = np.load(f"{ROOT}/{t['codes']}", mmap_mode="r")
    return feats, codes, t


def write_gff3(genes, seqid, path, source):
    kind = {"l": "lnc_RNA", "m": "mRNA"}
    with open(path, "w") as fh:
        fh.write("##gff-version 3\n")
        for i, g in enumerate(genes):
            ex, gid, c = g["exons"], f"{source}_{seqid}_{i}", f"{g['conf']:.6f}"
            a, b = ex[0][0], ex[-1][1]
            fh.write(f"{seqid}\t{source}\tgene\t{a}\t{b}\t{c}\t{g['strand']}\t.\tID={gid};confidence={c}\n")
            fh.write(f"{seqid}\t{source}\t{kind[g['type']]}\t{a}\t{b}\t{c}\t{g['strand']}\t.\tID={gid}.t1;Parent={gid};confidence={c}\n")
            for j, (x, y) in enumerate(ex):
                fh.write(f"{seqid}\t{source}\texon\t{x}\t{y}\t.\t{g['strand']}\t.\tID={gid}.t1.e{j};Parent={gid}.t1\n")


def posteriors(run, ck, feats, codes, device):
    """Exact posteriors of a G1B / R run: the emission and marginal code of screen_g1b / screen2.evaluate_chromosome."""
    import torch
    from ncrna_annot.decode import chrom
    from ncrna_annot.grammar import g1
    from ncrna_annot.heads.c0 import C0BHead, C0Head
    from ncrna_annot.train import screen2 as s2
    from ncrna_annot.train import screen_g1b as sg
    L = len(feats) * sg.TOK
    trans = g1.G1Transitions().to(device)
    trans.load_state_dict(ck["trans"])
    coef = None
    if run["arm"] == "G1B":
        head, coef = C0BHead(feats.shape[-1]).to(device), g1.BoundaryCoef().to(device)
        coef.load_state_dict(ck["coef"])
    else:
        head = C0Head(feats.shape[-1]).to(device)
    head.load_state_dict(ck["head"])
    head.eval()

    def emit_fn(s, e):
        f = torch.as_tensor(np.asarray(feats[s // sg.TOK:e // sg.TOK], np.float32), device=device)[None]
        c = torch.as_tensor(np.asarray(codes[s:e], np.int64), device=device)[None]
        with torch.no_grad():
            if run["arm"] == "G1B":
                p, b = head(f, c)
                return sg._mask_lnc_only(p, run["lnc_only"])[0], b[0]
            return s2._emit(head, f, c, run["lnc_only"])[0]

    T = trans().detach()
    return chrom.chromosome_marginals(emit_fn, L, T, codes[:L], device, coef=coef, exact=True), L


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--target", required=True, help="a target of the config (chr18, chr20, t2t_chr20; v2.0: chr9)")
    ap.add_argument("--runs", default="", help="comma list (default: the config's runs for this target)")
    a = ap.parse_args()
    cfg_path = a.config if os.path.isabs(a.config) else f"{ROOT}/{a.config}"
    cfg, cfg_sha = json.load(open(cfg_path)), sha256(cfg_path)
    P15 = f"{ROOT}/{cfg.get('out_dir', 'outputs/PHASE1/P15')}"
    if a.target not in cfg["targets"]:
        sys.exit(f"{a.target} is not a target of {a.config}")
    if a.target != "chr18":
        ok = f"{P15}/chr18/REPRO_OK"
        if not os.path.exists(ok) or json.load(open(ok)).get("config_sha256") != cfg_sha:
            sys.exit(f"{a.target} refused: no chr18 reproduction for this config ({ok})")
    check_code(cfg)
    runs = [r for r in a.runs.split(",") if r] or cfg["targets"][a.target]["runs"]
    import torch
    torch.cuda.set_per_process_memory_fraction(float(os.environ.get("NCRNA_GPU_FRAC", "0.30")), 0)
    from ncrna_annot.decode import chrom
    from ncrna_annot.grammar import g1
    from ncrna_annot.train import screen2 as s2
    feats, codes, tg = load_target(cfg, a.target)
    d = cfg["decode"]
    Tg = g1.G1Transitions()().detach().double().numpy()  # mea_duration uses only the allowed transitions
    for name in runs:
        run = cfg["runs"][name]
        out = f"{P15}/{a.target}/{name}"
        if os.path.exists(f"{out}/read_meta.json"):
            print(json.dumps({"skip": name, "done": True}), flush=True)
            continue
        if gate_busy():
            sys.exit(75)
        os.makedirs(out, exist_ok=True)
        cks = [f"{ROOT}/{p}" for p in run["ckpts"]]
        bad = [p for p, h in zip(run["ckpts"], run["ckpt_sha256"]) if sha256(f"{ROOT}/{p}") != h]
        if bad:
            sys.exit(f"{name}: checkpoint hash mismatch {bad}")
        ck = s2.load_ckpt(cks, "cuda")
        t0 = time.time()
        meta = {"name": name, "target": a.target, "config": a.config, "config_sha256": cfg_sha, "run": run,
                "host": os.uname().nodename, "torch_threads": torch.get_num_threads(),
                "started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
        if run["arm"] == "G2PH":  # AB6: G2-PH's own decoding and scoring (protocol v1.8 section 15, report only)
            from ncrna_annot.heads.c0 import C0BHead
            from ncrna_annot.train import screen_g2ph as sp
            space = sp.g2ph.Space(*sp.CONFIGS["G2PH"])
            head, prm = C0BHead(feats.shape[-1]).to("cuda"), sp.g2ph.G2PHParams(space).to("cuda")
            head.load_state_dict(ck["head"])
            prm.load_state_dict(ck["grammar"])
            ch = dict(feats=feats, codes=codes, truth_gtf=f"{ROOT}/{tg['truth']}", chrom=tg["truth_chrom"],
                      gff_out=f"{out}/pred.gff3")
            res = {"arm": "G2PH", "epochs": ck["epoch"], **sp.evaluate_chromosome(head, prm, space, ch, "cuda")}
            json.dump(res, open(f"{out}/result.json", "w"), indent=1)
            meta["read_s"] = round(time.time() - t0)
            json.dump(meta, open(f"{out}/read_meta.json", "w"), indent=1)
            print(json.dumps({"done": name, "read_s": meta["read_s"]}), flush=True)
            continue
        post, L = posteriors(run, ck, feats, codes, "cuda")
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
        write_gff3(kept, tg["seqid"], f"{out}/pred.gff3", "ncrna_annot_v19")
        meta.update(L=L, n_genes=len(genes), n_pred_lnc_all=sum(g["type"] == "l" for g in genes),
                    n_pred_lnc_kept=sum(g["type"] == "l" for g in kept), n_pred_mrna=sum(g["type"] == "m" for g in kept),
                    post_sha256=sha256(f"{out}/post.npy"), read_s=round(time.time() - t0),
                    mb_per_h=round(L / 1e6 / max(time.time() - t0, 1) * 3600, 1))
        json.dump(meta, open(f"{out}/read_meta.json", "w"), indent=1)
        print(json.dumps({"done": name, **{k: meta[k] for k in ("post_s", "decode_s", "n_pred_lnc_kept", "n_pred_mrna")}}),
              flush=True)
        del post
    print("RUN_DONE", flush=True)


if __name__ == "__main__":
    main()
