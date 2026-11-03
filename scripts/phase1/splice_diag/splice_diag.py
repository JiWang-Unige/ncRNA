#!/usr/bin/env python3
"""chr18 splice diagnostic, step 2 (exploratory, chr18 only; user decision 2026-11-03 on card p15-next).

A frozen v1.9 main run (configs/p15_v19b.json: weights, gamma*, tau*, MEAX) is decoded on chr18 with an extra
log-potential on the intron boundary states, no retraining:
  L (first intron base) at i: w * phi(site i - 1), column donor + / acceptor - ;
  R (last intron base) at j:  w * phi(site j + 1), column acceptor + / donor -   (sai_score.py site convention)
  phi(p) = log max(s[p], 1e-6) - log pi,  pi = mean s over the sites whose boundary state the motif gate allows
(SpliceAI as an independent expert with its prior divided out); the same bonus for the mRNA and the lnc states.
--oracle replaces s by the chr18 truth splice sites (1 at every lnc / mRNA truth site, 1e-6 elsewhere): the
decoder's ceiling with perfect splice information.  w = 0 must reproduce the P1.5 chr18 posteriors bit for bit.

  splice_diag.py --run g1b-42-ke --w 1 --sai outputs/PHASE1/SPLICE_DIAG/sai_chr18.npy --out outputs/PHASE1/SPLICE_DIAG/g1b-42-ke_w1

Writes <out>/genes.json (all MEAX genes with confidence), result.json (screen2.summarize core, ruler L1-L6 of the
final predictions at tau* without bootstrap, guardrails, timings, pi).  SpliceAI was trained on chr18 (and chr20):
diagnostic only, never a frozen model or a chr20 number.
"""
import argparse
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p15_read as R  # noqa: E402
import p15_score as P  # noqa: E402
import sai_score  # noqa: E402
from ncrna_annot.grammar import g1  # noqa: E402

# boundary state (part, strand) -> (SpliceAI column, site offset from the state's base)
BOUNDARY = {("L", "+"): (1, -1), ("R", "+"): (0, +1), ("L", "-"): (2, -1), ("R", "-"): (3, +1)}


def site_allowed(codes):
    """[L, 4] bool: site p of column c can be used, i.e. the motif gate allows its boundary state."""
    d = np.full(len(codes) + 2, -1, np.int16)
    c = np.asarray(codes).astype(np.int16)
    d[:len(c) - 1] = c[:-1] * 5 + c[1:]  # d[i] = dinucleotide (i, i + 1)
    ok = np.zeros((len(c), 4), bool)
    for (part, s), (col, off) in BOUNDARY.items():
        motifs = list(g1._dinuc_set(g1.MOTIF_L[s] if part == "L" else g1.MOTIF_R[s]))
        st = np.isin(d[:len(c)], motifs) if part == "L" else np.isin(np.concatenate([[-1], d[:len(c) - 1]]), motifs)
        # state at t allowed -> site t + off allowed
        if off < 0:
            ok[:-1, col] = st[1:]
        else:
            ok[1:, col] = st[:-1]
    return ok


def state_bonus(sc, codes, w):
    """[L, 4] float32 bonus for the boundary states (L+, R+, L-, R-) at every base, and pi per column."""
    ok = site_allowed(codes)
    L = len(codes)
    out, pis = np.zeros((L, 4), np.float32), {}
    for k, ((part, s), (col, off)) in enumerate(BOUNDARY.items()):
        v = np.asarray(sc[:L, col], np.float32)
        pi = float(v[ok[:, col]].mean())
        pis[f"{part}{s}"] = pi
        phi = np.where(ok[:, col], w * (np.log(np.maximum(v, 1e-6)) - np.log(pi)), 0.0).astype(np.float32)
        if off < 0:   # state at t uses site t - 1
            out[1:, k] = phi[:-1]
        else:         # state at t uses site t + 1
            out[:-1, k] = phi[1:]
    return out, pis


def oracle_scores(truth, L):
    sc = np.full((L, 4), 1e-6, np.float32)
    for cls in ("mrna", "lnc"):
        for c, pos in sai_score.truth_sites(truth, cls).items():
            pos = pos[(pos >= 0) & (pos < L)]
            sc[pos, c] = 1.0
    return sc


def posteriors(run, ck, feats, codes, bonus, device):
    """p15_read.posteriors with `bonus` added to the boundary-state emissions."""
    import torch
    from ncrna_annot.decode import chrom
    from ncrna_annot.heads.c0 import C0BHead
    from ncrna_annot.train import screen_g1b as sg
    assert run["arm"] == "G1B" and not run["lnc_only"]
    L = len(feats) * sg.TOK
    trans = g1.G1Transitions().to(device)
    trans.load_state_dict(ck["trans"])
    head, coef = C0BHead(feats.shape[-1]).to(device), g1.BoundaryCoef().to(device)
    coef.load_state_dict(ck["coef"])
    head.load_state_dict(ck["head"])
    head.eval()
    cols = [[g1.IDX[f"{k}{part}{s}"] for k in g1.TYPES] for (part, s) in BOUNDARY]

    def emit_fn(s, e):
        f = torch.as_tensor(np.asarray(feats[s // sg.TOK:e // sg.TOK], np.float32), device=device)[None]
        c = torch.as_tensor(np.asarray(codes[s:e], np.int64), device=device)[None]
        with torch.no_grad():
            p, b = head(f, c)
            p = sg._mask_lnc_only(p, False)[0]
            if bonus is not None:
                bt = torch.as_tensor(bonus[s:e], device=device, dtype=p.dtype)
                p = p.clone()
                for k, idx in enumerate(cols):
                    p[:, idx] += bt[:, k:k + 1]
        return p, b[0]

    T = trans().detach()
    return chrom.chromosome_marginals(emit_fn, L, T, codes[:L], device, coef=coef, exact=True), L


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/p15_v19b.json")
    ap.add_argument("--run", required=True)
    ap.add_argument("--w", type=float, required=True)
    ap.add_argument("--sai", default="outputs/PHASE1/SPLICE_DIAG/sai_chr18.npy")
    ap.add_argument("--oracle", action="store_true")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    import torch
    torch.set_num_threads(int(os.environ.get("NCRNA_THREADS", "2")))
    torch.cuda.set_per_process_memory_fraction(float(os.environ.get("NCRNA_GPU_FRAC", "0.15")), 0)
    from ncrna_annot.decode import chrom
    from ncrna_annot.eval import ruler
    from ncrna_annot.train import screen2 as s2
    cfg = json.load(open(f"{ROOT}/{a.config}"))
    run, d = cfg["runs"][a.run], cfg["decode"]
    out = f"{ROOT}/{a.out}"
    if os.path.exists(f"{out}/result.json"):
        print(json.dumps({"skip": a.out}))
        return
    if R.gate_busy():
        sys.exit(75)
    os.makedirs(out, exist_ok=True)
    bad = [p for p, h in zip(run["ckpts"], run["ckpt_sha256"]) if R.sha256(f"{ROOT}/{p}") != h]
    assert not bad, bad
    feats, codes, tg = R.load_target(cfg, "chr18")
    truth = ruler.read_gtf(f"{ROOT}/{tg['truth']}", tg["truth_chrom"])
    L = len(feats) * 6
    t0 = time.time()
    meta = {"run": a.run, "w": a.w, "oracle": a.oracle, "config": a.config, "started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "code_sha256": {p: R.sha256(f"{ROOT}/{p}") for p in ("src/ncrna_annot/grammar/g1.py", "src/ncrna_annot/decode/chrom.py",
                                                                  "src/ncrna_annot/eval/ruler.py", "scripts/phase1/splice_diag/splice_diag.py")}}
    bonus = None
    if a.w != 0 or a.oracle:
        sc = oracle_scores(truth, L) if a.oracle else np.load(f"{ROOT}/{a.sai}", mmap_mode="r")
        bonus, meta["pi"] = state_bonus(sc, np.asarray(codes[:L]), a.w)
    ck = s2.load_ckpt([f"{ROOT}/{p}" for p in run["ckpts"]], "cuda")
    post, L = posteriors(run, ck, feats, codes, bonus, "cuda")
    meta["post_s"] = round(time.time() - t0)
    if a.w == 0 and not a.oracle:  # bit-for-bit check against the P1.5 chr18 read of the same run
        np.save(f"{out}/post.npy", post)
        ref = json.load(open(f"{ROOT}/{cfg['out_dir']}/chr18/{a.run}/read_meta.json"))["post_sha256"]
        meta["post_sha256"] = R.sha256(f"{out}/post.npy")
        meta["post_matches_p15"] = meta["post_sha256"] == ref
        os.remove(f"{out}/post.npy")
    t1 = time.time()
    Tg = g1.G1Transitions()().detach().double().numpy()
    genes = chrom.decode_blocks(post, codes[:L], Tg, d["min_exon"], d["min_intron"], mode="mea", gamma=d["gamma"],
                                level=d["level"], log_gain=False)
    genes = [g for g in genes if g["type"] != "l" or sum(b - a_ + 1 for a_, b in g["exons"]) >= d["lnc_min_len"]]
    for g in genes:
        g["conf"] = chrom.gene_confidence(g, post)
    meta["decode_s"] = round(time.time() - t1)
    json.dump([{"type": g["type"], "strand": g["strand"], "exons": [list(x) for x in g["exons"]], "conf": g["conf"]}
               for g in genes], open(f"{out}/genes.json", "w"))
    genes = P.load_genes(f"{out}/genes.json")
    fin, _ = P.score_preds(truth, P.final_preds(genes, d["tau"]), boot=0)
    res = {**meta, "core": P.core(genes, post, truth, d["tau"]), "final": fin, "total_s": round(time.time() - t0)}
    json.dump(P.norm(res), open(f"{out}/result.json", "w"), indent=1)
    sl, sm = fin["lnc"], fin["mrna"]
    print(json.dumps({"done": a.out, "post_matches_p15": meta.get("post_matches_p15"), "auprc_lnc": res["core"]["auprc_lnc"],
                      "lnc_L1": [sl["L1_locus"]["sn"], sl["L1_locus"]["pr_strict"]], "lnc_L2_sn": sl["L2_intron_chain"]["sn"],
                      "lnc_L3": [sl["L3_junction"]["recall"], sl["L3_junction"]["precision"]],
                      "mrna_L2_sn": sm["L2_intron_chain"]["sn"], "mrna_L3": [sm["L3_junction"]["recall"], sm["L3_junction"]["precision"]],
                      "total_s": res["total_s"]}), flush=True)


if __name__ == "__main__":
    main()
