#!/usr/bin/env python3
"""chr18 splice diagnostic, round 2 (user decision 2026-11-03 on card splice-next, option 1): our own splice-site CNN,
trained on GENCODE v50 mRNA AND lnc splice sites, replacing SpliceAI as the external potential of splice_diag.py.

Design (written before training; docs/phase1_results/PLAN_SPLICE_OWN.md):
  chromosomes  train chr1-8, 10-17, 19, 21, X;  chr22 = the splice net's own validation (epoch and variant choice);
               chr9 (v2.0 test), chr18 (development), chr20 (read) are never used here except chr18 in `score`.
  labels       every splice site of every transcript of the v50 genes of class mrna or lnc (ruler.read_gtf classes),
               in transcript orientation: 1 acceptor (first exon base), 2 donor (last exon base), 0 elsewhere;
               N positions are ignored.  "-" strand sites are learnt on the reverse complement (SpliceAI convention).
  model        SpliceAI-style residual dilated CNN, 32 channels, 3-way softmax; s2k = RB(11,1)x4 RB(11,4)x4
               RB(21,10)x4 (+-1,000 nt context), s10k = s2k + RB(41,25)x4 (+-5,000 nt).
  training     per epoch 64,000 blocks of 5,000 output nt: 40% centred (uniform offset) on a random lnc site, 40% on
               a random mRNA site (block in the site's transcript orientation), 20% uniform (random orientation);
               cross-entropy, AdamW lr 1e-3 cosine to 1e-5, weight decay 1e-4, batch 16, 10 epochs, bf16 autocast,
               seed 0.
  selection    lnc average precision on chr22 (mean over acceptor / donor x both strands; positives = lnc sites,
               candidates = motif-allowed positions that are not mRNA-only sites), best epoch per variant, then the
               better variant.  Only the selected checkpoint scores chr18.

  own_splice.py prep   [--out data/splice_own]
  own_splice.py train  --variant s2k|s10k [--out outputs/PHASE1/SPLICE_OWN]
  own_splice.py select [--dir outputs/PHASE1/SPLICE_OWN]                     -> <dir>/SELECT.json
  own_splice.py score  --codes data/p1/chr18_codes.npy --truth outputs/PHASE1/DIAG_L2/chr18_v50.gtf --chrom chr18 \
                       --out outputs/PHASE1/SPLICE_OWN/own_chr18.npy           (sai_score.py output format)
  own_splice.py verdict [--dir outputs/PHASE1/SPLICE_OWN/diag]              -> <dir>/VERDICT.json (pre-written rule)

train exits 75 when the GPU gate closes (TE / other users / memory) after saving <out>/<variant>/resume.pt; rerunning
continues from the same step with the same data (each epoch's blocks are drawn from a seed fixed by the epoch).
"""
import argparse
import glob
import gzip
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

HG38 = f"{ROOT}/data/human/hg38.fa.gz"
GTF = f"{ROOT}/data/human/gencode.v50.annotation.gtf.gz"
TRAIN = [f"chr{i}" for i in list(range(1, 9)) + list(range(10, 18)) + [19, 21]] + ["chrX"]
VAL = "chr22"
NEVER = ("chr9", "chr18", "chr20")
BLOCK, BATCH, EPOCHS, PER_EPOCH, SEED = 5000, 16, 10, 64_000, 0
MIX = (("lnc", 0.4), ("mrna", 0.4), ("rand", 0.2))
SPECS = {"s2k": [(11, 1)] * 4 + [(11, 4)] * 4 + [(21, 10)] * 4}
SPECS["s10k"] = SPECS["s2k"] + [(41, 25)] * 4
COLS = ("acc+", "don+", "acc-", "don-")
assert not set(NEVER) & set(TRAIN + [VAL])


# ---------------------------------------------------------------- data
def read_gtf_multi(path, chroms):
    """ruler.read_gtf for several chromosomes in one pass -> {chrom: {gene_id: gene}} (same fields and classes)."""
    from ncrna_annot.eval import ruler
    genes = {c: {} for c in chroms}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            c = line.split("\t", 1)[0]
            if c not in genes:
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] != "exon":
                continue
            a = ruler._attrs_gtf(f[8])
            gtype = a.get("gene_type") or a.get("gene_biotype", "")
            g = genes[c].setdefault(a["gene_id"], {"strand": f[6], "gene_type": gtype, "cls": ruler.gene_class(gtype), "tx": {}})
            g["tx"].setdefault(a["transcript_id"], []).append((int(f[3]), int(f[4])))
    for gs in genes.values():
        for g in gs.values():
            for t in g["tx"].values():
                t.sort()
            if g["gene_type"] in ruler.LEN_SPLIT:
                longest = max(sum(e - s + 1 for s, e in ruler._merge(t)) for t in g["tx"].values())
                g["cls"] = ruler.gene_class(g["gene_type"], longest)
    return genes


def prep(out):
    import prep_p1
    import sai_score
    from ncrna_annot.grammar import g1
    os.makedirs(f"{ROOT}/{out}", exist_ok=True)
    chroms = TRAIN + [VAL]
    t0 = time.time()
    seqs = prep_p1.read_genome(set(chroms))
    genes = read_gtf_multi(GTF, set(chroms))
    sites, stats = {}, {"chroms": {}, "gtf_sha256": None}
    for c in chroms:
        codes = g1.encode(seqs.pop(c))
        np.save(f"{ROOT}/{out}/codes_{c}.npy", codes)
        st = {"len": int(len(codes)), "n": int((codes == 4).sum())}
        for cls in ("mrna", "lnc"):
            for col, pos in sai_score.truth_sites(genes[c], cls).items():
                pos = pos[(pos >= 0) & (pos < len(codes))]
                sites[f"{c}|{cls}|{col}"] = pos
                st[f"{cls}_{COLS[col]}"] = int(len(pos))
        stats["chroms"][c] = st
        print(c, json.dumps(st), flush=True)
    np.savez(f"{ROOT}/{out}/sites.npz", **sites)
    import hashlib
    stats["gtf_sha256"] = hashlib.sha256(open(GTF, "rb").read()).hexdigest()
    stats["prep_s"] = round(time.time() - t0)
    json.dump(stats, open(f"{ROOT}/{out}/prep.json", "w"), indent=1)


class Data:
    def __init__(self, d, chroms):
        z = np.load(f"{ROOT}/{d}/sites.npz")
        self.codes = {c: np.load(f"{ROOT}/{d}/codes_{c}.npy", mmap_mode="r") for c in chroms}
        self.sites = {k: z[k] for k in z.files if k.split("|")[0] in chroms}
        self.chroms = chroms
        # label lookup per chromosome and orientation: sorted positions + label (1 acc, 2 don), mrna and lnc merged
        self.lab = {}
        for c in chroms:
            for o, (ca, cd) in (("+", (0, 1)), ("-", (2, 3))):
                pos = np.concatenate([self.sites[f"{c}|{cls}|{col}"] for cls in ("mrna", "lnc") for col in (ca, cd)])
                lab = np.concatenate([np.full(len(self.sites[f"{c}|{cls}|{col}"]), 1 if col == ca else 2, np.int8)
                                      for cls in ("mrna", "lnc") for col in (ca, cd)])
                order = np.argsort(pos, kind="stable")
                self.lab[(c, o)] = (pos[order], lab[order])
        # sampling pools: (chrom, position, orientation) of every lnc / mRNA site
        self.pool = {}
        for cls in ("mrna", "lnc"):
            cs, ps, os_ = [], [], []
            for c in chroms:
                for col in range(4):
                    p = self.sites[f"{c}|{cls}|{col}"]
                    cs.append(np.full(len(p), chroms.index(c), np.int16))
                    ps.append(p)
                    os_.append(np.full(len(p), col >= 2, bool))
            self.pool[cls] = (np.concatenate(cs), np.concatenate(ps), np.concatenate(os_))
        self.lens = np.array([len(self.codes[c]) for c in chroms], np.float64)

    def epoch_blocks(self, epoch):
        rng = np.random.default_rng(SEED * 1000 + epoch)
        out = []
        for kind, frac in MIX:
            n = int(round(PER_EPOCH * frac))
            if kind == "rand":
                ci = rng.choice(len(self.chroms), n, p=self.lens / self.lens.sum())
                st = (rng.random(n) * (self.lens[ci] - BLOCK)).astype(np.int64)
                neg = rng.random(n) < 0.5
            else:
                cs, ps, os_ = self.pool[kind]
                k = rng.integers(0, len(ps), n)
                ci, neg = cs[k].astype(np.int64), os_[k]
                st = ps[k] - rng.integers(0, BLOCK, n)
            out += list(zip(ci.tolist(), st.tolist(), neg.tolist()))
        return [out[i] for i in rng.permutation(len(out))]

    def batch(self, blocks, flank):
        """-> x float32 [B, 4, BLOCK + 2 flank] (block orientation), y int64 [B, BLOCK] (-1 = ignore)."""
        import sai_torch
        B = len(blocks)
        x = np.zeros((B, 4, BLOCK + 2 * flank), np.float32)
        y = np.full((B, BLOCK), -1, np.int64)
        for i, (ci, st, neg) in enumerate(blocks):
            c = self.chroms[ci]
            codes = self.codes[c]
            a, b = st - flank, st + BLOCK + flank
            seg = np.full(b - a, 4, np.int8)
            lo, hi = max(a, 0), min(b, len(codes))
            seg[lo - a:hi - a] = codes[lo:hi]
            lab = np.where(seg[flank:flank + BLOCK] == 4, -1, 0).astype(np.int64)
            pos, l = self.lab[(c, "-" if neg else "+")]
            j0, j1 = np.searchsorted(pos, [st, st + BLOCK])
            lab[pos[j0:j1] - st] = np.where(lab[pos[j0:j1] - st] >= 0, l[j0:j1], -1)
            xi = sai_torch.one_hot(seg)
            if neg:
                xi, lab = xi[::-1, ::-1], lab[::-1]
            x[i], y[i] = xi, lab
        return x, y


# ---------------------------------------------------------------- model
def build(variant):
    import torch
    nn = torch.nn

    class RB(nn.Module):
        def __init__(self, C, k, d):
            super().__init__()
            self.f = nn.Sequential(nn.BatchNorm1d(C), nn.ReLU(), nn.Conv1d(C, C, k, dilation=d, padding=d * (k - 1) // 2),
                                   nn.BatchNorm1d(C), nn.ReLU(), nn.Conv1d(C, C, k, dilation=d, padding=d * (k - 1) // 2))

        def forward(self, x):
            return x + self.f(x)

    class Net(nn.Module):
        def __init__(self, spec, C=32):
            super().__init__()
            self.flank = sum(2 * (k - 1) * d for k, d in spec) // 2
            self.inp, self.skip0 = nn.Conv1d(4, C, 1), nn.Conv1d(C, C, 1)
            groups = [spec[i:i + 4] for i in range(0, len(spec), 4)]
            self.groups = nn.ModuleList(nn.Sequential(*[RB(C, k, d) for k, d in g]) for g in groups)
            self.skips = nn.ModuleList(nn.Conv1d(C, C, 1) for _ in groups)
            self.out = nn.Conv1d(C, 3, 1)

        def forward(self, x):  # [B, 4, L + 2 flank] -> logits [B, 3, L]
            h = self.inp(x)
            s = self.skip0(h)
            for g, sk in zip(self.groups, self.skips):
                h = g(h)
                s = s + sk(h)
            f = self.flank
            return self.out(s[..., f:s.shape[-1] - f])

    return Net(SPECS[variant])


@np.errstate(invalid="ignore")
def score_chrom(net, codes, device="cuda", block=100_000):
    """[L, 4] float32 (acc+, don+, acc-, don-) at genomic positions, like sai_score.py."""
    import torch
    import sai_torch
    L, f = len(codes), net.flank
    out = np.zeros((L, 4), np.float32)
    net.eval()
    for k, strand in enumerate("+-"):
        x = sai_torch.one_hot(codes)
        if strand == "-":
            x = x[::-1, ::-1].copy()
        xp = np.zeros((4, L + 2 * f), np.float32)
        xp[:, f:f + L] = x
        o = np.zeros((L, 2), np.float32)
        with torch.no_grad():
            for s in range(0, L, block):
                e = min(s + block, L)
                xt = torch.as_tensor(xp[None, :, s:e + 2 * f], device=device)
                with torch.autocast("cuda", dtype=torch.bfloat16):
                    lg = net(xt)
                o[s:e] = torch.softmax(lg.float(), 1)[0, 1:].T.cpu().numpy()
        out[:, 2 * k:2 * k + 2] = o[::-1] if strand == "-" else o
    return out


def average_precision(scores, is_pos):
    order = np.argsort(-scores, kind="stable")
    hit = is_pos[order]
    n = hit.sum()
    if n == 0:
        return None
    prec = np.cumsum(hit) / np.arange(1, len(hit) + 1)
    return float(prec[hit].sum() / n)


def site_ap(sc, codes, sites_of):
    """AP per column for lnc and mRNA sites among motif-allowed positions; for each class the other class's sites
    are removed from the candidates (they are real sites, neither positive nor negative)."""
    import splice_diag
    ok = splice_diag.site_allowed(np.asarray(codes))
    res = {}
    for col in range(4):
        cand = np.flatnonzero(ok[:, col])
        for cls, other in (("lnc", "mrna"), ("mrna", "lnc")):
            pos, oth = sites_of(cls, col), np.setdiff1d(sites_of(other, col), sites_of(cls, col))
            keep = cand[~np.isin(cand, oth)]
            res[f"{cls}_{COLS[col]}"] = average_precision(sc[keep, col], np.isin(keep, pos))
            res[f"{cls}_{COLS[col]}_n"] = int(np.isin(pos, cand).sum())
    for cls in ("lnc", "mrna"):
        v = [res[f"{cls}_{c}"] for c in COLS if res[f"{cls}_{c}"] is not None]
        res[f"{cls}_ap"] = float(np.mean(v))
    return res


# ---------------------------------------------------------------- train / select / score
def gate_busy():
    import p15_read as R
    return R.gate_busy()


def train(variant, out, data_dir):
    import torch
    torch.backends.cudnn.benchmark = True
    torch.set_num_threads(int(os.environ.get("NCRNA_THREADS", "2")))
    torch.cuda.set_per_process_memory_fraction(float(os.environ.get("NCRNA_GPU_FRAC", "0.30")), 0)
    od = f"{ROOT}/{out}/{variant}"
    os.makedirs(od, exist_ok=True)
    if os.path.exists(f"{od}/DONE"):
        print(json.dumps({"skip": od}))
        return
    if gate_busy():
        sys.exit(75)
    data = Data(data_dir, TRAIN)
    val = Data(data_dir, [VAL])
    torch.manual_seed(SEED)
    net = build(variant).cuda()
    opt = torch.optim.AdamW(net.parameters(), lr=1e-3, weight_decay=1e-4)
    steps_per_epoch = PER_EPOCH // BATCH
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, EPOCHS * steps_per_epoch, eta_min=1e-5)
    ep0, st0, hist = 0, 0, []
    if os.path.exists(f"{od}/resume.pt"):
        ck = torch.load(f"{od}/resume.pt", map_location="cpu", weights_only=False)
        net.load_state_dict(ck["net"])
        opt.load_state_dict(ck["opt"])
        sched.load_state_dict(ck["sched"])
        ep0, st0, hist = ck["epoch"], ck["step"], ck["hist"]
        print(f"resume epoch {ep0} step {st0}", flush=True)
    loss_fn = torch.nn.CrossEntropyLoss(ignore_index=-1)
    meta = {"variant": variant, "flank": net.flank, "params": sum(p.numel() for p in net.parameters())}
    print(json.dumps(meta), flush=True)

    def save_resume(ep, st):
        torch.save({"net": net.state_dict(), "opt": opt.state_dict(), "sched": sched.state_dict(), "epoch": ep,
                    "step": st, "hist": hist}, f"{od}/resume.pt.tmp")
        os.replace(f"{od}/resume.pt.tmp", f"{od}/resume.pt")

    for ep in range(ep0, EPOCHS):
        blocks = data.epoch_blocks(ep)
        net.train()
        t0, run = time.time(), []
        for st in range(st0 if ep == ep0 else 0, steps_per_epoch):
            if st % 200 == 0 and st and gate_busy():
                save_resume(ep, st)
                print(f"gate closed at epoch {ep} step {st}", flush=True)
                sys.exit(75)
            x, y = data.batch(blocks[st * BATCH:(st + 1) * BATCH], net.flank)
            x, y = torch.as_tensor(x, device="cuda"), torch.as_tensor(y, device="cuda")
            with torch.autocast("cuda", dtype=torch.bfloat16):
                lg = net(x)
            loss = loss_fn(lg.float(), y)
            opt.zero_grad(set_to_none=True)
            loss.backward()
            opt.step()
            sched.step()
            run.append(loss.item())
            if st % 500 == 0:
                print(f"ep {ep} step {st} loss {np.mean(run[-500:]):.5f} {time.time() - t0:.0f}s", flush=True)
        t1 = time.time()
        sc = score_chrom(net, val.codes[VAL])
        ap = site_ap(sc, val.codes[VAL], lambda cls, col: val.sites[f"{VAL}|{cls}|{col}"])
        hist.append({"epoch": ep, "train_loss": float(np.mean(run)) if run else None, "train_s": round(t1 - t0),
                     "val_s": round(time.time() - t1), **ap})
        torch.save(net.state_dict(), f"{od}/ep{ep:02d}.pt")
        save_resume(ep + 1, 0)
        json.dump({**meta, "hist": hist}, open(f"{od}/train.json", "w"), indent=1)
        print(json.dumps(hist[-1]), flush=True)
    open(f"{od}/DONE", "w").write(time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))


def select(d):
    rows = {}
    for f in sorted(glob.glob(f"{ROOT}/{d}/*/train.json")):
        t = json.load(open(f))
        assert len(t["hist"]) == EPOCHS, f
        best = max(t["hist"], key=lambda h: h["lnc_ap"])
        rows[t["variant"]] = {"epoch": best["epoch"], "lnc_ap": best["lnc_ap"], "mrna_ap": best["mrna_ap"],
                              "ckpt": f"{d}/{t['variant']}/ep{best['epoch']:02d}.pt"}
    assert set(rows) == set(SPECS), rows
    pick = max(rows, key=lambda v: rows[v]["lnc_ap"])
    import p15_read as R
    res = {"rule": "best chr22 lnc AP epoch per variant, then the better variant", "variants": rows, "selected": pick,
           "ckpt": rows[pick]["ckpt"], "ckpt_sha256": R.sha256(f"{ROOT}/{rows[pick]['ckpt']}")}
    json.dump(res, open(f"{ROOT}/{d}/SELECT.json", "w"), indent=1)
    print(json.dumps(res, indent=1))


def score(a):
    import torch
    import sai_score
    from ncrna_annot.eval import ruler
    torch.cuda.set_per_process_memory_fraction(float(os.environ.get("NCRNA_GPU_FRAC", "0.30")), 0)
    sel = json.load(open(f"{ROOT}/{a.dir}/SELECT.json"))
    import p15_read as R
    assert R.sha256(f"{ROOT}/{sel['ckpt']}") == sel["ckpt_sha256"]
    if gate_busy():
        sys.exit(75)
    net = build(sel["selected"]).cuda()
    net.load_state_dict(torch.load(f"{ROOT}/{sel['ckpt']}", map_location="cuda"))
    codes = np.load(f"{ROOT}/{a.codes}", mmap_mode="r")
    t0 = time.time()
    sc = score_chrom(net, np.asarray(codes))
    meta = {"select": sel, "score_s": round(time.time() - t0), "started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    out = f"{ROOT}/{a.out}"
    np.save(out, sc.astype(np.float16))
    truth = ruler.read_gtf(f"{ROOT}/{a.truth}", a.chrom)
    sc16 = np.load(out, mmap_mode="r")
    meta["sites"] = {cls: sai_score.site_stats(sc16, sai_score.truth_sites(truth, cls)) for cls in ("mrna", "lnc")}
    ts = {cls: sai_score.truth_sites(truth, cls) for cls in ("mrna", "lnc")}
    meta["ap"] = site_ap(np.asarray(sc16, np.float32), np.asarray(codes), lambda cls, col: ts[cls][col])
    json.dump(meta, open(out + ".json", "w"), indent=1)
    print(json.dumps({k: meta[k] for k in ("score_s", "ap")}), flush=True)


def verdict(d):
    """The rule written before the first diagnostic (2026-11-03, node L3.splicediag) applied to <d>/SUMMARY.json."""
    s = json.load(open(f"{ROOT}/{d}/SUMMARY.json"))
    base = s["mean"]["0"]
    runs = ("g1b-42-ke", "g1b-123-ke", "g1b-456-ke")
    out = {}
    for k, m in s["mean"].items():
        if k in ("0", "oracle"):
            continue
        rise = {met: all(s["rows"][f"{r}|{k}"][met] > s["rows"][f"{r}|0"][met] for r in runs) for met in ("lnc_L3_sn", "lnc_L2_sn")}
        struct = (m["lnc_L3_sn"] >= 0.10 and rise["lnc_L3_sn"]) or (m["lnc_L2_sn"] >= 0.05 and rise["lnc_L2_sn"])
        guard = (base["lnc_L1_sn"] - m["lnc_L1_sn"] <= 0.03) and m["lnc_L1_pr"] >= 0.40
        out[k] = {"lnc_L3_sn": m["lnc_L3_sn"], "lnc_L2_sn": m["lnc_L2_sn"], "lnc_L1_sn": m["lnc_L1_sn"],
                  "lnc_L1_pr": m["lnc_L1_pr"], "all_seeds_rise": rise, "structure": struct, "guard": guard,
                  "effective": bool(struct and guard)}
    res = {"baseline": base, "w": out, "effective": any(v["effective"] for v in out.values())}
    json.dump(res, open(f"{ROOT}/{d}/VERDICT.json", "w"), indent=1)
    print(json.dumps(res, indent=1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("prep", "train", "select", "score", "verdict"))
    ap.add_argument("--variant", choices=sorted(SPECS))
    ap.add_argument("--data", default="data/splice_own")
    ap.add_argument("--out")
    ap.add_argument("--dir", default="outputs/PHASE1/SPLICE_OWN")
    ap.add_argument("--codes")
    ap.add_argument("--truth")
    ap.add_argument("--chrom")
    a = ap.parse_args()
    if a.cmd == "prep":
        prep(a.out or a.data)
    elif a.cmd == "train":
        train(a.variant, a.out or a.dir, a.data)
    elif a.cmd == "select":
        select(a.dir)
    elif a.cmd == "score":
        score(a)
    else:
        verdict(a.dir)


if __name__ == "__main__":
    main()
