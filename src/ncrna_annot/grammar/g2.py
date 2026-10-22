"""G2 grammar (protocol v1.6 section 13): sparse-event semi-Markov model (GHMM) on candidate boundaries.

Labels (Y = 21): IG, and per strand s and type k the exon classes by GENOMIC position
  S single exon, A leftmost exon of a multi-exon gene, M middle exon, Z rightmost exon, and N intron.
  Transcript roles (mirror-tied parameters): + strand A = first, Z = last; - strand A = last, Z = first.
Paths, left to right:  IG -> ( kSs | kAs (kNs kMs)* kNs kZs ) -> IG ...   (a gene never changes k or s)

Segments are half-open [a, b) between NODES.  Nodes come from candidate events (decode/g2cand.py channels):
  gene-left base p of strand s  -> node p    (GL_s: may start S / A, may end IG)
  gene-right base p             -> node p+1  (GR_s: may end S / Z, may start IG)
  intron first base p (L)       -> node p    (IL_s: may end A / M, may start N)
  intron last base p (R)        -> node p+1  (IR_s: may end N, may start M / Z)
Node 0 and node L are open window ends: any label may start at 0 / end at L ("open" sides get no event,
transition, length or end score; minimum lengths are not enforced there).

Segment score (all levels): sum of the per-base G1 potentials of the segment's state (exon: kEs; intron: kLs at
the first base, kIs inside, kRs at the last base when those sides are closed; IG) + a per-base slope per
role (linear length = geometric, G2-1) + boundary-head bonuses at closed gene ends (coef x TSS / TES logit, as
in G1B) + transition scores, each attached to the closed side of the adjacent exon (IG->S, S->IG, IG->first,
last->IG, donor, acceptor; mirror-tied).  G2-2 adds bucketed log2-length scores (exon roles, intron) on segments
closed on both sides; G2-3 adds an exon segment MLP on prefix-pooled C0 hidden features (segment mean and
50-bp flanks).  IG has no length table (its lengths exceed windows and the look-back; logged in the protocol).

Inference: exact forward / backward over nodes in float64 (no autograd through the loop); log Z gradients
via the surrogate sum(mu * segment score) with detached posteriors mu (grad log Z = E[grad score]).
Non-IG segments look back over at most K valid start nodes of their group; IG uses an unbounded accumulator.
"""
import math

import numpy as np
import torch
from torch import nn

from . import g1

NEGM = -1e9
STRANDS, TYPES = g1.STRANDS, g1.TYPES
CLASSES = ("S", "A", "M", "Z", "N")
LABELS = ["IG"] + [f"{k}{c}{s}" for s in STRANDS for k in TYPES for c in CLASSES]
LIDX = {n: i for i, n in enumerate(LABELS)}
Y = len(LABELS)
J = Y - 1
ROLES = ("single", "first", "internal", "last", "intron")
NB = 20
KINDS = [(g, s) for s in STRANDS for g in ("GL", "GR", "IL", "IR")]
KIDX = {x: i for i, x in enumerate(KINDS)}
GROUPS = [(s, g) for s in STRANDS for g in ("GL", "IR", "IL")]
START_GROUP = {"S": "GL", "A": "GL", "M": "IR", "Z": "IR", "N": "IL"}
END_KIND = {"S": "GR", "Z": "GR", "A": "IL", "M": "IL", "N": "IR"}
CAND_NODE = {0: (0, "IL", "+"), 1: (1, "IR", "+"), 2: (0, "IL", "-"), 3: (1, "IR", "-"),
             4: (0, "GL", "+"), 5: (1, "GR", "+"), 6: (0, "GL", "-"), 7: (1, "GR", "-")}
EXON_LABELS = [n for n in LABELS[1:] if n[1] != "N"]          # 16, MLP output order
MAXLEN = {"exon": 20_000, "intron": 500_000}
TRANS = ("s_in", "s_out", "e5_in", "e3_out", "donor", "acceptor")


def role(c, s):
    return {"S": "single", "M": "internal", "N": "intron", "A": "first" if s == "+" else "last",
            "Z": "last" if s == "+" else "first"}[c]


def trans_names(c, s):
    """(transition param on the closed LEFT side, on the closed RIGHT side) of a segment of class c, strand s."""
    p = s == "+"
    return {"S": ("s_in" if p else "s_out", "s_out" if p else "s_in"),
            "A": ("e5_in" if p else "e3_out", "donor" if p else "acceptor"),
            "M": ("acceptor" if p else "donor", "donor" if p else "acceptor"),
            "Z": ("acceptor" if p else "donor", "e3_out" if p else "e5_in"),
            "N": (None, None)}[c]


def end_channels(s):
    """(boundary channel at a gene's left end, coef index), (right end, coef index); g1.BOUNDARY order, coef 0 = start."""
    return ((0, 0), (1, 1)) if s == "+" else ((3, 1), (2, 0))


def _transition_structure():
    T = np.full((Y, Y), NEGM)
    for s in STRANDS:
        for k in TYPES:
            i = lambda c: LIDX[f"{k}{c}{s}"]
            for a, b in ((0, i("S")), (0, i("A")), (i("S"), 0), (i("Z"), 0), (i("A"), i("N")), (i("M"), i("N")),
                         (i("N"), i("M")), (i("N"), i("Z"))):
                T[a, b] = 0.0
    return T


TS = _transition_structure()


class G2Params(nn.Module):
    """Grammar-level parameters of G2 (the per-base potentials and boundary logits come from a C0BHead)."""

    def __init__(self, level=3, hidden=256, mlp=16):
        super().__init__()
        self.level = level
        self.trans = nn.Parameter(torch.zeros(len(TYPES), len(TRANS)))
        self.slope = nn.Parameter(torch.zeros(len(TYPES), len(ROLES)))
        self.slope_ig = nn.Parameter(torch.zeros(()))
        self.coef = nn.Parameter(torch.ones(len(TYPES), 2))
        if level >= 2:
            self.psi = nn.Parameter(torch.zeros(len(TYPES), len(ROLES), NB))
        if level >= 3:
            self.mlp = mlp
            self.proj = nn.Linear(hidden, 3 * mlp, bias=False)
            self.u_len = nn.Parameter(torch.zeros(mlp))
            self.b1 = nn.Parameter(torch.zeros(mlp))
            self.out = nn.Linear(mlp, len(EXON_LABELS))
            nn.init.zeros_(self.out.weight)
            nn.init.zeros_(self.out.bias)


# ---- lattice -----------------------------------------------------------------------------------------------

def build_lattice(L, cands, K=256):
    """cands: 8 arrays of base positions (g2cand channel order).  -> dict(pos [N], kind [N, 8], SI [N, K, 6]).
    SI[n, :, g] = the (up to) K most recent valid start nodes < n of group g, right-aligned, -1 = none."""
    pts, flags = {0, L}, {}
    for ch, arr in enumerate(cands):
        off, kd, s = CAND_NODE[ch]
        for p in np.asarray(arr, np.int64).ravel():
            x = int(p) + off
            if 0 < x < L:
                pts.add(x)
                flags.setdefault(x, set()).add((kd, s))
    pos = np.array(sorted(pts), np.int64)
    N = len(pos)
    kind = np.zeros((N, len(KINDS)), bool)
    for i, x in enumerate(pos):
        for f in flags.get(int(x), ()):
            kind[i, KIDX[f]] = True
    SI = np.full((N, K, len(GROUPS)), -1, np.int64)
    for gi, (s, g) in enumerate(GROUPS):
        ok = kind[:, KIDX[(g, s)]].copy()
        ok[0] = True
        idx = np.flatnonzero(ok)
        cnt = np.searchsorted(idx, np.arange(N))
        for n in range(1, N):
            sel = idx[max(0, cnt[n] - K):cnt[n]]
            SI[n, K - len(sel):, gi] = sel
    return {"pos": pos, "kind": kind, "SI": SI, "L": L}


def collate(lats, device="cpu"):
    """pad a list of lattices (same L) to a batch of torch tensors."""
    B, N = len(lats), max(len(l["pos"]) for l in lats)
    K = lats[0]["SI"].shape[1]
    pos = np.zeros((B, N), np.int64)
    kind = np.zeros((B, N, len(KINDS)), bool)
    SI = np.full((B, N, K, len(GROUPS)), -1, np.int64)
    for b, l in enumerate(lats):
        n = len(l["pos"])
        pos[b, :n], pos[b, n:] = l["pos"], l["L"]
        kind[b, :n], SI[b, :n] = l["kind"], l["SI"]
    t = lambda a: torch.as_tensor(a, device=device)
    return {"pos": t(pos), "kind": t(kind), "SI": t(SI), "nlast": t(np.array([len(l["pos"]) - 1 for l in lats])),
            "L": lats[0]["L"]}


# ---- segment scores ------------------------------------------------------------------------------------------

def _at(T, p):
    """T [B, X, C], p [B, ...] long -> [B, ..., C]"""
    B = T.shape[0]
    flat = p.reshape(B, -1)
    return T.gather(1, flat[..., None].expand(-1, -1, T.shape[-1])).reshape(*p.shape, T.shape[-1])


def _bucket(ln):
    return torch.clamp(torch.floor(torch.log2(ln.clamp(min=1).to(torch.float64))), 0, NB - 1).long()


def segment_scores(x, gates, bnd, lat, prm, h=None, allowed=None, min_len=None):
    """x [B, L, 17] per-base potentials (ungated), gates [B, L, 17] (0 / g1.NEG), bnd [B, L, 4] boundary logits,
    h [B, L, H] C0 hidden (level 3), allowed bool [B, L, 17] (observed partition) or None,
    min_len = (min_exon, min_intron) for decoding (closed-closed segments only) or None.
    -> SEG [B, N, K, J] float64 (NEGM = invalid), IG dict(Q [B, N], end, start, reset [B, N] bool)."""
    dt = torch.float64
    B, L, _ = x.shape
    pos, SI, kind = lat["pos"], lat["SI"], lat["kind"]
    N, K = pos.shape[1], SI.shape[2]
    dev = x.device
    xd = x.to(dt)
    P = torch.cat([xd.new_zeros(B, 1, xd.shape[-1]), xd.cumsum(1)], 1)
    pa, pb = pos.clamp(max=L - 1), (pos - 1).clamp(min=0)
    Pn, xa, xb = _at(P, pos), _at(xd, pa), _at(xd, pb)
    gok_a, gok_b = _at(gates, pa) > g1.NEG / 2, _at(gates, pb) > g1.NEG / 2
    bd = bnd.to(dt)
    ba, bb = _at(bd, pa), _at(bd, pb)
    node = torch.arange(N, device=dev)[None].expand(B, -1)
    open_e = pos == L
    obs = allowed is not None
    if obs:
        al = allowed.bool()
        D = torch.cat([torch.zeros(B, 1, al.shape[-1], dtype=torch.long, device=dev), (~al).long().cumsum(1)], 1)
        anyI = torch.stack([al[..., g1.IDX[f"{k}L{s}"]] | al[..., g1.IDX[f"{k}I{s}"]] | al[..., g1.IDX[f"{k}R{s}"]]
                            for s in STRANDS for k in TYPES], -1)
        D2 = torch.cat([torch.zeros(B, 1, 4, dtype=torch.long, device=dev), (~anyI).long().cumsum(1)], 1)
        Dn, Dp1, Dm1 = _at(D, pos), _at(D, (pos + 1).clamp(max=L)), _at(D, (pos - 1).clamp(min=0))
        D2n, D2p1, D2m1 = _at(D2, pos), _at(D2, (pos + 1).clamp(max=L)), _at(D2, (pos - 1).clamp(min=0))
        ala, alb = _at(al.long(), pa) > 0, _at(al.long(), pb) > 0
    mlp = h is not None and prm.level >= 3
    if mlp:
        U = (h.to(prm.proj.weight.dtype) @ prm.proj.weight.T).to(dt)
        U = torch.cat([U.new_zeros(B, 1, U.shape[-1]), U.cumsum(1)], 1)
        Un, Um50, Up50 = _at(U, pos), _at(U, (pos - 50).clamp(min=0)), _at(U, (pos + 50).clamp(max=L))
    SEG = torch.full((B, N, K, J), NEGM, dtype=dt, device=dev)
    posf = pos.to(dt)
    for gi, (s, g) in enumerate(GROUPS):
        m = SI[..., gi]
        ok0 = m >= 0
        mc = m.clamp(min=0)
        a = pos.gather(1, mc.reshape(B, -1)).reshape(B, N, K)
        b = pos[:, :, None].expand(-1, -1, K)
        ln = b - a
        lnf = ln.to(dt)
        cl = mc != 0                                     # closed left side (start is not node 0)
        cr = (~open_e)[:, :, None].expand(-1, -1, K)     # closed right side
        both = cl & cr
        for name in LABELS[1:]:
            k, c, s2 = name[0], name[1], name[2]
            if s2 != s or START_GROUP[c] != g:
                continue
            j = LIDX[name] - 1
            ki = TYPES.index(k)
            r = ROLES.index(role(c, s))
            end_ok = open_e | kind[..., KIDX[(END_KIND[c], s)]]
            valid = ok0 & end_ok[:, :, None] & (ln >= 1)
            if c == "N":
                Ic, Lc, Rc = (g1.IDX[f"{k}{p}{s}"] for p in ("I", "L", "R"))
                sc = Pn[..., Ic][:, :, None] - _at(Pn[..., Ic:Ic + 1], mc)[..., 0]
                sc = sc + torch.where(cl, _at((xa[..., Lc] - xa[..., Ic])[..., None], mc)[..., 0], 0.0)
                sc = sc + torch.where(cr, (xb[..., Rc] - xb[..., Ic])[:, :, None], 0.0)
                valid &= (ln >= 2) | ~both
                valid &= ~cl | _at(gok_a[..., Lc:Lc + 1].long(), mc)[..., 0].bool()
                valid &= ~cr | gok_b[..., Rc][:, :, None]
                valid &= ln <= MAXLEN["intron"]
            else:
                Ec = g1.IDX[f"{k}E{s}"]
                sc = Pn[..., Ec][:, :, None] - _at(Pn[..., Ec:Ec + 1], mc)[..., 0]
                valid &= ln <= MAXLEN["exon"]
            sc = sc + prm.slope[ki, r].to(dt) * lnf
            # gene-end bonuses (closed sides)
            (chL, cfL), (chR, cfR) = end_channels(s)
            if c in ("S", "A"):
                sc = sc + torch.where(cl, prm.coef[ki, cfL].to(dt) * _at(ba[..., chL:chL + 1], mc)[..., 0], 0.0)
            if c in ("S", "Z"):
                sc = sc + torch.where(cr, prm.coef[ki, cfR].to(dt) * bb[..., chR][:, :, None], 0.0)
            tl, tr = trans_names(c, s)
            if tl:
                sc = sc + torch.where(cl, prm.trans[ki, TRANS.index(tl)].to(dt), 0.0)
            if tr:
                sc = sc + torch.where(cr, prm.trans[ki, TRANS.index(tr)].to(dt), 0.0)
            if prm.level >= 2:
                sc = sc + torch.where(both, prm.psi[ki, r].to(dt)[_bucket(ln)], 0.0)
            if mlp and c != "N":
                M = prm.mlp
                body = (Un[..., :M][:, :, None] - _at(Un[..., :M], mc)) / lnf.clamp(min=1)[..., None]
                fl_n = (a - (a - 50).clamp(min=0)).clamp(min=1).to(dt)
                left = (_at(Un[..., M:2 * M] - Um50[..., M:2 * M], mc)) / fl_n[..., None]
                fr_n = ((b + 50).clamp(max=L) - b).clamp(min=1).to(dt)
                right = (Up50[..., 2 * M:] - Un[..., 2 * M:])[:, :, None] / fr_n[..., None]
                z = body + left + right + prm.u_len.to(dt) * torch.log(lnf.clamp(min=1))[..., None] + prm.b1.to(dt)
                col = EXON_LABELS.index(name)
                wd = prm.out.weight.dtype
                sc = sc + (nn.functional.gelu(z.to(wd)) @ prm.out.weight[col] + prm.out.bias[col]).to(dt)
            if obs:
                if c == "N":
                    valid &= ~cl | _at(ala[..., Lc:Lc + 1].long(), mc)[..., 0].bool()
                    valid &= ~cr | alb[..., Rc][:, :, None]
                    lo_both = _at(Dp1[..., Ic:Ic + 1], mc)[..., 0]
                    hi_both = Dm1[..., Ic][:, :, None]
                    ksi = STRANDS.index(s) * 2 + ki
                    lo_any = torch.where(cl, _at(D2p1[..., ksi:ksi + 1], mc)[..., 0], _at(D2n[..., ksi:ksi + 1], mc)[..., 0])
                    hi_any = torch.where(cr, D2m1[..., ksi][:, :, None], D2n[..., ksi][:, :, None])
                    valid &= torch.where(both, (hi_both - lo_both) <= 0, (hi_any - lo_any) <= 0)
                else:
                    valid &= (Dn[..., Ec][:, :, None] - _at(Dn[..., Ec:Ec + 1], mc)[..., 0]) == 0
            if min_len is not None:
                mn = min_len[1] if c == "N" else min_len[0]
                valid &= ~both | (ln >= mn)
            SEG[..., j] = torch.where(valid, sc, torch.full_like(sc, NEGM))
    IGc = g1.IDX["IG"]
    Q = Pn[..., IGc] + prm.slope_ig.to(dt) * posf
    ig_end = open_e | kind[..., KIDX[("GL", "+")]] | kind[..., KIDX[("GL", "-")]]
    ig_start = (node == 0) | kind[..., KIDX[("GR", "+")]] | kind[..., KIDX[("GR", "-")]]
    if obs:
        prev = torch.cat([pos[:, :1], pos[:, :-1]], 1)
        ig_reset = (Dn[..., IGc] - _at(D[..., IGc:IGc + 1], prev)[..., 0]) > 0
    else:
        ig_reset = torch.zeros_like(ig_end)
    return SEG, {"Q": Q, "end": ig_end, "start": ig_start, "reset": ig_reset}


# ---- dynamic programming -------------------------------------------------------------------------------------

def _lab_si(SI):
    """[B, N, K, 6] group start indices -> [B, N, K, J] per label."""
    gidx = torch.tensor([GROUPS.index((n[2], START_GROUP[n[1]])) for n in LABELS[1:]], device=SI.device)
    return SI[..., gidx]


def _gather_beta(beta, m):
    """beta [B, N, Y], m [B, K, J] (or [B, N, K, J]) start nodes -> beta[b, m, label] with NEGM where m < 0."""
    B, N, _ = beta.shape
    lab = torch.arange(1, Y, device=beta.device)
    idx = m.clamp(min=0) * Y + lab
    v = beta.reshape(B, N * Y).gather(1, idx.reshape(B, -1)).reshape(m.shape)
    return torch.where(m >= 0, v, torch.full_like(v, NEGM))


def forward(SEG, IG, lat, mode="sum"):
    """-> dict(alpha, beta, logZ [B]) (+ back-pointers in max mode).  No autograd."""
    with torch.no_grad():
        B, N, K, _ = SEG.shape
        dev, dt = SEG.device, SEG.dtype
        TSt = torch.as_tensor(TS, dtype=dt, device=dev)
        SIl = _lab_si(lat["SI"])
        Q, ig_end, ig_start, ig_reset = IG["Q"], IG["end"], IG["start"], IG["reset"]
        alpha = torch.full((B, N, Y), NEGM, dtype=dt, device=dev)
        beta = torch.full((B, N, Y), NEGM, dtype=dt, device=dev)
        beta[:, 0] = 0.0
        acc = beta[:, 0, 0] - Q[:, 0]
        mx = mode == "max"
        if mx:
            bpK = torch.zeros((B, N, J), dtype=torch.long, device=dev)
            bpIG = torch.full((B, N), -1, dtype=torch.long, device=dev)
            bpPrev = torch.zeros((B, N, Y), dtype=torch.long, device=dev)
            accarg = torch.zeros(B, dtype=torch.long, device=dev)
        neg = torch.full((B,), NEGM, dtype=dt, device=dev)
        for n in range(1, N):
            tot = _gather_beta(beta, SIl[:, n]) + SEG[:, n]
            if mx:
                a_n, arg = tot.max(1)
                bpK[:, n] = arg
            else:
                a_n = torch.logsumexp(tot, 1)
            acc = torch.where(ig_reset[:, n], neg, acc)
            if mx:
                accarg = torch.where(ig_reset[:, n], torch.full_like(accarg, -1), accarg)
                bpIG[:, n] = accarg
            a_ig = torch.where(ig_end[:, n], acc + Q[:, n], neg)
            alpha[:, n, 0] = torch.maximum(a_ig, neg)
            alpha[:, n, 1:] = torch.maximum(a_n, torch.full_like(a_n, NEGM))
            t = alpha[:, n, :, None] + TSt
            if mx:
                b_n, prev = t.max(1)
                bpPrev[:, n] = prev
            else:
                b_n = torch.logsumexp(t, 1)
            b_n = torch.maximum(b_n, torch.full_like(b_n, NEGM))
            b_n[:, 0] = torch.where(ig_start[:, n], b_n[:, 0], neg)
            beta[:, n] = b_n
            cand = b_n[:, 0] - Q[:, n]
            if mx:
                better = cand > acc
                accarg = torch.where(better, torch.full_like(accarg, n), accarg)
                acc = torch.where(better, cand, acc)
            else:
                acc = torch.logaddexp(acc, cand)
        last = alpha[torch.arange(B, device=dev), lat["nlast"]]
        out = {"alpha": alpha, "beta": beta, "SIl": SIl}
        if mx:
            out.update(score=last.max(-1).values, ylast=last.argmax(-1), bpK=bpK, bpIG=bpIG, bpPrev=bpPrev)
        else:
            out["logZ"] = torch.logsumexp(last, -1)
        return out


def backward(SEG, IG, lat, fw):
    """-> dict(gamma, delta) given forward sums.  gamma[n, y]: completions after a y segment ends at node n."""
    with torch.no_grad():
        B, N, K, _ = SEG.shape
        dev, dt = SEG.device, SEG.dtype
        TSt = torch.as_tensor(TS, dtype=dt, device=dev)
        SIl, beta, logZ = fw["SIl"], fw["beta"], fw["logZ"]
        Q, ig_end, ig_reset = IG["Q"], IG["end"], IG["reset"]
        nlast = lat["nlast"]
        gamma = torch.full((B, N, Y), NEGM, dtype=dt, device=dev)
        delta = torch.full((B, N, Y), NEGM, dtype=dt, device=dev)
        Dt = torch.zeros((B, N * Y), dtype=dt, device=dev)
        lab = torch.arange(1, Y, device=dev)
        accR = torch.full((B,), NEGM, dtype=dt, device=dev)
        neg = torch.full((B,), NEGM, dtype=dt, device=dev)
        for n in range(N - 1, -1, -1):
            dl = Dt[:, n * Y:(n + 1) * Y]
            d = torch.where((dl > 0) & (beta[:, n] > NEGM / 2), torch.log(dl.clamp(min=1e-300)) - beta[:, n] + logZ[:, None],
                            torch.full_like(dl, NEGM))
            d[:, 0] = torch.where(IG["start"][:, n], torch.maximum(accR - Q[:, n], neg), neg)
            delta[:, n] = d
            if n == 0:
                break
            gm = torch.logsumexp(d[:, None, :] + TSt, -1)
            gm = torch.where((nlast == n)[:, None], torch.zeros_like(gm),
                             torch.where((n < nlast)[:, None], torch.maximum(gm, torch.full_like(gm, NEGM)), torch.full_like(gm, NEGM)))
            gamma[:, n] = gm
            m = SIl[:, n]
            mu = torch.exp(_gather_beta(beta, m) + SEG[:, n] + gm[:, None, 1:] - logZ[:, None, None])
            mu = torch.where(m >= 0, mu, torch.zeros_like(mu))
            Dt.scatter_add_(1, (m.clamp(min=0) * Y + lab).reshape(B, -1), mu.reshape(B, -1))
            accR = torch.where(ig_end[:, n], torch.logaddexp(accR, gm[:, 0] + Q[:, n]), accR)
            accR = torch.where(ig_reset[:, n], neg, accR)
        return {"gamma": gamma, "delta": delta}


def posteriors(SEG, IG, lat, fw=None, bw=None):
    """segment posteriors mu [B, N, K, J], expected IG segment ends / starts per node [B, N]; fw, bw."""
    fw = fw or forward(SEG, IG, lat)
    bw = bw or backward(SEG, IG, lat, fw)
    with torch.no_grad():
        logZ = fw["logZ"][:, None, None, None]
        mu = torch.exp(_gather_beta(fw["beta"], fw["SIl"]) + SEG + bw["gamma"][:, :, None, 1:] - logZ)
        mu = torch.where(fw["SIl"] >= 0, mu, torch.zeros_like(mu))
        e_end = torch.exp(fw["alpha"][..., 0] + bw["gamma"][..., 0] - fw["logZ"][:, None])
        e_start = torch.exp(fw["beta"][..., 0] + bw["delta"][..., 0] - fw["logZ"][:, None])
    return mu, e_end, e_start, fw, bw


def log_partition(SEG, IG, lat):
    """-> (logZ [B] float64 without grad, surrogate [B] whose gradient equals grad log Z)."""
    mu, e_end, e_start, fw, _ = posteriors(SEG, IG, lat)
    seg = torch.where(mu > 0, SEG, torch.zeros_like(SEG))
    sur = (mu * seg).sum((1, 2, 3)) + (e_end * IG["Q"]).sum(1) - (e_start * IG["Q"]).sum(1)
    return fw["logZ"], sur


def nll(SEG, IG, SEG_obs, IG_obs, lat):
    """partial-observation CRF loss per window: value = log Z - log Z_obs, gradient through the surrogates.
    Returns (loss [B] with grad, feasible [B] bool); infeasible windows (log Z_obs at the mask floor) get 0."""
    z, sz = log_partition(SEG, IG, lat)
    zo, so = log_partition(SEG_obs, IG_obs, lat)
    feas = zo > NEGM / 2
    val = (z - zo).detach()
    loss = (sz - so) + (val - (sz - so).detach())
    return torch.where(feas, loss, torch.zeros_like(loss)), feas


# ---- decoding ------------------------------------------------------------------------------------------------

def viterbi(SEG, IG, lat):
    """MAP segmentation per window -> list of [(a, b, label name)] in genomic order, and scores."""
    fw = forward(SEG, IG, lat, mode="max")
    out = []
    pos = lat["pos"].cpu().numpy()
    SIl = fw["SIl"].cpu().numpy()
    bpK, bpIG, bpPrev = fw["bpK"].cpu().numpy(), fw["bpIG"].cpu().numpy(), fw["bpPrev"].cpu().numpy()
    for b in range(SEG.shape[0]):
        n, y = int(lat["nlast"][b]), int(fw["ylast"][b])
        segs = []
        while n > 0:
            m = int(bpIG[b, n]) if y == 0 else int(SIl[b, n, bpK[b, n, y - 1], y - 1])
            segs.append((int(pos[b, m]), int(pos[b, n]), LABELS[y]))
            if m == 0:
                break
            y, n = int(bpPrev[b, m, y]), m
        out.append(segs[::-1])
    return out, fw["score"]


def segments_to_genes(segs, offset=0):
    """[(a, b, label)] -> genes as in g1.path_to_genes (type, strand, 1-based closed exons)."""
    genes, cur = [], None
    for a, b, lab in segs:
        if lab == "IG":
            if cur:
                genes.append(cur)
                cur = None
            continue
        k, c, s = lab[0], lab[1], lab[2]
        if cur is None or (cur["type"], cur["strand"]) != (k, s) or c in ("S", "A"):
            if cur:
                genes.append(cur)
            cur = {"type": k, "strand": s, "exons": []}
        if c != "N":
            cur["exons"].append([offset + a + 1, offset + b])
    if cur:
        genes.append(cur)
    return [g for g in genes if g["exons"]]
