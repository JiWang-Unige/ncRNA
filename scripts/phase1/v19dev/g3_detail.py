"""Per-seed exon-MEA operating points (gamma x tau) for three checkpoint sets of the main model (chr18, v1.9 development).
  g3_detail.py <dir with dev_decode outputs s<seed>_<ckpt>.json [+ g46_*.json]>"""
import json, sys, os
D = sys.argv[1]
def load(n):
    r = json.load(open(f"{D}/{n}.json"))
    if os.path.exists(f"{D}/g46_{n}.json"):
        r.update({k: v for k, v in json.load(open(f"{D}/g46_{n}.json")).items() if k.startswith("MEAX")})
    return r
for sel_name, sel in (("avg", ["s42_avg", "s123_avg", "s456_avg"]), ("C", ["s42_ep11", "s123_avg", "s456_avg"]), ("C2", ["s42_ep11", "s123_ep10", "s456_avg"])):
    rs = [load(n) for n in sel]
    gs = [k for k in rs[0] if k.startswith("MEAX") and all(k in r for r in rs)]
    print("==", sel_name, sel)
    for g in sorted(gs, key=lambda k: float(k.split("_g")[1])):
        ms = [r[g] for r in rs]
        print(f" {g:9s} AUPRC {sum(m['auprc_lnc'] for m in ms)/3:.3f} mRNA Sn {[round(m['mrna_L1']['sn'],3) for m in ms]} Pr {[round(m['mrna_L1']['pr'],3) for m in ms]}")
        for t in (0.2, 0.5, 0.6, 0.7, 0.8):
            v = [next(x for x in m["conf_grid"] if x["tau"] == t) for m in ms]
            print(f"    tau {t}: mean {sum(x['sn'] for x in v)/3:.3f}/{sum((x['pr'] or 0) for x in v)/3:.3f}  n {[x['n'] for x in v]} sn {[round(x['sn'],3) for x in v]} pr {[round(x['pr'] or 0,3) for x in v]}")
