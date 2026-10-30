"""Summarize the 12-checkpoint exon-MEA gamma grid (chr18, v1.9 development): rules x operating points + half split.
  g3_analyze.py <dir with dev_decode outputs s<seed>_<ckpt>.json and half_split outputs half_*.json>"""
import json, os, sys, itertools
D = sys.argv[1]
SEEDS, CKS, GS = ["s42", "s123", "s456"], ["ep10", "ep11", "ep12", "avg"], ["0.5", "1", "1.5", "2", "3"]
VALNLL = {"s42": "ep11", "s123": "ep12", "s456": "ep12"}
R = {}
for s in SEEDS:
    for c in CKS:
        r = json.load(open(f"{D}/{s}_{c}.json"))
        for g in GS:
            m = r[f"MEAX_g{g}"]
            R[s, c, g] = {"auprc": m["auprc_lnc"], "n": m["n_pred_lnc"], "sn": m["lnc_L1"]["sn"], "pr": m["lnc_L1"]["pr"],
                          "msn": m["mrna_L1"]["sn"], "mpr": m["mrna_L1"]["pr"],
                          "grid": {x["tau"]: (x["n"], x["sn"], x["pr"] or 0) for x in m["conf_grid"]}}
print("== AUPRC by gamma (rows seed/ckpt)")
for s in SEEDS:
    for c in CKS:
        print(f"{s:5s} {c:5s}", " ".join(f"{R[s,c,g]['auprc']:.3f}" for g in GS),
              "| g1 n %d sn %.3f pr %.3f msn %.3f mpr %.3f" % tuple(R[s, c, "1"][k] for k in ("n", "sn", "pr", "msn", "mpr")))
rules = {"valNLL": {s: VALNLL[s] for s in SEEDS}, "avg": {s: "avg" for s in SEEDS},
         "auprc_g1": {s: max(CKS, key=lambda c: R[s, c, "1"]["auprc"]) for s in SEEDS}}
taus = sorted(R["s42", "ep10", "1"]["grid"])
out = {"rules": rules}
for rn, sel in rules.items():
    print(f"\n== rule {rn}: {sel}")
    print("  mean AUPRC by gamma:", " ".join(f"g{g} {sum(R[s,sel[s],g]['auprc'] for s in SEEDS)/3:.3f}" for g in GS))
    rows = []
    for g in GS:
        for t in taus:
            v = [R[s, sel[s], g]["grid"][t] for s in SEEDS]
            rows.append((g, t, sum(x[1] for x in v) / 3, sum(x[2] for x in v) / 3, [round(x[1], 3) for x in v], [round(x[2], 3) for x in v], [x[0] for x in v]))
    for floor in (0.40, 0.42, 0.45):
        ok = [r for r in rows if r[3] >= floor]
        if ok:
            b = max(ok, key=lambda r: (r[2], r[3]))
            print(f"  Pr>={floor:.2f}: gamma {b[0]} tau {b[1]} meanSn {b[2]:.3f} meanPr {b[3]:.3f} per-seed Sn {b[4]} Pr {b[5]} n {b[6]}")
        else:
            print(f"  Pr>={floor:.2f}: none")
    out[rn] = rows
    print("  full (gamma tau meanSn meanPr):")
    for g in GS:
        print("   g" + g, " ".join(f"{t}:{r[2]:.3f}/{r[3]:.3f}" for r in rows if r[0] == g for t in [r[1]]))
# half split
print("\n== half split (gamma 1): AUPRC A / B and rank of checkpoints")
for s in SEEDS:
    a = {c: json.load(open(f"{D}/half_{s}_{c}.json")) for c in CKS if os.path.exists(f"{D}/half_{s}_{c}.json")}
    if len(a) < 4:
        print(s, "incomplete"); continue
    print(f"{s:5s}", " ".join(f"{c}:{a[c]['A']['auprc']:.3f}/{a[c]['B']['auprc']:.3f}" for c in CKS),
          "| best A", max(CKS, key=lambda c: a[c]["A"]["auprc"]), "best B", max(CKS, key=lambda c: a[c]["B"]["auprc"]),
          "| nTruth A/B", a["ep10"]["A"]["n_truth"], a["ep10"]["B"]["n_truth"])
    # cross-half: choose on A, report B auprc of choice vs mean of all on B (and vice versa)
    ca, cb = max(CKS, key=lambda c: a[c]["A"]["auprc"]), max(CKS, key=lambda c: a[c]["B"]["auprc"])
    mB = sum(a[c]["B"]["auprc"] for c in CKS) / 4; mA = sum(a[c]["A"]["auprc"] for c in CKS) / 4
    print(f"       choose-on-A -> B {a[ca]['B']['auprc']:.3f} (mean B {mB:.3f}, avg-ckpt B {a['avg']['B']['auprc']:.3f}); "
          f"choose-on-B -> A {a[cb]['A']['auprc']:.3f} (mean A {mA:.3f}, avg-ckpt A {a['avg']['A']['auprc']:.3f})")
json.dump(out, open(f"{D}/g3_summary.json", "w"), indent=1)
