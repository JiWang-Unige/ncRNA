"""Descriptive (not in the plan): chance level of the best-of-273-tracks 'oracle' -- the same rule scored at decoy
positions TSS +- 10 kb (inside the +-20 kb window) instead of the true start."""
import json, numpy as np
O = "outputs/M1/AG"
rows = {r["id"]: r for r in json.load(open(f"{O}/tss.json"))}
lnc = [r for r in rows.values() if r["chrom"] in ("chr8", "chr9") and r["k"] == "l"]
q = np.quantile([r["cpg"] for r in lnc], [0.25, 0.5, 0.75])
acc = {}
for m in ("fold-0", "fold-1", "fold-2", "fold-3"):
    z = np.load(f"{O}/colab_auto/raw/{m}.npz")
    for j, rid in enumerate(z["rid"]):
        r = rows[str(rid)]
        if r["chrom"] not in ("chr8", "chr9"):
            continue
        g = "mrna_cage" if r["k"] == "m" and r["TSS_ev"] else ("lnc_lowcpg" if r["k"] == "l" and r["cpg"] <= q[1] else None)
        if g is None:
            continue
        ca, am = z["cage_arg"][j].astype(int), int(z["agg_cage_mean_arg"][j])
        d = acc.setdefault(g, {"n": 0, "best_true": 0, "best_decoy": 0, "mean_true": 0, "mean_decoy": 0, "ntracks_hit_true": [], "ntracks_hit_decoy": []})
        d["n"] += 1
        d["best_true"] += np.abs(ca).min() <= 1000
        d["mean_true"] += abs(am) <= 1000
        d["ntracks_hit_true"].append(int((np.abs(ca) <= 1000).sum()))
        for off in (-10000, 10000):
            d["best_decoy"] += 0.5 * (np.abs(ca - off).min() <= 1000)
            d["mean_decoy"] += 0.5 * (abs(am - off) <= 1000)
            d["ntracks_hit_decoy"].append(int((np.abs(ca - off) <= 1000).sum()))
out = {}
for g, d in acc.items():
    n = d["n"]
    out[g] = {"n": n, "best_track_true": round(d["best_true"] / n, 3), "best_track_decoy_pm10kb": round(d["best_decoy"] / n, 3),
              "mean_track_true": round(d["mean_true"] / n, 3), "mean_track_decoy_pm10kb": round(d["mean_decoy"] / n, 3),
              "median_tracks_hitting_true": float(np.median(d["ntracks_hit_true"])),
              "frac_starts_ge_half_tracks_hit_true": round(float(np.mean(np.array(d["ntracks_hit_true"]) >= 137)), 3)}
json.dump(out, open(f"{O}/DECOY_best_track.json", "w"), indent=1)
print(json.dumps(out, indent=1))
