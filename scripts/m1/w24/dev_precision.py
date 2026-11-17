#!/usr/bin/env python3
"""chr18 development: S50 precision with 'unknown' predictions (on WTC11-transcribed loci below the truth bar) left out,
for the DNA-only calls, snapping (w 1 kb), pairing (d 50 kb) and the two long-read tools.  -> outputs/M1/W24/dev/PRECISION_chr18.json"""
import json
import numpy as np
import w24_lib as L
from dev_baselines import pair, snap
from dev_longread_score import loci_of, tib_hit
import simple_det as D

C = "chr18"
truth = L.load_truth(C)
tib = D.read_tiberius(C)
ours = D.load_ours(C)
model = [m for m in ours if m != "v1.9"][0]
tau = ours[model]["tau"]
dna = {r: [p for p in v if p["conf"] >= tau] for r, v in ours[model]["runs"].items()}
ev = L.evidence(C, 1.0, 0)
O = f"{L.W}/longread/{C}"
sets = {f"dna_{r}": p for r, p in dna.items()}
sets.update({f"snap1k_{r}": snap(p, ev, 1000) for r, p in dna.items()})
sets["pair50k"] = pair(ev, tib, 50_000)
for tool, gtf in (("stringtie", f"{O}/stringtie.gtf"), ("isoquant", f"{O}/isoquant/rep1/rep1.transcript_models.gtf")):
    sets[tool] = [p for p in loci_of(gtf, C) if not tib_hit(tib, p["strand"], p["exons"][0][0], p["exons"][-1][1])]
out = {}
for k, p in sets.items():
    r = L.score_unknown(truth, p, C)
    r.pop("S50_ids")
    out[k] = {kk: (round(v, 3) if isinstance(v, float) else v) for kk, v in r.items()}
    print(k, out[k])
json.dump(out, open(f"{L.W}/dev/PRECISION_{C}.json", "w"), indent=1)
