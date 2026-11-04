#!/usr/bin/env python3
"""v2.0 tau (decision card v20-tau, user 2026-11-05 option 1; protocol section 17 amendment 1), from the chr18
development curve outputs/PHASE1/V20/dev/TAU_CURVE.json (tau 0.20-0.90 in steps of 0.01, three-run means):
  tau_base = section 16 rule on the 0.01 grid: among tau with mean strict Pr >= 0.43, the largest mean lnc L1 Sn;
  tau_v20  = among tau with mean Sn >= Sn(tau_base) - 0.02, the largest mean strict Pr (ties: larger Sn).
-> outputs/PHASE1/V20/dev/TAU_FINAL.json.  chr18 only."""
import hashlib
import json
import os

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
DEV = f"{ROOT}/outputs/PHASE1/V20/dev"


def main():
    p = f"{DEV}/TAU_CURVE.json"
    c = json.load(open(p))
    ok = [r for r in c["baseline"] if r["pr_strict"] >= 0.43]
    b = max(ok, key=lambda r: r["sn"])
    cand = [r for r in c["v20"] if r["sn"] >= b["sn"] - 0.02]
    v = max(cand, key=lambda r: (round(r["pr_strict"], 12), r["sn"]))
    out = {"rule": __doc__.split("\n", 1)[1].split("->")[0].strip(), "tau_curve_sha256": hashlib.sha256(open(p, "rb").read()).hexdigest(),
           "tau_baseline": b["tau"], "baseline_chr18": {k: b[k] for k in ("sn", "pr_strict", "sn_runs", "pr_runs")},
           "tau_v20": v["tau"], "v20_chr18": {k: v[k] for k in ("sn", "pr_strict", "sn_runs", "pr_runs")}}
    json.dump(out, open(f"{DEV}/TAU_FINAL.json", "w"), indent=1)
    print(json.dumps({k: out[k] for k in ("tau_baseline", "tau_v20")} | {"base": [round(b["sn"], 3), round(b["pr_strict"], 3)], "v20": [round(v["sn"], 3), round(v["pr_strict"], 3)]}))


if __name__ == "__main__":
    main()
