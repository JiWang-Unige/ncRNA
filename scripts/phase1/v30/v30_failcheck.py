#!/usr/bin/env python3
"""v3.0 failed-run rule (protocol section 18.1 item 1): a run fails when its epoch-12 training NLL exceeds 3 x the
largest epoch-12 training NLL of the H arm (v20-g1b-{42,123,456}-ke).  Training losses only; nothing is read.

  v30_failcheck.py   -> outputs/PHASE1/V30/FAILCHECK.json  {threshold, runs: {name: {train_nll, failed}}, arms}
Arm run sets (section 18.1 item 2): seeds 42 / 123 / 456 that pass; a failed 123 is replaced by 789 when 789 is
trained and passes; otherwise the arm keeps its passing seeds.
"""
import json
import os

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
F = f"{ROOT}/outputs/PHASE1/FORMAL"
H = [f"v20-g1b-{s}-ke" for s in (42, 123, 456)]


def ep12(name):
    p = f"{F}/{name}/metrics.jsonl"
    if not os.path.exists(f"{F}/{name}/ep12.pt"):
        return None
    rows = [json.loads(l) for l in open(p)]
    return [r for r in rows if r["epoch"] == 12][-1]["train_nll"]


def main():
    thr = 3 * max(ep12(n) for n in H)
    runs, arms = {}, {}
    for arm in ("h6", "hm"):
        keep = []
        for s in (42, 123, 456, 789):
            n = f"v30{arm}-g1b-{s}-ke"
            v = ep12(n)
            if v is None:
                continue
            runs[n] = {"train_nll": v, "failed": v > thr}
            if s != 789 and v <= thr:
                keep.append(n)
        if f"v30{arm}-g1b-123-ke" in runs and runs[f"v30{arm}-g1b-123-ke"]["failed"]:
            r789 = runs.get(f"v30{arm}-g1b-789-ke")
            if r789 is not None and not r789["failed"]:
                keep.append(f"v30{arm}-g1b-789-ke")
        arms[arm] = {"runs": keep, "replacement_pending": f"v30{arm}-g1b-123-ke" in runs and runs[f"v30{arm}-g1b-123-ke"]["failed"]
                     and f"v30{arm}-g1b-789-ke" not in runs}
    out = {"rule": "epoch-12 train NLL > 3 x max over H arm", "H": {n: ep12(n) for n in H}, "threshold": thr, "runs": runs, "arms": arms}
    json.dump(out, open(f"{ROOT}/outputs/PHASE1/V30/FAILCHECK.json", "w"), indent=1)
    print(json.dumps(out), flush=True)


if __name__ == "__main__":
    main()
