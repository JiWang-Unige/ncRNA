#!/usr/bin/env python3
"""Protocol v1.9 section 16: write the P1.5 read configuration (configs/p15_v19.json) on ont, from chr18 results only.

  p15_config.py --commit <git sha of the code> --out configs/p15_v19.json [--draft | --partial]
  p15_config.py --commit <git sha of the code> --out configs/p15_v19b.json --previous configs/p15_v19.json

Contents: decoding settings (gamma* = 3, tau* = 0.6, exon-level MEA, exact posteriors); every run's epoch 10-12
weights with sha256; the sha256 of every src/ncrna_annot/**/*.py and scripts/phase1/p15_*.py file (p15_read.py
refuses other code); the chr18 v1.9 evaluation of every run (MEAX gamma*: FORMAL/v19-<run>/result.json,
V19/v19-ab2-f-123.json for ab2-f-123, DEV/g3/s<seed>_avg.json for the main runs); the main-model selection
(section 15 rule, section 16 item 5: an arm replaces G1B only if mean Delta AUPRC_lnc >= +0.01 and all three
seeds' Delta > 0; AB1, AB4, AB6 do not take part); the chr18 reference that the read must reproduce (main runs:
sha256 of the development posterior dump + the development MEAX gamma* summary); targets and the runs read on each.
--draft: allow missing runs / evaluations (selection "pending"); a draft is for testing p15_read.py on chr18 only.
Two read batches (user decision 2026-11-01, dashboard card p15-read-order; protocol deviation log):
--partial  batch 1: runs whose epoch 10-12 weights do not exist yet are deferred (not read); every run with weights
           needs its chr18 evaluation, and the main model must already be fixed by the selection rule: an arm with
           any seed's Delta <= 0 cannot have "Delta > 0 in all three seeds", so it is ruled out before its other
           seeds finish.  The main runs are read on chr18 / chr20 / t2t_chr20, the finished ablations on chr20.
--previous batch 2 (complete): chr20 lists only the runs the previous config did not read (every run is read
           once), t2t_chr20 none; the decoding settings, the main model, the reading code (src/ncrna_annot,
           p15_read.py) and the weights of runs already read must equal the previous config.  Outputs go to
           outputs/PHASE1/P15b (own chr18 reproduction and REPRO_OK).
"""
import argparse
import glob
import hashlib
import json
import os
import sys
import time

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
FORMAL, DEV, V19 = (f"{ROOT}/outputs/PHASE1/{d}" for d in ("FORMAL", "DEV/g3", "V19"))
GAMMA, TAU, KEY = 3.0, 0.6, "MEAX_g3"
SEEDS = ("42", "123", "456")
MAIN = {"42": "g1b-42-ke", "123": "g1b-123-ke", "456": "g1b-456-ke"}
ARMS = {  # arm: (train_formal_ont arm, lnc_only, labels, runs by seed)
    "G1B": ("G1B", False, "v50uev", MAIN),
    "AB0": ("R", False, "v50uev", {"42": "ab0-r-42-ke", "123": "ab0-r-123", "456": "ab0-r-456"}),
    "AB2": ("G1B", False, "v50f", {"42": "ab2-f-42-ke", "123": "ab2-f-123", "456": "ab2-f-456"}),
    "AB3": ("G1B", True, "v50uevlnc", {"42": "ab3-lnc-42", "123": "ab3-lnc-123", "456": "ab3-lnc-456"}),
    "AB4": ("G1B", False, "v41uev", {"42": "ab4-v41-42-ke", "123": "ab4-v41-123", "456": "ab4-v41-456"}),
    "AB6": ("G2PH", False, "v50uev", {"42": "ab6-g2ph-42"}),
}
SELECTABLE = ("AB0", "AB2", "AB3")
DECODE = {"posteriors": "exact (float64 CPU forward-backward per tile, chrom.chromosome_marginals exact=True)",
          "checkpoint": "uniform weight average of epochs 10-12", "mode": "mea", "level": "exon", "delta": 1e-4,
          "gamma": GAMMA, "tau": TAU, "min_exon": 20, "min_intron": 50, "lnc_min_len": 200,
          "transitions": "allowed transitions of g1.G1Transitions() (mea_duration scores none)"}
DEV_DUMP = {"42": "dev-g1b42ke-avg", "123": "dev-g1b123ke-avg", "456": "dev-g1b456ke-avg"}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def summary(path):
    r = json.load(open(path))[KEY]
    return {k: v for k, v in r.items() if k not in ("curve", "decode_s")}


def chr18_eval(arm, seed, name):
    if arm == "G1B":
        p = f"{DEV}/s{seed}_avg.json"
    elif name == "ab2-f-123":
        p = f"{V19}/v19-ab2-f-123.json"
    else:
        p = f"{FORMAL}/v19-{name}/result.json"
    return (p, summary(p)) if os.path.exists(p) else (p, None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--commit", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--draft", action="store_true")
    ap.add_argument("--partial", action="store_true")
    ap.add_argument("--previous", default="")
    a = ap.parse_args()
    assert not (a.partial and a.previous), "--partial is batch 1, --previous batch 2"
    missing, deferred, runs, ev = [], [], {}, {}
    for arm, (tarm, lnc_only, labels, by_seed) in ARMS.items():
        for seed, name in by_seed.items():
            cks = [f"outputs/PHASE1/FORMAL/{name}/ep{e:02d}.pt" for e in (10, 11, 12)]
            if not all(os.path.exists(f"{ROOT}/{c}") for c in cks):
                if a.partial and arm != "G1B":
                    deferred.append(name)
                else:
                    missing.append(f"{name}: weights")
                continue
            runs[name] = {"arm": tarm, "family": arm, "seed": int(seed), "labels": labels, "lnc_only": lnc_only,
                          "ckpts": cks, "ckpt_sha256": [sha256(f"{ROOT}/{c}") for c in cks]}
            if arm != "AB6":
                p, s = chr18_eval(arm, seed, name)
                if s is None:
                    missing.append(f"{name}: chr18 v1.9 evaluation ({os.path.relpath(p, ROOT)})")
                else:
                    ev[name] = {"source": os.path.relpath(p, ROOT), "auprc_lnc": s["auprc_lnc"], "n_pred_lnc": s["n_pred_lnc"],
                                "tau": [x for x in s["conf_grid"] if x["tau"] == TAU][0], "mrna_L1": s["mrna_L1"]}
    if missing and not a.draft:
        sys.exit("incomplete:\n  " + "\n  ".join(missing))
    sel = {"rule": "replace G1B by an arm only if mean Delta AUPRC_lnc (chr18, MEAX gamma*, no threshold) >= +0.01 "
                   "and Delta > 0 in all three seeds; AB1 / AB4 / AB6 excluded; an arm with any seed's Delta <= 0 "
                   "is ruled out before its other seeds finish", "arms": {}}
    winners = []
    for arm in SELECTABLE:
        by_seed = ARMS[arm][3]
        dl = {s: ev[by_seed[s]]["auprc_lnc"] - ev[MAIN[s]]["auprc_lnc"] for s in SEEDS if by_seed[s] in ev and MAIN[s] in ev}
        if len(dl) < 3 and all(v > 0 for v in dl.values()):
            sel["arms"][arm] = "pending"
            continue
        mdl = sum(dl.values()) / 3 if len(dl) == 3 else None
        rep = len(dl) == 3 and mdl >= 0.01 and all(v > 0 for v in dl.values())
        sel["arms"][arm] = {"delta": dl, "mean_delta": mdl, "replaces": rep}
        if rep:
            winners.append((mdl, arm))
    pending = any(v == "pending" for v in sel["arms"].values())
    main_arm = "pending" if pending else (max(winners)[1] if winners else "G1B")
    if a.partial and (pending or main_arm != "G1B" or missing):
        sys.exit(f"batch 1 needs the main model fixed and every trained run evaluated: main {main_arm}, missing {missing}")
    sel["main_family"] = main_arm
    main_runs = ARMS[main_arm][3] if main_arm != "pending" else MAIN
    arms = {k: v[3] for k, v in ARMS.items() if k not in (main_arm, "AB6")}
    code = sorted(glob.glob(f"{ROOT}/src/ncrna_annot/**/*.py", recursive=True) + glob.glob(f"{ROOT}/scripts/phase1/p15_*.py"))
    ref = {}
    for s in SEEDS:
        dump = f"{FORMAL}/{DEV_DUMP[s]}/chr18_post.npy"
        if os.path.exists(dump) and os.path.exists(f"{DEV}/s{s}_avg.json"):
            ref[MAIN[s]] = {"post": os.path.relpath(dump, ROOT), "post_sha256": sha256(dump), "summary": summary(f"{DEV}/s{s}_avg.json")}
    all_runs = [n for n in runs]
    t2t_runs, out_dir, prev_info = list(main_runs.values()), "outputs/PHASE1/P15", None
    if a.previous:
        pp = a.previous if os.path.isabs(a.previous) else f"{ROOT}/{a.previous}"
        prev = json.load(open(pp))
        reading = lambda files: {k: v for k, v in files.items() if k.startswith("src/") or k.endswith("p15_read.py")}
        code_now = {os.path.relpath(p, ROOT): sha256(p) for p in code}
        bad = [k for k in ("decode",) if prev[k] != json.loads(json.dumps(DECODE))]
        bad += ["main"] if prev["main"]["runs"] != main_runs else []
        bad += ["reading code"] if reading(prev["code_files"]) != reading(code_now) else []
        bad += [f"{n}: weights" for n in prev["runs"] if n in runs and prev["runs"][n]["ckpt_sha256"] != runs[n]["ckpt_sha256"]]
        if bad:
            sys.exit(f"batch 2 differs from {a.previous}: {bad}")
        done = set(prev["targets"]["chr20"]["runs"])
        all_runs = [n for n in all_runs if n not in done]
        t2t_runs = [n for n in t2t_runs if n not in prev["targets"]["t2t_chr20"]["runs"]]
        out_dir, prev_info = "outputs/PHASE1/P15b", {"config": a.previous, "config_sha256": sha256(pp), "chr20_read": sorted(done)}
    cfg = {
        "protocol": "PHASE1_PROTOCOL.md section 16 (v1.9, frozen 2026-10-30)", "draft": bool(missing) or a.draft,
        "missing": missing, "batch": 2 if a.previous else (1 if a.partial else None), "deferred": deferred,
        "previous": prev_info, "out_dir": out_dir, "code_commit": a.commit, "created": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "decode": DECODE,
        "main": {"family": main_arm, "runs": main_runs}, "selection": sel, "arms": arms, "runs": runs, "chr18_v19": ev,
        "chr18_reference": ref,
        "targets": {
            "chr18": {"feats": "data/p1_feat/chr18.npy", "cover": "data/p1_feat/chr18.cover.json", "tiles": 3926,
                      "codes": "data/p1/chr18_codes.npy", "truth": "outputs/PHASE1/DIAG_L2/chr18_v50.gtf",
                      "truth_chrom": "chr18", "seqid": "chr18", "runs": sorted(ref)},
            "chr20": {"feats": "data/p1_feat/chr20.npy", "cover": "data/p1_feat/chr20.cover.json", "tiles": 3148,
                      "codes": "data/p1/chr20_codes.npy", "truth": "data/human/gencode.v50.annotation.gtf.gz",
                      "truth_chrom": "chr20", "seqid": "chr20", "runs": all_runs},
            "t2t_chr20": {"feats": "data/p1_feat/t2t_chr20.npy", "cover": "data/p1_feat/t2t_chr20.cover.json", "tiles": 3234,
                          "codes": "data/p1/t2t_chr20_codes.npy",
                          "truth": "data/baselines/t2t_chr20/leaderboard_data/ground_truth/chr20.gff",
                          "truth_chrom": "NC_060944.1", "seqid": "NC_060944.1", "runs": t2t_runs}},
        "code_files": {os.path.relpath(p, ROOT): sha256(p) for p in code},
    }
    out = a.out if os.path.isabs(a.out) else f"{ROOT}/{a.out}"
    os.makedirs(os.path.dirname(out), exist_ok=True)
    json.dump(cfg, open(out, "w"), indent=1)
    print(json.dumps({"out": out, "draft": cfg["draft"], "batch": cfg["batch"], "missing": missing, "deferred": deferred,
                      "chr20": all_runs, "t2t_chr20": t2t_runs, "selection": sel}, indent=1))


if __name__ == "__main__":
    main()
