#!/usr/bin/env python3
"""v3.0 chr18 development (protocol section 18, execution detail 8): for each new arm, the v2.0 calibration and tau
procedure unchanged (v20_calib.calibrate -> v20_tau_curve -> v20_tau_final, i.e. the section 17.1 rule), on the
human chr18 reads of that arm's three runs.  The H arm keeps configs/v20.json and is not refitted.

  v30_dev.py config --commit <sha>   -> configs/v30_dev.json (chr18 target, the six new runs, out_dir V30/dev)
  v30_dev.py arm h6|hm                -> outputs/PHASE1/V30/dev/<arm>/{CALIB,TAU_CURVE,TAU_FINAL}.json
Each arm works in its own directory whose chr18/ links to the shared chr18 reads, so the v2.0 scripts run as they are.
"""
import glob
import json
import os
import sys
import time

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1/v20"))
DEV = f"{ROOT}/outputs/PHASE1/V30/dev"
ARMS = {"h6": [f"v30h6-g1b-{s}-ke" for s in (42, 123, 456)], "hm": [f"v30hm-g1b-{s}-ke" for s in (42, 123, 456)]}


def config(commit):
    import v20_config as VC
    v20 = json.load(open(f"{ROOT}/configs/v20.json"))
    runs = {}
    for arm, names in ARMS.items():
        for name in names:
            ck = [f"outputs/PHASE1/FORMAL/{name}/ep{e}.pt" for e in (10, 11, 12)]
            runs[name] = {"arm": "G1B", "family": "G1B", "seed": int(name.split("-")[2]), "labels": "v50uev", "lnc_only": False,
                          "data": f"v30{arm}", "ckpts": ck, "ckpt_sha256": [VC.sha256(f"{ROOT}/{p}") for p in ck]}
    targets = {"chr18": dict(v20["targets"]["chr18"], runs=sorted(runs))}
    code = sorted(glob.glob(f"{ROOT}/src/ncrna_annot/**/*.py", recursive=True)
                  + [f"{ROOT}/scripts/phase1/{p}" for p in ("p15_read.py", "p15_score.py")]
                  + glob.glob(f"{ROOT}/scripts/phase1/v20/*.py") + glob.glob(f"{ROOT}/scripts/phase1/v30/*.py"))
    cfg = {"protocol": "PHASE1_PROTOCOL.md section 18 (v3.0)", "stage": "dev", "code_commit": commit,
           "created": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "main": {"family": "G1B", "runs": {n: n for n in sorted(runs)}}, "runs": runs, "decode": dict(v20["decode"]),
           "out_dir": "outputs/PHASE1/V30/dev", "targets": targets,
           "code_files": {os.path.relpath(p, ROOT): VC.sha256(p) for p in code}}
    json.dump(cfg, open(f"{ROOT}/configs/v30_dev.json", "w"), indent=1)
    print(json.dumps({"out": "configs/v30_dev.json", "runs": sorted(runs)}), flush=True)


def arm(name):
    import v20_calib as C
    import v20_tau_curve as TC
    import v20_tau_final as TF
    d = f"{DEV}/{name}"
    os.makedirs(d, exist_ok=True)
    if not os.path.exists(f"{d}/chr18"):
        os.symlink("../chr18", f"{d}/chr18")
    for r in ARMS[name]:
        assert os.path.exists(f"{d}/chr18/{r}/genes.json"), f"no chr18 read of {r}"
    C.DEV, C.RUNS = d, tuple(ARMS[name])
    TF.DEV = d
    C.calibrate()
    TC.main()
    TF.main()
    print(json.dumps({"arm": name, **{k: v for k, v in json.load(open(f"{d}/TAU_FINAL.json")).items() if k in ("tau_baseline", "tau_v20")},
                      "keep": json.load(open(f"{d}/CALIB.json"))["keep"]}), flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "config":
        config(sys.argv[sys.argv.index("--commit") + 1])
    else:
        arm(sys.argv[2])
