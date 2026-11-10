#!/usr/bin/env python3
"""v3.0 test read configuration (protocol section 18, execution details 8-9), written on ont after the chr18
development and before any test set is read.

  v30_config.py --commit <sha> [--stage final|preview]  -> configs/v30.json | configs/v30_preview.json

Section 18.1: preview = H plus seeds 42 / 456 of H6 and HM (calibration and tau from V30/dev_preview); final = H plus
the FAILCHECK.json run sets (calibration and tau from V30/dev; the seed-789 runs come from configs/v30_dev2.json).
Both read into V30/read; the final config must carry the preview's code_files unchanged and records the preview
config's sha256, so reads made with the preview config count as the final reads of those runs.

Same structure as the P1.6 configs, so p16_read.py reads it unchanged:
- runs: the nine runs (H = v20-g1b-*-ke with the sha256 recorded in configs/v20.json; H6 / HM = v30h6 / v30hm-g1b-*-ke,
  epoch 10-12 weights, sha256 now);
- arms: per arm its runs, calibration (keep, model, hexamer table and sha256) and tau_v20 (H: configs/v20.json;
  H6 / HM: outputs/PHASE1/V30/dev/<arm>/{CALIB,TAU_FINAL}.json, with their sha256);
- decode: the v2.0 decoding (decode.tau = H's tau_base, used only for pred.gff3);
- chr18_reference: sha256 of each run's development chr18 posterior and genes (H: V20/dev, H6 / HM: V30/dev), which
  the chr18 reproduction must match before any test set is read;
- targets: chr18 + A mouse_chr19 + B chicken_5, zebrafish_22, fly_2L, worm_I + C (data/v30/C_CHROMS.json);
- code_files: sha256 of src/ncrna_annot/**/*.py, p15_read.py and p16_read.py.
"""
import glob
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1/v20"))
from v20_config import sha256  # noqa: E402

V30 = "outputs/PHASE1/V30"
C_NAMES = {"rat": "rat", "cow": "cow", "xtrop": "xtrop"}


def c_targets():
    cc = json.load(open(f"{ROOT}/data/v30/C_CHROMS.json"))
    out = {}
    for short in C_NAMES:
        c = cc[short]
        out[f"{short}_{c['chrom']}"] = {"truth": c["gtf"], "truth_chrom": c["chrom"], "seqid": c["chrom"], "in_verdict": c["in_verdict"]}
    return out


def main():
    commit = sys.argv[sys.argv.index("--commit") + 1]
    stage = sys.argv[sys.argv.index("--stage") + 1] if "--stage" in sys.argv else "final"
    from ncrna_annot.features import generanno as G
    v20 = json.load(open(f"{ROOT}/configs/v20.json"))
    dev_cfgs = ["configs/v30_dev.json"] + (["configs/v30_dev2.json"] if stage == "final" and os.path.exists(f"{ROOT}/configs/v30_dev2.json") else [])
    dev = {"runs": {}}
    for dc in dev_cfgs:
        dev["runs"].update(json.load(open(f"{ROOT}/{dc}"))["runs"])
    if stage == "preview":
        use = {arm: [f"v30{arm}-g1b-{s}-ke" for s in (42, 456)] for arm in ("h6", "hm")}
    else:
        fc = json.load(open(f"{ROOT}/{V30}/FAILCHECK.json"))["arms"]
        assert not any(fc[a]["replacement_pending"] for a in fc), "seed-789 replacement pending"
        use = {arm: fc[arm]["runs"] for arm in ("h6", "hm")}
    dev["runs"] = {n: r for n, r in dev["runs"].items() if n in use["h6"] + use["hm"]}
    devdir = f"{V30}/dev" if stage == "final" else f"{V30}/dev_preview"
    runs = {}
    for name, r in v20["runs"].items():
        bad = [c for c, h in zip(r["ckpts"], r["ckpt_sha256"]) if sha256(f"{ROOT}/{c}") != h]
        assert not bad, f"{name}: weights changed {bad}"
        runs[name] = r
    for name, r in dev["runs"].items():
        bad = [c for c, h in zip(r["ckpts"], r["ckpt_sha256"]) if sha256(f"{ROOT}/{c}") != h]
        assert not bad, f"{name}: weights changed since the dev read {bad}"
        runs[name] = r
    arms = {"H": {"runs": sorted(v20["runs"]), "calibration": v20["calibration"], "tau_baseline": v20["decode"]["tau"],
                  "tau_v20": v20["calibration"]["tau_v20"], "source": "configs/v20.json", "source_sha256": sha256(f"{ROOT}/configs/v20.json")}}
    for arm in ("h6", "hm"):
        d = f"{ROOT}/{devdir}/{arm}"
        cal, tf = json.load(open(f"{d}/CALIB.json")), json.load(open(f"{d}/TAU_FINAL.json"))
        assert sorted(cal["by_run"]) == sorted(use[arm]), f"{arm}: calibration fitted on {sorted(cal['by_run'])}"
        arms[arm.upper()] = {"runs": sorted(use[arm]),
                             "calibration": {"keep": cal["keep"], "model": cal["model"], "tau_v20": tf["tau_v20"],
                                             "hexamer_table": "data/p1v20/hexamer_table.json", "hexamer_table_sha256": cal["hexamer_table_sha256"],
                                             "dev_result": f"{devdir}/{arm}/CALIB.json", "dev_result_sha256": sha256(f"{d}/CALIB.json")},
                             "tau_baseline": tf["tau_baseline"], "tau_v20": tf["tau_v20"],
                             "tau_final": f"{devdir}/{arm}/TAU_FINAL.json", "tau_final_sha256": sha256(f"{d}/TAU_FINAL.json"),
                             "chr18_dev": {"baseline": tf["baseline_chr18"], "v20": tf["v20_chr18"]}}
    arms["H"]["chr18_dev"] = {"v20": json.load(open(f"{ROOT}/outputs/PHASE1/V20/dev/TAU_FINAL.json"))["v20_chr18"]}
    ref = {}
    for name in runs:
        rd = f"{ROOT}/outputs/PHASE1/{'V20' if name.startswith('v20-') else 'V30'}/dev/chr18/{name}"
        m = json.load(open(f"{rd}/read_meta.json"))
        want = ["configs/v20_dev.json"] if name.startswith("v20-") else dev_cfgs
        assert m["config_sha256"] in [sha256(f"{ROOT}/{w}") for w in want], f"{name}: chr18 dev read of another config"
        ref[name] = {"read": os.path.relpath(rd, ROOT), "post_sha256": m["post_sha256"], "genes_sha256": sha256(f"{rd}/genes.json")}
    allr = sorted(runs)
    targets = {"chr18": dict(v20["targets"]["chr18"], runs=allr)}
    p16, p17 = json.load(open(f"{ROOT}/configs/p16_v19.json")), json.load(open(f"{ROOT}/configs/p17_mouse.json"))
    targets["mouse_chr19"] = dict(p17["targets"]["mouse_chr19"], runs=allr, group="A", in_verdict=False)
    for t in ("chicken_5", "zebrafish_22", "fly_2L", "worm_I"):
        targets[t] = dict(p16["targets"][t], runs=allr, group="B", in_verdict=True)
    for t, c in c_targets().items():
        codes = f"data/p1/{t}_codes.npy"
        L = len(np.load(f"{ROOT}/{codes}", mmap_mode="r"))
        targets[t] = {"feats": f"data/p1_feat/{t}.npy", "cover": f"data/p1_feat/{t}.cover.json", "tiles": len(G.chromosome_tiles(L)),
                      "codes": codes, "codes_meta": json.load(open(f"{ROOT}/{codes[:-4]}.json")), "truth": c["truth"],
                      "truth_chrom": c["truth_chrom"], "seqid": c["seqid"], "runs": allr, "group": "C", "in_verdict": c["in_verdict"]}
    code = sorted(glob.glob(f"{ROOT}/src/ncrna_annot/**/*.py", recursive=True)
                  + [f"{ROOT}/scripts/phase1/p15_read.py", f"{ROOT}/scripts/phase1/p16_read.py"])
    cfg = {"protocol": "PHASE1_PROTOCOL.md section 18 (v3.0), 18.1", "stage": stage, "code_commit": commit,
           "created": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "out_dir": f"{V30}/read",
           "decode": dict(v20["decode"]), "main": {"family": "G1B", "runs": {n: n for n in allr}}, "runs": runs, "arms": arms,
           "chr18_reference": ref, "targets": targets,
           "c_chroms": {"file": "data/v30/C_CHROMS.json", "sha256": sha256(f"{ROOT}/data/v30/C_CHROMS.json")},
           "code_files": {os.path.relpath(p, ROOT): sha256(p) for p in code}}
    out = "configs/v30.json" if stage == "final" else "configs/v30_preview.json"
    if stage == "final" and os.path.exists(f"{ROOT}/configs/v30_preview.json"):
        pv = json.load(open(f"{ROOT}/configs/v30_preview.json"))
        assert pv["code_files"] == cfg["code_files"], "read code changed since the preview reads"
        for n, r in pv["runs"].items():
            assert n not in runs or runs[n]["ckpt_sha256"] == r["ckpt_sha256"], n
        cfg["preview_config"] = {"path": "configs/v30_preview.json", "sha256": sha256(f"{ROOT}/configs/v30_preview.json")}
    json.dump(cfg, open(f"{ROOT}/{out}", "w"), indent=1)
    print(json.dumps({"out": out, "runs": allr, "arms": {k: (v["tau_v20"], v["calibration"]["keep"]) for k, v in arms.items()},
                      "targets": {k: v.get("group", "dev") for k, v in targets.items()}}), flush=True)


if __name__ == "__main__":
    main()
