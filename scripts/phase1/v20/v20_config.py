#!/usr/bin/env python3
"""v2.0 read configurations (protocol section 17, detail 5), written on ont.

  v20_config.py dev   --commit <sha>  -> configs/v20_dev.json  (chr18 only; the development reads in outputs/PHASE1/V20/dev)
  v20_config.py final --commit <sha>  -> configs/v20.json      (chr18 reproduction + chr9; outputs/PHASE1/V20/read)

Same structure as the P1.5 configs, so p15_read.py reads them unchanged: the three retrained main runs
(v20-g1b-<seed>-ke, epoch 10-12 weights with sha256), the v1.9 decoding (exon-level MEA gamma* = 3, exact
posteriors), and the sha256 of the reading code.  "decode.tau" is the baseline tau (raw confidence, section 16 rule
re-applied on chr18); the v2.0 calibration (coding-potential logistic, its keep decision, tau_v20, hexamer table
hash) is under "calibration" and is applied by the scoring step, never by the decoder.  final adds the chr18
reference (sha256 of each development posterior + the p15_score.core summary at the baseline tau), which the chr18
reproduction must match before chr9 is read.
"""
import glob
import hashlib
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts/phase1"))
RUNS = {"42": "v20-g1b-42-ke", "123": "v20-g1b-123-ke", "456": "v20-g1b-456-ke"}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def main():
    stage = sys.argv[1]
    commit = sys.argv[sys.argv.index("--commit") + 1]
    prev = json.load(open(f"{ROOT}/configs/p15_v19b.json"))
    runs = {}
    for seed, name in RUNS.items():
        ck = [f"outputs/PHASE1/FORMAL/{name}/ep{e}.pt" for e in (10, 11, 12)]
        runs[name] = {"arm": "G1B", "family": "G1B", "seed": int(seed), "labels": "v50uev", "lnc_only": False, "data": "v20",
                      "ckpts": ck, "ckpt_sha256": [sha256(f"{ROOT}/{p}") for p in ck]}
    decode = dict(prev["decode"])
    targets = {"chr18": dict(prev["targets"]["chr18"], runs=sorted(runs))}
    code = sorted(glob.glob(f"{ROOT}/src/ncrna_annot/**/*.py", recursive=True)
                  + [f"{ROOT}/scripts/phase1/{p}" for p in ("p15_read.py", "p15_score.py")]
                  + glob.glob(f"{ROOT}/scripts/phase1/v20/*.py"))
    cfg = {"protocol": "PHASE1_PROTOCOL.md section 17 (v2.0)", "stage": stage, "code_commit": commit,
           "created": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "main": {"family": "G1B", "runs": RUNS},
           "runs": runs, "decode": decode}
    if stage == "dev":
        cfg["out_dir"] = "outputs/PHASE1/V20/dev"
        out = f"{ROOT}/configs/v20_dev.json"
    else:
        import p15_score as P
        from ncrna_annot.eval import ruler
        from ncrna_annot.features import generanno as G
        cal = json.load(open(f"{ROOT}/outputs/PHASE1/V20/dev/CALIB.json"))
        tf = json.load(open(f"{ROOT}/outputs/PHASE1/V20/dev/TAU_FINAL.json"))   # section 17 amendment 1 (card v20-tau)
        decode["tau"] = tf["tau_baseline"]
        cfg["tau_final"] = {"file": "outputs/PHASE1/V20/dev/TAU_FINAL.json", "sha256": sha256(f"{ROOT}/outputs/PHASE1/V20/dev/TAU_FINAL.json"),
                            "rule": tf["rule"]}
        cfg["calibration"] = {"keep": cal["keep"], "model": cal["model"], "tau_v20": tf["tau_v20"],
                              "hexamer_table": "data/p1v20/hexamer_table.json", "hexamer_table_sha256": cal["hexamer_table_sha256"],
                              "dev_result": "outputs/PHASE1/V20/dev/CALIB.json", "dev_result_sha256": sha256(f"{ROOT}/outputs/PHASE1/V20/dev/CALIB.json")}
        truth = ruler.read_gtf(f"{ROOT}/{targets['chr18']['truth']}", "chr18")
        ref = {}
        for name in runs:
            rd = f"{ROOT}/outputs/PHASE1/V20/dev/chr18/{name}"
            genes = P.load_genes(f"{rd}/genes.json")
            post = np.load(f"{rd}/post.npy", mmap_mode="r")
            ref[name] = {"post": os.path.relpath(f"{rd}/post.npy", ROOT), "post_sha256": sha256(f"{rd}/post.npy"),
                         "summary": P.norm(P.core(genes, post, truth, decode["tau"]))}
        cfg["chr18_reference"] = ref
        L9 = len(np.load(f"{ROOT}/data/p1/chr9_codes.npy", mmap_mode="r"))
        targets["chr9"] = {"feats": "data/p1_feat/chr9.npy", "cover": "data/p1_feat/chr9.cover.json", "tiles": len(G.chromosome_tiles(L9)),
                           "codes": "data/p1/chr9_codes.npy", "truth": "data/human/gencode.v50.annotation.gtf.gz",
                           "truth_chrom": "chr9", "seqid": "chr9", "runs": sorted(runs)}
        cfg["out_dir"] = "outputs/PHASE1/V20/read"
        out = f"{ROOT}/configs/v20.json"
    cfg["targets"] = targets
    cfg["code_files"] = {os.path.relpath(p, ROOT): sha256(p) for p in code}
    json.dump(cfg, open(out, "w"), indent=1)
    print(json.dumps({"out": out, "runs": sorted(runs), "tau": decode["tau"], "targets": sorted(targets)}), flush=True)


if __name__ == "__main__":
    main()
