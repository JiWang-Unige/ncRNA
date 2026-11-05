#!/usr/bin/env python3
"""v3.0 C-group test chromosomes (protocol section 18, execution detail 2), chosen from the annotation and assembly
only, before any model output on these species exists: among the autosomes of 40-80 Mb (Ensembl REST info/assembly,
karyotype minus X / Y / Z / W / MT), the one with the most genes the ruler classes as lnc in the Ensembl 116 GTF;
a species whose pick has < 200 lnc genes is descriptive only.

  select_c_chroms.py --out data/v30/C_CHROMS.json
"""
import argparse
import glob
import json
import os
import sys
import urllib.request

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
from ncrna_annot.eval import ruler  # noqa: E402

SPECIES = {"rat": "rattus_norvegicus", "cow": "bos_taurus", "xtrop": "xenopus_tropicalis"}
LO, HI, MIN_LNC = 40_000_000, 80_000_000, 200


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="data/v30/C_CHROMS.json")
    a = ap.parse_args()
    res = {"rule": "protocol 18 detail 2: autosome 40-80 Mb with the most ruler-lnc genes; < 200 lnc genes = descriptive only"}
    for short, sp in SPECIES.items():
        url = f"https://rest.ensembl.org/info/assembly/{sp}?content-type=application/json"
        asm = json.load(urllib.request.urlopen(url, timeout=120))
        lens = {r["name"]: r["length"] for r in asm["top_level_region"] if r.get("coord_system") == "chromosome"}
        auto = [c for c in asm.get("karyotype", []) if c.upper() not in ("X", "Y", "Z", "W", "MT") and c in lens]
        gtf = glob.glob(f"{ROOT}/data/v30/ensembl/{short}_*.116.gtf.gz")
        assert len(gtf) == 1, gtf
        cand = {}
        for c in auto:
            if LO <= lens[c] <= HI:
                g = ruler.read_gtf(gtf[0], c)
                cand[c] = {"length": lens[c], "lnc_genes": sum(v["cls"] == "lnc" for v in g.values()),
                           "mrna_genes": sum(v["cls"] == "mrna" for v in g.values())}
        assert cand, f"{sp}: no autosome of 40-80 Mb"
        pick = max(cand, key=lambda c: (cand[c]["lnc_genes"], -int(c) if c.isdigit() else 0))
        res[short] = {"species": sp, "assembly": asm.get("assembly_name"), "gtf": os.path.relpath(gtf[0], ROOT),
                      "chrom": pick, **cand[pick], "in_verdict": cand[pick]["lnc_genes"] >= MIN_LNC, "candidates": cand}
        print(short, pick, cand[pick], flush=True)
    json.dump(res, open(f"{ROOT}/{a.out}", "w"), indent=1)


if __name__ == "__main__":
    main()
