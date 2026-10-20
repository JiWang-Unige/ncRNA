"""Phase 0.6: small-ncRNA scan (Infernal/Rfam + tRNAscan-SE) on Modal CPU containers.

Volume `ncrna-annot`: /vol/data/rfam/{Rfam.cm*, Rfam.clanin}, /vol/data/<species>/<chrom>.fa
Each chromosome is cut into 1 Mb chunks (10 kb overlap) scanned in parallel with .map.
Run from ont:
  MODAL_PROFILE=ji-wang-unige ~/envs/modal/bin/modal run -d smallrna_modal.py \
      --fasta data/worm/I.fa --tag PHASE0/0.6_smallrna/worm_I [--max-mb 2]   # --max-mb = smoke
"""
import json
import os
import subprocess
import tempfile
import time

import modal

image = modal.Image.micromamba(python_version="3.11").micromamba_install(
    "infernal=1.1.5", "trnascan-se=2.0.12", channels=["conda-forge", "bioconda"])
app = modal.App("ncrna-smallrna", image=image)
vol = modal.Volume.from_name("ncrna-annot", create_if_missing=True)
CHUNK, OVERLAP = 1_000_000, 10_000
COLS = ["chrom", "start", "end", "strand", "model", "rfam_acc", "clan", "score", "evalue", "olp", "trunc"]


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
            else:
                seq.append(line.strip())
    return name, "".join(seq)


@app.function(cpu=4, memory=8192, timeout=6 * 3600, volumes={"/vol": vol})
def cmscan_chunk(fasta: str, start: int, end: int):
    chrom, seq = read_fasta(f"/vol/{fasta}")
    sub = seq[start:end]
    t0 = time.time()
    with tempfile.TemporaryDirectory() as td:
        open(f"{td}/q.fa", "w").write(f">{chrom}\n{sub}\n")
        subprocess.run(["cmscan", "--rfam", "--cut_ga", "--nohmmonly", "--oclan", "--fmt", "2",
                        "--clanin", "/vol/data/rfam/Rfam.clanin", "--cpu", "4",
                        "--tblout", f"{td}/t.tbl", "-o", "/dev/null",
                        "/vol/data/rfam/Rfam.cm", f"{td}/q.fa"], check=True)
        rows = []
        for line in open(f"{td}/t.tbl"):
            if line.startswith("#"):
                continue
            f = line.split(None, 26)
            # fmt 2: idx tname tacc qname qacc clan mdl mdlfrom mdlto seqfrom seqto strand trunc
            #        pass gc bias score evalue inc olp ...
            s, e = int(f[9]), int(f[10])
            lo, hi = min(s, e) + start, max(s, e) + start
            rows.append([chrom, lo, hi, f[11], f[1], f[2], f[5], float(f[16]), f[17], f[19], f[12]])
    return {"start": start, "end": end, "rows": rows, "secs": round(time.time() - t0, 1)}


@app.function(cpu=4, memory=8192, timeout=6 * 3600, volumes={"/vol": vol})
def trnascan(fasta: str, max_bp: int = 0):
    chrom, seq = read_fasta(f"/vol/{fasta}")
    if max_bp:
        seq = seq[:max_bp]
    t0 = time.time()
    with tempfile.TemporaryDirectory() as td:
        open(f"{td}/q.fa", "w").write(f">{chrom}\n{seq}\n")
        subprocess.run(["tRNAscan-SE", "-E", "--thread", "4", "-q", "-o", f"{td}/out.txt",
                        f"{td}/q.fa"], check=True)
        txt = open(f"{td}/out.txt").read()
    return {"txt": txt, "secs": round(time.time() - t0, 1)}


@app.function(timeout=20 * 3600, volumes={"/vol": vol})
def run(fasta: str, tag: str, max_mb: int = 0):
    out = f"/vol/outputs/{tag}"
    os.makedirs(out, exist_ok=True)
    _, seq = read_fasta(f"/vol/{fasta}")
    n = min(len(seq), max_mb * 1_000_000) if max_mb else len(seq)
    spans = [(s, min(s + CHUNK + OVERLAP, n)) for s in range(0, n, CHUNK)]
    t0 = time.time()
    trna = trnascan.spawn(fasta, n if max_mb else 0)
    res = list(cmscan_chunk.map([fasta] * len(spans), [s for s, _ in spans], [e for _, e in spans]))
    # drop overlap duplicates: same model + same coordinates seen from two chunks
    seen, rows = set(), []
    for r in sorted((x for c in res for x in c["rows"]), key=lambda x: (x[1], -x[7])):
        k = (r[4], r[1], r[2], r[3])
        if k not in seen:
            seen.add(k)
            rows.append(r)
    with open(f"{out}/cmscan.tsv", "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    t = trna.get()
    open(f"{out}/trnascan.txt", "w").write(t["txt"])
    json.dump({"fasta": fasta, "bp_scanned": n, "n_chunks": len(spans),
               "chunk_secs": [c["secs"] for c in res], "cmscan_cpu_h": round(sum(c["secs"] for c in res) * 4 / 3600, 2),
               "trnascan_secs": t["secs"], "wall_s": round(time.time() - t0, 1), "n_cmscan_hits": len(rows)},
              open(f"{out}/run.json", "w"), indent=1)
    vol.commit()
    return len(rows)


@app.local_entrypoint()
def main(fasta: str, tag: str, max_mb: int = 0):
    print("hits", run.remote(fasta, tag, max_mb))
