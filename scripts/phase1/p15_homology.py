#!/usr/bin/env python3
"""Protocol section 1, homology isolation for the chr20 read (truth and genome only; no prediction is read).

Every merged exon (>= 20 nt) of the GENCODE v50 lncRNA genes on chr20 is searched with BLASTn (blastn 2.16+,
-task blastn, e < 1e-10) against the merged exons of the v50 lncRNA genes on the training chromosomes
(chr1-17, 19, 21, 22, X).  An exon hits when one HSP has identity >= 80 % and covers >= 50 % of the exon; a chr20
gene with a hitting exon is "homologous".  Classes, strands and merged exons come from ruler.read_gtf.
Two query variants, both fixed before the read:
  masked (primary)     hg38 soft-masked repeats (lower case) cannot seed a hit (-lcase_masking; DUST on), so a
                       shared transposon fragment alone does not make a gene homologous
  unmasked (sensitivity)  plain sequence, DUST off
  p15_homology.py --out outputs/PHASE1/P15/homology_chr20.json [--threads 6] [--work <dir>]
"""
import argparse
import gzip
import json
import os
import subprocess
import sys
from collections import defaultdict

ROOT = os.path.expanduser(os.environ.get("NCRNA_ROOT", "~/ncrna-annot"))
sys.path.insert(0, os.path.join(ROOT, "src"))
GTF = f"{ROOT}/data/human/gencode.v50.annotation.gtf.gz"
FASTA = f"{ROOT}/data/human/hg38.fa.gz"
BLAST = os.path.expanduser("~/envs/earlgrey/bin")
TRAIN = [f"chr{i}" for i in list(range(1, 18)) + [19, 21, 22]] + ["chrX"]
TEST = "chr20"
EVALUE, PIDENT, COVER, MIN_EXON = 1e-10, 80.0, 0.5, 20


def split_gtf(work, chroms):
    """One pass over the GENCODE GTF -> per-chromosome GTFs (so ruler.read_gtf reads each quickly)."""
    fh = {c: open(f"{work}/{c}.gtf", "w") for c in chroms}
    with gzip.open(GTF, "rt") as g:
        for line in g:
            c = line.split("\t", 1)[0]
            if c in fh:
                fh[c].write(line)
    for f in fh.values():
        f.close()


def read_chroms(chroms):
    seqs, cur, buf = {}, None, []
    with gzip.open(FASTA, "rt") as g:
        for line in g:
            if line.startswith(">"):
                if cur in chroms:
                    seqs[cur] = "".join(buf)
                cur, buf = line[1:].split()[0], []
            elif cur in chroms:
                buf.append(line.rstrip("\n"))
    if cur in chroms:
        seqs[cur] = "".join(buf)
    return seqs


def lnc_exons(work, chrom):
    from ncrna_annot.eval import ruler
    genes = ruler.read_gtf(f"{work}/{chrom}.gtf", chrom)
    return {gid: (g["strand"], ruler._merge(sorted(x for t in g["tx"].values() for x in t)))
            for gid, g in genes.items() if g["cls"] == "lnc"}


def write_fasta(path, records, seq, upper):
    with open(path, "a") as fh:
        for rid, a, b in records:
            s = seq[a - 1:b]
            fh.write(f">{rid}\n{s.upper() if upper else s}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--work", default=f"{ROOT}/outputs/PHASE1/P15/homology_work")
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)
    split_gtf(a.work, TRAIN + [TEST])
    seqs = read_chroms(set(TRAIN + [TEST]))
    test = lnc_exons(a.work, TEST)
    q = [(f"{gid}|{i}|{e - s + 1}", s, e) for gid, (_, ex) in test.items() for i, (s, e) in enumerate(ex)
         if e - s + 1 >= MIN_EXON]
    for f in ("query_masked.fa", "query_unmasked.fa", "db.fa"):
        open(f"{a.work}/{f}", "w").close()
    write_fasta(f"{a.work}/query_masked.fa", q, seqs[TEST], upper=False)
    write_fasta(f"{a.work}/query_unmasked.fa", q, seqs[TEST], upper=True)
    n_db = 0
    for c in TRAIN:
        rec = [(f"{c}|{gid}|{i}", s, e) for gid, (_, ex) in lnc_exons(a.work, c).items() for i, (s, e) in enumerate(ex)]
        write_fasta(f"{a.work}/db.fa", rec, seqs[c], upper=True)
        n_db += len(rec)
    subprocess.run([f"{BLAST}/makeblastdb", "-in", f"{a.work}/db.fa", "-dbtype", "nucl", "-out", f"{a.work}/db"],
                   check=True, stdout=subprocess.DEVNULL)
    res = {"params": {"blastn": "2.16.0+", "task": "blastn", "evalue": EVALUE, "pident": PIDENT, "cover": COVER,
                      "min_exon": MIN_EXON, "train": TRAIN, "test": TEST, "truth": "GENCODE v50 (ruler.read_gtf, lnc)",
                      "genome": "hg38.fa.gz (UCSC, soft-masked)"},
           "n_lnc_genes": len(test), "n_query_exons": len(q), "n_db_exons": n_db, "primary": "masked", "variants": {}}
    for var, extra in (("masked", ["-lcase_masking"]), ("unmasked", ["-dust", "no"])):
        tsv = f"{a.work}/hits_{var}.tsv"
        subprocess.run([f"{BLAST}/blastn", "-task", "blastn", "-query", f"{a.work}/query_{var}.fa", "-db", f"{a.work}/db",
                        "-evalue", str(EVALUE), "-outfmt", "6 qseqid sseqid pident length qstart qend qlen evalue",
                        "-max_target_seqs", "50", "-max_hsps", "5", "-num_threads", str(a.threads), "-out", tsv] + extra,
                       check=True)
        hom, n_hsp = defaultdict(int), 0
        for line in open(tsv):
            qid, _, pid, _, qs, qe, ql, ev = line.rstrip("\n").split("\t")
            n_hsp += 1
            if float(ev) < EVALUE and float(pid) >= PIDENT and (int(qe) - int(qs) + 1) >= COVER * int(ql):
                hom[qid.split("|")[0]] += 1
        res["variants"][var] = {"n_hsp": n_hsp, "n_homologous": len(hom), "homologous": sorted(hom)}
        print(json.dumps({var: {"n_hsp": n_hsp, "n_homologous": len(hom), "of": len(test)}}), flush=True)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    json.dump(res, open(a.out, "w"), indent=1)


if __name__ == "__main__":
    main()
