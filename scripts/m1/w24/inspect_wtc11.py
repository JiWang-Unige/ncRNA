"""Describe the downloaded WTC11 files (structure only; no prediction is read)."""
import gzip, re, sys, collections, os
D = sys.argv[1]
def op(p): return gzip.open(p, "rt") if p.endswith(".gz") else open(p)
for g in ("ENCFF370SDY", "ENCFF798DWV", "ENCFF422CGK"):
    p = f"{D}/encode/{g}.gtf.gz"
    if not os.path.exists(p): continue
    feat, keys, tx_attr = collections.Counter(), collections.Counter(), collections.Counter()
    chrom_tx = collections.Counter(); ex = []
    with op(p) as fh:
        for line in fh:
            if line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")
            feat[f[2]] += 1
            at = dict(re.findall(r'(\S+) "([^"]*)"', f[8]))
            if f[2] == "transcript":
                chrom_tx[f[0]] += 1
                for k in at: keys[k] += 1
                for k in ("transcript_status", "gene_status", "gene_type", "transcript_type", "ISM_transcript", "NIC_transcript",
                          "NNC_transcript", "antisense_gene", "intergenic_novel", "genomic_transcript", "talon_transcript"):
                    if k in at and k != "talon_transcript": tx_attr[(k, at[k])] += 1
                if len(ex) < 2: ex.append(f[8][:500])
    print("==", g, dict(feat))
    print("  tx attr keys:", keys.most_common(40))
    print("  values:", sorted(tx_attr.items(), key=lambda x: -x[1])[:40])
    print("  chr12/chr18 tx:", chrom_tx["chr12"], chrom_tx["chr18"], "total", sum(chrom_tx.values()))
    print("  example:", ex)
for q in ("ENCFF803XYM",):
    p = f"{D}/encode/{q}.tsv"
    if os.path.exists(p):
        with open(p) as fh:
            print("== quant", q, [next(fh).rstrip()[:400] for _ in range(3)])
for b in ("ENCFF603KSV",):
    with op(f"{D}/encode/{b}.bed.gz") as fh:
        print("== quantseq", [next(fh).rstrip() for _ in range(3)])
for c in ("GSM5625727_WTC11-1_rep1.count.bed.gz", "GSM5625727_WTC11-1_rep1.ctss.bed.gz"):
    with op(f"{D}/geo/{c}") as fh:
        print("==", c, [next(fh).rstrip() for _ in range(3)])
