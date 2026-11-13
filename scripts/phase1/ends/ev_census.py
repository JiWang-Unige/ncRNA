#!/usr/bin/env python3
"""PLAN_END.md section 0: census of end labels in the v2.0 training windows and their CAGE / polyA support (CPU, ont)."""
import gzip, json, sys, os, numpy as np
ROOT=os.path.expanduser("~/ncrna-annot"); sys.path.insert(0,ROOT+"/scripts/phase1/s3"); sys.path.insert(0,ROOT+"/src")
import s3_features as S3
from collections import defaultdict
WIN=24576
wins=[l.split() for l in open(ROOT+"/data/p1v20/train_windows.tsv")]
wins=[(c,int(s)) for c,s,*_ in wins]
print("windows",len(wins), "chroms", len({c for c,_ in wins}))
# canonical transcripts v50
tx=defaultdict(lambda: {"ex":[]})
with gzip.open(ROOT+"/data/human/gencode.v50.annotation.gtf.gz","rt") as fh:
    for line in fh:
        if line[0]=="#": continue
        f=line.split("\t")
        if f[2]!="exon" or "Ensembl_canonical" not in f[8]: continue
        gt=f[8].split('gene_type "')[1].split('"')[0]
        if gt not in ("lncRNA","protein_coding"): continue
        tid=f[8].split('transcript_id "')[1].split('"')[0]
        t=tx[tid]; t["c"]=f[0]; t["s"]=f[6]; t["k"]="lnc" if gt=="lncRNA" else "mrna"; t["ex"].append((int(f[3]),int(f[4])))
bychr=defaultdict(list)
for t in tx.values():
    lo=min(a for a,_ in t["ex"]); hi=max(b for _,b in t["ex"])
    tss,tes=(lo-1,hi-1) if t["s"]=="+" else (hi-1,lo-1)
    bychr[t["c"]].append((t["k"],t["s"],tss,tes,hi-lo+1))
cage={}; pa={}
def peaks(c):
    if c not in cage:
        cage[c]=S3.load_peaks(ROOT+"/data/bioval/raw/fantom5_cage_peaks.bed.gz",c,5,6,7)
        pa[c]=S3.load_peaks(ROOT+"/data/bioval/raw/polyasite2_GRCh38.96.tsv.gz",c,5,1,2,skip_header=True,chrom_name=c.replace("chr",""))
    return cage[c],pa[c]
cnt=defaultdict(int); spans=defaultdict(list)
wbychr=defaultdict(list)
for c,s in wins: wbychr[c].append(s)
for c,ws in wbychr.items():
    cg,pp=peaks(c)
    for k,st,tss,tes,L in bychr.get(c,[]):
        intss=any(s<=tss<s+WIN for s in ws); intes=any(s<=tes<s+WIN for s in ws)
        both=any(s<=min(tss,tes) and max(tss,tes)<s+WIN for s in ws)
        if intss:
            cnt[(k,"tss")]+=1; cnt[(k,"tss_cage")]+=S3.near(cg[st],tss,200)
        if intes:
            cnt[(k,"tes")]+=1; cnt[(k,"tes_polya")]+=S3.near(pp[st],tes,100)
        if intss or intes:
            cnt[(k,"any_end")]+=1; cnt[(k,"both_in_one_window")]+=both; spans[k].append(L)
for k in ("lnc","mrna"):
    print(k,{x:cnt[(k,x)] for x in ("tss","tss_cage","tes","tes_polya","any_end","both_in_one_window")}, "median span", int(np.median(spans[k])), "frac span>24.6kb", round(float(np.mean(np.array(spans[k])>WIN)),3))
ev=np.load(ROOT+"/data/p1v20/train_events.npy",mmap_mode="r")
print("events shape",ev.shape)
for ch,n in ((2,"TSS+"),(3,"TES+"),(6,"TSS-"),(7,"TES-")):
    x=np.asarray(ev[:,:,ch]); print(n,"pos bases",int((x==1).sum()),"ign",int((x==-1).sum()),"neg",int((x==0).sum()))
# whole chromosome-wide fraction of canonical lnc TSS with CAGE (for reference)
tot=defaultdict(int)
for c in ["chr%d"%i for i in range(1,23)]+["chrX"]:
    if c in ("chr18","chr20","chr9"): continue
    cg,pp=peaks(c)
    for k,st,tss,tes,L in bychr.get(c,[]):
        tot[(k,"n")]+=1; tot[(k,"cage")]+=S3.near(cg[st],tss,200); tot[(k,"pa")]+=S3.near(pp[st],tes,100); tot[(k,"both")]+=S3.near(cg[st],tss,200) and S3.near(pp[st],tes,100)
print({f"{a}_{b}":v for (a,b),v in tot.items()})
