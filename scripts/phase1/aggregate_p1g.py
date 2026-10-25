"""Aggregate R (P1-S2) vs G1B vs G2-PH on chr18 (protocol v1.6 section 13 rule 1, v1.7 section 14 rules 2-5).
Run from a directory holding outputs/PHASE1/{P1S2,P1G} result.json files; writes p1g_summary.json."""
import json, statistics as st
S=[42,123,456]
def load(arm, s):
    if arm=="R":   d=json.load(open(f"outputs/PHASE1/P1S2/r-{s}/p1s2_r_{s}/result.json")); return d["S3"], d
    if arm=="G1B": d=json.load(open(f"outputs/PHASE1/P1G/ncrna-p1g-g1b-{s}-eval/p1g_g1b_{s}/result.json")); return d["S3"], d
    if arm=="G2PH":d=json.load(open(f"outputs/PHASE1/P1G/ncrna-p1g-g2ph-{s}-eval/p1g_g2ph_{s}/result.json")); return d["PH"], d
f1=lambda x: 0 if x["sn"]+x["pr"]==0 else 2*x["sn"]*x["pr"]/(x["sn"]+x["pr"])
rows={}
for arm in ["R","G1B","G2PH"]:
    rows[arm]={}
    for s in S:
        m,d=load(arm,s)
        rows[arm][s]=dict(auprc=m["auprc_lnc"], f1=m["best_locus_f1"], lncL1sn=m["lnc_L1"]["sn"], lncL1pr=m["lnc_L1"]["pr"],
            L2sn=m["lnc_L2"]["sn"], L2pr=m["lnc_L2"]["pr"], mL1sn=m["mrna_L1"]["sn"], mL1pr=m["mrna_L1"]["pr"], mF1=f1(m["mrna_L1"]),
            npl=m["n_pred_lnc"], npm=m["n_pred_mrna"], ep=d.get("best_epoch"))
print("per seed")
for arm in rows:
    for s in S:
        r=rows[arm][s]; print(f"{arm:5s} {s:4d} AUPRC {r['auprc']:.4f} F1 {r['f1']:.3f} lncL1 Sn {r['lncL1sn']:.3f} Pr {r['lncL1pr']:.3f} L2 Sn {r['L2sn']:.4f} Pr {r['L2pr']:.4f} | mRNA L1 Sn {r['mL1sn']:.3f} Pr {r['mL1pr']:.3f} F1 {r['mF1']:.3f} | nlnc {r['npl']} nmrna {r['npm']} ep {r['ep']}")
mean=lambda arm,k: st.mean(rows[arm][s][k] for s in S)
print("\nmeans"); 
for arm in rows: print(arm, {k: round(mean(arm,k),4) for k in ["auprc","f1","lncL1sn","lncL1pr","L2sn","L2pr","mF1","npl","npm"]})
def delta(a,b,k="auprc"): return [rows[a][s][k]-rows[b][s][k] for s in S]
def same_sign(ds): return all(x>0 for x in ds) or all(x<0 for x in ds)
out={}
d1=delta("G1B","R"); out["rule1"]=dict(delta=d1, mean=st.mean(d1), same_sign=same_sign(d1), pass_=st.mean(d1)>=0.01 and same_sign(d1) and all(x>0 for x in d1))
base="G1B" if out["rule1"]["pass_"] else "R"
d2=delta("G2PH",base); L2ok = mean("G2PH","L2sn")>=0.01 and all(rows["G2PH"][s]["L2sn"]>rows[base][s]["L2sn"] for s in S)
out["rule2"]=dict(baseline=base, delta=d2, mean=st.mean(d2), same_sign=same_sign(d2), L2sn_mean=mean("G2PH","L2sn"), L2_each_above=[rows["G2PH"][s]["L2sn"]>rows[base][s]["L2sn"] for s in S], pass_=st.mean(d2)>=0.02 and same_sign(d2) and all(x>0 for x in d2) and L2ok)
out["rule3_mrna"]=dict(G2PH_mF1=mean("G2PH","mF1"), R_mF1=mean("R","mF1"), G1B_mF1=mean("G1B","mF1"), drop_vs_R=mean("R","mF1")-mean("G2PH","mF1"), G1B_drop_vs_R=mean("R","mF1")-mean("G1B","mF1"))
dPR=delta("G2PH","R"); out["rule4_g2ph_vs_R"]=dict(delta=dPR, mean=st.mean(dPR), same_sign=same_sign(dPR))
print("\nrules"); print(json.dumps(out, indent=1))
json.dump(dict(rows=rows, rules=out), open("p1g_summary.json","w"), indent=1)
