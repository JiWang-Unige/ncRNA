#!/usr/bin/env python
"""POST-HOC (not part of the frozen gates): how much of HGB_E's gain over ref_only comes from a non-linear
function of the reference reactivity and reference RPKM alone (regression to the mean + coverage), versus
sequence / folding features?  Same GroupKFold(5) and HGB settings as b0_lite.py."""
import json, os, sys
import numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import GroupKFold
sys.path.insert(0, os.path.dirname(__file__))
import b0_lite as B
task = sys.argv[1]
df = pd.read_pickle(f"{B.OUT}/{task}.pkl"); y = df.delta.values
XE = B.featurize(df, True)
sets = {"ref_rpkm_only": np.c_[df.r_ref.values, np.log1p(df.rpkm_ref.values)],
        "ref_rpkm_tlen_type": np.c_[df.r_ref.values, np.log1p(df.rpkm_ref.values), XE[:, -5 - 2:-2][:, -4:]],
        "ref_plus_unpaired_gc": np.c_[df.r_ref.values, np.log1p(df.rpkm_ref.values), df.unpaired.fillna(0.5).values, df.gc_cnt.values]}
res = {}
for name, X in sets.items():
    p = np.zeros(len(df))
    for tr, te in GroupKFold(5).split(df, groups=df.gid.values):
        p[te] = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.05, random_state=0).fit(X[tr], y[tr]).predict(X[te])
    res[name] = float(np.corrcoef(y, p)[0, 1])
d = pd.read_pickle(f"{B.OUT}/{task}_preds.pkl")
res["ref_only_linear"] = float(np.corrcoef(d.y, d.ref_only)[0, 1]); res["HGB_E_full"] = float(np.corrcoef(d.y, d.HGB_E)[0, 1])
json.dump(res, open(f"{B.OUT}/{task}_posthoc_ablation.json", "w"), indent=1); print(task, json.dumps(res))
