#!/usr/bin/env python
"""D1 step 0: inspect TSP33 total-RNA h5ad structure (metadata only; no localisation statistics)."""
import anndata as ad, numpy as np, sys
a = ad.read_h5ad(sys.argv[1], backed="r")
print(a)
print("obs cols:", list(a.obs.columns))
print("var cols:", list(a.var.columns)); print(a.var.head())
for c in a.obs.columns:
    if a.obs[c].dtype.name in ("category", "object") and a.obs[c].nunique() < 200:
        print("\n==", c, a.obs[c].nunique()); print(a.obs[c].value_counts().head(30).to_string())
print("layers:", list(a.layers.keys()) if a.layers is not None else None)
x = a.X[:5]; print(type(x), getattr(x, "dtype", None))
try:
    print(x.toarray()[:, :20] if hasattr(x, "toarray") else x[:, :20])
except Exception as e: print(e)
