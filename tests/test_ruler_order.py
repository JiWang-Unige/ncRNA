import json
import os
import subprocess
import sys

# L2 bootstrap must not depend on Python's string-hash seed (the gene vector used to follow set order)
_SCRIPT = """
import json
from ncrna_annot.eval import ruler
truth, preds = {}, []
for i in range(60):
    s = 10000 * i
    ex = [(s, s + 100), (s + 300, s + 400), (s + 700, s + 800)]
    truth["ENSG%011d" % i] = {"strand": "+", "cls": "lnc", "gene_type": "lncRNA", "tx": {"t%d" % i: ex}}
    if i % 3 == 0:
        preds.append({"strand": "+", "cls": "lnc", "exons": ex})
r = ruler.score_class("lnc", truth, preds, ruler.ExonIndex(truth), 200)
print(json.dumps(r["L2_intron_chain"]["sn_ci95"]))
"""


def test_l2_bootstrap_independent_of_hash_seed():
    src = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src")
    out = set()
    for seed in ("0", "1", "2", "3"):
        env = dict(os.environ, PYTHONHASHSEED=seed, PYTHONPATH=src + os.pathsep + os.environ.get("PYTHONPATH", ""))
        res = subprocess.run([sys.executable, "-c", _SCRIPT], env=env, capture_output=True, text=True, check=True)
        out.add(json.dumps(json.loads(res.stdout)))
    assert len(out) == 1, out
