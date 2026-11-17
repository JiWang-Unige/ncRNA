"""Stdlib readers for the FANTOM5 sample table (xlsx) and the FF ontology (obo)."""
import re
import zipfile
import xml.etree.ElementTree as ET

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def col_index(ref):
    letters = re.match(r"[A-Z]+", ref).group(0)
    n = 0
    for ch in letters:
        n = n * 26 + ord(ch) - 64
    return n - 1


def read_xlsx(path, sheet=1):
    z = zipfile.ZipFile(path)
    ss = []
    if "xl/sharedStrings.xml" in z.namelist():
        for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", NS):
            ss.append("".join(t.text or "" for t in si.iter(f"{{{NS['m']}}}t")))
    root = ET.fromstring(z.read(f"xl/worksheets/sheet{sheet}.xml"))
    rows = []
    for row in root.iter(f"{{{NS['m']}}}row"):
        vals = {}
        for c in row.findall("m:c", NS):
            v = c.find("m:v", NS)
            if v is None:
                isel = c.find("m:is", NS)
                txt = "".join(t.text or "" for t in isel.iter(f"{{{NS['m']}}}t")) if isel is not None else None
            else:
                txt = ss[int(v.text)] if c.get("t") == "s" else v.text
            vals[col_index(c.get("r"))] = txt
        if vals:
            rows.append([vals.get(i) for i in range(max(vals) + 1)])
    return rows


def read_obo(path):
    terms = {}
    for block in open(path, encoding="utf-8").read().split("\n[")[1:]:
        if not block.startswith("Term]"):
            continue
        t = {"is_a": [], "rel": [], "xref": [], "isect": []}
        for line in block.split("\n")[1:]:
            if ": " not in line:
                continue
            k, v = line.split(": ", 1)
            v = v.split(" ! ")[0].strip()
            if k == "id":
                t["id"] = v
            elif k == "name":
                t["name"] = v
            elif k == "is_a":
                t["is_a"].append(v.split()[0])
            elif k == "relationship":
                t["rel"].append(tuple(v.split()[:2]))
            elif k == "intersection_of":
                t["isect"].append(v)
            elif k == "xref":
                t["xref"].append(v.split()[0])
        if "id" in t:
            terms[t["id"]] = t
    return terms
