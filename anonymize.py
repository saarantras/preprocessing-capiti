#!/usr/bin/env python3
"""Build an anonymized target list for downstream ML.

Keeps only rows with at least one PDB entry. Replaces all CARD/UniProt
identity fields with a target ID of the form T-C<N>. Emits one row per
(target, PDB) pair, with chain info for structure retrieval.

Output columns:
    target_id, PDB, chains, method, resolution
"""

import csv, json, urllib.parse, urllib.request
from collections import OrderedDict
from pathlib import Path

ROOT = Path(__file__).parent
IN = ROOT / "capiti_C_targets_list.tsv"
OUT = ROOT / "capiti_C_targets_anonymized.tsv"
UNIPROT = "https://rest.uniprot.org/uniprotkb"

def _get(url):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    return urllib.request.urlopen(req, timeout=60)

def fetch_pdb_properties(uniprot_accs):
    """Return {uniprot: {pdb_id: {chains, method, resolution}}}."""
    result = {}
    accs = sorted(set(a for a in uniprot_accs if a))
    CHUNK = 50
    import re
    NEXT = re.compile(r'<([^>]+)>\s*;\s*rel="next"')
    for i in range(0, len(accs), CHUNK):
        chunk = accs[i:i+CHUNK]
        q = " OR ".join(f"accession:{a}" for a in chunk)
        url = (f"{UNIPROT}/search?query={urllib.parse.quote(q)}"
               f"&fields=accession,xref_pdb&format=json&size=500")
        while url:
            with _get(url) as r:
                p = json.load(r)
                for e in p.get("results", []):
                    up = e["primaryAccession"]
                    per_up = result.setdefault(up, {})
                    for x in e.get("uniProtKBCrossReferences", []):
                        if x.get("database") != "PDB":
                            continue
                        props = {kv["key"]: kv["value"]
                                 for kv in x.get("properties", [])}
                        per_up[x["id"]] = {
                            "chains":     props.get("Chains", ""),
                            "method":     props.get("Method", ""),
                            "resolution": props.get("Resolution", ""),
                        }
                m = NEXT.search(r.headers.get("Link", "") or "")
                url = m.group(1) if m else None
    return result

# ---- Load and filter rows with PDB
rows = list(csv.DictReader(IN.open(), delimiter="\t"))
rows_with_pdb = [r for r in rows if r.get("PDB")]
rows_with_pdb.sort(key=lambda r: r["ARO Accession"])  # deterministic ordering
print(f"{len(rows_with_pdb)} rows have at least one PDB")

# ---- Gather UniProts we need chain info for
uniprots = set()
for r in rows_with_pdb:
    uniprots.update(u for u in r["UniProt"].split(";") if u)
print(f"fetching PDB properties for {len(uniprots)} UniProt accessions")
props_by_up = fetch_pdb_properties(uniprots)

def chains_for(r, pdb):
    """Look up chain info across any UniProt linked to this row."""
    for up in r["UniProt"].split(";"):
        info = props_by_up.get(up, {}).get(pdb)
        if info:
            return info
    return {"chains": "", "method": "", "resolution": ""}

# ---- Emit anonymized long table
fields = ["target_id", "PDB", "chains", "method", "resolution"]
n_out = 0
with OUT.open("w") as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    w.writeheader()
    for i, r in enumerate(rows_with_pdb, start=1):
        tid = f"T-C{i}"
        for pdb in r["PDB"].split(";"):
            if not pdb:
                continue
            info = chains_for(r, pdb)
            w.writerow({"target_id": tid, "PDB": pdb, **info})
            n_out += 1

print(f"targets: {len(rows_with_pdb)}")
print(f"rows emitted: {n_out}")
print(f"wrote {OUT}")
