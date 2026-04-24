#!/usr/bin/env python3
"""Build an anonymized CAPITI-E target list for downstream ML.

Input:  capiti_E_splitdb_pdb.tsv (PDB + AlphaFold rows with identity columns).
Output: capiti_E_targets_anonymized.tsv with T-E<N> IDs, only structural ref.

Columns: target_id, db, id, chain, method, resolution, confidence
  - db = PDB or AF
  - id = PDB code or AlphaFold UniProt accession
  - method/resolution fetched from RCSB for PDB rows; empty for AF rows
  - confidence = mean pLDDT for AF rows; empty for PDB rows

Targets are assigned in sorted order of (db, id, chain) for reproducibility.
"""

import csv, json, urllib.request
from pathlib import Path

ROOT = Path(__file__).parent
IN = ROOT / "capiti_E_splitdb_pdb.tsv"
OUT = ROOT / "capiti_E_targets_anonymized.tsv"
RCSB = "https://data.rcsb.org/rest/v1/core/entry"

def fetch_rcsb(pdb_id):
    """Return (method, resolution) from RCSB for a PDB entry."""
    url = f"{RCSB}/{pdb_id.lower()}"
    with urllib.request.urlopen(url, timeout=60) as r:
        data = json.load(r)
    methods = [e.get("method", "") for e in data.get("exptl", []) if e.get("method")]
    method = "; ".join(methods)
    res = data.get("rcsb_entry_info", {}).get("resolution_combined") or []
    resolution = f"{res[0]:.2f}" if res else ""
    return method, resolution

rows = list(csv.DictReader(IN.open(), delimiter="\t"))
rows.sort(key=lambda r: (r["db"], r["id"], r["chain"]))

pdb_ids = sorted({r["id"] for r in rows if r["db"] == "PDB"})
print(f"fetching RCSB metadata for {len(pdb_ids)} PDB entries")
rcsb_cache = {pid: fetch_rcsb(pid) for pid in pdb_ids}

fields = ["target_id", "db", "id", "chain", "method", "resolution", "confidence"]
with OUT.open("w") as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    w.writeheader()
    for i, r in enumerate(rows, start=1):
        method, resolution = rcsb_cache.get(r["id"], ("", "")) if r["db"] == "PDB" else ("", "")
        w.writerow({
            "target_id":  f"T-E{i}",
            "db":         r["db"],
            "id":         r["id"],
            "chain":      r["chain"],
            "method":     method,
            "resolution": resolution,
            "confidence": r.get("confidence", ""),
        })

print(f"targets: {len(rows)}")
print(f"wrote {OUT}")
