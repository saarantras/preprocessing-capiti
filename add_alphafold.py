#!/usr/bin/env python3
"""Rescue unresolved capiti_E rows via AlphaFold.

For each unresolved id, resolve to a UniProt accession (if needed), query the
AlphaFold EBI API, and read `globalMetricValue` (mean pLDDT). AlphaFold
confidence bands: <50 very low, 50-70 low, 70-90 confident, >90 very high.
Anything with mean pLDDT < 50 ("very low") is dropped.

Effects:
  - Appends kept rows (db=AF) to capiti_E_splitdb_pdb.tsv.
  - Rewrites capiti_E_unresolved.tsv to contain only rows that are still
    unresolvable (no UniProt mapping, no AF entry, or mean pLDDT < 50).
  - Adds a `confidence` column to the main output (blank for PDB rows).
"""

import csv, json, re, time, urllib.parse, urllib.request
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).parent
MAIN  = ROOT / "capiti_E_splitdb_pdb.tsv"
UNRES = ROOT / "capiti_E_unresolved.tsv"

PLDDT_DROP = 50.0  # drop if mean pLDDT < 50 ("very low confidence")

API  = "https://rest.uniprot.org/idmapping"
AFDB = "https://alphafold.ebi.ac.uk/api/prediction"

UNIPROT_RE = re.compile(
    r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$")

def _get(url, timeout=60):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    return urllib.request.urlopen(req, timeout=timeout)

def _post(url, data):
    req = urllib.request.Request(url,
        data=urllib.parse.urlencode(data).encode(), method="POST")
    return urllib.request.urlopen(req, timeout=60)

def classify(acc):
    bare = acc.split(".")[0]
    if UNIPROT_RE.match(bare): return "uniprot", bare
    if bare.startswith(("WP_", "NP_")): return "refseq", acc
    return "genbank", acc

def idmap(from_db, ids):
    if not ids: return {}
    job = json.load(_post(f"{API}/run",
        {"from": from_db, "to": "UniProtKB", "ids": ",".join(ids)})
    )["jobId"]
    while True:
        s = json.load(_get(f"{API}/status/{job}"))
        st = s.get("jobStatus")
        if st in (None, "FINISHED"): break
        if st == "ERROR": raise RuntimeError(f"job {job} failed")
        time.sleep(2)
    url = f"{API}/uniprotkb/results/{job}?format=json&size=500&fields=accession"
    out = defaultdict(list)
    NEXT = re.compile(r'<([^>]+)>\s*;\s*rel="next"')
    while url:
        with _get(url) as r:
            payload = json.load(r)
            for rec in payload.get("results", []):
                out[rec["from"]].append(rec["to"].get("primaryAccession"))
            m = NEXT.search(r.headers.get("Link","") or "")
            url = m.group(1) if m else None
    return dict(out)

def best_af(uniprot):
    """Return (plddt, entry_id) for highest-confidence AF model, or None."""
    try:
        with _get(f"{AFDB}/{uniprot}") as r:
            preds = json.load(r)
    except Exception:
        return None
    if not preds: return None
    best = max(preds, key=lambda p: p.get("globalMetricValue") or 0)
    pl = best.get("globalMetricValue")
    if pl is None: return None
    return (float(pl), best.get("modelEntityId", f"AF-{uniprot}-F1"))

# ---- Load current files
unres_rows = list(csv.DictReader(UNRES.open(), delimiter="\t"))
main_rows  = list(csv.DictReader(MAIN.open(),  delimiter="\t"))
print(f"unresolved rows: {len(unres_rows)}")

# ---- Resolve to UniProt
kinds = defaultdict(list)
for r in unres_rows:
    r["id"] = r["id"].strip()
    k, norm = classify(r["id"])
    r["_kind"] = k; r["_norm"] = norm
    kinds[k].append(r)
refseq_ids  = sorted({r["_norm"] for r in kinds["refseq"]})
genbank_ids = sorted({r["_norm"] for r in kinds["genbank"]})
r_to_up = idmap("RefSeq_Protein", refseq_ids) if refseq_ids else {}
g_to_up = idmap("EMBL-GenBank-DDBJ_CDS", genbank_ids) if genbank_ids else {}

def uniprots_for(row):
    if row["_kind"] == "uniprot": return [row["_norm"]]
    if row["_kind"] == "refseq":  return r_to_up.get(row["_norm"], [])
    return g_to_up.get(row["_norm"], [])

# ---- Query AlphaFold per unresolved row
rescued_af = []           # rows to append with db=AF
still_unresolved = []     # rows that remain unresolvable
for r in unres_rows:
    ups = uniprots_for(r)
    best = None
    best_up = None
    for up in ups:
        res = best_af(up)
        if res is None: continue
        if best is None or res[0] > best[0]:
            best = res; best_up = up
    if best is None:
        still_unresolved.append({**{k: r[k] for k in ("","chain","db","id")},
                                 "note": "no PDB, no AlphaFold entry"})
        continue
    plddt, model_id = best
    if plddt < PLDDT_DROP:
        still_unresolved.append({**{k: r[k] for k in ("","chain","db","id")},
                                 "note": f"AF mean pLDDT {plddt:.1f} < {PLDDT_DROP} "
                                         f"(very low); dropped"})
        continue
    rescued_af.append({
        "":            r[""],
        "chain":       "A",      # AlphaFold monomers are single-chain A
        "db":          "AF",
        "id":          best_up,
        "source":      "alphafold",
        "original_id": r["id"],
        "confidence":  f"{plddt:.2f}",
        "note":        f"AF model {model_id}",
    })

# ---- Rewrite main output with new schema (adds `confidence`) + AF rows
fields = ["", "chain", "db", "id", "source", "original_id", "confidence", "note"]
with MAIN.open("w") as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    w.writeheader()
    for r in main_rows:
        r.setdefault("confidence", "")
        w.writerow({k: r.get(k, "") for k in fields})
    for r in rescued_af:
        w.writerow(r)

# ---- Rewrite unresolved file
if still_unresolved:
    with UNRES.open("w") as f:
        w = csv.DictWriter(f, fieldnames=["","chain","db","id","note"],
                           delimiter="\t")
        w.writeheader()
        for r in still_unresolved:
            w.writerow(r)
else:
    UNRES.unlink(missing_ok=True)

# ---- Report
from collections import Counter
bands = Counter()
for r in rescued_af:
    pl = float(r["confidence"])
    if   pl >= 90: bands["very high (>=90)"] += 1
    elif pl >= 70: bands["confident (70-90)"] += 1
    else:          bands["low (50-70)"] += 1
print(f"\nrescued via AlphaFold: {len(rescued_af)}")
for b, n in bands.most_common():
    print(f"  {b}: {n}")
print(f"still unresolved: {len(still_unresolved)}")
print(f"wrote {MAIN}")
if still_unresolved:
    print(f"wrote {UNRES}")
else:
    print("(all unresolved entries resolved; capiti_E_unresolved.tsv removed)")
