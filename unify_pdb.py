#!/usr/bin/env python3
"""Convert rows in capiti_E_splitdb.tsv where db != PDB into uniform PDB rows.

Strategy: for each non-PDB id, resolve via UniProt ID-mapping (when needed) to
get PDB cross-references. Prefer a PDB whose Chains property contains the
row's existing chain letter; otherwise fall back to the best-resolution PDB.
Rows with no resolvable PDB are written to capiti_E_unresolved.tsv.
"""

import csv, json, re, time, urllib.parse, urllib.request
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).parent
IN  = ROOT / "capiti_E_splitdb.tsv"
OUT = ROOT / "capiti_E_splitdb_pdb.tsv"
UNRES = ROOT / "capiti_E_unresolved.tsv"

API = "https://rest.uniprot.org/idmapping"
UNIPROT = "https://rest.uniprot.org/uniprotkb"

def _get(url):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    return urllib.request.urlopen(req, timeout=60)

def _post(url, data):
    req = urllib.request.Request(url,
        data=urllib.parse.urlencode(data).encode(), method="POST")
    return urllib.request.urlopen(req, timeout=60)

UNIPROT_ACC_RE = re.compile(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$"
                            r"|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$")

def classify(acc):
    bare = acc.split(".")[0]
    if UNIPROT_ACC_RE.match(bare):
        return "uniprot", bare
    if bare.startswith(("WP_", "NP_")):
        return "refseq", acc
    return "genbank", acc

def idmap(from_db, ids):
    """Map NCBI-style ids to UniProt primary accessions."""
    if not ids: return {}
    job = json.load(_post(f"{API}/run",
        {"from": from_db, "to": "UniProtKB", "ids": ",".join(ids)})
    )["jobId"]
    while True:
        s = json.load(_get(f"{API}/status/{job}"))
        if s.get("jobStatus") in (None, "FINISHED"): break
        if s.get("jobStatus") == "ERROR":
            raise RuntimeError(f"job {job} failed")
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

def fetch_pdb_xrefs(uniprots):
    """UniProt -> list of {pdb, chains_str, method, resolution}."""
    out = defaultdict(list)
    accs = sorted(set(a for a in uniprots if a))
    if not accs: return out
    CHUNK = 50
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
                    for x in e.get("uniProtKBCrossReferences", []):
                        if x.get("database") != "PDB": continue
                        props = {kv["key"]: kv["value"]
                                 for kv in x.get("properties", [])}
                        out[up].append({
                            "pdb":        x["id"],
                            "chains_str": props.get("Chains", ""),
                            "method":     props.get("Method", ""),
                            "resolution": props.get("Resolution", ""),
                        })
                m = NEXT.search(r.headers.get("Link","") or "")
                url = m.group(1) if m else None
    return out

def chains_of(chains_str):
    """'A/B/C=1-178,D=10-200' -> {'A','B','C','D'}."""
    if not chains_str: return set()
    letters = set()
    for seg in chains_str.split(","):
        prefix = seg.split("=")[0]
        letters.update(c.strip() for c in prefix.split("/") if c.strip())
    return letters

def resolution_float(s):
    if not s: return float("inf")
    try: return float(s.split()[0])
    except Exception: return float("inf")

def pick_pdb(pdbs, wanted_chain):
    """Prefer PDB containing wanted_chain; tie-break by best resolution."""
    if not pdbs: return None
    with_chain = [p for p in pdbs if wanted_chain in chains_of(p["chains_str"])]
    pool = with_chain or pdbs
    pool = sorted(pool, key=lambda p: resolution_float(p["resolution"]))
    return pool[0]

# ---- Load and classify
rows = list(csv.DictReader(IN.open(), delimiter="\t"))
# Strip surrounding whitespace in id and chain
for r in rows:
    r["id"] = (r.get("id") or "").strip()
    r["chain"] = (r.get("chain") or "").strip()
    r["db"] = (r.get("db") or "").strip()

non_pdb = [r for r in rows if r["db"] != "PDB"]
print(f"total rows: {len(rows)}; non-PDB rows: {len(non_pdb)}")

by_kind = defaultdict(list)
for r in non_pdb:
    kind, norm = classify(r["id"])
    r["_kind"] = kind
    r["_norm"] = norm
    by_kind[kind].append(r)
for k, v in by_kind.items():
    print(f"  {k}: {len(v)}")

# ---- Resolve to UniProt
uniprot_accs = set()
refseq_ids   = sorted({r["_norm"] for r in by_kind["refseq"]})
genbank_ids  = sorted({r["_norm"] for r in by_kind["genbank"]})
uniprot_ids  = sorted({r["_norm"] for r in by_kind["uniprot"]})

print(f"[idmap] refseq: {len(refseq_ids)}  genbank: {len(genbank_ids)}  "
      f"direct uniprot: {len(uniprot_ids)}")
r_to_up = idmap("RefSeq_Protein", refseq_ids) if refseq_ids else {}
g_to_up = idmap("EMBL-GenBank-DDBJ_CDS", genbank_ids) if genbank_ids else {}

for r in by_kind["uniprot"]:
    uniprot_accs.add(r["_norm"])
for r in by_kind["refseq"]:
    uniprot_accs.update(r_to_up.get(r["_norm"], []))
for r in by_kind["genbank"]:
    uniprot_accs.update(g_to_up.get(r["_norm"], []))

print(f"[fetch] PDB xrefs for {len(uniprot_accs)} UniProt accessions")
pdb_by_up = fetch_pdb_xrefs(uniprot_accs)

def uniprots_for(row):
    if row["_kind"] == "uniprot": return [row["_norm"]]
    if row["_kind"] == "refseq":  return r_to_up.get(row["_norm"], [])
    return g_to_up.get(row["_norm"], [])

# ---- Build output
resolved, unresolved = [], []
for r in rows:
    if r["db"] == "PDB":
        resolved.append({"": r[""], "chain": r["chain"],
                         "db": "PDB", "id": r["id"],
                         "source": "original",
                         "original_id": "",
                         "note": ""})
        continue
    ups = uniprots_for(r)
    candidates = []
    for up in ups:
        candidates.extend(pdb_by_up.get(up, []))
    # De-dup by PDB id (keep best-resolution)
    best_by_id = {}
    for c in candidates:
        prev = best_by_id.get(c["pdb"])
        if prev is None or resolution_float(c["resolution"]) < resolution_float(prev["resolution"]):
            best_by_id[c["pdb"]] = c
    picked = pick_pdb(list(best_by_id.values()), r["chain"])
    if picked is None:
        unresolved.append({**{k: r[k] for k in ("", "chain", "db", "id")},
                           "note": "no PDB found"})
        continue
    new_chain = r["chain"] if r["chain"] in chains_of(picked["chains_str"]) \
                           else sorted(chains_of(picked["chains_str"]))[0]
    note = ""
    if r["chain"] and r["chain"] != new_chain:
        note = f"chain {r['chain']}->{new_chain} (not in PDB)"
    resolved.append({"": r[""], "chain": new_chain, "db": "PDB",
                     "id": picked["pdb"],
                     "source": "resolved",
                     "original_id": r["id"],
                     "note": note})

# ---- Write
fields = ["", "chain", "db", "id", "source", "original_id", "note"]
with OUT.open("w") as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    w.writeheader()
    for r in resolved:
        w.writerow(r)

if unresolved:
    with UNRES.open("w") as f:
        w = csv.DictWriter(f, fieldnames=["", "chain", "db", "id", "note"],
                           delimiter="\t")
        w.writeheader()
        for r in unresolved:
            w.writerow(r)

print(f"\nresolved:   {len(resolved)}")
print(f"unresolved: {len(unresolved)}")
print(f"wrote {OUT}")
if unresolved:
    print(f"wrote {UNRES}")
