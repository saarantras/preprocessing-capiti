#!/usr/bin/env python3
"""Annotate capiti_C_targets_list.tsv with UniProt, UniParc, and PDB IDs.

Pass 1: UniProt ID-mapping from RefSeq_Protein / EMBL-GenBank-DDBJ_CDS -> UniProtKB.
Pass 2 (rescue): for accessions not found in pass 1, search UniParc by accession,
  record the UPI, and harvest any UniProtKB cross-references inside the UniParc
  entry. Then fetch PDB xrefs for those additional UniProt accessions.
"""

import csv, json, re, sys, time, urllib.parse, urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).parent
IO_PATH = ROOT / "capiti_C_targets_list.tsv"

API = "https://rest.uniprot.org/idmapping"
UNIPROT = "https://rest.uniprot.org/uniprotkb"
UNIPARC = "https://rest.uniprot.org/uniparc"

_NEXT_RE = re.compile(r'<([^>]+)>\s*;\s*rel="next"')

def _get(url):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    return urllib.request.urlopen(req, timeout=60)

def _next_link(link_header):
    m = _NEXT_RE.search(link_header or "")
    return m.group(1) if m else None

# ---------- Pass 1: idmapping ----------

def submit(from_db, ids):
    req = urllib.request.Request(f"{API}/run",
        data=urllib.parse.urlencode({"from": from_db, "to": "UniProtKB",
                                     "ids": ",".join(ids)}).encode(),
        method="POST")
    return json.load(urllib.request.urlopen(req, timeout=60))["jobId"]

def wait(job):
    while True:
        s = json.load(_get(f"{API}/status/{job}"))
        st = s.get("jobStatus")
        if st in (None, "FINISHED"): return
        if st == "ERROR": raise RuntimeError(f"job {job} failed: {s}")
        time.sleep(3)

def extract_pdbs(entry):
    return [x["id"] for x in entry.get("uniProtKBCrossReferences", [])
            if x.get("database") == "PDB"]

def idmap_batch(from_db, ids, mapping):
    if not ids: return
    print(f"  idmapping {from_db}: {len(ids)} ids")
    job = submit(from_db, ids)
    wait(job)
    url = (f"{API}/uniprotkb/results/{job}"
           f"?fields=accession,xref_pdb&format=json&size=500")
    seen = set()
    while url:
        with _get(url) as r:
            payload = json.load(r)
            for rec in payload.get("results", []):
                src = rec["from"]
                entry = rec["to"]
                seen.add(src)
                up = entry.get("primaryAccession", "")
                pdbs = extract_pdbs(entry)
                cur = mapping.setdefault(src, {"uniprot": set(),
                                               "uniparc": set(),
                                               "pdb": set()})
                if up: cur["uniprot"].add(up)
                cur["pdb"].update(pdbs)
            url = _next_link(r.headers.get("Link", ""))
    print(f"    mapped {len(seen)}/{len(ids)}")

# ---------- Pass 2: UniParc rescue ----------

def uniparc_lookup(acc):
    """Search UniParc for an accession. Return (upi, [uniprot_accs], [pdb_ids])
    or None if no hit."""
    # Strip version for search; UniParc indexes unversioned IDs too.
    bare = acc.split(".")[0]
    q = urllib.parse.quote(bare)
    try:
        with _get(f"{UNIPARC}/search?query={q}&format=json&size=1") as r:
            p = json.load(r)
    except Exception:
        return None
    results = p.get("results") or []
    if not results: return None
    upi = results[0].get("uniParcId")
    if not upi: return None
    # Fetch full entry for xrefs
    try:
        with _get(f"{UNIPARC}/{upi}?format=json") as r:
            entry = json.load(r)
    except Exception:
        return (upi, [], [])
    uniprots = []
    pdbs = []
    for x in entry.get("uniParcCrossReferences", []):
        db = x.get("database", "")
        if db.startswith("UniProtKB"):
            uniprots.append(x["id"])
        elif db == "PDB":
            pdbs.append(x["id"])
    return (upi, uniprots, pdbs)

def uniprot_pdbs_batch(accs):
    """Fetch PDB xrefs for a list of UniProt accessions via uniprotkb/search."""
    out = {}
    accs = [a for a in accs if a]
    CHUNK = 50
    for i in range(0, len(accs), CHUNK):
        chunk = accs[i:i+CHUNK]
        q = " OR ".join(f"accession:{a}" for a in chunk)
        url = (f"{UNIPROT}/search?query={urllib.parse.quote(q)}"
               f"&fields=accession,xref_pdb&format=json&size=500")
        while url:
            with _get(url) as r:
                p = json.load(r)
                for e in p.get("results", []):
                    out[e["primaryAccession"]] = extract_pdbs(e)
                url = _next_link(r.headers.get("Link", ""))
    return out

# ---------- Main ----------

def main():
    rows = list(csv.DictReader(IO_PATH.open(), delimiter="\t"))
    all_acc = sorted({r["Protein Accession"] for r in rows if r["Protein Accession"]})
    refseq = [a for a in all_acc if a.split(".")[0].startswith(("WP_", "NP_"))]
    genbank = [a for a in all_acc if a not in set(refseq)]
    print(f"total: {len(all_acc)}  refseq: {len(refseq)}  genbank: {len(genbank)}")

    mapping = {}  # acc -> {uniprot: set, uniparc: set, pdb: set}

    print("[pass 1] UniProt idmapping")
    idmap_batch("RefSeq_Protein", refseq, mapping)
    idmap_batch("EMBL-GenBank-DDBJ_CDS", genbank, mapping)

    unmapped = [a for a in all_acc if a not in mapping]
    print(f"[pass 2] UniParc rescue: {len(unmapped)} unmapped accessions")

    extra_uniprots = set()
    done = 0
    with ThreadPoolExecutor(max_workers=12) as pool:
        futs = {pool.submit(uniparc_lookup, a): a for a in unmapped}
        for fut in as_completed(futs):
            acc = futs[fut]
            done += 1
            if done % 200 == 0:
                print(f"    {done}/{len(unmapped)}")
            res = fut.result()
            if res is None:
                continue
            upi, ups, pdbs = res
            cur = mapping.setdefault(acc, {"uniprot": set(),
                                           "uniparc": set(),
                                           "pdb": set()})
            cur["uniparc"].add(upi)
            cur["uniprot"].update(ups)
            cur["pdb"].update(pdbs)
            extra_uniprots.update(ups)

    if extra_uniprots:
        print(f"[pass 3] fetch PDB xrefs for {len(extra_uniprots)} "
              "UniProt accessions discovered via UniParc")
        pdb_by_up = uniprot_pdbs_batch(sorted(extra_uniprots))
        # Propagate discovered PDBs back to source accessions.
        for acc, info in mapping.items():
            for up in info["uniprot"]:
                info["pdb"].update(pdb_by_up.get(up, []))

    # ---- Write output
    fields = [k for k in rows[0].keys() if k not in
              ("UniProt", "UniParc", "PDB", "PDB_Links")]
    fields += ["UniProt", "UniParc", "PDB", "PDB_Links"]
    with IO_PATH.open("w") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in rows:
            info = mapping.get(r["Protein Accession"],
                               {"uniprot": set(), "uniparc": set(), "pdb": set()})
            r["UniProt"] = ";".join(sorted(info["uniprot"]))
            r["UniParc"] = ";".join(sorted(info["uniparc"]))
            pdbs = sorted(info["pdb"])
            r["PDB"] = ";".join(pdbs)
            r["PDB_Links"] = ";".join(
                f"https://www.rcsb.org/structure/{p}" for p in pdbs)
            w.writerow(r)

    # ---- Report
    n_up = sum(1 for r in rows if mapping.get(r["Protein Accession"],
                                              {}).get("uniprot"))
    n_uc = sum(1 for r in rows if mapping.get(r["Protein Accession"],
                                              {}).get("uniparc"))
    n_pdb = sum(1 for r in rows if mapping.get(r["Protein Accession"],
                                               {}).get("pdb"))
    print(f"\nrows with UniProt: {n_up}/{len(rows)}")
    print(f"rows with UniParc: {n_uc}/{len(rows)}")
    print(f"rows with PDB:     {n_pdb}/{len(rows)}")
    print(f"wrote {IO_PATH}")

if __name__ == "__main__":
    main()
