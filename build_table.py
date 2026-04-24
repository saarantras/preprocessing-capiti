#!/usr/bin/env python3
"""Filter CARD aro_index.tsv for proteins conferring resistance to drugs in
better_list.tsv. Emits:
  - capiti_C_targets_list.tsv: filtered CARD rows
  - class_map.tsv: which drugs (from better_list) map to which CARD class token

Matching rule: each better_list Class -> one CARD "target". A target is either
a drug-class regex (e.g. 'cephalosporin') matched against the CARD Drug Class
column, or a name-level regex (e.g. 'daptomycin') matched against AMR Gene
Family + ARO Name (for classes that aren't their own CARD drug class).
"""

import csv, re
from pathlib import Path
from collections import defaultdict

ROOT = Path(__file__).parent
INDEX = ROOT / "card" / "aro_index.tsv"
BETTER = ROOT / "better_list.tsv"
OUT_ROWS = ROOT / "capiti_C_targets_list.tsv"
OUT_MAP  = ROOT / "class_map.tsv"

# better_list Class value -> (match_mode, CARD regex, CARD token label)
# match_mode: "class" searches Drug Class column; "name" searches
# AMR Gene Family + ARO Name for genes linked to a specific antibiotic.
# None means intentionally unmapped (documented in class_map.tsv).
CLASS_MAP = {
    "Aminoglycosides":                        ("class", r"\baminoglycoside antibiotic\b",   "aminoglycoside antibiotic"),
    "Carbapenems":                            ("class", r"\bcarbapenem\b",                  "carbapenem"),
    "Third-generation-cephalosporins":        ("class", r"\bcephalosporin\b",               "cephalosporin"),
    "Fifth-generation-cephalosporins":        ("class", r"\bcephalosporin\b",               "cephalosporin"),
    "Other-cephalosporins":                   ("class", r"\bcephalosporin\b",               "cephalosporin"),
    "Glycopeptides":                          ("class", r"\bglycopeptide antibiotic\b",     "glycopeptide antibiotic"),
    "Glycylcyclines":                         ("class", r"\bglycylcycline\b",               "glycylcycline"),
    "Lipopeptides":                           ("name",  r"\bdaptomycin\b",                  "daptomycin (lipopeptide)"),
    "Monobactams":                            ("class", r"\bmonobactam\b",                  "monobactam"),
    "Monobactam/beta-lactamase-inhibitor":    ("class", r"\bmonobactam\b",                  "monobactam"),
    "Oxazolidinones":                         ("class", r"\boxazolidinone antibiotic\b",    "oxazolidinone antibiotic"),
    "Penems":                                 (None,    None,                               "(no CARD class; faropenem absent from CARD)"),
    "Phosphonics":                            ("class", r"\bphosphonic acid antibiotic\b",  "phosphonic acid antibiotic"),
    "Pleuromutilin":                          ("class", r"\bpleuromutilin antibiotic\b",    "pleuromutilin antibiotic"),
    "Polymyxins":                             ("name",  r"\b(polymyxin|colistin)\b",        "polymyxin/colistin (peptide antibiotic)"),
    "Streptogramins":                         ("class", r"\bstreptogramin( [AB])? antibiotic\b", "streptogramin antibiotic"),
    "Tetracyclines":                          ("class", r"\btetracycline antibiotic\b",     "tetracycline antibiotic"),
    "Trimethoprim-derivatives":               ("class", r"\bdiaminopyrimidine antibiotic\b","diaminopyrimidine antibiotic"),
}

# ---- Load better_list, group drugs by class
drugs_by_class = defaultdict(list)
with BETTER.open() as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        ab = (row.get("Antibiotic") or "").strip()
        cl = (row.get("Class") or "").strip()
        if ab and cl:
            drugs_by_class[cl].append(ab)

unknown = [c for c in drugs_by_class if c not in CLASS_MAP]
if unknown:
    raise SystemExit(f"Unmapped class(es) in better_list: {unknown}")

# ---- Write class_map.tsv (for manuscript)
with OUT_MAP.open("w") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["Better_List_Class", "Match_Mode", "CARD_Token", "Drugs"])
    for cl in sorted(drugs_by_class):
        mode, _rx, tok = CLASS_MAP[cl]
        w.writerow([cl, mode or "unmapped", tok, ";".join(drugs_by_class[cl])])

# ---- Compile one combined matcher per (match_mode, regex) pair, remembering
# which better_list classes it represents.
compiled = []  # (mode, pattern, [better_list_classes], CARD_token)
by_target = defaultdict(list)
for cl, (mode, rx, tok) in CLASS_MAP.items():
    if mode is None or cl not in drugs_by_class:
        continue
    by_target[(mode, rx, tok)].append(cl)
for (mode, rx, tok), cls in by_target.items():
    compiled.append((mode, re.compile(rx, re.I), cls, tok))

# ---- Scan CARD index
rows_out = []
with INDEX.open() as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        class_hay = row["Drug Class"]
        name_hay  = row["AMR Gene Family"] + " | " + row["ARO Name"]
        hits_tokens = []
        hits_classes = []
        for mode, pat, cls, tok in compiled:
            hay = class_hay if mode == "class" else name_hay
            if pat.search(hay):
                hits_tokens.append(tok)
                hits_classes.extend(cls)
        if hits_tokens:
            row["Matched_CARD_Tokens"] = ";".join(sorted(set(hits_tokens)))
            row["Matched_Better_Classes"] = ";".join(sorted(set(hits_classes)))
            rows_out.append(row)

fields = ["ARO Accession", "ARO Name", "CARD Short Name",
          "Protein Accession", "DNA Accession",
          "AMR Gene Family", "Drug Class", "Resistance Mechanism",
          "Matched_CARD_Tokens", "Matched_Better_Classes"]

with OUT_ROWS.open("w") as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for r in rows_out:
        w.writerow(r)

# ---- Report
from collections import Counter
print(f"matched rows: {len(rows_out)}")
print(f"unique protein accessions: "
      f"{len({r['Protein Accession'] for r in rows_out if r['Protein Accession']})}")
c = Counter()
for r in rows_out:
    for t in r["Matched_CARD_Tokens"].split(";"):
        c[t] += 1
print("rows per CARD token:")
for t, n in c.most_common():
    print(f"  {n:6d}  {t}")
print(f"\nwrote {OUT_ROWS}")
print(f"wrote {OUT_MAP}")
