# cardpreprocessing

Preprocessing pipelines for the CAPITI ML training targets. Two independent
scopes live here:

- **CAPITI-C**: antibiotic-resistance protein targets derived from CARD,
  filtered by the WHO Reserve-tier drug classes in `better_list.tsv`.
- **CAPITI-E**: a separate curated list (`capiti_E_splitdb.tsv`) where mixed
  PDB / NCBI / UniProt identifiers are homogenized to a single structural
  reference per row (PDB, else AlphaFold).

## Data sources

- CARD: retrieved from <https://card.mcmaster.ca/latest/data> on
  2026-04-24 10:25 AM EST. The downloaded `card/` directory is gitignored;
  re-download it before running the CAPITI-C pipeline.
- UniProt ID-mapping and cross-references:
  <https://rest.uniprot.org/idmapping>, <https://rest.uniprot.org/uniprotkb>.
- UniParc (rescue pass for obsolete/unmapped accessions):
  <https://rest.uniprot.org/uniparc>.
- AlphaFold DB predictions (CAPITI-E fallback):
  <https://alphafold.ebi.ac.uk/api/prediction>.

## CAPITI-C pipeline

1. `build_table.py` -- reads `better_list.tsv` and `card/aro_index.tsv`,
   emits `capiti_C_targets_list.tsv` (prefilter) and `class_map.tsv`
   (drug -> CARD-class mapping, for manuscript methods).
2. `fetch_pdb.py` -- annotates `capiti_C_targets_list.tsv` in place with
   `UniProt`, `UniParc`, `PDB`, and `PDB_Links` columns. Pass 1 uses
   UniProt ID-mapping (`RefSeq_Protein` / `EMBL-GenBank-DDBJ_CDS` ->
   `UniProtKB`). Pass 2 rescues unmapped accessions via a UniParc
   accession search and harvests any `UniProtKB` cross-references found
   inside the UniParc record, then pass 3 resolves PDB xrefs for those.
3. `anonymize.py` -- emits `capiti_C_targets_anonymized.tsv`, a leakage-free
   target list for downstream ML: rows become `T-C<N>` IDs with only
   `PDB`, `chains`, `method`, `resolution`. Targets are assigned in sorted
   order of ARO Accession for reproducibility. One row per (target, PDB).

### CAPITI-C numbers

- 5,218 CARD rows pass the class filter (5,122 unique protein accessions).
- UniProt mapping: 4,104 / 5,218 rows (78.6%).
- UniParc handles for otherwise-unmapped accessions: 1,136 additional rows.
- PDB coverage: 246 / 5,218 rows (4.7%). These become the 246 anonymized
  targets in `capiti_C_targets_anonymized.tsv` (2,106 target-PDB rows).

## CAPITI-E pipeline

`capiti_E_splitdb.tsv` mixes PDB codes with UniProt and NCBI protein
accessions. Goal: homogenize to a structural reference per row.

1. `unify_pdb.py` -- for each non-PDB row, resolves the identifier to one
   or more UniProt accessions (idmapping when needed) and picks the best
   PDB cross-reference. Chain selection prefers PDBs whose Chains property
   contains the row's existing chain letter; ties broken by resolution.
   Rows with no PDB xref go to `capiti_E_unresolved.tsv`. Output:
   `capiti_E_splitdb_pdb.tsv` with `source`, `original_id`, `note` columns.
2. `add_alphafold.py` -- rescues unresolved rows via the AlphaFold DB.
   Uses `globalMetricValue` (mean pLDDT) as the model confidence score.
   Rows with pLDDT below `PLDDT_DROP` (default 50, AlphaFold's "very low"
   band) are kept in `capiti_E_unresolved.tsv`; everything else is
   appended to `capiti_E_splitdb_pdb.tsv` with `db=AF`, `chain=A`,
   `id=<UniProt>`, and `confidence=<mean pLDDT>`.

### CAPITI-E numbers

- 72 input rows. 53 already had `db=PDB`.
- 19 non-PDB rows: 1 resolved to an experimental PDB (Q57398 -> 1UYJ),
  7 rescued via AlphaFold (1 very high, 6 confident), 11 remain unresolved
  (GenBank / RefSeq WP accessions with no UniProt mapping and no AF entry).
- Final `capiti_E_splitdb_pdb.tsv`: 61 rows (54 PDB + 7 AF).

## Files

- `better_list.tsv` -- CAPITI-C input drug list (WHO Reserve tier).
- `class_map.tsv` -- CAPITI-C drug -> CARD drug-class token mapping.
- `capiti_C_targets_list.tsv` -- CAPITI-C annotated target list
  (CARD fields + UniProt/UniParc/PDB cross-references).
- `capiti_C_targets_anonymized.tsv` -- CAPITI-C anonymized ML target list.
- `capiti_E_splitdb.tsv` -- CAPITI-E input list (mixed db identifiers).
- `capiti_E_splitdb_pdb.tsv` -- CAPITI-E homogenized output (PDB + AF).
- `capiti_E_unresolved.tsv` -- CAPITI-E rows with no available structure.
