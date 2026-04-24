# cardpreprocessing

Preprocessing pipeline that builds `capiti_C_targets_list.tsv`, the target
list used to train CAPITI (version C). Filters the CARD resistance database
down to proteins conferring resistance to the antibiotic classes in
`better_list.tsv` (WHO Reserve tier), and annotates each protein with
UniProt, UniParc, and PDB cross-references.

## Data source

Retrieved CARD from <https://card.mcmaster.ca/latest/data> on
2026-04-24 10:25 AM EST. The downloaded `card/` directory is gitignored;
re-download it into `card/` before running the pipeline.

## Pipeline

1. `build_table.py` -- reads `better_list.tsv` and `card/aro_index.tsv`,
   emits `capiti_C_targets_list.tsv` (prefilter) and `class_map.tsv`
   (drug -> CARD-class mapping, for manuscript methods).
2. `fetch_pdb.py` -- annotates `capiti_C_targets_list.tsv` in place with
   `UniProt`, `UniParc`, `PDB`, and `PDB_Links` columns via the UniProt
   REST API (idmapping + UniParc rescue pass).

## Files

- `better_list.tsv` -- input drug list (WHO Reserve tier).
- `class_map.tsv` -- drug -> CARD drug-class token mapping.
- `capiti_C_targets_list.tsv` -- final annotated target list.
