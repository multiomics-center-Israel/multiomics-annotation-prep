# R GMT builders

R utilities for turning organism annotations into GMT gene-set files for
`multiomics-core`. These complement the Python KEGG/GO downloaders in `src/` and
`scripts/`; they are used when a project ships its own annotation (a GAF, an
eggNOG-mapper table, a KEGG dump) that needs to be re-keyed onto the DE feature
IDs of a specific dataset.

All of them source `gmt_utils.R` (vendored `read_gmt`/`write_gmt`/`validate_gmt`/
`filter_gmt_by_size`/`generate_gmt_from_*` from `multiomics-core`) — no
dependency on the pipeline source tree.

| Script | What it does |
|--------|--------------|
| `build_ankri_go_gmt.R` | Serge Ankri *E. histolytica*: bridges the NCBI GAF (GeneID→GO) to the DIA-NN `XP_` Protein.Group via the genome GTF, writes `GO_ehistolytica.gmt`. |
| `build_elad_go_gmt.R` | Elad Chiel *Spalangia*+*Sodalis*: eggNOG (host) + RefSeq GAF (symbiont) → per-species `GO_spalangia.gmt` / `GO_sodalis.gmt`. |
| `generate_gmt.R` | Generic GMT builder (GO/KEGG via biomaRt/KEGGREST or a custom table). |
| `extract_kegg_compounds.R` | Pull KEGG compound→pathway sets for metabolomics. |
| `gmt_utils.R` | Vendored shared helpers (keep in sync with `multiomics-core/R/core/13_gmt_utils.R`). |

## Data location

Project inputs (GAF/GTF/eggNOG/pg_matrix) are **not** stored in this repo. Point
the builders at the analysis project's data with:

```bash
export ANNPREP_DATA_ROOT=/path/to/multiomics-core/data
conda activate mofa2_env
Rscript R/build_ankri_go_gmt.R
```

Without `ANNPREP_DATA_ROOT` the scripts look under this repo's own `data/`.
