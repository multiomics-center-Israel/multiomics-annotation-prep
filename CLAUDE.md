# CLAUDE.md

Guidance for Claude Code (and other agents) working in this repository.

## What this project is

`multiomic-annotation-prep` prepares **organism-specific annotation files for
enrichment analysis** (KEGG & GO) that are consumed by the `multiomic-core`
pipeline. Unlike the manual workflow it is based on
([Neat_Annotation](https://github.com/veredcc/Neat_Annotation)), this repo
**downloads all reference files automatically** from the original resources
(KEGG REST API, Gene Ontology, Ensembl/BioMart, UniProt) and writes ready-to-use
`.tab` files for `clusterProfiler`.

Language: **R** (>= 4.0). Primary use case: **non-model organisms**
(transcriptome + KAAS). Ensembl/UniProt modules are secondary (model organisms).

## Output contract (do not break without a reason)

These file names and column headers are what `multiomic-core` / `clusterProfiler`
expect. All writers are centralized in `R/utils.R` (`write_*` functions) — change
the format in ONE place if the pipeline's expected layout differs.

| File | Header row | Meaning |
|------|-----------|---------|
| `KEGG_pathway2gene.tab` | `v1` \t `index` | pathway id, gene id |
| `KEGG_pathway2name.tab` | `pathway` \t `info` | pathway id, pathway name |
| `GO2gene_{BP,MF,CC}.tab` | `GO` \t `Gene` | GO id, gene id |
| `GO2name_{BP,MF,CC}.tab` | `GO` \t `Term` | GO id, term name |
| `KEGG_annot_genes.txt` | multi-col | descriptive per-gene annotation |
| `Annotation.tab` | multi-col | descriptive annot (Ensembl/UniProt modules) |

> **OPEN QUESTION (highest priority):** the exact format `multiomic-core`
> expects has not been confirmed. The current format follows the
> Neat_RNA-Seq / Neat_Proteomics (clusterProfiler) convention. Before relying on
> outputs, verify against a real `multiomic-core` run and adjust the `write_*`
> helpers if needed.

## Repository map

```
R/
  bootstrap.R              # finds repo root, sources all modules (used by scripts/)
  install_deps.R           # installs CRAN + Bioconductor deps
  utils.R                  # logging, ensure_dir, cached_download, %||%, write_* (OUTPUT FORMAT lives here)
  download_kegg.R          # KEGG REST downloads + parsers (ko->name, ko->path, path->name)
  download_kegg_org.R      # KEGG for a specific organism code (model organisms)
  prepare_kegg_nonmodel.R  # KAAS query.ko.txt -> KEGG enrichment + annot files  [PRIMARY]
  download_go.R            # GO term names/namespaces (GO.db preferred, else go-basic.obo)
  prepare_go.R             # per-gene GO table -> expanded GO enrichment files    [PRIMARY]
  prepare_ensembl.R        # biomaRt path (model organisms)                       [secondary]
  prepare_uniprot.R        # UniProt REST path (proteomics)                       [secondary]
scripts/
  run_kegg_nonmodel.R      # CLI wrapper (optparse)
  run_go.R                 # CLI wrapper
  run_all.R                # config-driven driver (reads config/config.yml)
config/config.yml          # which modules run + their inputs
examples/                  # tiny inputs so the tool runs end-to-end
data/                      # download cache (git-ignored)
results/                   # output .tab files (git-ignored)
```

## How to run

```bash
Rscript -e 'source("R/install_deps.R")'                 # once
Rscript scripts/run_kegg_nonmodel.R --kaas examples/query.ko.txt --out results --cache data
Rscript scripts/run_go.R --go-table examples/trinotate_go.txt --out results --cache data
Rscript scripts/run_all.R --config config/config.yml    # config-driven
```

Requires network access (downloads from rest.kegg.jp, purl.obolibrary.org,
ensembl.org, rest.uniprot.org). Downloads are cached under `data/`; pass
`--refresh` (or `refresh: true` in config) to force re-download.

## Conventions

- **Every remote fetch goes through `cached_download()`** in `utils.R` so runs
  are reproducible and offline after the first fetch. Do not add ad-hoc
  `download.file`/`httr::GET` calls elsewhere.
- **All output writing goes through the `write_*` helpers** in `utils.R`. Keep
  the output format decisions there, not scattered across modules.
- **Logging** via `log_msg()` (timestamped). Avoid bare `print`/`cat`.
- Modules are plain sourced R files (no package namespace). `bootstrap.R` sources
  them in dependency order; add new files to its `files` vector.
- Use base R + the listed Bioconductor/CRAN deps; avoid adding heavy new deps
  without reason.

## Gotchas / things to verify when editing

- `prepare_kegg_nonmodel.R` strips isoform suffixes with `sub("_i\\d+$", "", ...)`
  to collapse transcripts to gene level (matches the original Perl). Disable with
  `--no-strip-isoform` if IDs are already gene-level.
- KEGG `ko->path` keeps only reference pathways (`^map`). Organism-code pathways
  (`mmu00010` etc.) come from `download_kegg_org.R` instead.
- `clusterProfiler::buildGOmap`'s returned column names vary by version;
  `expand_go()` defaults to the deterministic GO.db ancestor maps and only falls
  back to `buildGOmap` (with defensive column detection) if GO.db is missing.
- Trinotate GO strings carry `^namespace^description` suffixes and mixed
  separators; the parser extracts `GO:#######` tokens by regex — keep that robust.

## Open tasks / roadmap

1. **Confirm & lock the `multiomic-core` output format** (see OPEN QUESTION above).
2. Add a KEGG Mapper module (colored pathway maps + per-pathway gene counts),
   like Neat_Annotation's `create_colored_KEGG_maps_and_summarize_per_path`.
3. Add automated tests (e.g. `testthat`) for the parsers and writers using the
   `examples/` inputs with a small mocked/cached reference set.
4. Consider packaging as a proper R package (DESCRIPTION/NAMESPACE) if it grows.
5. Descriptive-annotation file for the non-model GO path is not yet produced
   (only enrichment files) — add if the pipeline needs it.

## Verification (no CI yet)

There is no test suite yet. To sanity-check a change, run both `run_*` scripts on
the `examples/` inputs and confirm the expected `.tab` files appear in `results/`
with the correct headers. The example inputs are deliberately tiny and use real
KO/GO ids so the downloads resolve.
