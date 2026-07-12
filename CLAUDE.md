# CLAUDE.md

Guidance for Claude Code (and other agents) working in this repository.

## What this project is

`multiomic-annotation-prep` prepares **organism-specific annotation files for
enrichment analysis** (KEGG & GO) that are consumed by the `multiomic-core`
pipeline and the `Neat_RNA-Seq` / `Neat_Proteomics` workflows. Unlike the
manual workflow it is based on
([Neat_Annotation](https://github.com/veredcc/Neat_Annotation)), this repo
**downloads all reference files automatically** from the original resources
(KEGG REST API, Gene Ontology, Ensembl/BioMart, UniProt) and writes ready-to-use
output files for `clusterProfiler`.

Language: **Python** (>= 3.8). Dependencies: `requests`, `pyyaml`.
Primary use case: **non-model organisms** (transcriptome + KAAS).
Ensembl/UniProt modules are secondary (model organisms).

## Output contract (do not break without a reason)

### .tab files (for Neat_RNA-Seq / Neat_Proteomics / clusterProfiler::enricher)

These file names and column headers are what `clusterProfiler::enricher()`
expects positionally (column 1 = term ID, column 2 = gene or name).
All writers are centralized in `src/utils.py` (`write_*` functions).

| File | Header row | Meaning |
|------|-----------|---------|
| `KEGG_pathway2gene.tab` | `v1` \t `index` | pathway id, gene id |
| `KEGG_pathway2name.tab` | `pathway` \t `info` | pathway id, pathway name |
| `GO2gene_{BP,MF,CC}.tab` | `GO` \t `Gene` | GO id, gene id |
| `GO2name_{BP,MF,CC}.tab` | `GO` \t `Term` | GO id, term name |
| `KEGG_annot_genes.txt` | multi-col | descriptive per-gene annotation |
| `Annotation.tab` | multi-col | descriptive annot (Ensembl/UniProt modules) |

### .gmt files (for multiomic-core)

`multiomic-core` reads **GMT format** for non-model organisms via its
`read_gmt()` / `load_gene_sets()` functions. Each line is:
`TERM_ID<tab>Description<tab>GENE1<tab>GENE2<tab>...`

| File | Content |
|------|---------|
| `KEGG_pathway.gmt` | KEGG pathway gene sets |
| `GO_{BP,MF,CC}.gmt` | GO gene sets per namespace |

> **FORMAT RESOLVED:** The `.tab` files follow the Neat_RNA-Seq / Neat_Proteomics
> (clusterProfiler) convention and are correct for that workflow. `multiomic-core`
> uses GMT format for custom (non-model) gene sets — both formats are now produced.

## Repository map

```
src/
  __init__.py
  utils.py                 # logging, ensure_dir, cached_download, write_* (OUTPUT FORMAT lives here)
  download_kegg.py         # KEGG REST downloads + parsers (ko/path gene-sets AND cpd/rn/pathway entities)
  download_kegg_org.py     # KEGG for a specific organism code (model organisms)
  prepare_kegg_nonmodel.py # KAAS query.ko.txt -> KEGG enrichment + annot files  [PRIMARY]
  download_go.py           # GO term names/namespaces (go-basic.obo parser + ancestor graph)
  prepare_go.py            # per-gene GO table -> expanded GO enrichment files    [PRIMARY]
  prepare_ensembl.py       # pybiomart path (model organisms)                     [secondary]
  prepare_uniprot.py       # UniProt REST path (proteomics)                       [secondary]
  masses.py                # neutral monoisotopic mass (validates formula, delegates to mass2chem)
  prepare_mummichog_model.py # KEGG org code -> mummichog metabolic model JSON + manifest  [NEW output type]
  model_registry.py        # MODELS.md registry helpers (append_model_row, idempotent on tag)
scripts/
  run_kegg_nonmodel.py     # CLI wrapper (argparse)
  run_go.py                # CLI wrapper
  run_mummichog_model.py   # CLI wrapper for the mummichog model
  publish_model.py         # release a model + auto-append a row to MODELS.md, then commit it
  run_all.py               # config-driven driver (reads config/config.yml)
config/config.yml          # which modules run + their inputs
requirements-mummichog.txt # scoped, pinned optional deps for the mummichog model only
MODELS.md                  # registry/index of published models (tracked source, not an artifact)
tests/                     # pytest: masses, KEGG parsers, model assembly, mummichog smoke, registry
examples/                  # tiny inputs so the tool runs end-to-end
data/                      # download cache (git-ignored)
results/                   # output files (git-ignored)
```

The **mummichog metabolic model** is a separate, COMPOUND-centric output type
(compounds w/ neutral formula+mass, real substrate/product reactions,
pathways-of-reactions) consumed by `mummichog -n <model>.json`. It is NOT derived
from the gene sets; it is pulled from KEGG cpd/rn/pathway entities and serialized
to the metDataModel shape. Both inputs (a KEGG organism code, or a KAAS KO list)
reduce to the same pipeline: `KO list -> reactions -> compounds -> pathways`.
Contract: `MODEL_CONTRACT.md`. Its deps (`metDataModel`, `mass2chem`,
`mummichog`) are optional imports scoped to this module
(`requirements-mummichog.txt`) so the gene-set modules stay `requests`+
`pyyaml`-only.

## How to run

```bash
pip install requests pyyaml                             # once (only deps)

python scripts/run_kegg_nonmodel.py --kaas examples/query.ko.txt --out results --cache data
python scripts/run_go.py --go-table examples/trinotate_go.txt --out results --cache data
python scripts/run_all.py --config config/config.yml    # config-driven
```

Requires network access (downloads from rest.kegg.jp, purl.obolibrary.org,
ensembl.org, rest.uniprot.org). Downloads are cached under `data/`; pass
`--refresh` to force re-download.

## Conventions

- **Every remote fetch goes through `cached_download()`** in `src/utils.py` so runs
  are reproducible and offline after the first fetch.
- **All output writing goes through the `write_*` helpers** in `src/utils.py`. Keep
  the output format decisions there, not scattered across modules.
- **Logging** via `log_msg()` (timestamped to stderr). Avoid bare `print`.
- Use only stdlib + `requests` + `pyyaml`; the Ensembl module optionally needs
  `pybiomart`. Avoid adding heavy new deps without reason.

## Gotchas / things to verify when editing

- `prepare_kegg_nonmodel.py` strips isoform suffixes with `re.sub(r"_i\d+$", "", ...)`
  to collapse transcripts to gene level (matches the original Perl). Disable with
  `--no-strip-isoform` if IDs are already gene-level.
- KEGG `ko->path` keeps only reference pathways (`^map`). Organism-code pathways
  (`mmu00010` etc.) come from `download_kegg_org.py` instead.
- KEGG has **no gene->reaction link** for an organism (`link/reaction/<org>` is an
  invalid query → HTTP 400). The mummichog model resolves reactions via KO:
  `link/ko/<org>` (gene→KO) intersected with `link/reaction/ko` (KO→reaction,
  KEGG-wide). The KAAS path supplies the same KO list from a file. Do not
  "restore" a direct gene→reaction lookup.
- GO hierarchy expansion uses the `is_a` relationships from `go-basic.obo` to build
  a transitive ancestor closure — every direct annotation propagates to all parent
  terms.
- Trinotate GO strings carry `^namespace^description` suffixes and mixed
  separators; the parser extracts `GO:\d{7}` tokens by regex — keep that robust.

## Open tasks / roadmap

1. Add a KEGG Mapper module (colored pathway maps + per-pathway gene counts),
   like Neat_Annotation's `create_colored_KEGG_maps_and_summarize_per_path`.
2. Add automated tests (e.g. `pytest`) for the parsers and writers using the
   `examples/` inputs with a small mocked/cached reference set.
3. Descriptive-annotation file for the non-model GO path is not yet produced
   (only enrichment files) — add if the pipeline needs it.

## Verification (no CI yet)

To sanity-check a change, run both `run_*` scripts on the `examples/` inputs
and confirm the expected output files appear in `results/` with the correct
headers:

```bash
rm -rf results
python scripts/run_kegg_nonmodel.py --kaas examples/query.ko.txt --out results --cache data
python scripts/run_go.py --go-table examples/trinotate_go.txt --out results --cache data
```
