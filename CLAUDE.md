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
  kegg_entities.py         # KEGG source loading (load_source: KO->reactions->compounds->pathways); shared, light
  prepare_mummichog_model.py # KEGG org code -> mummichog metabolic model JSON + manifest  [NEW output type]
  prepare_kegg_compound_sets.py # KEGG -> pathway compound-set GMT + table (ID-based enrichment)  [NEW output type]
scripts/
  run_kegg_nonmodel.py     # CLI wrapper (argparse)
  run_go.py                # CLI wrapper
  run_mummichog_model.py   # CLI wrapper for the mummichog model
  run_kegg_compound_sets.py # CLI wrapper for the compound-set GMT + table
  run_ensembl.py           # CLI for the Ensembl/BioMart module (model organisms)
  run_uniprot.py           # CLI for the UniProt module (by taxon id)
  run_all.py               # config-driven driver (reads config/config.yml)
.github/workflows/         # workflow_dispatch pipelines: build on a runner + publish a dated GitHub Release
  publish-organism-artifacts.yml # KEGG code -> mummichog model + compound-set GMT (compound-based)
  build-rna-annotation.yml       # rna_inputs/<project>/ (KAAS+Trinotate) -> annotation_dir tables (gene-based)
  build-model-annotation.yml     # Ensembl/UniProt -> annotation_dir tables (gene-based, model organisms)
rna_inputs/                # committed KAAS/Trinotate inputs per project, read by build-rna-annotation.yml
config/config.yml          # which modules run + their inputs
requirements-mummichog.txt # scoped, pinned optional deps for the mummichog model only
docs/team-guide-he.md      # Hebrew operator guide (build / run / add-organism / input-prep)
tests/                     # pytest: masses, KEGG parsers, model assembly, mummichog smoke, compound sets
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

That shared pipeline lives in **`kegg_entities.py`** (`load_source`,
`is_metabolic_map`, `SOURCE_TOKENS`), a dependency-light module (`requests` +
`pyyaml` only) so both compound-centric outputs reuse it without importing each
other. The **`prepare_kegg_compound_sets.py`** module builds the ID-based
enrichment inputs (`<stem>.compound_pathway.gmt` + `<stem>.pathway2compound.tab`,
same stem as the model) from that same `load_source`. It applies NO mass filter,
so its compound set is a **superset** of the model's compounds (ID-based ORA /
GSEA / QEA need no mass). It has no heavy deps. Set
`mummichog_model.emit_compound_sets: true` (or `--emit-compound-sets`) to emit the
model + these companions from a single `load_source` (one KEGG snapshot);
companions are then listed in the model manifest's `companion_files`.
`multiomic-core` loads the GMT as a plain data file (no sha256 pinning, unlike the
model), so the sidecar manifest is provenance-only.

**Two consumption paths (know which you're building for).** `multiomic-core`
consumes two distinct kinds of artifact:

1. **COMPOUND-based** metabolomics inputs — the mummichog model (`model_ref`:
   URL+sha256, fetched+verified) and the compound-set GMT (`gmt_file`: a local
   path). Universal per KEGG organism code; built + released by
   `publish-organism-artifacts.yml`.
2. **GENE-based** enrichment tables — `KEGG_pathway2gene.tab` +
   `KEGG_pathway2name.tab` and `GO2gene_{BP,MF,CC}.tab` + `GO2name_{BP,MF,CC}.tab`,
   read from a local directory `enrichment.annotation_dir` (multiomic-core v2 /
   PR #128; the loader takes the first two columns, GO must be hierarchy-expanded).
   These are keyed on the project's **gene IDs**, so they are NOT universal per
   organism — they come from either a non-model transcriptome (KAAS +
   Trinotate → `prepare_kegg_nonmodel` + `prepare_go`) or a model organism
   (Ensembl/BioMart or UniProt → `prepare_ensembl` / `prepare_uniprot`, which
   emit the same `GO2gene`/`GO2name` files, hierarchy-expanded, and KEGG tables
   too when `kegg_org` is set).

All three `workflow_dispatch` workflows publish dated, immutable Releases:
`build-rna-annotation.yml` (transcriptome inputs from `rna_inputs/<project>/`) and
`build-model-annotation.yml` (`source: ensembl|uniprot`) build the gene-based
`annotation_dir` bundle; `run_ensembl.py` / `run_uniprot.py` are their local CLIs.
The recurring gotcha for the gene-based path: the table gene IDs must match the
RNA counts' `gene_id` (the pipeline warns on <5% overlap); `--no-strip-isoform`
controls `_iN` collapsing.

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
- KEGG's `/conv/` has **no `ensembl` gene database** (`conv/ensembl/<org>` → HTTP
  400); its gene outside-DBs are only `ncbi-geneid`, `ncbi-proteinid`, `uniprot`.
  So `prepare_ensembl` keys `KEGG_pathway2gene.tab` on Ensembl gene ids by
  *bridging*: KEGG gene → NCBI gene id (`conv/ncbi-geneid/<org>`) → Ensembl gene
  id (BioMart NCBI cross-ref, attr name varies by division — see
  `_fetch_ncbi_xref`). Organisms whose KEGG gene ids already ARE Ensembl locus
  codes (e.g. Arabidopsis `ath`) also resolve via the direct `ext_id_universe`
  fallback in `prepare_kegg_by_org`. `prepare_uniprot` uses `conv/uniprot/<org>`
  directly (no bridge). Do not "restore" a bare `conv/ensembl/<org>` call.
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

## Verification (no test CI)

The GitHub Actions here are **build-and-publish** `workflow_dispatch` pipelines
(they release artifacts), NOT a test CI — nothing runs the suite automatically,
so run it yourself.

Offline: run the test suite. The parser / writer / compound-set tests need no
network and no heavy deps; the model-assembly + mummichog-smoke tests skip
automatically unless the scoped optional deps are installed.

```bash
python -m pytest tests/ -q
```

To sanity-check a change end-to-end, run the `run_*` scripts on the `examples/`
inputs and confirm the expected output files appear in `results/` with the
correct headers:

```bash
rm -rf results
python scripts/run_kegg_nonmodel.py --kaas examples/query.ko.txt --out results --cache data
python scripts/run_go.py --go-table examples/trinotate_go.txt --out results --cache data
```
