# multiomic-annotation-prep

Prepare **organism-specific annotation files for enrichment analysis** (KEGG & GO)
so they can be consumed directly by the `multiomic-core` pipeline.

Unlike the manual approach in
[Neat_Annotation](https://github.com/veredcc/Neat_Annotation), this repo
**downloads the required source files automatically** from the original
resources (KEGG REST API, Gene Ontology, Ensembl/BioMart, UniProt) and produces
ready-to-use `.tab` files for `clusterProfiler`.

---

## What it produces

For **enrichment analysis** (over-representation / GSEA with `clusterProfiler`):

| File | Columns | Source module |
|------|---------|---------------|
| `KEGG_pathway2gene.tab` | `v1`, `index` (pathway, gene) | KEGG |
| `KEGG_pathway2name.tab` | `pathway`, `info` | KEGG |
| `GO2gene_BP.tab` / `_MF.tab` / `_CC.tab` | `GO`, `Gene` | GO |
| `GO2name_BP.tab` / `_MF.tab` / `_CC.tab` | `GO`, `Term` | GO |

For **descriptive annotation** (result tables / Excel):

| File | Content |
|------|---------|
| `KEGG_annot_genes.txt` | per-gene KO, names, EC, pathways |

> **Format note.** These column names/formats match the
> `Neat_RNA-Seq` / `Neat_Proteomics` (clusterProfiler) convention. If
> `multiomic-core` expects a different layout, adjust the writers in
> `R/*.R` (they are isolated in small `write_*` helpers) — nothing else changes.

---

## Primary use case: non-model organism (transcriptome / KAAS)

For an organism **without a reference genome** (e.g. a Trinity transcriptome),
KEGG Orthology (KO) assignments come from the
[KEGG-KAAS server](https://www.genome.jp/kegg/kaas/), and GO assignments
typically come from [Trinotate](https://github.com/Trinotate/Trinotate).

This repo **auto-downloads everything else**:

* KEGG: `ko -> name`, `ko -> pathway`, `pathway -> name` (from the KEGG REST API)
* GO:   term names + hierarchy (`go-basic.obo` / `GO.db`), and expands direct
  annotations to parental terms with `clusterProfiler::buildGOmap`

You only provide two organism-specific inputs:

1. `query.ko.txt` — KAAS output (transcript/gene `<TAB>` K number)
2. a Trinotate-style GO table — gene `<TAB>` comma-separated GO IDs *(optional)*

---

## Quick start

```bash
# 1. install dependencies (once)
Rscript -e 'source("R/install_deps.R")'

# 2. non-model KEGG (auto-downloads KEGG REST files, then processes KAAS output)
Rscript scripts/run_kegg_nonmodel.R \
    --kaas       examples/query.ko.txt \
    --out        results \
    --cache      data

# 3. non-model GO (auto-downloads GO term names, expands hierarchy)
Rscript scripts/run_go.R \
    --go-table   examples/trinotate_go.txt \
    --out        results \
    --cache      data

# or run everything driven by config/config.yml
Rscript scripts/run_all.R --config config/config.yml
```

All network downloads are **cached** under `data/` and re-used on the next run.
Delete the cache (or pass `--refresh`) to force a fresh download — recommended
periodically, since KEGG/GO are updated frequently.

---

## Source modules

| Module | File | Use when | Auto-downloads |
|--------|------|----------|----------------|
| KEGG (non-model) | `R/prepare_kegg_nonmodel.R` | transcriptome + KAAS | KEGG REST tables |
| GO (non-model)   | `R/prepare_go.R`           | Trinotate GO         | GO term names + hierarchy |
| Ensembl/BioMart  | `R/prepare_ensembl.R`      | genome in Ensembl    | BioMart attributes, KEGG map |
| UniProt          | `R/prepare_uniprot.R`      | proteome in UniProt  | UniProt REST, KEGG map |

The non-model KEGG + GO modules are the primary, fully-featured path. Ensembl
and UniProt modules cover model organisms and share the same writers/output
formats.

---

## Layout

```
multiomic-annotation-prep/
├── R/
│   ├── install_deps.R          # install Bioconductor/CRAN deps
│   ├── utils.R                 # logging, caching, download helpers, writers
│   ├── download_kegg.R         # KEGG REST API downloads (cached)
│   ├── prepare_kegg_nonmodel.R # KAAS -> KEGG enrichment + annot files
│   ├── download_go.R           # GO term/hierarchy downloads (cached)
│   ├── prepare_go.R            # GO table -> expanded GO enrichment files
│   ├── prepare_ensembl.R       # Ensembl/BioMart path (model organisms)
│   └── prepare_uniprot.R       # UniProt path (proteomics)
├── scripts/
│   ├── run_kegg_nonmodel.R
│   ├── run_go.R
│   └── run_all.R               # config-driven, runs selected modules
├── config/config.yml
├── examples/                   # tiny inputs so it runs out of the box
├── data/                       # download cache (git-ignored)
└── results/                    # output .tab files (git-ignored)
```

## Dependencies

R (>= 4.0) with: `KEGGREST`, `clusterProfiler`, `GO.db`, `AnnotationDbi`
(Bioconductor); `biomaRt` (for the Ensembl module); `httr`, `yaml`, `optparse`
(CRAN). Run `R/install_deps.R` to install them.

## Credit

Formats and workflow follow Vered Chalifa-Caspi's
[Neat_Annotation](https://github.com/veredcc/Neat_Annotation)
(Bioinformatics Core Facility, Ben-Gurion University).
=======
