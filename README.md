# multiomic-annotation-prep

Prepare **organism-specific annotation files for enrichment analysis** (KEGG & GO)
so they can be consumed directly by the `multiomic-core` pipeline and the
`Neat_RNA-Seq` / `Neat_Proteomics` workflows.

Unlike the manual approach in
[Neat_Annotation](https://github.com/veredcc/Neat_Annotation), this repo
**downloads the required source files automatically** from the original
resources (KEGG REST API, Gene Ontology, Ensembl/BioMart, UniProt) and produces
ready-to-use output files.

> 🌿 **מדריך צוות (Hebrew team guide):** [`docs/team-guide-he.md`](docs/team-guide-he.md)
> — how to build artifacts, run the pipeline, and add a new organism, end to end.

---

## What it produces

For **enrichment analysis** with `clusterProfiler::enricher()` (Neat workflows):

| File | Columns | Source module |
|------|---------|---------------|
| `KEGG_pathway2gene.tab` | `v1`, `index` (pathway, gene) | KEGG |
| `KEGG_pathway2name.tab` | `pathway`, `info` | KEGG |
| `GO2gene_BP.tab` / `_MF.tab` / `_CC.tab` | `GO`, `Gene` | GO |
| `GO2name_BP.tab` / `_MF.tab` / `_CC.tab` | `GO`, `Term` | GO |

For **multiomic-core** (GMT format for non-model organisms):

| File | Content |
|------|---------|
| `KEGG_pathway.gmt` | KEGG pathway gene sets |
| `GO_BP.gmt` / `GO_MF.gmt` / `GO_CC.gmt` | GO gene sets per namespace |

For **descriptive annotation** (result tables / Excel):

| File | Content |
|------|---------|
| `KEGG_annot_genes.txt` | per-gene KO, names, EC, pathways |

For **mummichog** metabolomics pathway analysis (a **compound-centric**
metabolic model, distinct from the gene-set files above):

| File | Content |
|------|---------|
| `<org>_<source>_<date>.json` | metabolic model (compounds + reactions + pathways) for `mummichog -n` |
| `<org>_<source>_<date>.manifest.json` | provenance, counts, sha256, validation |

For **multiomic-core** metabolomics **ID-based** enrichment (ORA / GSEA / QEA over
compound IDs, distinct from mummichog's m/z-based analysis):

| File | Content |
|------|---------|
| `<org>_<source>_<date>.compound_pathway.gmt` | KEGG pathway **compound** sets (GMT) |
| `<org>_<source>_<date>.pathway2compound.tab` | readable table: pathway id/name, compound id/name |

---

## Primary use case: non-model organism (transcriptome / KAAS)

For an organism **without a reference genome** (e.g. a Trinity transcriptome),
KEGG Orthology (KO) assignments come from the
[KEGG-KAAS server](https://www.genome.jp/kegg/kaas/), and GO assignments
typically come from [Trinotate](https://github.com/Trinotate/Trinotate).

This repo **auto-downloads everything else**:

* KEGG: `ko -> name`, `ko -> pathway`, `pathway -> name` (from the KEGG REST API)
* GO:   term names + hierarchy (`go-basic.obo`), and expands direct annotations
  to parental terms using the `is_a` relationship graph

You only provide two organism-specific inputs:

1. `query.ko.txt` — KAAS output (transcript/gene `<TAB>` K number)
2. a Trinotate-style GO table — gene `<TAB>` comma-separated GO IDs *(optional)*

---

## Quick start

```bash
# 1. install dependencies (once)
pip install requests pyyaml

# 2. non-model KEGG (auto-downloads KEGG REST files, then processes KAAS output)
python scripts/run_kegg_nonmodel.py \
    --kaas       examples/query.ko.txt \
    --out        results \
    --cache      data

# 3. non-model GO (auto-downloads GO term names, expands hierarchy)
python scripts/run_go.py \
    --go-table   examples/trinotate_go.txt \
    --out        results \
    --cache      data

# or run everything driven by config/config.yml
python scripts/run_all.py --config config/config.yml
```

All network downloads are **cached** under `data/` and re-used on the next run.
Delete the cache (or pass `--refresh`) to force a fresh download — recommended
periodically, since KEGG/GO are updated frequently.

---

## Metabolic model for mummichog (compound-centric)

A separate, self-contained module builds an organism-specific **metabolic model**
for the [mummichog](http://mummichog.org) metabolomics pathway tool. Unlike the
gene-set files above, this model is **compound-centric**: compounds carry a
neutral formula + neutral monoisotopic mass, reactions carry real
substrate/product links, and pathways carry lists of reactions. It is pulled
directly from KEGG compound/reaction/pathway entities and serialized to the
[metDataModel](https://github.com/shuzhao-li/metDataModel) shape that
`mummichog -n <model>.json` consumes (verified against mummichog 2.7.0).

Both supported inputs reduce to the **same pipeline** —
`KO list -> reactions -> compounds -> pathways` — differing only in where the KO
list comes from. (KEGG does not link organism genes directly to reactions, so
the organism's reactions are resolved via KO: `link/ko/<code>` intersected with
`link/reaction/ko`.)

> **Run this in a clean virtualenv, not `conda base`.** A broken Intel MKL in an
> Anaconda `base` environment can crash the `mummichog -n` validation step (numpy
> import fault). Use `python -m venv .venv && source .venv/bin/activate` (or a
> dedicated conda env) before installing.

```bash
# scoped optional deps (do NOT affect the gene-set modules)
pip install -r requirements-mummichog.txt

# (A) from a KEGG organism code (e.g. cre = Chlamydomonas reinhardtii,
#     used here as a surrogate for Coelastrella, which is not in KEGG)
python scripts/run_mummichog_model.py \
    --kegg-code       cre \
    --model-organism  "Chlamydomonas reinhardtii" \
    --target-organism "Coelastrella sp." \
    --out results --cache data --validate

# (B) from a KAAS KO list, for a non-model organism not in KEGG
python scripts/run_mummichog_model.py \
    --source kaas --kaas examples/query.ko.txt \
    --model-organism "Coelastrella sp." --model-kegg-code coel \
    --out results --cache data --validate
```

Path (A) writes `cre_kegg_<date>.json` + `cre_kegg_<date>.manifest.json`. The
manifest records that the model organism is **cre (a relative)**, not
Coelastrella, so the surrogate is transparent in results/publication, along with
KO->reaction coverage. `--validate` additionally runs `mummichog -n` on a
synthetic feature table and records the outcome. The input source (`kegg_org` vs
`kaas`) is also a config choice under `mummichog_model:` in `config/config.yml`.
Path (B) uses KEGG reference pathways (`map#####`), since a non-model organism
has none of its own.

---

## Compound sets for ID-based metabolomics enrichment (GMT)

The mummichog model (above) drives **m/z-based** enrichment. For **ID-based**
enrichment (ORA / GSEA / QEA over compound IDs) `multiomic-core` reads a **GMT of
pathway compound sets**. This is built from the **same organism source** as the
model (`KO list -> reactions -> compounds -> pathways`), so both enrichment paths
describe the same pathways and the same biology. Because ID-based enrichment
needs no mass, the compound set is a **superset** of the model's compounds: no
mass filter is applied, so compounds the model drops (no computable neutral mass)
are still included. Light deps only — `requests` + `pyyaml`, no `metDataModel` /
`mass2chem`.

```bash
# same organism source as the model, same file stem
python scripts/run_kegg_compound_sets.py \
    --kegg-code       cre \
    --model-organism  "Chlamydomonas reinhardtii" \
    --target-organism "Coelastrella sp." \
    --out results --cache data
```

This writes `cre_kegg_<date>.compound_pathway.gmt` (consumed by multiomic-core's
`read_gmt` / `load_gene_sets`), `cre_kegg_<date>.pathway2compound.tab` (a readable
table), and a small `cre_kegg_<date>.compound_sets.manifest.json` sidecar
(files + sha256 + counts + provenance).

To produce the model **and** these companions in one pass — one download, one
KEGG snapshot — add `--emit-compound-sets` to `run_mummichog_model.py` (or set
`mummichog_model.emit_compound_sets: true` in `config/config.yml`); the companions
are then recorded in the model manifest under `companion_files`. Both source
inputs (`kegg_org` / `kaas`) work here too, under `kegg_compound_sets:` in the
config.

---

## Runbook — build & publish organism artifacts

The commands below **download live from KEGG** (`rest.kegg.jp`), so run them where
that host is reachable: your **local machine** (normal internet — nothing to
configure), or a cloud/CI environment whose network policy allows `rest.kegg.jp`.
A restricted sandbox that blocks KEGG egress will fail with a proxy `403`.

### 1. Set up (once)

```bash
git clone https://github.com/multiomics-center-Israel/multiomics-annotation-prep.git
cd multiomics-annotation-prep
python -m venv .venv && source .venv/bin/activate     # Windows: .venv\Scripts\activate
pip install requests pyyaml                            # gene-set + compound-set modules
pip install -r requirements-mummichog.txt              # model: metDataModel / mass2chem / mummichog
```

### 2. Build (model + compound-set companions, one KEGG snapshot each)

`--emit-compound-sets` writes the model **and** the ID-based `.gmt` + `.tab` from a
single `load_source`. `--validate` runs a `mummichog -n` smoke test. These green
algae stand in for *Coelastrella sp.* (not in KEGG); drop `--target-organism` to
record an organism as itself (`model_is_surrogate: false`).

```bash
# cre — Chlamydomonas reinhardtii
python scripts/run_mummichog_model.py --kegg-code cre \
    --model-organism "Chlamydomonas reinhardtii" --target-organism "Coelastrella sp." \
    --emit-compound-sets --validate --out results --cache data

# cvr — Chlorella variabilis
python scripts/run_mummichog_model.py --kegg-code cvr \
    --model-organism "Chlorella variabilis" --target-organism "Coelastrella sp." \
    --emit-compound-sets --validate --out results --cache data

# mng — Monoraphidium neglectum
python scripts/run_mummichog_model.py --kegg-code mng \
    --model-organism "Monoraphidium neglectum" --target-organism "Coelastrella sp." \
    --emit-compound-sets --validate --out results --cache data
```

Each run writes to `results/` (per organism, `<code>_kegg_<YYYYMMDD>`):
`.json`, `.manifest.json`, `.compound_pathway.gmt`, `.pathway2compound.tab`.

### 3. Publish to a GitHub Release (do **not** commit the artifacts)

`results/` and `data/` are git-ignored on purpose — per `MODEL_CONTRACT.md`
artifacts are **dated, immutable, and published as GitHub Releases** (one tag per
organism/version), not checked into the repo. Upload via the web UI (Releases →
*Draft a new release* → drag the files in) or the CLI:

```bash
# example for one organism (repeat per organism/date)
gh release create cre_kegg_$(date -u +%Y%m%d) results/cre_kegg_*.json \
    results/cre_kegg_*.manifest.json results/cre_kegg_*.compound_pathway.gmt \
    results/cre_kegg_*.pathway2compound.tab \
    --title "cre_kegg_$(date -u +%Y%m%d)" --notes "KEGG-derived model + compound sets"
```

`multiomic-core` then pins the Release download URLs in its config
(`mummichog.model_ref` for the model, `enrichment.gmt_file` for the GMT).

### Automated: one-click build + publish (recommended)

Steps 1–3 are wrapped in a manually-dispatched GitHub Actions workflow,
`.github/workflows/publish-organism-artifacts.yml`. Trigger it from the repo's
**Actions** tab (*Build & publish organism artifacts* → **Run workflow**) with
the organism codes (default `cre,cvr,mng`), a surrogate `target_organism`
(default `Coelastrella sp.`), and an optional `date`. Each organism builds on a
GitHub-hosted runner (open internet → KEGG reachable) and is published as a
Release via `GITHUB_TOKEN`. Releases are immutable: a re-run on the same date
fails fast rather than overwriting — bump the `date` for a rebuild. The workflow
tolerates the known flaky `mummichog` colored-map crash (`rescale_color`) as long
as the pathway + module analysis tables were produced.

---

## Source modules

| Module | File | Use when | Auto-downloads |
|--------|------|----------|----------------|
| KEGG (non-model) | `src/prepare_kegg_nonmodel.py` | transcriptome + KAAS | KEGG REST tables |
| GO (non-model)   | `src/prepare_go.py`           | Trinotate GO         | GO term names + hierarchy |
| Ensembl/BioMart  | `src/prepare_ensembl.py`      | genome in Ensembl    | BioMart attributes, KEGG map |
| UniProt          | `src/prepare_uniprot.py`      | proteome in UniProt  | UniProt REST, KEGG map |
| mummichog model  | `src/prepare_mummichog_model.py` | metabolomics, m/z-based (`mummichog -n`) | KEGG compounds/reactions/pathways |
| KEGG compound sets | `src/prepare_kegg_compound_sets.py` | metabolomics, ID-based (GMT) | (reuses the model's KEGG source) |

The non-model KEGG + GO modules are the primary, fully-featured path. Ensembl
and UniProt modules cover model organisms and share the same writers/output
formats.

---

## Layout

```
multiomic-annotation-prep/
├── src/
│   ├── utils.py                 # logging, caching, download helpers, writers
│   ├── download_kegg.py         # KEGG REST API downloads (cached)
│   ├── prepare_kegg_nonmodel.py # KAAS -> KEGG enrichment + annot files
│   ├── download_go.py           # GO term/hierarchy downloads (cached)
│   ├── prepare_go.py            # GO table -> expanded GO enrichment files
│   ├── prepare_ensembl.py       # Ensembl/BioMart path (model organisms)
│   ├── prepare_uniprot.py       # UniProt path (proteomics)
│   ├── masses.py                # neutral monoisotopic mass (via mass2chem)
│   ├── kegg_entities.py         # KEGG source loading (KO -> reactions -> compounds -> pathways)
│   ├── prepare_mummichog_model.py     # KEGG -> mummichog metabolic model (JSON)
│   └── prepare_kegg_compound_sets.py  # KEGG -> compound-set GMT + table (ID-based)
├── scripts/
│   ├── run_kegg_nonmodel.py
│   ├── run_go.py
│   ├── run_mummichog_model.py
│   ├── run_kegg_compound_sets.py
│   ├── run_ensembl.py           # Ensembl/BioMart GO (+ KEGG) tables (model organisms)
│   ├── run_uniprot.py           # UniProt GO (+ KEGG) tables (by taxon id)
│   └── run_all.py               # config-driven, runs selected modules
├── config/config.yml
├── examples/                    # tiny inputs so it runs out of the box
├── tests/                       # pytest (parsers, masses, model, mummichog smoke)
├── data/                        # download cache (git-ignored)
└── results/                     # output files (git-ignored)
```

## Dependencies

Python (>= 3.8) with: `requests`, `pyyaml` (install via pip).
The Ensembl module additionally requires `pybiomart`.
The mummichog metabolic-model module has its own **scoped, pinned** optional
deps in `requirements-mummichog.txt` (`metDataModel`, `mass2chem`, and
`mummichog` for validation) — the gene-set modules do not need them.

## Credit

Formats and workflow follow Vered Chalifa-Caspi's
[Neat_Annotation](https://github.com/veredcc/Neat_Annotation)
(Bioinformatics Core Facility, Ben-Gurion University).
