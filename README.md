# multiomics-annotation-prep

This repo makes the **annotation files that `multiomics-core` needs for
enrichment analysis** — the files that say which genes, proteins or metabolites
belong to which GO term or KEGG pathway.

You don't need to run any code. Everything is done with **GitHub Actions**: you
fill in a short form, a GitHub machine downloads the reference data and builds
the files, and the result is published as a **Release** (a dated, permanent
download page) that you point your pipeline at.

> 🌿 Hebrew team guide: [`docs/team-guide-he.md`](docs/team-guide-he.md)

---

## Which Action do I need?

| I am analysing… | My organism | Use this Action |
|---|---|---|
| **Proteomics** | any organism, when I want GO / KEGG files keyed on my own protein ids (needed for non-model organisms) | [Build gene annotation (RNA / proteomics)](#1-build-gene-annotation-rna--proteomics) |
| **RNA-seq** | non-model (de novo transcriptome, or an organism database such as GiardiaDB) | [Build gene annotation (RNA / proteomics)](#1-build-gene-annotation-rna--proteomics) |
| **RNA-seq** | model organism in Ensembl or UniProt | [Build model-organism annotation](#2-build-model-organism-annotation) |
| **Metabolomics** | any organism KEGG knows (or a close relative) | [Build & publish organism artifacts](#3-build--publish-organism-artifacts-metabolomics) |
| — | I published something by mistake | [Delete releases](#4-delete-releases-maintenance) |

> Proteomics of a **model organism** (human, mouse, …) usually needs nothing from
> this repo: `multiomics-core` loads GO / KEGG for it by itself from
> `annotation.organism`.

**How to run any of them:** open the repo's **Actions** tab → pick the Action on
the left → **Run workflow** (top right) → fill in the fields → **Run workflow**.
When the run turns green, the files are on the repo's **Releases** page.

---

## 1. Build gene annotation (RNA / proteomics)

**What it does:** turns your project's own annotation into GO and KEGG gene sets
keyed on **the ids in your data** (gene ids for RNA, protein ids for proteomics).

### Step 1 — put the inputs in the repo

Create a folder `annotation_inputs/<project>/` (on GitHub: go into
`annotation_inputs`, **Add file → Upload files**, and drag in a *folder* named
after your project). Put in it whichever of these you have:

| File | What it is | Gives you |
|---|---|---|
| `query.ko.txt` | KAAS result: one line per id, `id <TAB> KO` | KEGG |
| `trinotate_go.txt` | GO table: `id <TAB> GO ids` (from Trinotate, eggNOG, …) | GO |
| a `.gaf` file | GO annotation file from an organism database (e.g. GiardiaDB). **Keep the downloaded file name** — it contains the database release, which then appears in the Release record. | GO |
| `id_map.txt` | *optional*, translates the ids of the **`.gaf` and `kegg_org` sources only**: `id_in_file <TAB> id_in_your_data`. It does **not** touch `query.ko.txt` or `trinotate_go.txt` — those must already use the ids in your data. | — |

Rules: use `trinotate_go.txt` **or** a `.gaf`, not both; use `query.ko.txt`
**or** the `kegg_org` field below, not both. You can skip the folder entirely
if all you need is KEGG from `kegg_org`.

### Step 2 — fill in the form

| Field | What to enter |
|---|---|
| `project` | the folder name from step 1 |
| `omics` | `rna` or `proteomics` |
| `kegg_org` | *optional.* A KEGG organism code (e.g. `gla` for Giardia, `cre` for Chlamydomonas). Builds KEGG straight from KEGG's own annotation — **no KAAS needed**. Check the code exists at `https://rest.kegg.jp/list/<code>`. |
| `ncbi_assembly` | *optional.* An NCBI assembly accession (e.g. `GCA_000002435.2`). Use it when your data uses NCBI protein ids (like `KAE8301220.1`) but the GAF / KEGG use gene names (like `GL50803_0017558`): the Action downloads NCBI's table that links them and translates for you. |
| `id_label` | *optional.* A short word for the id type, added to the Release name (e.g. `kae`). Lets you keep two builds of the same project on the same day. Defaults to `ncbi` when `ncbi_assembly` is set. |
| `date` | leave empty (today) |
| `strip_isoform` | RNA only, and only for `query.ko.txt`: ticked turns transcript ids (`..._g1_i1`) into gene ids (`..._g1`). It does **not** change `trinotate_go.txt` — make that file use the same level as your counts yourself. |

### Step 3 — what you get

A Release named `prot_annot_<project>[_<label>]_<date>` (proteomics) or
`rna_annot_<project>[_<label>]_<date>` (RNA). What is in it depends on what you
gave it:

| You gave | You get |
|---|---|
| a GO source (`trinotate_go.txt` or `.gaf`) | `GO_BP.gmt`, `GO_MF.gmt`, `GO_CC.gmt` + `GO2gene_*.tab` / `GO2name_*.tab` — one per GO branch that has annotations (a branch with none is skipped) |
| a KEGG source (`query.ko.txt` or `kegg_org`) | `KEGG_pathway.gmt` + `KEGG_pathway2gene.tab` / `KEGG_pathway2name.tab` |
| a `.gaf` | also `GO_provenance.json` and `Annotation.tab` |
| `kegg_org` | also `KEGG_provenance.json` |

Every Release also has `SHA256SUMS.txt`, a zip of everything, and notes listing
the inputs (with checksums) and the commit and run that built it. Check the
Release's file list before writing your config.

### Step 4 — use it in `multiomics-core`

**Proteomics** reads the `.gmt` files. Download them into `annotation/` inside
your project's **raw-data folder** (`paths.raw`; relative paths below are
relative to it), then:

```yaml
modes:
  proteomics:
    annotation:
      organism: "Giardia lamblia"
      skip_annotation: true        # the GMT uses the raw Protein.Group ids
    pathway:
      gmt_file:
        - annotation/GO_BP.gmt
        - annotation/GO_MF.gmt
        - annotation/GO_CC.gmt
        - annotation/KEGG_pathway.gmt
```

**RNA** reads the `.tab` files. Unzip the Release zip into a folder in your
project and point at it (a relative path here is relative to the **project
folder**, `project.dir`):

```yaml
modes:
  rna:
    enrichment:
      annotation_dir: "data/annotation"   # the folder with the .tab files
```

> ⚠️ **The ids must match your data.** If the GMT has `KAE…` ids and your data
> has `GL50803_…` ids, nothing matches and enrichment finds nothing. When in
> doubt, open one GMT, look at a few ids, and compare with your Protein.Group /
> gene_id column.

**Example — Giardia proteomics:** the GiardiaDB GAF was uploaded to
`annotation_inputs/giardia/`, and the Action was run with `omics: proteomics`,
`kegg_org: gla`, `ncbi_assembly: GCA_000002435.2`. Result: GO + KEGG gene sets
keyed on the `KAE…` protein ids the proteomics data uses.

---

## 2. Build model-organism annotation

**What it does:** for an organism with a reference genome in **Ensembl** or a
proteome in **UniProt**, downloads GO (and optionally KEGG) directly — no input
files needed.

| Field | What to enter |
|---|---|
| `source` | `ensembl` or `uniprot` |
| `label` | a short name for the Release (e.g. `athaliana`) |
| `dataset`, `mart`, `host` | Ensembl only. The BioMart dataset (e.g. `athaliana_eg_gene`); for plants/algae/protists also the matching `mart` (e.g. `plants_mart`) and `host` (e.g. `https://plants.ensembl.org`) |
| `taxon_id` | UniProt only. The NCBI taxon id (e.g. `3702`) |
| `reviewed` | UniProt only. Ticked = reviewed (Swiss-Prot) entries only |
| `kegg_org` | *optional.* KEGG organism code to also build KEGG tables |
| `date` | leave empty |

**You get:** a Release `<label>_<source>_annot_<date>` with the `GO2gene_*` /
`GO2name_*` (and KEGG) **`.tab` tables**.

**Use it in:** the **RNA** pipeline — unzip and set
`modes.rna.enrichment.annotation_dir` to that folder (as in Action 1). Gene ids
are Ensembl gene ids or UniProt accessions, so they must match your counts.

---

## 3. Build & publish organism artifacts (metabolomics)

**What it does:** builds metabolomics pathway files from KEGG for one or more
KEGG organism codes. For an organism KEGG doesn't have (e.g. *Coelastrella*), use
close relatives and name the real organism as the target.

| Field | What to enter |
|---|---|
| `organisms` | comma-separated KEGG codes (default `cre,cvr,mng`) |
| `target_organism` | the organism you actually study, if the codes above are stand-ins (default `Coelastrella sp.`); leave empty if they are the organism itself |

> The Action knows the organism **names** only for `cre`, `cvr` and `mng`. For any other code the files are still built, but the manifest records no organism name. To get it recorded, add the code to the `NAMES` list in `.github/workflows/publish-organism-artifacts.yml` first.
| `date` | leave empty |

**You get:** one Release per organism, `<code>_kegg_<date>`, containing:

- `<code>_kegg_<date>.json` — metabolic model for **mummichog** (m/z-based analysis)
- `<code>_kegg_<date>.compound_pathway.gmt` — pathway compound sets for **ID-based** enrichment
- a manifest with versions and checksums

**Use it in** the **metabolomics** pipeline:

```yaml
modes:
  metabolomics:
    enrichment:
      run_enrichment: true
      gmt_file: annotation/cre_kegg_20260711.compound_pathway.gmt   # downloaded; relative to paths.raw
      mummichog:
        enabled: true
        ionization_mode: positive                                   # pos_default | positive | negative
        model_ref:                                                  # fetched + verified by the pipeline
          url: https://github.com/multiomics-center-Israel/multiomics-annotation-prep/releases/download/cre_kegg_20260711/cre_kegg_20260711.json
          sha256: <the sha256 from the Release>
```

---

## 4. Delete releases (maintenance)

Releases can't be overwritten, so to rebuild something under the same name you
first delete it. Enter the Release name(s) in `tags` (comma-separated) and type
`DELETE` in `confirm`. The Release and its tag are removed.

---

## Good to know

- **Releases are permanent.** Running the same Action twice on the same day
  with the same name fails on purpose. Delete the old Release first, or use a
  different `id_label` / `date`.
- **Only one thing is manual:** files that need a login to download (e.g.
  GiardiaDB) or that come from a web form (KAAS). Everything public — KEGG,
  Gene Ontology, NCBI, Ensembl, UniProt — is downloaded by the Action itself.
- **What is recorded differs per Action.** *Build gene annotation* lists its
  inputs with checksums plus the commit and run, and GAF / `kegg_org` builds
  add database versions in their provenance files. *Build & publish organism
  artifacts* writes a manifest with versions and checksums. *Build
  model-organism annotation* records only the output checksums and the source
  name — note the dataset / taxon you used yourself.

---

# Developer reference

Everything below is for running the scripts locally or changing the code.
`CLAUDE.md` has the conventions and the output contract.

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
| GO from GAF      | `src/prepare_gaf.py`          | a downloaded GAF (e.g. GiardiaDB), optional NCBI id map | GO term names + hierarchy |
| KEGG by organism | `src/prepare_kegg_org_sets.py` | organism KEGG annotates (e.g. `gla`), no KAAS; optional NCBI id map | KEGG REST links + names, writes GMT |
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
│   ├── prepare_gaf.py           # GAF (+ id map) -> expanded GO enrichment files
│   ├── prepare_kegg_org_sets.py # KEGG org code (+ id map) -> KEGG tables + GMT
│   ├── id_maps.py               # shared id maps (NCBI feature table / TSV)
│   ├── prepare_ensembl.py       # Ensembl/BioMart path (model organisms)
│   ├── prepare_uniprot.py       # UniProt path (proteomics)
│   ├── masses.py                # neutral monoisotopic mass (via mass2chem)
│   ├── kegg_entities.py         # KEGG source loading (KO -> reactions -> compounds -> pathways)
│   ├── prepare_mummichog_model.py     # KEGG -> mummichog metabolic model (JSON)
│   └── prepare_kegg_compound_sets.py  # KEGG -> compound-set GMT + table (ID-based)
├── scripts/
│   ├── run_kegg_nonmodel.py
│   ├── run_go.py
│   ├── run_gaf.py
│   ├── run_kegg_org_sets.py
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
The mummichog metabolic-model module has its own **scoped, pinned** optional
deps in `requirements-mummichog.txt` (`metDataModel`, `mass2chem`, and
`mummichog` for validation) — the gene-set modules do not need them.

## Credit

Formats and workflow follow Vered Chalifa-Caspi's
[Neat_Annotation](https://github.com/veredcc/Neat_Annotation)
(Bioinformatics Core Facility, Ben-Gurion University).
