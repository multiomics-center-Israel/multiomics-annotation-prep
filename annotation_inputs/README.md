# Gene-based annotation inputs

Project-specific inputs for the **Build gene annotation (RNA / proteomics)**
workflow (`.github/workflows/build-gene-annotation.yml`). Unlike the
compound-centric organism artifacts (which come from a KEGG organism code),
these enrichment files are keyed on **your project's IDs** (genes, transcripts
or proteins), so they are built from your project's own annotation files.

## Layout

Put each project's inputs in its own folder, using these **exact** filenames:

```
annotation_inputs/
└── <project>/
    ├── query.ko.txt        # KAAS output (id <TAB> KO)              -> KEGG
    ├── trinotate_go.txt    # per-id GO table (id <TAB> GO ids)      -> GO
    ├── <name>.gaf[.gz]     # one GO Annotation File, e.g. GiardiaDB-71_..._GO.gaf -> GO
    └── id_map.txt[.gz]     # optional: re-keys GAF ids (see below)
```

- At least one of `query.ko.txt`, `trinotate_go.txt`, the GAF is needed; you
  get the databases whose input is present.
- `trinotate_go.txt` and the GAF are alternatives -- keep only one.
- The GAF: VEuPathDB sites (GiardiaDB, PlasmoDB, ToxoDB, ...) need you to log in
  to download, so download the organism's GO GAF in the browser and upload it.
  **Keep the downloaded file name** (e.g. `GiardiaDB-71_GintestinalisAssemblageAWB_GO.gaf`):
  it carries the source release, and the Release records it. Only one GAF per
  folder.
- Re-keying the GAF ids (when the GAF's gene ids differ from the ids in your
  data): **prefer the workflow's `ncbi_assembly` input** (e.g.
  `GCA_000002435.2`) -- the runner downloads that assembly's public NCBI feature
  table (`locus_tag` -> protein `product_accession`) and records its URL +
  sha256, so nothing extra is committed by hand. For any other mapping, commit
  `id_map.txt[.gz]`: a two-column TSV `gaf_id <TAB> data_id` or an NCBI
  `*_feature_table.txt(.gz)`. Use one or the other, not both.

The IDs that end up in the output must match the IDs in your data (RNA counts
`gene_id`, or the proteomics `Protein.Group`).

## Run

GitHub **Actions** -> **Build gene annotation (RNA / proteomics)** -> **Run workflow**:

- `project` -- the folder name above (letters/digits/underscore/hyphen)
- `omics` -- `rna` or `proteomics` (sets the release tag and notes)
- `ncbi_assembly` -- optional, for the GAF only (see above)
- `date` -- empty = today (UTC)
- `strip_isoform` -- RNA only: `true` strips `_iN` isoform suffixes from KAAS
  ids to gene level (default); ignored for proteomics

Output: a Release `rna_annot_<project>_<date>` or `prot_annot_<project>_<date>`
with every built `.tab` and `.gmt` file, a convenience zip and `SHA256SUMS.txt`.
The release notes record the commit, the run, each input's sha256 and the
GAF's name and header (source database + release); GAF builds also attach
`GO_provenance.json` (GAF, id map and GO ontology versions + checksums).

- RNA: unzip into your project and point `modes.rna.enrichment.annotation_dir`
  at that folder.
- Proteomics: list the `GO_*.gmt` (and `KEGG_pathway.gmt`) files under
  `modes.proteomics.pathway.gmt_file`.

> עברית: מעלים ל‑`annotation_inputs/<project>/` את `query.ko.txt` ו/או
> `trinotate_go.txt` או the GAF (+ `id_map.txt` אופציונלי), מריצים את ה‑workflow
> עם שם התיקייה ו‑`omics`, ומקבלים Release עם קבצי ה‑`.tab` וה‑`.gmt`.
