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
    ├── go.gaf              # GO Annotation File, e.g. from GiardiaDB -> GO
    └── id_map.txt[.gz]     # optional: re-keys go.gaf ids (see below)
```

- At least one of `query.ko.txt`, `trinotate_go.txt`, `go.gaf` is needed; you
  get the databases whose input is present.
- `trinotate_go.txt` and `go.gaf` are alternatives -- keep only one.
- `go.gaf`: VEuPathDB sites (GiardiaDB, PlasmoDB, ToxoDB, ...) need you to log in
  to download, so download the organism's GO GAF in the browser and upload it.
- `id_map.txt[.gz]` applies to `go.gaf` only. Use it when the GAF's gene ids
  differ from the ids in your data: an NCBI `*_feature_table.txt.gz` (maps
  `locus_tag` -> protein `product_accession`), or a two-column TSV
  `gaf_id <TAB> data_id`.

The IDs that end up in the output must match the IDs in your data (RNA counts
`gene_id`, or the proteomics `Protein.Group`).

## Run

GitHub **Actions** -> **Build gene annotation (RNA / proteomics)** -> **Run workflow**:

- `project` -- the folder name above (letters/digits/underscore/hyphen)
- `omics` -- `rna` or `proteomics` (sets the release tag and notes)
- `date` -- empty = today (UTC)
- `strip_isoform` -- RNA only: `true` strips `_iN` isoform suffixes from KAAS
  ids to gene level (default); ignored for proteomics

Output: a Release `rna_annot_<project>_<date>` or `prot_annot_<project>_<date>`
with every built `.tab` and `.gmt` file, a convenience zip and `SHA256SUMS.txt`.

- RNA: unzip into your project and point `modes.rna.enrichment.annotation_dir`
  at that folder.
- Proteomics: list the `GO_*.gmt` (and `KEGG_pathway.gmt`) files under
  `modes.proteomics.pathway.gmt_file`.

> עברית: מעלים ל‑`annotation_inputs/<project>/` את `query.ko.txt` ו/או
> `trinotate_go.txt` או `go.gaf` (+ `id_map.txt` אופציונלי), מריצים את ה‑workflow
> עם שם התיקייה ו‑`omics`, ומקבלים Release עם קבצי ה‑`.tab` וה‑`.gmt`.
