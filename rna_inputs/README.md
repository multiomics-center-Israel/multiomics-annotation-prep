# RNA enrichment inputs

Project-specific inputs for the **Build RNA annotation tables** workflow
(`.github/workflows/build-rna-annotation.yml`). Unlike the compound-centric
organism artifacts (which come from a KEGG organism code), RNA enrichment tables
are **gene-ID based**, so they are built from *your project's* KAAS + Trinotate
outputs.

## Layout

Put each project's inputs in its own folder, using these **exact** filenames:

```
rna_inputs/
└── <project>/
    ├── query.ko.txt        # KAAS output (gene <TAB> KO)       -> KEGG tables
    └── trinotate_go.txt    # Trinotate GO (gene <TAB> GO ids)  -> GO tables
```

Either file may be omitted — you get the databases whose input is present. The
gene IDs in these files must match your RNA counts matrix `gene_id` column.

## Run

GitHub **Actions** → **Build RNA annotation tables** → **Run workflow**:

- `project` — the folder name above (letters/digits/underscore/hyphen)
- `date` — empty = today (UTC)
- `strip_isoform` — `true` strips `_iN` isoform suffixes to gene level (default);
  set `false` to keep transcript-level IDs if that is what your counts use

Output: a Release `rna_annot_<project>_<date>` with the `annotation_dir` files
(`KEGG_pathway2gene.tab`, `KEGG_pathway2name.tab`, `GO2gene_{BP,MF,CC}.tab`,
`GO2name_{BP,MF,CC}.tab`, …) plus a convenience zip. Download, unzip into your
RNA project, and point `modes.rna.enrichment.annotation_dir` at that folder.

> עברית: מעלים ל‑`rna_inputs/<project>/` את `query.ko.txt` ו/או `trinotate_go.txt`,
> מריצים את ה‑workflow עם שם התיקייה, ומקבלים Release עם קבצי ה‑annotation לעבוד איתם.
