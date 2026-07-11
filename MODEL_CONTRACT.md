# Model Contract — the seam between the builder repo and the pipeline

This is the **only** thing the two repos share. The builder produces artifacts
that satisfy this contract; the pipeline consumes them. Neither needs to know the
other's internals.

- **Builder repo (A):** builds + validates + publishes organism metabolic models.
  Heavy deps (KEGG / GEM / Bioconductor or Python metDataModel). Own test suite.
- **Pipeline repo (B):** pins a model version, verifies its checksum, feeds it to
  `mummichog -n <model>.json`. Light deps only. No build logic.

Anything not written here is a private implementation detail of one side.

---

## Artifact 1 — the model file (`<org>_<source>_<YYYYMMDD>.json`)

A mummichog-compatible metabolic model in **metDataModel** `MetabolicModel`
serialization. Target shape (verify exact field names against the installed
`metDataModel` version and a known-good `mummichog -n` example — see Acceptance):

```jsonc
{
  "id": "cre_kegg_20260711",           // model id (matches filename stem)
  "version": "20260711",
  "meta_data": {
    "model_organism": "Chlamydomonas reinhardtii",  // what the model IS
    "model_kegg_code": "cre",
    "target_organism": "Coelastrella sp.",          // the biology it stands in for
    "model_is_surrogate": true,                     // true when model != target
    "source": "KEGG REST",             // or "GEM:<id>@<version>"
    "source_version": "KEGG release 110.0",
    "builder_version": "<git SHA of builder repo>"
  },
  "list_of_compounds": [
    {
      "id": "C00031",
      "name": "D-Glucose",
      "neutral_formula": "C6H12O6",     // REQUIRED — neutral, not a salt
      "neutral_mono_mass": 180.063388   // REQUIRED — monoisotopic mass of the neutral form
    }
  ],
  "list_of_reactions": [
    {
      "id": "R00299",
      "reactants": ["C00031"],          // compound ids
      "products":  ["C00668"],
      "enzymes":   ["2.7.1.1"]          // optional (EC numbers)
    }
  ],
  "list_of_pathways": [
    {
      "id": "cre00010",
      "name": "Glycolysis / Gluconeogenesis",
      "list_of_reactions": ["R00299"]   // reaction ids
    }
  ]
}
```

### Hard requirements (these are what break the science if wrong)
1. Every compound carries a **neutral_formula** and a **neutral_mono_mass**
   (monoisotopic mass of the neutral molecule — mummichog adds adducts itself).
   No salt formulas.
2. Reactions reference compound ids that exist in `list_of_compounds`.
3. Pathways reference reaction ids that exist in `list_of_reactions`.
4. Reactions are built from **real substrate/product links** (from KEGG/GEM),
   NOT "compounds that co-occur in a pathway" heuristics.
5. Compound ids are stable and consistent across the three lists.

---

## Artifact 2 — the sidecar manifest (`<same-stem>.manifest.json`)

Fully under our control; the pipeline reads this to pin + verify.

```jsonc
{
  "model_file": "cre_kegg_20260711.json",
  "sha256": "<64-hex checksum of the model file>",
  "model_organism": "Chlamydomonas reinhardtii",
  "model_kegg_code": "cre",
  "target_organism": "Coelastrella sp.",   // biology the model stands in for
  "model_is_surrogate": true,              // true when model_organism != target_organism
  "source": "KEGG REST",
  "source_version": "KEGG release 110.0",
  "build_timestamp_utc": "2026-07-11T09:00:00Z",
  "builder_git_sha": "<commit>",
  "counts": { "compounds": 0, "reactions": 0, "pathways": 0 },
  "validation": {
    "mummichog_version_tested": "2.7.0",
    "loads_via_-n": true,
    "smoke_run_exit_0": true,
    "mass_spotcheck_passed": true
  }
}
```

---

## Publishing & versioning

- Artifacts are **immutable** once published. A rebuild = a new dated version.
- Publish to **GitHub Releases** (tag per model version) and/or **Zenodo**
  (gives a citable DOI per version — recommended for anything used in a paper).
- Never a "latest" pointer in the pipeline.

## Pipeline consumption (repo B)

Pipeline config pins an exact model, e.g.:

```yaml
modes:
  metabolomics:
    organism: "Coelastrella sp."          # the real sample organism
mummichog:
  model_ref:
    model_organism: "Chlamydomonas reinhardtii (cre)"   # surrogate model actually used
    url: "https://github.com/<org>/<builder-repo>/releases/download/cre_kegg_20260711/cre_kegg_20260711.json"
    sha256: "<expected checksum>"
```

Pipeline behavior:
1. If no `model_ref` → use mummichog's built-in `human_mfn` (default).
2. If `model_ref` set → fetch (cached), **verify sha256** (fail loudly on
   mismatch), pass the local path to `mummichog.main -n <path>`.

## Acceptance criteria (the real contract)

A model artifact is valid **iff**:
1. It loads cleanly via `mummichog.main -n <model>.json` on the pinned mummichog
   version (2.7.0) and completes a run on a small synthetic feature table
   (exit 0, produces pathway + module tables).
2. A spot-check of ≥5 known compounds' `neutral_mono_mass` matches an independent
   source (e.g. PubChem/KEGG) within 1 mDa.
3. Compound/reaction/pathway counts are within a sane range of KEGG for that
   organism (documented in the manifest).
4. The sidecar manifest is present, well-formed, and its sha256 matches the file.

The loader in criterion #1 is authoritative: if mummichog's `-n` parses it and
runs, the field names are right; if not, fix the builder to match metDataModel.
