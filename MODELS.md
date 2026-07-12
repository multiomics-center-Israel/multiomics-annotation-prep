# Published models

Registry of organism-specific metabolic models built by this repo for
`mummichog -n`. **Before building a new model, check this table** -- if your
organism (or a suitable surrogate) is already listed, reuse its **Model URL** and
**sha256** in your pipeline config instead of rebuilding.

Rows are appended automatically by `scripts/publish_model.py` on each release
(idempotent on the release tag). Published artifacts are immutable -- a rebuild
is a new dated row, never an edit of an existing one.

| Target organism | Model organism (KEGG code) | Surrogate? | Release tag | Model URL | sha256 | KEGG snapshot | Build date |
|---|---|---|---|---|---|---|---|
| Coelastrella sp. | Chlamydomonas reinhardtii (cre) | yes | cre_kegg_20260711 | https://github.com/multiomics-center-Israel/multiomics-annotation-prep/releases/download/cre_kegg_20260711/cre_kegg_20260711.json | c403c96fbec8df9ae34b828fec01270c8ea3940acc36e4e5ff770868dc8b912b | 2026-07-10 | 2026-07-11 |
