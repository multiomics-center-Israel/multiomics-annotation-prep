"""Model registry helpers: maintain MODELS.md, the human-friendly index of
published mummichog metabolic models.

``scripts/publish_model.py`` calls :func:`append_model_row` after a successful
release so colleagues can discover an already-built model (and its URL + sha256)
instead of rebuilding. The append is idempotent on the release tag. All fields
come from the sidecar manifest produced by ``prepare_mummichog_model``.

No network, no heavy deps -- pure text/JSON, so this is safe to import and test
anywhere.
"""

import os
import re

DEFAULT_REPO = "multiomics-center-Israel/multiomics-annotation-prep"

MODELS_COLUMNS = [
    "Target organism",
    "Model organism (KEGG code)",
    "Surrogate?",
    "Release tag",
    "Model URL",
    "sha256",
    "KEGG snapshot",
    "Build date",
]
_TAG_COL = MODELS_COLUMNS.index("Release tag")

MODELS_HEADER = "| " + " | ".join(MODELS_COLUMNS) + " |"
MODELS_SEP = "|" + "|".join(["---"] * len(MODELS_COLUMNS)) + "|"

MODELS_INTRO = """# Published models

Registry of organism-specific metabolic models built by this repo for
`mummichog -n`. **Before building a new model, check this table** -- if your
organism (or a suitable surrogate) is already listed, reuse its **Model URL** and
**sha256** in your pipeline config instead of rebuilding.

Rows are appended automatically by `scripts/publish_model.py` on each release
(idempotent on the release tag). Published artifacts are immutable -- a rebuild
is a new dated row, never an edit of an existing one.
"""

_DATE = re.compile(r"(\d{4}-\d{2}-\d{2})")
_YYYYMMDD = re.compile(r"^(\d{4})(\d{2})(\d{2})$")


def tag_from_manifest(manifest):
    """Release tag = the model filename stem, e.g. ``cre_kegg_20260711``."""
    return os.path.splitext(manifest["model_file"])[0]


def _snapshot_date(manifest):
    """KEGG snapshot date (YYYY-MM-DD) from source_version, else the newest of
    the per-db snapshot dates, else ''."""
    m = _DATE.search(manifest.get("source_version") or "")
    if m:
        return m.group(1)
    dates = [d for d in (manifest.get("build_details", {})
                         .get("kegg_db_dates", {}) or {}).values() if d]
    return max(dates) if dates else ""


def _build_date(manifest, tag):
    """Build date (YYYY-MM-DD): from the tag's trailing YYYYMMDD (ties it to the
    release), else the manifest build timestamp's date, else ''."""
    m = _YYYYMMDD.match(tag.rsplit("_", 1)[-1])
    if m:
        return f"{m.group(1)}-{m.group(2)}-{m.group(3)}"
    return (manifest.get("build_timestamp_utc") or "")[:10]


def model_url(manifest, repo=DEFAULT_REPO):
    tag = tag_from_manifest(manifest)
    return (f"https://github.com/{repo}/releases/download/"
            f"{tag}/{manifest['model_file']}")


def build_row(manifest, repo=DEFAULT_REPO):
    """Return the ordered cell values for a MODELS.md row from a manifest."""
    model_org = manifest.get("model_organism", "") or ""
    code = manifest.get("model_kegg_code", "") or ""
    return [
        manifest.get("target_organism") or model_org,        # falls back to model
        f"{model_org} ({code})" if code else model_org,
        "yes" if manifest.get("model_is_surrogate") else "no",
        tag_from_manifest(manifest),
        model_url(manifest, repo),
        manifest.get("sha256", ""),
        _snapshot_date(manifest),
        _build_date(manifest, tag_from_manifest(manifest)),
    ]


def format_row(cells):
    return "| " + " | ".join(str(c) for c in cells) + " |"


def _iter_table_rows(text):
    """Yield (line_index, cells) for each data row of the models table."""
    for i, line in enumerate(text.splitlines()):
        s = line.strip()
        if not s.startswith("|"):
            continue
        cells = [c.strip() for c in s.strip("|").split("|")]
        if cells and cells[0] in ("Target organism", ""):
            continue                      # header row
        if set("".join(cells)) <= set("-: "):
            continue                      # separator row
        yield i, cells


def existing_tags(text):
    """Release tags already present in the table."""
    tags = set()
    for _i, cells in _iter_table_rows(text):
        if len(cells) > _TAG_COL and cells[_TAG_COL]:
            tags.add(cells[_TAG_COL])
    return tags


def ensure_models_file(models_path):
    """Create MODELS.md (intro + header + separator) if it does not exist."""
    if os.path.exists(models_path):
        return
    with open(models_path, "w") as f:
        f.write(MODELS_INTRO + "\n" + MODELS_HEADER + "\n" + MODELS_SEP + "\n")


def append_model_row(models_path, manifest, repo=DEFAULT_REPO):
    """Append a row for *manifest* to MODELS.md. Idempotent on the release tag.

    Returns True if a row was added, False if the tag was already present.
    """
    ensure_models_file(models_path)
    tag = tag_from_manifest(manifest)
    with open(models_path) as f:
        text = f.read()
    if tag in existing_tags(text):
        return False

    row = format_row(build_row(manifest, repo))
    lines = text.rstrip("\n").split("\n")
    table_rows = [i for i, line in enumerate(lines) if line.lstrip().startswith("|")]
    if table_rows:
        lines.insert(table_rows[-1] + 1, row)   # after the last table line
    else:
        lines += ["", MODELS_HEADER, MODELS_SEP, row]
    with open(models_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return True
