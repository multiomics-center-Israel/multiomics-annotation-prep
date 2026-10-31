"""manifest.json for a gene-based annotation_dir release (RNA / proteomics / model).

One machine-readable record per release, written into the output folder before
it is zipped, so it travels with the tables. multiomic-core reads it to state
which release, source and version each table came from.

Versions are read, never inferred:

* ``go_ontology`` -- from the ``go-basic.obo`` the build itself downloaded into
  the cache (exact). Same shape as ``go_ontology`` in ``GO_provenance.json``.
* ``gene_annotation`` -- where the gene -> GO links (and ``Annotation.tab``)
  came from. This is a different resource from the ontology and has its own
  version, if any: for Ensembl, the release name the host's BioMart registry
  reports at build time; otherwise ``null``. GAF details stay in
  ``GO_provenance.json`` (named in ``details_file``, not copied).
* ``kegg`` -- the per-database dates rest.kegg.jp/info/kegg reports during the
  build job. They are not a single KEGG release. Organism-code details stay in
  ``KEGG_provenance.json``.

The KEGG dates and the Ensembl release name come from live queries, so they
describe the build only when the build fetched its own data in the same run
(the release workflows: a fresh cache on every runner). They are therefore
queried only when the caller asks (``query_versions=True``), always re-fetched
into manifest-only cache files (never the builders' caches, never a stale or
other-host copy), and best-effort: if they fail, the value is ``null`` and the
build still succeeds. Without ``query_versions`` they are ``null``. Detailed input
provenance is not repeated here: inputs are identified by ``builder_git_sha``
(they are committed) and by the companion ``*_provenance.json`` files, which
are listed in ``files`` with their sha256.
"""

import json
import os
import re
from datetime import datetime, timezone

from .download_go import GO_OBO_URL
from .download_kegg import parse_kegg_info_dates
from .prepare_ensembl import _virtual_schema_for
from .prepare_gaf import obo_data_version
from .utils import cached_download, git_head_sha, log_msg, sha256_file

MANIFEST_NAME = "manifest.json"

GENE_ANNOTATION_SOURCES = ("trinotate", "gaf", "ensembl", "uniprot")

# prepare_ensembl's default when no host is given; the registry must be queried
# on the host the tables came from.
_ENSEMBL_DEFAULT_HOST = "https://www.ensembl.org"

# Which manifest blocks each output file is built from (by file name).
_FILE_SOURCES = (
    (re.compile(r"^GO2gene_.*\.tab$|^GO_.*\.gmt$"), ["go_ontology", "gene_annotation"]),
    (re.compile(r"^GO2name_.*\.tab$"), ["go_ontology"]),
    (re.compile(r"^Annotation\.tab$"), ["gene_annotation"]),
    (re.compile(r"^KEGG_.*\.(tab|gmt|txt)$"), ["kegg"]),
)


def file_sources(name):
    """Return the manifest blocks a built file derives from ([] if none)."""
    for pattern, sources in _FILE_SOURCES:
        if pattern.match(name):
            return list(sources)
    return []


def _best_effort(what, fn):
    """Run an optional version query; log and return None on any failure."""
    try:
        return fn()
    except Exception as exc:  # noqa: BLE001 - optional metadata must not fail a build
        log_msg("manifest: ", what, " unavailable (", exc, ") -- recorded as null")
        return None


def _go_ontology(cache_dir):
    """go_ontology block from the cached go-basic.obo, or None if not used."""
    obo = os.path.join(cache_dir, "go-basic.obo")
    if not os.path.exists(obo):
        return None
    version = obo_data_version(obo)
    return {"url": GO_OBO_URL,
            "data_version": None if version == "unknown" else version,
            "sha256": sha256_file(obo)}


def kegg_db_dates(cache_dir):
    """Per-database dates from info/kegg, fetched now, or None if unavailable.

    Re-fetched into a manifest-only file: the builders' ``kegg_info_kegg.txt``
    may be from an earlier run and must not be overwritten here.
    """
    def query():
        path = cached_download("https://rest.kegg.jp/info/kegg",
                               "manifest_kegg_info.txt", cache_dir, refresh=True)
        return parse_kegg_info_dates(path) or None
    return _best_effort("KEGG info/kegg", query)


def parse_biomart_registry(text, mart):
    """Return the ``displayName`` of *mart* in a BioMart registry XML, or None.

    The main Ensembl mart is registered as ``ENSEMBL_MART_ENSEMBL``; Ensembl
    Genomes marts under their own name (e.g. ``plants_mart``).
    """
    schema = _virtual_schema_for(mart)
    wanted = "ENSEMBL_MART_ENSEMBL" if schema == "default" else schema
    for tag in re.findall(r"<MartURLLocation\b[^>]*>", text):
        attrs = dict(re.findall(r'(\w+)="([^"]*)"', tag))
        if attrs.get("name") == wanted:
            return attrs.get("displayName") or None
    return None


def ensembl_registry_name(host, mart, cache_dir):
    """Release name *host*'s BioMart registry reports now (best-effort, or None).

    Always re-fetched, so the name belongs to this host, not a cached copy.
    """
    def query():
        url = host.rstrip("/") + "/biomart/martservice?type=registry"
        path = cached_download(url, "manifest_biomart_registry.xml", cache_dir,
                               refresh=True)
        with open(path, encoding="utf-8", errors="replace") as f:
            return parse_biomart_registry(f.read(), mart)
    return _best_effort("BioMart registry", query)


def build_annotation_manifest(out_dir, cache_dir, release_tag=None, run_url=None,
                              inputs=None, gene_annotation_source=None,
                              ensembl_host=None, ensembl_mart="ensembl",
                              kegg_ids_from_gene_annotation=False,
                              query_versions=False):
    """Build the manifest dict for the files in *out_dir* (see module docstring).

    *inputs* is recorded verbatim (the workflow inputs: project / label, ids,
    organism code, dataset, taxon ...). *kegg_ids_from_gene_annotation* is set
    by the caller when the build re-keyed the KEGG gene ids through the gene
    annotation source (prepare_ensembl's Ensembl bridge); the KEGG files then
    derive from both. Version queries run only with *query_versions* and only
    for sources the built files actually use.
    """
    if gene_annotation_source not in (None,) + GENE_ANNOTATION_SOURCES:
        raise ValueError(f"unknown gene_annotation_source: {gene_annotation_source!r}")

    names = sorted(n for n in os.listdir(out_dir)
                   if n != MANIFEST_NAME and os.path.isfile(os.path.join(out_dir, n)))
    files = [{"file": n, "sha256": sha256_file(os.path.join(out_dir, n)),
              "sources": file_sources(n)} for n in names]
    if kegg_ids_from_gene_annotation and gene_annotation_source:
        for f in files:
            if f["sources"] == ["kegg"]:
                f["sources"] = ["kegg", "gene_annotation"]
    used = {s for f in files for s in f["sources"]}

    go_ontology = _go_ontology(cache_dir) if "go_ontology" in used else None

    gene_annotation = None
    if gene_annotation_source and "gene_annotation" in used:
        gene_annotation = {"source": gene_annotation_source, "details_file": None}
        if gene_annotation_source == "ensembl":
            host = ensembl_host or _ENSEMBL_DEFAULT_HOST
            gene_annotation["registry_host"] = host
            gene_annotation["registry_display_name"] = (
                ensembl_registry_name(host, ensembl_mart, cache_dir)
                if query_versions else None)
        if "GO_provenance.json" in names:
            gene_annotation["details_file"] = "GO_provenance.json"

    kegg = None
    if "kegg" in used:
        kegg = {"kegg_db_dates": kegg_db_dates(cache_dir) if query_versions else None,
                "details_file": "KEGG_provenance.json"
                if "KEGG_provenance.json" in names else None}

    return {
        "release_tag": release_tag,
        "run_url": run_url,
        "builder_git_sha": git_head_sha(),
        "build_timestamp_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "inputs": dict(inputs or {}),
        "go_ontology": go_ontology,
        "gene_annotation": gene_annotation,
        "kegg": kegg,
        "files": files,
    }


def write_annotation_manifest(out_dir, cache_dir, **kwargs):
    """Write ``<out_dir>/manifest.json``; return its path."""
    manifest = build_annotation_manifest(out_dir, cache_dir, **kwargs)
    path = os.path.join(out_dir, MANIFEST_NAME)
    with open(path, "w") as f:
        json.dump(manifest, f, indent=2)
        f.write("\n")
    log_msg("wrote manifest: ", path, f"  ({len(manifest['files'])} files)")
    return path
