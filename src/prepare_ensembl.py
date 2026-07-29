"""Ensembl/BioMart path (model organisms) — secondary module.

Queries BioMart directly over its REST martservice (build the ``<Query>`` XML,
GET the TSV result) instead of going through ``pybiomart``. pybiomart depends on
the ``?type=registry`` endpoint, which Ensembl intermittently serves as
unparseable XML/HTML (crashing with ElementTree "mismatched tag"), and it cannot
talk HTTPS at all (its URL builder mangles an ``https://`` host). The direct REST
path needs no registry/configuration fetch, works over HTTPS, and reuses the
repo's ``cached_download`` + ``requests`` — no pybiomart/pandas dependency.
"""

import os
import time
import urllib.parse
from collections import defaultdict

from .download_go import go_term_table, build_ancestor_cache, get_ancestors
from .download_kegg_org import kegg_to_ensembl_map, prepare_kegg_by_org
from .utils import cached_download, log_msg, write_go2gene, write_go2name, write_tab


# Friendly ``--mart`` values -> BioMart virtual-schema names. Main Ensembl
# (vertebrates + fungi/yeast) answers in the "default" schema; each Ensembl
# Genomes division answers in a schema named after its mart.
_MART_ALIASES = {
    "plants": "plants_mart",
    "protists": "protists_mart",
    "fungi": "fungi_mart",
    "metazoa": "metazoa_mart",
}

# BioMart's attribute for the NCBI/Entrez gene cross-reference is named
# differently across divisions (main vs Ensembl Genomes). Try the common
# spellings and use the first one the server accepts.
_NCBI_XREF_ATTRS = ("entrezgene_id", "entrezgene", "entrezgene_trans_name")


def _virtual_schema_for(mart):
    """Map a ``--mart`` value to BioMart's ``virtualSchemaName``."""
    m = (mart or "ensembl").strip().lower()
    if m in ("", "ensembl", "ensembl_mart_ensembl", "default"):
        return "default"
    return _MART_ALIASES.get(m, mart)


def _biomart_with_retry(fn, what, attempts=3):
    """Call *fn* with retries for transient BioMart network failures."""
    import requests  # lazy, matches the repo's requests-on-demand style
    last = None
    for i in range(attempts):
        try:
            return fn()
        except requests.exceptions.RequestException as exc:  # noqa: BLE001
            last = exc
            log_msg("BioMart ", what, " request failed (attempt ", i + 1, "/",
                    attempts, "): ", exc)
            if i < attempts - 1:
                time.sleep(3 * (i + 1))
    raise RuntimeError(
        f"BioMart {what} failed after {attempts} attempts ({last!r}). Likely a "
        f"transient Ensembl outage -- re-run the workflow in a few minutes.")


def _biomart_query(host, dataset, attributes, virtual_schema,
                   cache_dir, dest, refresh=False):
    """Run one BioMart query via direct REST; return rows (incl. header).

    Builds the martservice ``<Query>`` XML for *attributes* on *dataset* and
    fetches the TSV through ``cached_download``. Raises a clear error (and drops
    the cache file) if BioMart answers with an HTML/error page or an in-band
    ``Query ERROR`` instead of a TSV table.
    """
    attr_xml = "".join(f'<Attribute name="{a}" />' for a in attributes)
    query_xml = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<!DOCTYPE Query>'
        f'<Query virtualSchemaName="{virtual_schema}" formatter="TSV" '
        'header="1" uniqueRows="1" datasetConfigVersion="0.6">'
        f'<Dataset name="{dataset}" interface="default">{attr_xml}</Dataset>'
        '</Query>'
    )
    url = (host.rstrip("/") + "/biomart/martservice?"
           + urllib.parse.urlencode({"query": query_xml}))
    path = _biomart_with_retry(
        lambda: cached_download(url, dest, cache_dir, refresh),
        f"{dest} query")

    with open(path, encoding="utf-8", errors="replace") as f:
        text = f.read()
    if text.lstrip().startswith("<") or "Query ERROR" in text[:2000]:
        snippet = " ".join(text[:300].split())
        try:
            os.remove(path)  # never keep a bad response in the cache
        except OSError:
            pass
        raise RuntimeError(
            f"BioMart returned a non-TSV response for '{dest}' "
            f"(host={host}, dataset={dataset}, schema={virtual_schema}). "
            f"First bytes: {snippet!r}")

    ncol = len(attributes)
    rows = []
    for line in text.splitlines():
        if line == "":
            continue
        parts = line.split("\t")
        if len(parts) < ncol:            # pad missing trailing empty fields
            parts += [""] * (ncol - len(parts))
        rows.append(parts[:ncol])
    return rows


def _fetch_ncbi_xref(host, dataset, virtual_schema, cache_dir, refresh=False):
    """Return ``{ncbi_gene_id: ensembl_gene_id}`` from BioMart, or ``{}``.

    Bridges KEGG gene ids (which resolve to NCBI gene ids via KEGG
    ``/conv/ncbi-geneid``) to Ensembl gene ids. Empty if no NCBI cross-reference
    attribute is available for the dataset.
    """
    for attr in _NCBI_XREF_ATTRS:
        try:
            rows = _biomart_query(host, dataset, ["ensembl_gene_id", attr],
                                  virtual_schema, cache_dir,
                                  f"biomart_{dataset}_{attr}.tsv", refresh)
        except RuntimeError:  # attribute not offered by this mart / release
            continue
        mapping = {}
        for r in rows[1:]:
            ens, ncbi = r[0].strip(), r[1].strip()
            if ens and ncbi:
                mapping[ncbi] = ens
        if mapping:
            log_msg("NCBI<->Ensembl cross-reference via BioMart attribute '",
                    attr, "': ", len(mapping), " ids")
            return mapping
    log_msg("No NCBI/Entrez cross-reference attribute available; KEGG mapping "
            "will rely on direct KEGG-id == Ensembl-id matches only")
    return {}


def prepare_ensembl(dataset, out_dir, cache_dir, mart="ensembl", host=None,
                    kegg_org=None, id_source="ensembl", expand=True,
                    refresh=False):
    vs = _virtual_schema_for(mart)
    server_host = host or "https://www.ensembl.org"
    log_msg("Connecting to BioMart (REST): dataset=", dataset, " mart=", mart,
            " host=", server_host, " virtual_schema=", vs)

    log_msg("Fetching gene annotations...")
    annot_rows = _biomart_query(
        server_host, dataset,
        ["ensembl_gene_id", "external_gene_name", "gene_biotype", "description"],
        vs, cache_dir, f"biomart_{dataset}_annot.tsv", refresh)[1:]
    ensembl_ids = {r[0] for r in annot_rows if r[0]}
    write_tab(annot_rows, ["Gene", "gene_name", "gene_biotype", "description"],
              os.path.join(out_dir, "Annotation.tab"))

    log_msg("Fetching GO annotations...")
    # Namespace (BP/MF/CC) comes from the OBO downstream, not BioMart, so we ask
    # only for (gene, go_id) -- fewer attributes = a more robust query.
    go_rows = _biomart_query(server_host, dataset, ["ensembl_gene_id", "go_id"],
                             vs, cache_dir, f"biomart_{dataset}_go.tsv", refresh)
    pairs = [(r[0], r[1]) for r in go_rows[1:] if r[0] and r[1]]

    if expand:
        terms, parents = go_term_table(cache_dir, refresh)
        pairs = _expand_go(pairs, parents, terms)
    else:
        terms, _ = go_term_table(cache_dir, refresh)

    by_ns = defaultdict(lambda: {"go2gene": [], "go2name": {}})
    for gene, go_id in pairs:
        info = terms.get(go_id)
        if not info:
            continue
        ns = info["namespace"]
        by_ns[ns]["go2gene"].append((go_id, gene))
        by_ns[ns]["go2name"][go_id] = info["name"]

    for ns in ["BP", "MF", "CC"]:
        data = by_ns.get(ns)
        if not data:
            continue
        write_go2gene(data["go2gene"],
                      os.path.join(out_dir, f"GO2gene_{ns}.tab"))
        name_pairs = [(gid, name) for gid, name in sorted(data["go2name"].items())]
        write_go2name(name_pairs,
                      os.path.join(out_dir, f"GO2name_{ns}.tab"))

    if kegg_org:
        # The KEGG companion is optional and additive; a failure here must NOT
        # discard the GO tables already written above.
        try:
            if id_source == "ensembl":
                # KEGG's /conv/ has no 'ensembl' database. Bridge KEGG gene ids
                # to Ensembl gene ids through ncbi-geneid, using the NCBI
                # cross-reference from this same BioMart dataset. Organisms
                # whose KEGG gene ids ARE Ensembl locus codes (e.g. Arabidopsis
                # 'ath') still resolve via the direct ext_id_universe fallback.
                ncbi_to_ensembl = _fetch_ncbi_xref(server_host, dataset, vs,
                                                   cache_dir, refresh)
                ext_id_map = (kegg_to_ensembl_map(kegg_org, ncbi_to_ensembl,
                                                  cache_dir, refresh)
                              if ncbi_to_ensembl else None)
                prepare_kegg_by_org(kegg_org, out_dir, cache_dir, refresh,
                                    id_source, ext_id_map=ext_id_map,
                                    ext_id_universe=ensembl_ids)
            else:
                prepare_kegg_by_org(kegg_org, out_dir, cache_dir, refresh,
                                    id_source)
        except Exception as exc:  # noqa: BLE001 - KEGG tables are optional
            log_msg("KEGG mapping for ", kegg_org, " (id_source=", id_source,
                    ") failed; skipping KEGG tables, GO tables are unaffected: ",
                    exc)


def _expand_go(pairs, parents, terms):
    ancestor_cache = build_ancestor_cache(parents)
    expanded = set()
    for gene, go_id in pairs:
        expanded.add((gene, go_id))
        for anc in get_ancestors(go_id, parents, ancestor_cache):
            expanded.add((gene, anc))
    return list(expanded)
