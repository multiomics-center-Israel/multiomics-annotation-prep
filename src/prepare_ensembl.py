"""Ensembl/BioMart path (model organisms) — secondary module."""

import os
import re
from collections import defaultdict

from .download_go import go_term_table, build_ancestor_cache, get_ancestors
from .download_kegg_org import kegg_to_ensembl_map, prepare_kegg_by_org
from .utils import log_msg, write_go2gene, write_go2name, write_tab


# pybiomart keys ``server.marts`` by the mart's internal name (e.g.
# "ENSEMBL_MART_ENSEMBL"), NOT "ensembl". Accept friendly aliases + a
# case-insensitive match on name/display_name so values like "ensembl" or
# "plants" work regardless of the exact internal key.
_MART_ALIASES = {
    "ensembl": "ENSEMBL_MART_ENSEMBL",
    "plants": "plants_mart",
    "protists": "protists_mart",
    "fungi": "fungi_mart",
    "metazoa": "metazoa_mart",
}


def _resolve_mart(marts, mart):
    """Resolve a mart name tolerantly against pybiomart's ``server.marts``."""
    if mart in marts:
        return marts[mart]
    alias = _MART_ALIASES.get(mart.lower())
    if alias and alias in marts:
        return marts[alias]
    low = mart.lower()
    for m in marts.values():
        names = {getattr(m, "name", "").lower(),
                 (getattr(m, "display_name", "") or "").lower()}
        if low in names:
            return m
    raise ValueError(
        f"BioMart mart {mart!r} not found. Available: {', '.join(marts)}. "
        f"For a non-vertebrate genome set --host to the right division "
        f"(e.g. https://plants.ensembl.org with --mart plants_mart).")


# BioMart's attribute for the NCBI/Entrez gene cross-reference is named
# differently across divisions (main vs Ensembl Genomes). Try the common
# spellings and use the first one the server accepts.
_NCBI_XREF_ATTRS = ("entrezgene_id", "entrezgene", "entrezgene_trans_name")


def _fetch_ncbi_xref(ds):
    """Return ``{ncbi_gene_id: ensembl_gene_id}`` from BioMart, or ``{}``.

    Used to bridge KEGG gene ids (which resolve to NCBI gene ids via KEGG
    ``/conv/ncbi-geneid``) to Ensembl gene ids. Returns an empty dict if no
    NCBI cross-reference attribute is available for the dataset.
    """
    for attr in _NCBI_XREF_ATTRS:
        try:
            df = ds.query(attributes=["ensembl_gene_id", attr])
        except Exception:  # noqa: BLE001 - attribute not offered by this mart
            continue
        df.columns = ["ensembl", "ncbi"]
        mapping = {}
        for ens, ncbi in zip(df["ensembl"], df["ncbi"]):
            ncbi = str(ncbi).strip()
            # pandas may read integer Entrez ids as floats ("839580.0").
            if ncbi.endswith(".0"):
                ncbi = ncbi[:-2]
            if ncbi and ncbi.lower() != "nan":
                mapping[ncbi] = ens
        if mapping:
            log_msg("NCBI<->Ensembl cross-reference via BioMart attribute '",
                    attr, "': ", len(mapping), " ids")
            return mapping
    log_msg("No NCBI/Entrez cross-reference attribute available; KEGG mapping "
            "will rely on direct KEGG-id == Ensembl-id matches only")
    return {}


def _query_go_pairs(ds, dataset):
    """Query ``(ensembl_gene_id, go_id)`` from BioMart, with retry + shape guard.

    The GO namespace (BP/MF/CC) is taken from the OBO downstream, NOT from
    BioMart, so we do not request ``namespace_1003``. Fewer attributes = a more
    robust query: BioMart returns a single-column error frame when any requested
    attribute is unavailable for a dataset/release, which previously crashed the
    column rename with a cryptic pandas length mismatch. Retry once (BioMart
    hiccups are common) and fail with a clear message if the shape is still off.
    """
    last_shape = None
    for _ in range(2):
        df = ds.query(attributes=["ensembl_gene_id", "go_id"])
        if df.shape[1] == 2:
            df.columns = ["gene", "go_id"]
            return df
        last_shape = df.shape
    raise RuntimeError(
        f"BioMart GO query for dataset {dataset!r} returned an unexpected shape "
        f"{last_shape} (expected 2 columns: ensembl_gene_id, go_id). BioMart "
        f"usually does this when it returns an error page -- retry the run.")


def prepare_ensembl(dataset, out_dir, cache_dir, mart="ensembl", host=None,
                    kegg_org=None, id_source="ensembl", expand=True,
                    refresh=False):
    try:
        import pybiomart
    except ImportError:
        raise ImportError(
            "pybiomart is required for the Ensembl module. "
            "Install with: pip install pybiomart"
        )

    # pybiomart's Server talks HTTP on port 80; an ``https://`` host makes it
    # build a malformed URL (it ends up treating "https" as the hostname ->
    # "Failed to resolve 'https'"). Normalize to http:// -- BioMart martservice
    # is served over http, matching the default www.ensembl.org host.
    if host and host.startswith("https://"):
        host = "http://" + host[len("https://"):]
    log_msg("Connecting to BioMart: dataset=", dataset, " mart=", mart,
            " host=", host or "http://www.ensembl.org")
    server = pybiomart.Server(host=host or "http://www.ensembl.org")

    mart_obj = _resolve_mart(server.marts, mart)
    datasets = mart_obj.datasets
    if dataset not in datasets:
        raise ValueError(
            f"BioMart dataset {dataset!r} not found in mart {mart_obj.name!r}. "
            f"For a non-vertebrate genome point --host at the right division "
            f"(e.g. https://plants.ensembl.org with --mart plants_mart). "
            f"Available datasets start with: {', '.join(list(datasets)[:6])} ...")
    ds = datasets[dataset]

    # Ensembl Genomes divisions (plants_mart, protists_mart, fungi_mart, ...)
    # answer queries in a virtual schema named after the mart, NOT "default".
    # pybiomart defaults the query's virtualSchemaName to "default", so BioMart
    # replies "Dataset <x> NOT FOUND". Force the dataset's schema to the mart
    # name for any non-main division. (Best-effort: pybiomart stores it on the
    # private _virtual_schema that Dataset.query reads.)
    if mart_obj.name != "ENSEMBL_MART_ENSEMBL":
        for _attr in ("_virtual_schema", "virtual_schema"):
            if hasattr(ds, _attr):
                try:
                    setattr(ds, _attr, mart_obj.name)
                except Exception:  # noqa: BLE001 - never fatal
                    pass

    log_msg("Fetching gene annotations...")
    annot = ds.query(attributes=[
        "ensembl_gene_id", "external_gene_name", "gene_biotype", "description"
    ])
    annot.columns = ["Gene", "gene_name", "gene_biotype", "description"]
    annot_rows = [tuple(r) for r in annot.values]
    ensembl_ids = {str(r[0]) for r in annot_rows}
    write_tab(annot_rows, list(annot.columns),
              os.path.join(out_dir, "Annotation.tab"))

    log_msg("Fetching GO annotations...")
    go_df = _query_go_pairs(ds, dataset)
    go_df = go_df.dropna(subset=["go_id"])
    go_df = go_df[go_df["go_id"].astype(str).str.len() > 0]

    pairs = list(zip(go_df["gene"], go_df["go_id"]))
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
                ncbi_to_ensembl = _fetch_ncbi_xref(ds)
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
