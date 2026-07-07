"""Ensembl/BioMart path (model organisms) — secondary module."""

import os
import re
from collections import defaultdict

from .download_go import go_term_table, build_ancestor_cache, get_ancestors
from .download_kegg_org import prepare_kegg_by_org
from .utils import log_msg, write_go2gene, write_go2name, write_tab


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

    log_msg("Connecting to BioMart: dataset=", dataset)
    if host:
        server = pybiomart.Server(host=host)
    else:
        server = pybiomart.Server(host="http://www.ensembl.org")
    ds = server.marts[mart].datasets[dataset]

    log_msg("Fetching gene annotations...")
    annot = ds.query(attributes=[
        "ensembl_gene_id", "external_gene_name", "gene_biotype", "description"
    ])
    annot.columns = ["Gene", "gene_name", "gene_biotype", "description"]
    annot_rows = [tuple(r) for r in annot.values]
    write_tab(annot_rows, list(annot.columns),
              os.path.join(out_dir, "Annotation.tab"))

    log_msg("Fetching GO annotations...")
    go_df = ds.query(attributes=[
        "ensembl_gene_id", "go_id", "namespace_1003"
    ])
    go_df.columns = ["gene", "go_id", "namespace"]
    go_df = go_df[go_df["go_id"].str.len() > 0]

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
        prepare_kegg_by_org(kegg_org, out_dir, cache_dir, refresh, id_source)


def _expand_go(pairs, parents, terms):
    ancestor_cache = build_ancestor_cache(parents)
    expanded = set()
    for gene, go_id in pairs:
        expanded.add((gene, go_id))
        for anc in get_ancestors(go_id, parents, ancestor_cache):
            expanded.add((gene, anc))
    return list(expanded)
