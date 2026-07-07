"""UniProt REST path (proteomics) — secondary module."""

import os
import re
from collections import defaultdict

from .download_go import go_term_table, build_ancestor_cache, get_ancestors
from .download_kegg_org import prepare_kegg_by_org
from .utils import cached_download, log_msg, write_go2gene, write_go2name, write_tab


def prepare_uniprot(taxon_id, out_dir, cache_dir, reviewed=True,
                    kegg_org=None, expand=True, refresh=False):
    query = f"organism_id:{taxon_id}"
    if reviewed:
        query += "+AND+reviewed:true"
    url = (f"https://rest.uniprot.org/uniprotkb/stream?"
           f"query={query}&format=tsv"
           f"&fields=accession,gene_names,protein_name,go_id")

    tsv_file = cached_download(url, f"uniprot_{taxon_id}.tsv",
                               cache_dir, refresh)

    rows = []
    go_re = re.compile(r"GO:\d{7}")
    with open(tsv_file) as f:
        header = f.readline()
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            accession = parts[0]
            gene_names = parts[1]
            protein_name = parts[2]
            go_ids_str = parts[3]
            go_ids = go_re.findall(go_ids_str)
            rows.append({
                "accession": accession,
                "gene_names": gene_names,
                "protein_name": protein_name,
                "go_ids": go_ids,
            })

    annot_rows = [(r["accession"], r["gene_names"], r["protein_name"])
                  for r in rows]
    write_tab(annot_rows, ["Gene", "gene_name", "protein_name"],
              os.path.join(out_dir, "Annotation.tab"))

    pairs = []
    for r in rows:
        for go_id in r["go_ids"]:
            pairs.append((r["accession"], go_id))

    terms, parents = go_term_table(cache_dir, refresh)
    if expand:
        pairs = _expand_go(pairs, parents)

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
        prepare_kegg_by_org(kegg_org, out_dir, cache_dir, refresh, "uniprot")

    log_msg("UniProt: ", len(rows), " proteins processed")


def _expand_go(pairs, parents):
    ancestor_cache = build_ancestor_cache(parents)
    expanded = set()
    for gene, go_id in pairs:
        expanded.add((gene, go_id))
        for anc in get_ancestors(go_id, parents, ancestor_cache):
            expanded.add((gene, anc))
    return list(expanded)
