"""GO enrichment file preparation: per-gene GO table -> expanded GO files."""

import os
import re
from collections import defaultdict

from .download_go import go_term_table, build_ancestor_cache, get_ancestors
from .utils import log_msg, write_gmt, write_go2gene, write_go2name


def prepare_go(go_table_file, out_dir, cache_dir, refresh=False, expand=True):
    gene2go = _parse_go_table(go_table_file)
    log_msg("GO input: ", len(gene2go), " genes with GO annotations")

    terms, parents = go_term_table(cache_dir, refresh)

    pairs = []
    for gene, go_ids in gene2go.items():
        for go_id in go_ids:
            pairs.append((gene, go_id))

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
            log_msg("GO ", ns, ": no annotations found")
            continue
        write_go2gene(data["go2gene"],
                      os.path.join(out_dir, f"GO2gene_{ns}.tab"))
        name_pairs = [(go_id, name) for go_id, name in sorted(data["go2name"].items())]
        write_go2name(name_pairs,
                      os.path.join(out_dir, f"GO2name_{ns}.tab"))
        log_msg("GO ", ns, ": ", len(set(g for _, g in data["go2gene"])),
                " genes, ", len(data["go2name"]), " terms")

        term2genes = defaultdict(list)
        for go_id, gene in data["go2gene"]:
            term2genes[go_id].append(gene)
        write_gmt(dict(term2genes), data["go2name"],
                  os.path.join(out_dir, f"GO_{ns}.gmt"))


def _parse_go_table(path):
    """Parse per-gene GO table. Handles Trinotate suffixes and mixed separators."""
    gene2go = defaultdict(set)
    go_re = re.compile(r"GO:\d{7}")
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t", 1)
            if len(parts) < 2:
                continue
            gene = parts[0]
            go_ids = go_re.findall(parts[1])
            for go_id in go_ids:
                gene2go[gene].add(go_id)
    return dict(gene2go)


def _expand_go(pairs, parents):
    """Expand GO annotations to include all ancestor terms."""
    ancestor_cache = build_ancestor_cache(parents)
    expanded = set()
    for gene, go_id in pairs:
        expanded.add((gene, go_id))
        for anc in get_ancestors(go_id, parents, ancestor_cache):
            expanded.add((gene, anc))
    return list(expanded)
