"""Non-model organism KEGG: KAAS query.ko.txt -> enrichment + annotation files."""

import os
import re
from collections import defaultdict

from .download_kegg import download_kegg_rest, parse_ko_to_name, parse_ko_to_path, parse_pathway_names
from .utils import log_msg, write_gmt, write_pathway2gene, write_pathway2name, write_tab


def prepare_kegg_nonmodel(kaas_file, out_dir, cache_dir, refresh=False,
                          strip_isoform=True, gene_list=None):
    gene2ko = defaultdict(set)
    with open(kaas_file) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t", 1)
            contig = parts[0]
            ko = parts[1].strip() if len(parts) > 1 else ""
            if not ko:
                continue
            if strip_isoform:
                contig = re.sub(r"_i\d+$", "", contig)
            ko = ko.replace("ko:", "")
            gene2ko[contig].add(ko)

    log_msg("KAAS input: ", len(gene2ko), " genes with KO assignments")

    ko_name_file, ko_path_file, pw_name_file = download_kegg_rest(cache_dir, refresh)
    ko2info = parse_ko_to_name(ko_name_file)
    ko2path = parse_ko_to_path(ko_path_file)
    path_names = parse_pathway_names(pw_name_file)

    pw_name_pairs = [(pid, name) for pid, name in sorted(path_names.items())]
    write_pathway2name(pw_name_pairs,
                       os.path.join(out_dir, "KEGG_pathway2name.tab"))

    path2gene_pairs = []
    for gene, kos in gene2ko.items():
        for ko in kos:
            for pth in ko2path.get(ko, []):
                path2gene_pairs.append((pth, gene))

    write_pathway2gene(path2gene_pairs,
                       os.path.join(out_dir, "KEGG_pathway2gene.tab"))
    log_msg("KEGG pathway2gene: ", len(set(p for p, _ in path2gene_pairs)),
            " pathways, ", len(set(g for _, g in path2gene_pairs)), " genes")

    term2genes = defaultdict(list)
    for pth, gene in path2gene_pairs:
        term2genes[pth].append(gene)
    write_gmt(dict(term2genes), path_names,
              os.path.join(out_dir, "KEGG_pathway.gmt"))

    annot_rows = _build_kegg_annot(gene2ko, ko2info, ko2path, path_names,
                                   gene_list)
    annot_cols = ["Gene", "KEGG_ID", "KEGG_names", "KEGG_description",
                  "EC_number", "Pathway_IDs", "Pathway_names"]
    write_tab(annot_rows, annot_cols,
              os.path.join(out_dir, "KEGG_annot_genes.txt"))


def _build_kegg_annot(gene2ko, ko2info, ko2path, path_names, gene_list=None):
    genes = gene_list if gene_list else sorted(gene2ko.keys())
    rows = []
    for gene in genes:
        kos = sorted(gene2ko.get(gene, set()))
        if not kos:
            rows.append((gene, "", "", "", "", "", ""))
            continue
        all_names = []
        all_titles = []
        all_ec = []
        all_pids = set()
        for ko in kos:
            info = ko2info.get(ko, {})
            all_names.append(info.get("names", ""))
            all_titles.append(info.get("title", ""))
            ec = info.get("ec", "")
            if ec:
                all_ec.append(ec)
            for p in ko2path.get(ko, []):
                all_pids.add(p)
        pid_list = sorted(all_pids)
        pnames = [path_names.get(p, "") for p in pid_list]
        rows.append((
            gene,
            " | ".join(kos),
            " | ".join(all_names),
            " | ".join(all_titles),
            " | ".join(all_ec),
            " | ".join(pid_list),
            " | ".join(pnames),
        ))
    return rows
