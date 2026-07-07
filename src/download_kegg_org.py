"""Organism-specific KEGG downloads (model organisms)."""

import os
import re

from .utils import cached_download, log_msg, write_pathway2gene, write_pathway2name


def prepare_kegg_by_org(kegg_org, out_dir, cache_dir, refresh=False,
                        id_source="kegg"):
    gene2path_url = f"https://rest.kegg.jp/link/pathway/{kegg_org}"
    gene2path_file = cached_download(gene2path_url,
                                     f"kegg_{kegg_org}_gene2path.txt",
                                     cache_dir, refresh)

    pw_names_url = f"https://rest.kegg.jp/list/pathway/{kegg_org}"
    pw_names_file = cached_download(pw_names_url,
                                    f"kegg_{kegg_org}_pathway_names.txt",
                                    cache_dir, refresh)

    genes = []
    pathways = []
    with open(gene2path_file) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t", 1)
            if len(parts) < 2:
                continue
            genes.append(parts[0])
            pathways.append(parts[1].replace("path:", ""))

    if id_source != "kegg":
        conv_url = f"https://rest.kegg.jp/conv/{id_source}/{kegg_org}"
        conv_file = cached_download(conv_url,
                                    f"kegg_{kegg_org}_to_{id_source}.txt",
                                    cache_dir, refresh)
        id_map = {}
        with open(conv_file) as f:
            for line in f:
                line = line.rstrip("\n")
                if not line:
                    continue
                parts = line.split("\t", 1)
                if len(parts) < 2:
                    continue
                kegg_id = parts[0]
                ext_id = re.sub(r"^[^:]+:", "", parts[1])
                id_map[kegg_id] = ext_id

        mapped_genes = []
        mapped_paths = []
        for g, p in zip(genes, pathways):
            if g in id_map:
                mapped_genes.append(id_map[g])
                mapped_paths.append(p)
        genes = mapped_genes
        pathways = mapped_paths
    else:
        genes = [re.sub(r"^[^:]+:", "", g) for g in genes]

    pw_names = {}
    with open(pw_names_file) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t", 1)
            if len(parts) < 2:
                continue
            pw_names[parts[0].replace("path:", "")] = parts[1]

    path2gene_pairs = list(zip(pathways, genes))
    write_pathway2gene(path2gene_pairs,
                       os.path.join(out_dir, "KEGG_pathway2gene.tab"))

    pw_list = [(pid, name) for pid, name in pw_names.items()]
    write_pathway2name(pw_list,
                       os.path.join(out_dir, "KEGG_pathway2name.tab"))

    log_msg("KEGG (org=", kegg_org, "): ",
            len(set(pathways)), " pathways, ", len(set(genes)), " genes")
