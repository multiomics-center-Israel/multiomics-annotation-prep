"""Organism-specific KEGG downloads (model organisms)."""

import os
import re

from .utils import cached_download, log_msg, write_pathway2gene, write_pathway2name


# KEGG's ``/conv/`` operation maps KEGG gene ids to a handful of *outside*
# databases. For genes those are exactly these three. There is NO ``ensembl``
# conv database (``conv/ensembl/<org>`` -> HTTP 400): to key KEGG pathways on
# Ensembl gene ids we bridge through ``ncbi-geneid`` plus an Ensembl-provided
# NCBI<->Ensembl cross-reference (see ``kegg_to_ensembl_map`` and
# ``prepare_ensembl``).
_KEGG_GENE_CONV_DBS = {"ncbi-geneid", "ncbi-proteinid", "uniprot"}


def kegg_conv_map(kegg_org, conv_db, cache_dir, refresh=False):
    """Return ``{kegg_gene_id: external_id}`` from ``/conv/<conv_db>/<org>``.

    ``conv_db`` must be a KEGG outside database valid for genes
    (``ncbi-geneid``, ``ncbi-proteinid``, ``uniprot``). The KEGG gene id keeps
    its organism prefix (e.g. ``mmu:11298``); the external id has its own
    ``<db>:`` prefix stripped (e.g. ``ncbi-geneid:839580`` -> ``839580``).
    """
    conv_url = f"https://rest.kegg.jp/conv/{conv_db}/{kegg_org}"
    conv_file = cached_download(conv_url,
                                f"kegg_{kegg_org}_to_{conv_db}.txt",
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
            id_map[parts[0]] = re.sub(r"^[^:]+:", "", parts[1])
    return id_map


def kegg_to_ensembl_map(kegg_org, ncbi_to_ensembl, cache_dir, refresh=False):
    """Compose ``{kegg_gene_id: ensembl_gene_id}`` via the ncbi-geneid bridge.

    KEGG gene -> NCBI gene id (``/conv/ncbi-geneid/<org>``) -> Ensembl gene id
    (``ncbi_to_ensembl``, built from Ensembl BioMart's NCBI cross-reference).
    Genes whose NCBI id has no Ensembl cross-reference are simply absent from
    the result; callers may still resolve them by a direct id match (see
    ``prepare_kegg_by_org``'s ``ext_id_universe`` fallback).
    """
    kegg_to_ncbi = kegg_conv_map(kegg_org, "ncbi-geneid", cache_dir, refresh)
    mapping = {}
    for kegg_id, ncbi in kegg_to_ncbi.items():
        ens = ncbi_to_ensembl.get(str(ncbi))
        if ens:
            mapping[kegg_id] = ens
    return mapping


def prepare_kegg_by_org(kegg_org, out_dir, cache_dir, refresh=False,
                        id_source="kegg", ext_id_map=None, ext_id_universe=None):
    """Write ``KEGG_pathway2gene.tab`` / ``KEGG_pathway2name.tab`` for an organism.

    Gene ids are re-keyed onto the id space of the companion GO/annotation
    tables so ``clusterProfiler`` can match them:

    * ``id_source="kegg"`` (default) -> bare KEGG gene id (org prefix stripped).
    * ``id_source`` in {ncbi-geneid, ncbi-proteinid, uniprot} -> KEGG ``/conv/``.
    * ``ext_id_map`` / ``ext_id_universe`` supplied -> caller-built mapping,
      used for id spaces KEGG's ``/conv/`` can't produce (Ensembl gene ids).
      For each KEGG gene, ``ext_id_map`` wins; otherwise, if the prefix-stripped
      KEGG id is in ``ext_id_universe`` it maps to itself (covers organisms like
      Arabidopsis whose KEGG gene ids ARE the Ensembl locus codes). Genes
      resolved by neither are dropped.
    """
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

    if ext_id_map is not None or ext_id_universe:
        # Caller-supplied mapping (e.g. the Ensembl bridge). Prefer an explicit
        # cross-reference, fall back to a direct KEGG-id == final-id match.
        ext_id_map = ext_id_map or {}
        universe = set(ext_id_universe or ())
        mapped_genes = []
        mapped_paths = []
        for g, p in zip(genes, pathways):
            if g in ext_id_map:
                mapped_genes.append(ext_id_map[g])
                mapped_paths.append(p)
            else:
                stripped = re.sub(r"^[^:]+:", "", g)
                if stripped in universe:
                    mapped_genes.append(stripped)
                    mapped_paths.append(p)
        genes = mapped_genes
        pathways = mapped_paths
    elif id_source != "kegg":
        if id_source not in _KEGG_GENE_CONV_DBS:
            raise ValueError(
                f"KEGG /conv has no gene database {id_source!r}; supported: "
                f"{', '.join(sorted(_KEGG_GENE_CONV_DBS))}. For Ensembl gene "
                f"ids call via prepare_ensembl, which bridges through "
                f"ncbi-geneid.")
        id_map = kegg_conv_map(kegg_org, id_source, cache_dir, refresh)
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
