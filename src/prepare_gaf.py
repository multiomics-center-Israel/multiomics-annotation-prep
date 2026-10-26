"""GO from a GAF file (GO Annotation File, e.g. downloaded from GiardiaDB/VEuPathDB).

Organism databases such as the VEuPathDB sites publish per-organism GO
annotation as GAF, but downloads need a logged-in browser session, so the file
is fetched by hand and committed as a project input. This module parses it,
optionally re-keys its gene ids onto the project's id space (e.g. the NCBI
protein accessions a proteomics FASTA carries), and writes the same GO tables
and GMTs as the other GO paths via ``prepare_go.write_go_sets``.
"""

import gzip
import os
from collections import defaultdict

from .prepare_go import write_go_sets
from .utils import log_msg, write_tab


def _open_text(path):
    """Open a plain or gzip-compressed text file (detected by magic bytes)."""
    with open(path, "rb") as f:
        is_gz = f.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if is_gz else open(path)


def parse_gaf(path):
    """Parse a GAF 2.x file -> (gene2go, gene_info).

    gene2go:   gene id -> set of GO ids (rows qualified ``NOT`` are dropped --
               they assert the gene does NOT have the term).
    gene_info: gene id -> (symbol, name) from GAF columns 3 and 10.
    """
    gene2go = defaultdict(set)
    gene_info = {}
    with _open_text(path) as f:
        for line in f:
            if line.startswith("!"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5 or not cols[4].startswith("GO:"):
                continue
            if "NOT" in cols[3].split("|"):
                continue
            gene = cols[1]
            gene2go[gene].add(cols[4])
            if gene not in gene_info:
                name = cols[9] if len(cols) > 9 else ""
                gene_info[gene] = (cols[2], name)
    return dict(gene2go), gene_info


def parse_id_map(path):
    """Parse an id map -> dict: source id -> set of target ids.

    Two formats are accepted:

    * an NCBI ``*_feature_table.txt(.gz)`` -- CDS rows map ``locus_tag`` to
      ``product_accession`` (protein accession, e.g. the ids in an NCBI FASTA);
    * a two-column TSV ``source_id<TAB>target_id`` (``#`` lines ignored).
    """
    mapping = defaultdict(set)
    with _open_text(path) as f:
        first = f.readline()
        if first.startswith("# feature"):
            header = first.lstrip("# ").rstrip("\n").split("\t")
            try:
                i_feat = header.index("feature")
                i_loc = header.index("locus_tag")
                i_acc = header.index("product_accession")
            except ValueError as exc:
                raise ValueError(f"NCBI feature table {path} lacks a column: {exc}")
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if len(cols) <= max(i_feat, i_loc, i_acc) or cols[i_feat] != "CDS":
                    continue
                if cols[i_loc] and cols[i_acc]:
                    mapping[cols[i_loc]].add(cols[i_acc])
        else:
            for line in [first] + list(f):
                if not line.strip() or line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) >= 2 and cols[0] and cols[1]:
                    mapping[cols[0]].add(cols[1])
    return dict(mapping)


def remap_genes(gene2go, gene_info, id_map):
    """Re-key gene2go / gene_info onto the target ids of *id_map*.

    A source id with several targets (e.g. one locus, several protein
    isoforms) gives each target the same GO set. Unmapped genes are dropped,
    since an id the project data never carries can never be matched.
    """
    new_go = defaultdict(set)
    new_info = {}
    unmapped = 0
    for gene, go_ids in gene2go.items():
        targets = id_map.get(gene)
        if not targets:
            unmapped += 1
            continue
        for t in targets:
            new_go[t].update(go_ids)
            new_info.setdefault(t, gene_info.get(gene, ("", "")))
    if unmapped:
        log_msg("id map: ", unmapped, " of ", len(gene2go),
                " annotated genes have no mapping and were dropped")
    if not new_go:
        raise RuntimeError("id map matched none of the GAF gene ids -- check that "
                           "its source column uses the same id format as the GAF")
    return dict(new_go), new_info


def prepare_gaf(gaf, out_dir, cache_dir, id_map=None, expand=True, refresh=False):
    """Write GO tables + GMTs (and Annotation.tab) from a GAF file.

    gaf:    path to the GAF (plain or .gz).
    id_map: optional path of an id map (see ``parse_id_map``) that re-keys the
            GAF gene ids onto the project's id space.
    cache_dir / refresh: where go-basic.obo is cached, and whether to refetch it.
    """
    gene2go, gene_info = parse_gaf(gaf)
    log_msg("GAF: ", len(gene2go), " genes with GO annotations")
    if not gene2go:
        raise RuntimeError(f"No GO annotations parsed from {gaf} -- is it a GAF file?")

    if id_map:
        gene2go, gene_info = remap_genes(gene2go, gene_info, parse_id_map(id_map))
        log_msg("after id mapping: ", len(gene2go), " ids with GO annotations")

    annot_rows = [(g, *gene_info.get(g, ("", ""))) for g in sorted(gene2go)]
    write_tab(annot_rows, ["Gene", "gene_name", "protein_name"],
              os.path.join(out_dir, "Annotation.tab"))
    write_go_sets(gene2go, out_dir, cache_dir, refresh, expand)
