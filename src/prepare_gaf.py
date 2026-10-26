"""GO from a GAF file (GO Annotation File, e.g. downloaded from GiardiaDB/VEuPathDB).

Organism databases such as the VEuPathDB sites publish per-organism GO
annotation as GAF, but downloads need a logged-in browser session, so the file
is fetched by hand and committed as a project input. This module parses it,
optionally re-keys its gene ids onto the project's id space (e.g. the NCBI
protein accessions a proteomics FASTA carries -- the NCBI feature table that
maps them is public, so it can be fetched by assembly accession), and writes
the same GO tables and GMTs as the other GO paths via
``prepare_go.write_go_sets``, plus a ``GO_provenance.json`` recording where
every input came from.
"""

import json
import os
from collections import defaultdict
from datetime import datetime, timezone

from .download_go import GO_OBO_URL
from .id_maps import _open_text, resolve_id_map
from .prepare_go import write_go_sets
from .utils import git_head_sha, log_msg, sha256_file, write_tab


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


def gaf_header(path, limit=50):
    """Return the GAF's ``!`` header lines (source, release, date generated)."""
    lines = []
    with _open_text(path) as f:
        for line in f:
            if not line.startswith("!"):
                break
            if len(lines) < limit:
                lines.append(line.rstrip("\n"))
    return lines


def obo_data_version(path):
    """Return the ``data-version`` of a GO OBO file, or ``"unknown"``."""
    with open(path) as f:
        for line in f:
            if line.startswith("data-version:"):
                return line.split(":", 1)[1].strip()
            if line.startswith("[Term]"):
                break
    return "unknown"


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


def prepare_gaf(gaf, out_dir, cache_dir, id_map=None, ncbi_assembly=None,
                expand=True, refresh=False):
    """Write GO tables + GMTs, Annotation.tab and GO_provenance.json from a GAF.

    gaf:           path to the GAF (plain or .gz).
    id_map:        optional path of an id map (see ``id_maps.parse_id_map``) that
                   re-keys the GAF gene ids onto the project's id space.
    ncbi_assembly: alternative to ``id_map``: a versioned NCBI assembly
                   accession whose public feature table is downloaded and used
                   as the id map (locus_tag -> protein accession).
    cache_dir / refresh: where downloads are cached, and whether to refetch.
    """
    if id_map and ncbi_assembly:
        raise ValueError("Give either id_map or ncbi_assembly, not both")

    gene2go, gene_info = parse_gaf(gaf)
    log_msg("GAF: ", len(gene2go), " genes with GO annotations")
    if not gene2go:
        raise RuntimeError(f"No GO annotations parsed from {gaf} -- is it a GAF file?")

    mapping, map_info = resolve_id_map(id_map, ncbi_assembly, cache_dir, refresh)
    if mapping is not None:
        gene2go, gene_info = remap_genes(gene2go, gene_info, mapping)
        log_msg("after id mapping: ", len(gene2go), " ids with GO annotations")

    annot_rows = [(g, *gene_info.get(g, ("", ""))) for g in sorted(gene2go)]
    write_tab(annot_rows, ["Gene", "gene_name", "protein_name"],
              os.path.join(out_dir, "Annotation.tab"))
    write_go_sets(gene2go, out_dir, cache_dir, refresh, expand)

    obo = os.path.join(cache_dir, "go-basic.obo")
    provenance = {
        "go_source": "GAF",
        "gaf": {"file": os.path.basename(gaf), "sha256": sha256_file(gaf),
                "header": gaf_header(gaf)},
        "id_map": map_info,
        "go_ontology": {"url": GO_OBO_URL, "data_version": obo_data_version(obo),
                        "sha256": sha256_file(obo)},
        "hierarchy_expanded": expand,
        "annotated_ids": len(gene2go),
        "builder_git_sha": git_head_sha() or "unknown",
        "build_timestamp_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }
    path = os.path.join(out_dir, "GO_provenance.json")
    with open(path, "w") as f:
        json.dump(provenance, f, indent=2)
    log_msg("wrote provenance: ", path)
