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

import gzip
import json
import os
import re
from collections import defaultdict
from datetime import datetime, timezone

from .download_go import GO_OBO_URL
from .prepare_go import write_go_sets
from .utils import cached_download, git_head_sha, log_msg, sha256_file, write_tab

NCBI_GENOMES = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
_ASSEMBLY = re.compile(r"^(GC[AF])_(\d{3})(\d{3})(\d{3})\.\d+$")


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


def pick_assembly_dir(listing_html, accession):
    """Pick ``<accession>_<asm_name>`` from an NCBI genomes directory listing."""
    names = sorted(set(re.findall(r'href="(GC[AF]_[^"/]+)/"', listing_html)))
    hits = [n for n in names if n.startswith(accession + "_")]
    if not hits:
        raise RuntimeError(f"Assembly {accession} not found on NCBI; folders listed: "
                           f"{', '.join(names) or 'none'}")
    return hits[0]


def ncbi_feature_table_url(accession, cache_dir, refresh=False):
    """Resolve the ``*_feature_table.txt.gz`` URL of an NCBI assembly accession.

    The assembly name is part of the path (``GCA_000002435.2_<name>``), so it is
    read from the parent directory listing rather than guessed.
    """
    m = _ASSEMBLY.match(accession)
    if not m:
        raise ValueError(f"'{accession}' is not a versioned NCBI assembly accession "
                         f"(expected e.g. GCA_000002435.2)")
    parent = f"{NCBI_GENOMES}/{m.group(1)}/{m.group(2)}/{m.group(3)}/{m.group(4)}/"
    index = cached_download(parent, f"ncbi_{accession}_index.html", cache_dir, refresh)
    with open(index) as f:
        folder = pick_assembly_dir(f.read(), accession)
    return f"{parent}{folder}/{folder}_feature_table.txt.gz"


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


def prepare_gaf(gaf, out_dir, cache_dir, id_map=None, ncbi_assembly=None,
                expand=True, refresh=False):
    """Write GO tables + GMTs, Annotation.tab and GO_provenance.json from a GAF.

    gaf:           path to the GAF (plain or .gz).
    id_map:        optional path of an id map (see ``parse_id_map``) that
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

    map_info = None
    if ncbi_assembly:
        url = ncbi_feature_table_url(ncbi_assembly, cache_dir, refresh)
        id_map = cached_download(url, os.path.basename(url), cache_dir, refresh)
        map_info = {"source": "ncbi_assembly", "assembly": ncbi_assembly, "url": url}
    elif id_map:
        map_info = {"source": "file"}
    if id_map:
        map_info.update(file=os.path.basename(id_map), sha256=sha256_file(id_map))
        gene2go, gene_info = remap_genes(gene2go, gene_info, parse_id_map(id_map))
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
