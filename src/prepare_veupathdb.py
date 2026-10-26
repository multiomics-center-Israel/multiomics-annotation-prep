"""VEuPathDB path (GiardiaDB, PlasmoDB, ToxoDB, ...) -- GO from the site's GAF.

Organisms hosted by VEuPathDB are usually absent from Ensembl/UniProt in a form
whose IDs match a project (e.g. Giardia's GL50803_* locus tags), but each site
publishes a per-organism GO annotation file (GAF). This module downloads it,
optionally re-keys its gene ids onto the project's id space (e.g. the NCBI
protein accessions a DIA-NN FASTA carries), and writes the same GO tables and
GMTs as the other GO paths via ``prepare_go.write_go_sets``.
"""

import gzip
import os
import re
from collections import defaultdict
from urllib.parse import urlparse

from .prepare_go import write_go_sets
from .utils import cached_download, log_msg, write_tab

# lower-case host label -> the capitalised project name used in file names
SITES = {
    "amoebadb": "AmoebaDB",
    "cryptodb": "CryptoDB",
    "fungidb": "FungiDB",
    "giardiadb": "GiardiaDB",
    "hostdb": "HostDB",
    "microsporidiadb": "MicrosporidiaDB",
    "piroplasmadb": "PiroplasmaDB",
    "plasmodb": "PlasmoDB",
    "toxodb": "ToxoDB",
    "trichdb": "TrichDB",
    "tritrypdb": "TriTrypDB",
    "vectorbase": "VectorBase",
}

_SAFE = re.compile(r"^[A-Za-z0-9_][A-Za-z0-9_.-]*$")


def gaf_url(site, organism, release="CURRENT"):
    """Build the VEuPathDB download URL of an organism's GO GAF.

    ``release`` is ``CURRENT`` (moving) or a release number such as ``68``
    (pinned, reproducible).
    """
    key = site.lower()
    if key not in SITES:
        raise ValueError(f"Unknown VEuPathDB site '{site}'. "
                         f"Known: {', '.join(sorted(SITES))}")
    release = str(release)
    for label, value in (("organism", organism), ("release", release)):
        if not _SAFE.match(value):
            raise ValueError(f"VEuPathDB {label} has unexpected characters: '{value}'")
    if release.upper() == "CURRENT":
        rel_dir, rel_tag = "Current_Release", "CURRENT"
    else:
        rel_dir, rel_tag = f"release-{release}", release
    return (f"https://{key}.org/common/downloads/{rel_dir}/{organism}/gaf/"
            f"{SITES[key]}-{rel_tag}_{organism}_GO.gaf")


def _fetch(src, cache_dir, refresh):
    """Return a local path for *src*: downloaded (cached) if a URL, else as-is."""
    if re.match(r"^https?://", src):
        name = os.path.basename(urlparse(src).path) or "download"
        return cached_download(src, name, cache_dir, refresh)
    return src


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


def prepare_veupathdb(site, organism, out_dir, cache_dir, release="CURRENT",
                      gaf=None, id_map=None, expand=True, refresh=False):
    """Write GO tables + GMTs (and Annotation.tab) for a VEuPathDB organism.

    gaf:    optional path or URL overriding the GAF built from site/organism.
    id_map: optional path or URL of an id map (see ``parse_id_map``) that
            re-keys the GAF gene ids onto the project's id space.
    """
    src = gaf or gaf_url(site, organism, release)
    gaf_path = _fetch(src, cache_dir, refresh)
    gene2go, gene_info = parse_gaf(gaf_path)
    log_msg("GAF: ", len(gene2go), " genes with GO annotations")
    if not gene2go:
        raise RuntimeError(f"No GO annotations parsed from {src}")

    if id_map:
        gene2go, gene_info = remap_genes(
            gene2go, gene_info, parse_id_map(_fetch(id_map, cache_dir, refresh)))
        log_msg("after id mapping: ", len(gene2go), " ids with GO annotations")

    annot_rows = [(g, *gene_info.get(g, ("", ""))) for g in sorted(gene2go)]
    write_tab(annot_rows, ["Gene", "gene_name", "protein_name"],
              os.path.join(out_dir, "Annotation.tab"))
    write_go_sets(gene2go, out_dir, cache_dir, refresh, expand)
