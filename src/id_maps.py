"""Id maps that re-key annotation onto a project's id space.

Shared by every gene-based source that may need its ids translated before
they can match the data (a GAF's locus tags, KEGG's organism gene ids, ...),
so all of them resolve, parse and record the map the same way. The NCBI
feature table is public, so it can be fetched by assembly accession; any
other mapping is a committed file.
"""

import gzip
import os
import re
from collections import defaultdict

from .utils import cached_download, sha256_file

NCBI_GENOMES = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
_ASSEMBLY = re.compile(r"^(GC[AF])_(\d{3})(\d{3})(\d{3})\.\d+$")


def _open_text(path):
    """Open a plain or gzip-compressed text file (detected by magic bytes)."""
    with open(path, "rb") as f:
        is_gz = f.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if is_gz else open(path)


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


def resolve_id_map(id_map, ncbi_assembly, cache_dir, refresh=False):
    """Return ``(mapping, map_info)`` for an id-map file or an NCBI assembly.

    mapping:  source id -> set of target ids, or ``None`` when neither is given.
    map_info: provenance for the map (source, file, sha256, and the assembly
              accession + URL when fetched from NCBI), or ``None``.
    """
    if id_map and ncbi_assembly:
        raise ValueError("Give either id_map or ncbi_assembly, not both")
    if ncbi_assembly:
        url = ncbi_feature_table_url(ncbi_assembly, cache_dir, refresh)
        id_map = cached_download(url, os.path.basename(url), cache_dir, refresh)
        map_info = {"source": "ncbi_assembly", "assembly": ncbi_assembly, "url": url}
    elif id_map:
        map_info = {"source": "file"}
    else:
        return None, None
    map_info.update(file=os.path.basename(id_map), sha256=sha256_file(id_map))
    return parse_id_map(id_map), map_info
