#!/usr/bin/env python3
"""CLI wrapper: KEGG pathway gene sets for a KEGG organism code (no KAAS needed).

Writes KEGG_pathway2gene.tab + KEGG_pathway2name.tab + KEGG_pathway.gmt and
KEGG_provenance.json. With --id-map or --ncbi-assembly the KEGG gene ids are
re-keyed (e.g. Giardia locus tags -> the NCBI protein accessions a proteomics
FASTA uses). Light deps (requests + pyyaml).
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.prepare_kegg_org_sets import prepare_kegg_org_sets
from src.utils import log_msg


def main():
    parser = argparse.ArgumentParser(
        description="Prepare KEGG pathway gene sets (tables + GMT) for a KEGG organism")
    parser.add_argument("--kegg-org", required=True,
                        help="KEGG organism code, e.g. gla (Giardia lamblia)")
    mapping = parser.add_mutually_exclusive_group()
    mapping.add_argument("--id-map", default=None,
                         help="Optional id map: an NCBI *_feature_table.txt(.gz) "
                              "(locus_tag -> protein accession) or a two-column TSV "
                              "(KEGG gene id -> project id)")
    mapping.add_argument("--ncbi-assembly", default=None,
                         help="Optional NCBI assembly accession (e.g. GCA_000002435.2): "
                              "download its feature table and use it as the id map")
    parser.add_argument("--out", default="results",
                        help="Output directory (default: results)")
    parser.add_argument("--cache", default="data",
                        help="Cache directory (default: data)")
    parser.add_argument("--refresh", action="store_true",
                        help="Force re-download of cached files")
    args = parser.parse_args()

    log_msg("=== KEGG pathway sets (organism) ===")
    prepare_kegg_org_sets(
        kegg_org=args.kegg_org,
        out_dir=args.out,
        cache_dir=args.cache,
        id_map=args.id_map,
        ncbi_assembly=args.ncbi_assembly,
        refresh=args.refresh,
    )
    log_msg("=== Done ===")


if __name__ == "__main__":
    main()
