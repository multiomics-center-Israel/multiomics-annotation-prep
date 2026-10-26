#!/usr/bin/env python3
"""CLI wrapper: GO enrichment files from a GAF (e.g. downloaded from GiardiaDB).

Writes GO2gene_{BP,MF,CC}.tab + GO2name_{BP,MF,CC}.tab + GO_{BP,MF,CC}.gmt
(hierarchy-expanded), a descriptive Annotation.tab and GO_provenance.json.
With --id-map or --ncbi-assembly the gene ids are re-keyed (e.g. Giardia locus
tags -> the NCBI protein accessions a proteomics FASTA uses). Light deps
(requests + pyyaml).
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.prepare_gaf import prepare_gaf
from src.utils import log_msg


def main():
    parser = argparse.ArgumentParser(
        description="Prepare GO enrichment files (tables + GMTs) from a GAF")
    parser.add_argument("--gaf", required=True,
                        help="Path to the GAF file (plain or .gz)")
    mapping = parser.add_mutually_exclusive_group()
    mapping.add_argument("--id-map", default=None,
                         help="Optional id map: an NCBI *_feature_table.txt(.gz) "
                              "(locus_tag -> protein accession) or a two-column TSV "
                              "(GAF gene id -> project id)")
    mapping.add_argument("--ncbi-assembly", default=None,
                         help="Optional NCBI assembly accession (e.g. GCA_000002435.2): "
                              "download its feature table and use it as the id map")
    parser.add_argument("--out", default="results",
                        help="Output directory (default: results)")
    parser.add_argument("--cache", default="data",
                        help="Cache directory (default: data)")
    parser.add_argument("--no-expand", action="store_true",
                        help="Skip GO hierarchy expansion (default: expand)")
    parser.add_argument("--refresh", action="store_true",
                        help="Force re-download of cached files")
    args = parser.parse_args()

    log_msg("=== GO annotation (GAF) ===")
    prepare_gaf(
        gaf=args.gaf,
        out_dir=args.out,
        cache_dir=args.cache,
        id_map=args.id_map,
        ncbi_assembly=args.ncbi_assembly,
        expand=not args.no_expand,
        refresh=args.refresh,
    )
    log_msg("=== Done ===")


if __name__ == "__main__":
    main()
