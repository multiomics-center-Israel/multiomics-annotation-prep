#!/usr/bin/env python3
"""CLI wrapper: non-model organism KEGG annotation from KAAS output."""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.prepare_kegg_nonmodel import prepare_kegg_nonmodel
from src.utils import log_msg


def main():
    parser = argparse.ArgumentParser(
        description="Prepare KEGG annotation files from KAAS output")
    parser.add_argument("--kaas", required=True,
                        help="Path to KAAS query.ko.txt file")
    parser.add_argument("--out", default="results",
                        help="Output directory (default: results)")
    parser.add_argument("--cache", default="data",
                        help="Cache directory (default: data)")
    parser.add_argument("--genes", default=None,
                        help="Optional file with gene IDs (one per line) "
                             "for descriptive annotation")
    parser.add_argument("--refresh", action="store_true",
                        help="Force re-download of cached files")
    parser.add_argument("--no-strip-isoform", action="store_true",
                        help="Keep _iN isoform suffixes (default: strip them)")
    args = parser.parse_args()

    gene_list = None
    if args.genes:
        with open(args.genes) as f:
            gene_list = [line.strip() for line in f if line.strip()]

    log_msg("=== KEGG non-model annotation ===")
    prepare_kegg_nonmodel(
        kaas_file=args.kaas,
        out_dir=args.out,
        cache_dir=args.cache,
        refresh=args.refresh,
        strip_isoform=not args.no_strip_isoform,
        gene_list=gene_list,
    )
    log_msg("=== Done ===")


if __name__ == "__main__":
    main()
