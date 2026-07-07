#!/usr/bin/env python3
"""CLI wrapper: GO enrichment file preparation."""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.prepare_go import prepare_go
from src.utils import log_msg


def main():
    parser = argparse.ArgumentParser(
        description="Prepare GO enrichment files from a per-gene GO table")
    parser.add_argument("--go-table", required=True,
                        help="Path to per-gene GO table (Trinotate style)")
    parser.add_argument("--out", default="results",
                        help="Output directory (default: results)")
    parser.add_argument("--cache", default="data",
                        help="Cache directory (default: data)")
    parser.add_argument("--no-expand", action="store_true",
                        help="Skip GO hierarchy expansion")
    parser.add_argument("--refresh", action="store_true",
                        help="Force re-download of cached files")
    args = parser.parse_args()

    log_msg("=== GO annotation ===")
    prepare_go(
        go_table_file=args.go_table,
        out_dir=args.out,
        cache_dir=args.cache,
        refresh=args.refresh,
        expand=not args.no_expand,
    )
    log_msg("=== Done ===")


if __name__ == "__main__":
    main()
