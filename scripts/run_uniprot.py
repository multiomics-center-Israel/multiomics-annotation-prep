#!/usr/bin/env python3
"""CLI wrapper: UniProt GO (+ optional KEGG) annotation by NCBI taxon id.

Writes the same enrichment tables the pipeline reads from annotation_dir:
GO2gene_{BP,MF,CC}.tab + GO2name_{BP,MF,CC}.tab (hierarchy-expanded), keyed on
UniProt accessions, a descriptive Annotation.tab, and -- when --kegg-org is
given -- the KEGG_pathway2gene/2name tables. Light deps (requests + pyyaml).
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.prepare_uniprot import prepare_uniprot
from src.utils import log_msg


def main():
    parser = argparse.ArgumentParser(
        description="Prepare GO (+ optional KEGG) annotation tables from UniProt")
    parser.add_argument("--taxon-id", required=True,
                        help="NCBI taxon id (e.g. 3702 = Arabidopsis thaliana)")
    parser.add_argument("--no-reviewed", action="store_true",
                        help="Include unreviewed (TrEMBL) too; default: Swiss-Prot only")
    parser.add_argument("--kegg-org", default=None,
                        help="Optional KEGG organism code -> also build the "
                             "KEGG_pathway2gene/2name tables (e.g. ath)")
    parser.add_argument("--out", default="results",
                        help="Output directory (default: results)")
    parser.add_argument("--cache", default="data",
                        help="Cache directory (default: data)")
    parser.add_argument("--no-expand", action="store_true",
                        help="Skip GO hierarchy expansion (default: expand)")
    parser.add_argument("--refresh", action="store_true",
                        help="Force re-download of cached files")
    args = parser.parse_args()

    log_msg("=== UniProt annotation ===")
    prepare_uniprot(
        taxon_id=args.taxon_id,
        out_dir=args.out,
        cache_dir=args.cache,
        reviewed=not args.no_reviewed,
        kegg_org=args.kegg_org,
        expand=not args.no_expand,
        refresh=args.refresh,
    )
    log_msg("=== Done ===")


if __name__ == "__main__":
    main()
