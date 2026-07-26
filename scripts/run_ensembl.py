#!/usr/bin/env python3
"""CLI wrapper: Ensembl/BioMart GO (+ optional KEGG) annotation for a model organism.

Writes the same enrichment tables the pipeline reads from annotation_dir:
GO2gene_{BP,MF,CC}.tab + GO2name_{BP,MF,CC}.tab (hierarchy-expanded), a
descriptive Annotation.tab, and -- when --kegg-org is given -- the
KEGG_pathway2gene/2name tables. Needs `pybiomart`.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.prepare_ensembl import prepare_ensembl
from src.utils import log_msg


def main():
    parser = argparse.ArgumentParser(
        description="Prepare GO (+ optional KEGG) annotation tables from Ensembl/BioMart")
    parser.add_argument("--dataset", required=True,
                        help="BioMart dataset (e.g. athaliana_eg_gene, mmusculus_gene_ensembl)")
    parser.add_argument("--mart", default="ensembl",
                        help="BioMart mart (ensembl | plants_mart | protists_mart | ...)")
    parser.add_argument("--host", default=None,
                        help="BioMart host (default: www.ensembl.org; "
                             "e.g. https://plants.ensembl.org for Ensembl Plants)")
    parser.add_argument("--kegg-org", default=None,
                        help="Optional KEGG organism code -> also build the "
                             "KEGG_pathway2gene/2name tables (e.g. ath, cre)")
    parser.add_argument("--id-source", default="ensembl",
                        help="id space for KEGG mapping (default: ensembl)")
    parser.add_argument("--out", default="results",
                        help="Output directory (default: results)")
    parser.add_argument("--cache", default="data",
                        help="Cache directory (default: data)")
    parser.add_argument("--no-expand", action="store_true",
                        help="Skip GO hierarchy expansion (default: expand)")
    parser.add_argument("--refresh", action="store_true",
                        help="Force re-download of cached files")
    args = parser.parse_args()

    log_msg("=== Ensembl/BioMart annotation ===")
    prepare_ensembl(
        dataset=args.dataset,
        out_dir=args.out,
        cache_dir=args.cache,
        mart=args.mart,
        host=args.host,
        kegg_org=args.kegg_org,
        id_source=args.id_source,
        expand=not args.no_expand,
        refresh=args.refresh,
    )
    log_msg("=== Done ===")


if __name__ == "__main__":
    main()
