#!/usr/bin/env python3
"""Generate a GMT gene-set file for an organism — Python port of R/generate_gmt.R.

KEGG pathways come from the KEGG REST API (any organism code). GO gene sets come
from Ensembl/biomaRt (model organisms only, needs `pybiomart`); for non-model
organisms build GO from a GAF/eggNOG table with build_ankri_go_gmt.py /
build_elad_go_gmt.py instead.

Usage:
  python scripts/generate_gmt.py --kegg-org ehi --output results/ehi_kegg.gmt
  python scripts/generate_gmt.py --organism "Danio rerio" --database GO --ontology BP
  python scripts/generate_gmt.py --list-kegg --search elegans
"""

import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.gmt_utils import (KEGG_REST, filter_gmt_by_size, generate_gmt_from_go,
                           generate_gmt_from_kegg, generate_gmt_from_reactome,
                           write_gmt)
from src.utils import cached_download, log_msg


def list_kegg_organisms(pattern=None, cache_dir="data"):
    """Print KEGG organism codes and names, optionally filtered by a pattern."""
    path = cached_download(f"{KEGG_REST}/list/organism", "kegg_organisms.txt",
                           cache_dir, refresh=False)
    rows = []
    with open(path) as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 3:
                rows.append((cols[1], cols[2]))  # (code, name)
    if pattern:
        rx = re.compile(pattern, re.IGNORECASE)
        rows = [r for r in rows if rx.search(r[1])]

    print(f"\n{'Code':<6} {'Organism':<60}")
    print("-" * 68)
    for code, name in rows[:50]:
        print(f"{code:<6} {name[:60]:<60}")
    if len(rows) > 50:
        print(f"\n... and {len(rows) - 50} more. Use --search to filter.")


def main():
    parser = argparse.ArgumentParser(
        description="Generate a GMT gene-set file for an organism")
    parser.add_argument("-o", "--organism",
                        help="Organism scientific name (e.g. 'Danio rerio')")
    parser.add_argument("-k", "--kegg-org", dest="kegg_org",
                        help="KEGG organism code (e.g. 'dre', 'ehi')")
    parser.add_argument("-d", "--database", default="GO",
                        help="Database: GO, Reactome, KEGG [default: GO]")
    parser.add_argument("--ontology", default="BP",
                        help="GO ontology: BP, MF, CC, all [default: BP]")
    parser.add_argument("-i", "--id-type", dest="id_type", default="symbol",
                        help="Gene ID type: symbol, ensembl, entrez [default: symbol]")
    parser.add_argument("--output", help="Output GMT file path")
    parser.add_argument("--min-genes", dest="min_genes", type=int, default=10,
                        help="Minimum genes per pathway [default: 10]")
    parser.add_argument("--max-genes", dest="max_genes", type=int, default=500,
                        help="Maximum genes per pathway [default: 500]")
    parser.add_argument("--cache", default="data", help="Download cache dir")
    parser.add_argument("--refresh", action="store_true",
                        help="Force re-download of cached files")
    parser.add_argument("--list-kegg", action="store_true",
                        help="List available KEGG organisms")
    parser.add_argument("--search", help="Search pattern for organism lists")
    args = parser.parse_args()

    if args.list_kegg:
        list_kegg_organisms(args.search, args.cache)
        return

    if args.kegg_org:
        log_msg(f"=== Generating KEGG GMT for organism: {args.kegg_org} ===")
        output = args.output or f"kegg_{args.kegg_org}.gmt"
        pathways = generate_gmt_from_kegg(
            args.kegg_org, id_type=args.id_type, cache_dir=args.cache,
            refresh=args.refresh, name_strip=rf" - {re.escape(args.kegg_org)}$")
        write_gmt(pathways, output)
        log_msg(f"Unique genes: {len({g for genes in pathways.values() for g in genes})}")
    elif args.organism:
        log_msg(f"=== Generating {args.database} GMT for {args.organism} ===")
        base = args.organism.lower().replace(" ", "_")
        if args.database.upper() == "REACTOME":
            output = args.output or f"{base}_Reactome.gmt"
            generate_gmt_from_reactome(
                args.organism, id_type=args.id_type, output_file=output,
                min_genes=args.min_genes, max_genes=args.max_genes,
                cache_dir=args.cache, refresh=args.refresh)
        else:
            output = args.output or f"{base}_{args.database}_{args.ontology}.gmt"
            generate_gmt_from_go(
                args.organism, id_type=args.id_type, ontology=args.ontology,
                output_file=output, min_genes=args.min_genes,
                max_genes=args.max_genes, cache_dir=args.cache,
                refresh=args.refresh)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
