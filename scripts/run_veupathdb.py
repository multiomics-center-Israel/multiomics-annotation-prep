#!/usr/bin/env python3
"""CLI wrapper: GO annotation for a VEuPathDB organism (GiardiaDB, PlasmoDB, ...).

Downloads the organism's GO GAF from its VEuPathDB site and writes
GO2gene_{BP,MF,CC}.tab + GO2name_{BP,MF,CC}.tab + GO_{BP,MF,CC}.gmt
(hierarchy-expanded) and a descriptive Annotation.tab. With --id-map the gene
ids are re-keyed (e.g. Giardia locus tags -> the NCBI protein accessions a
proteomics FASTA uses). Light deps (requests + pyyaml).
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.prepare_veupathdb import SITES, prepare_veupathdb
from src.utils import log_msg


def main():
    parser = argparse.ArgumentParser(
        description="Prepare GO annotation tables + GMTs from a VEuPathDB GAF")
    parser.add_argument("--site", required=True, choices=sorted(SITES),
                        help="VEuPathDB site, e.g. giardiadb")
    parser.add_argument("--organism", required=True,
                        help="Organism as named in the site's download tree, "
                             "e.g. GintestinalisAssemblageAWB")
    parser.add_argument("--release", default="CURRENT",
                        help="VEuPathDB release number to pin (e.g. 68); "
                             "default: CURRENT")
    parser.add_argument("--gaf", default=None,
                        help="Path or URL of a GAF to use instead of the one "
                             "built from --site/--organism/--release")
    parser.add_argument("--id-map", default=None,
                        help="Path or URL of an id map: an NCBI *_feature_table.txt(.gz) "
                             "(locus_tag -> protein accession) or a two-column TSV "
                             "(GAF gene id -> project id)")
    parser.add_argument("--out", default="results",
                        help="Output directory (default: results)")
    parser.add_argument("--cache", default="data",
                        help="Cache directory (default: data)")
    parser.add_argument("--no-expand", action="store_true",
                        help="Skip GO hierarchy expansion (default: expand)")
    parser.add_argument("--refresh", action="store_true",
                        help="Force re-download of cached files")
    args = parser.parse_args()

    log_msg("=== VEuPathDB annotation ===")
    prepare_veupathdb(
        site=args.site,
        organism=args.organism,
        out_dir=args.out,
        cache_dir=args.cache,
        release=args.release,
        gaf=args.gaf,
        id_map=args.id_map,
        expand=not args.no_expand,
        refresh=args.refresh,
    )
    log_msg("=== Done ===")


if __name__ == "__main__":
    main()
