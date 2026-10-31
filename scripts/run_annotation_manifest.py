#!/usr/bin/env python3
"""CLI wrapper: write manifest.json into a built annotation_dir (before zipping)."""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.annotation_manifest import GENE_ANNOTATION_SOURCES, write_annotation_manifest
from src.utils import log_msg


def _key_value(text):
    key, sep, value = text.partition("=")
    if not sep or not key:
        raise argparse.ArgumentTypeError(f"expected KEY=VALUE, got {text!r}")
    return key, value


def main():
    parser = argparse.ArgumentParser(
        description="Write manifest.json (files, sha256, source versions) into "
                    "a built annotation_dir")
    parser.add_argument("--out", default="results",
                        help="Built output directory (default: results)")
    parser.add_argument("--cache", default="data",
                        help="Cache directory the build used (default: data)")
    parser.add_argument("--release-tag", default=None, help="Release tag")
    parser.add_argument("--run-url", default=None, help="URL of the build run")
    parser.add_argument("--gene-annotation-source", default=None,
                        choices=GENE_ANNOTATION_SOURCES,
                        help="Where the gene -> GO links came from (omit if no GO)")
    parser.add_argument("--ensembl-host", default=None,
                        help="[ensembl] BioMart host the build used (default: www.ensembl.org)")
    parser.add_argument("--ensembl-mart", default="ensembl",
                        help="[ensembl] BioMart mart the build used")
    parser.add_argument("--kegg-ids-from-gene-annotation", action="store_true",
                        help="The build re-keyed KEGG gene ids through the gene "
                             "annotation source (Ensembl bridge)")
    parser.add_argument("--query-source-versions", action="store_true",
                        help="Query the KEGG database dates and the Ensembl release "
                             "name now. Use only when the build fetched its data in "
                             "this same run (fresh cache); otherwise they stay null")
    parser.add_argument("--input", action="append", type=_key_value, default=[],
                        metavar="KEY=VALUE",
                        help="Build input to record verbatim (repeatable)")
    args = parser.parse_args()

    log_msg("=== annotation manifest ===")
    write_annotation_manifest(
        args.out, args.cache,
        release_tag=args.release_tag,
        run_url=args.run_url,
        inputs=dict(args.input),
        gene_annotation_source=args.gene_annotation_source,
        ensembl_host=args.ensembl_host or None,
        ensembl_mart=args.ensembl_mart or "ensembl",
        kegg_ids_from_gene_annotation=args.kegg_ids_from_gene_annotation,
        query_versions=args.query_source_versions,
    )
    log_msg("=== Done ===")


if __name__ == "__main__":
    main()
