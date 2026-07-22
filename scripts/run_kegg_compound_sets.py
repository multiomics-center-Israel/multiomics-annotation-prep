#!/usr/bin/env python3
"""CLI wrapper: build a KEGG compound-set GMT + readable table for ID-based
metabolomics enrichment (ORA / GSEA / QEA in multiomic-core).

Built from the same organism source as the mummichog model, with the same file
stem, so `<stem>.compound_pathway.gmt` sits alongside `<stem>.json` in a Release.
Mirrors the style of run_mummichog_model.py.

    python scripts/run_kegg_compound_sets.py \\
        --kegg-code cre \\
        --model-organism "Chlamydomonas reinhardtii" \\
        --target-organism "Coelastrella sp." \\
        --out results --cache data
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.prepare_kegg_compound_sets import prepare_kegg_compound_sets
from src.utils import log_msg


def main():
    parser = argparse.ArgumentParser(
        description="Build a KEGG compound-set GMT + readable table")
    parser.add_argument("--source", choices=["kegg_org", "kaas"],
                        default="kegg_org",
                        help="Where the KO list comes from (default: kegg_org). "
                             "kegg_org: link/ko/<code>; kaas: a KAAS KO file.")
    parser.add_argument("--kegg-code",
                        help="KEGG organism code to build from (e.g. cre). "
                             "Required for --source kegg_org.")
    parser.add_argument("--model-organism", default=None,
                        help='Species the sets describe (e.g. "Chlamydomonas '
                             'reinhardtii")')
    parser.add_argument("--model-kegg-code", default=None,
                        help="Short label used in the output filename stem "
                             "(defaults to --kegg-code; required for "
                             "--source kaas, e.g. 'coel')")
    parser.add_argument("--target-organism", default=None,
                        help='Biology the sets stand in for (e.g. '
                             '"Coelastrella sp."); leave unset if same')
    parser.add_argument("--source-version", default=None,
                        help="Override the recorded source version "
                             "(else auto-detected from KEGG)")
    parser.add_argument("--date", default=None,
                        help="Build date stamp YYYYMMDD (default: today, UTC)")
    parser.add_argument("--out", default="results",
                        help="Output directory (default: results)")
    parser.add_argument("--cache", default="data",
                        help="Cache directory (default: data)")
    parser.add_argument("--refresh", action="store_true",
                        help="Force re-download of cached KEGG files")
    parser.add_argument("--kaas", default=None,
                        help="KAAS KO list (gene<TAB>KO) for --source kaas")
    args = parser.parse_args()

    if args.source == "kegg_org" and not args.kegg_code:
        parser.error("--kegg-code is required for --source kegg_org")
    if args.source == "kaas":
        if not args.kaas:
            parser.error("--kaas <query.ko.txt> is required for --source kaas")
        if not (args.model_kegg_code or args.kegg_code):
            parser.error("--model-kegg-code (filename label, e.g. 'coel') is "
                         "required for --source kaas")

    log_msg("=== KEGG compound sets (source=", args.source, ") ===")
    gmt_path, tab_path, manifest_path = prepare_kegg_compound_sets(
        kegg_code=args.kegg_code,
        out_dir=args.out,
        cache_dir=args.cache,
        source=args.source,
        model_kegg_code=args.model_kegg_code,
        model_organism=args.model_organism,
        target_organism=args.target_organism,
        source_version=args.source_version,
        date=args.date,
        refresh=args.refresh,
        kaas_file=args.kaas,
    )
    log_msg("=== Done: ", os.path.basename(gmt_path), " + ",
            os.path.basename(tab_path), " + ",
            os.path.basename(manifest_path), " ===")


if __name__ == "__main__":
    main()
