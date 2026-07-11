#!/usr/bin/env python3
"""CLI wrapper: build an organism metabolic model for mummichog (`mummichog -n`).

Primary path is a KEGG organism code (e.g. cre for Chlamydomonas reinhardtii,
standing in for Coelastrella). Writes <stem>.json + <stem>.manifest.json per
MODEL_CONTRACT.md. Mirrors the style of run_kegg_nonmodel.py.

    python scripts/run_mummichog_model.py \\
        --kegg-code cre \\
        --model-organism "Chlamydomonas reinhardtii" \\
        --target-organism "Coelastrella sp." \\
        --out results --cache data --validate
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.prepare_mummichog_model import prepare_mummichog_model
from src.utils import log_msg


def main():
    parser = argparse.ArgumentParser(
        description="Build an organism metabolic model for mummichog (-n)")
    parser.add_argument("--source", choices=["kegg_org", "kaas"],
                        default="kegg_org",
                        help="Input source (default: kegg_org). "
                             "kaas is a planned seam, not yet implemented.")
    parser.add_argument("--kegg-code",
                        help="KEGG organism code to build from (e.g. cre). "
                             "Required for --source kegg_org.")
    parser.add_argument("--model-organism", default=None,
                        help='Species the model IS (e.g. "Chlamydomonas '
                             'reinhardtii")')
    parser.add_argument("--model-kegg-code", default=None,
                        help="KEGG code recorded for the model "
                             "(defaults to --kegg-code)")
    parser.add_argument("--target-organism", default=None,
                        help='Biology the model stands in for (e.g. '
                             '"Coelastrella sp."); leave unset if same as model')
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
    parser.add_argument("--validate", action="store_true",
                        help="Run `mummichog -n` on a synthetic feature table "
                             "and record the result (needs mummichog installed)")
    parser.add_argument("--kaas", default=None,
                        help="(future) KAAS KO list for --source kaas")
    args = parser.parse_args()

    if args.source == "kegg_org" and not args.kegg_code:
        parser.error("--kegg-code is required for --source kegg_org")

    log_msg("=== mummichog metabolic model (source=", args.source, ") ===")
    model_path, manifest_path = prepare_mummichog_model(
        kegg_code=args.kegg_code,
        out_dir=args.out,
        cache_dir=args.cache,
        source=args.source,
        model_organism=args.model_organism,
        model_kegg_code=args.model_kegg_code,
        target_organism=args.target_organism,
        source_version=args.source_version,
        date=args.date,
        refresh=args.refresh,
        validate=args.validate,
        kaas_file=args.kaas,
    )
    log_msg("=== Done: ", os.path.basename(model_path), " + ",
            os.path.basename(manifest_path), " ===")


if __name__ == "__main__":
    main()
