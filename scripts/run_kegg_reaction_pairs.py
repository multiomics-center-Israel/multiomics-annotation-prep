#!/usr/bin/env python3
"""CLI wrapper: build the KEGG reaction-pair reference for metabolite networks.

Fetches every KEGG reaction, generates cross-side compound pairs
(equation_side_cartesian_product), and writes a dated, checksummed .tsv.gz plus
a sidecar .manifest.json and an .excluded.tsv report. Mirrors the style of
run_mummichog_model.py. This script does NOT publish or attach a release asset.

    python scripts/run_kegg_reaction_pairs.py --out results --cache data
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.prepare_kegg_reaction_pairs import PAIR_METHOD, prepare_kegg_reaction_pairs
from src.utils import log_msg


def main():
    parser = argparse.ArgumentParser(
        description="Build the KEGG reaction-pair reference (cross-side pairs)")
    parser.add_argument("--out", default="results",
                        help="Output directory (default: results)")
    parser.add_argument("--cache", default="data",
                        help="Cache directory (default: data)")
    parser.add_argument("--date", default=None,
                        help="Build date stamp YYYYMMDD (default: today, UTC)")
    parser.add_argument("--refresh", action="store_true",
                        help="Force re-download of cached KEGG files")
    parser.add_argument("--retries", type=int, default=3,
                        help="Retry rounds for unresolved reaction ids (default: 3)")
    parser.add_argument("--rate-limit", type=float, default=0.34,
                        help="Seconds to sleep after each KEGG download (default: 0.34)")
    parser.add_argument("--allow-incomplete", action="store_true",
                        help="DEV ONLY: proceed despite fetch_failed ids. Produces "
                             "NO publishable asset (writes only the failure report).")
    args = parser.parse_args()

    log_msg("=== KEGG reaction-pair reference (method=", PAIR_METHOD, ") ===")
    data_path, manifest_path, excluded_path = prepare_kegg_reaction_pairs(
        out_dir=args.out, cache_dir=args.cache, date=args.date,
        refresh=args.refresh, retries=args.retries, rate_limit_s=args.rate_limit,
        allow_incomplete=args.allow_incomplete)

    if data_path is None:
        # allow_incomplete path: no asset was produced by design.
        log_msg("=== Incomplete build: NO asset produced; see ",
                os.path.basename(excluded_path), " ===")
        sys.exit(1)

    log_msg("=== Done: ", os.path.basename(data_path), " + ",
            os.path.basename(manifest_path), " + ",
            os.path.basename(excluded_path), " ===")


if __name__ == "__main__":
    main()
