#!/usr/bin/env python3
"""Config-driven runner: reads config/config.yml and runs selected modules."""

import argparse
import os
import sys

import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.utils import log_msg


def main():
    parser = argparse.ArgumentParser(
        description="Run annotation modules based on config file")
    parser.add_argument("--config", default="config/config.yml",
                        help="Path to config YAML file")
    args = parser.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    out_dir = cfg.get("out_dir", "results")
    cache_dir = cfg.get("cache_dir", "data")
    refresh = cfg.get("refresh", False)
    expand_go = cfg.get("expand_go", True)
    modules = cfg.get("modules", {})

    # KEGG non-model
    kegg_nm = modules.get("kegg_nonmodel", {})
    if kegg_nm.get("enabled", False):
        from src.prepare_kegg_nonmodel import prepare_kegg_nonmodel
        log_msg("=== Module: kegg_nonmodel ===")
        gene_list = None
        genes_file = kegg_nm.get("genes")
        if genes_file:
            with open(genes_file) as f:
                gene_list = [line.strip() for line in f if line.strip()]
        prepare_kegg_nonmodel(
            kaas_file=kegg_nm["kaas"],
            out_dir=out_dir,
            cache_dir=cache_dir,
            refresh=refresh,
            strip_isoform=kegg_nm.get("strip_isoform", True),
            gene_list=gene_list,
        )

    # GO non-model
    go_nm = modules.get("go_nonmodel", {})
    if go_nm.get("enabled", False):
        from src.prepare_go import prepare_go
        log_msg("=== Module: go_nonmodel ===")
        prepare_go(
            go_table_file=go_nm["go_table"],
            out_dir=out_dir,
            cache_dir=cache_dir,
            refresh=refresh,
            expand=expand_go,
        )

    # Ensembl
    ens = modules.get("ensembl", {})
    if ens.get("enabled", False):
        from src.prepare_ensembl import prepare_ensembl
        log_msg("=== Module: ensembl ===")
        prepare_ensembl(
            dataset=ens["dataset"],
            out_dir=out_dir,
            cache_dir=cache_dir,
            mart=ens.get("mart", "ensembl"),
            host=ens.get("host"),
            kegg_org=ens.get("kegg_org"),
            id_source=ens.get("id_source", "ensembl"),
            expand=expand_go,
            refresh=refresh,
        )

    # UniProt
    uni = modules.get("uniprot", {})
    if uni.get("enabled", False):
        from src.prepare_uniprot import prepare_uniprot
        log_msg("=== Module: uniprot ===")
        prepare_uniprot(
            taxon_id=uni["taxon_id"],
            out_dir=out_dir,
            cache_dir=cache_dir,
            reviewed=uni.get("reviewed", True),
            kegg_org=uni.get("kegg_org"),
            expand=expand_go,
            refresh=refresh,
        )

    log_msg("=== All modules done ===")


if __name__ == "__main__":
    main()
