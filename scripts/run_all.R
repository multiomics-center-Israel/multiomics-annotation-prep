#!/usr/bin/env Rscript
# Config-driven driver: runs whichever modules are enabled in config.yml.
# Usage:  Rscript scripts/run_all.R --config config/config.yml

suppressWarnings(suppressMessages({ library(optparse); library(yaml) }))

.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(.this), "..", "R", "bootstrap.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config", type = "character", default = "config/config.yml",
              help = "path to config.yml [%default]")
)))
cfg <- yaml::read_yaml(opt$config)

out_dir   <- cfg$out_dir   %||% "results"
cache_dir <- cfg$cache_dir %||% "data"
refresh   <- isTRUE(cfg$refresh)
expand_go <- if (is.null(cfg$expand_go)) TRUE else isTRUE(cfg$expand_go)
m         <- cfg$modules

enabled <- function(x) !is.null(x) && isTRUE(x$enabled)

if (enabled(m$kegg_nonmodel)) {
  log_msg("== module: kegg_nonmodel ==")
  prepare_kegg_nonmodel(
    kaas_file     = m$kegg_nonmodel$kaas,
    out_dir       = out_dir, cache_dir = cache_dir, refresh = refresh,
    strip_isoform = if (is.null(m$kegg_nonmodel$strip_isoform)) TRUE
                    else isTRUE(m$kegg_nonmodel$strip_isoform))
}

if (enabled(m$go_nonmodel)) {
  log_msg("== module: go_nonmodel ==")
  prepare_go(go_table = m$go_nonmodel$go_table,
             out_dir = out_dir, cache_dir = cache_dir,
             refresh = refresh, expand = expand_go)
}

if (enabled(m$ensembl)) {
  log_msg("== module: ensembl ==")
  prepare_ensembl(
    dataset = m$ensembl$dataset, out_dir = out_dir, cache_dir = cache_dir,
    mart = m$ensembl$mart %||% "ensembl", host = m$ensembl$host,
    kegg_org = m$ensembl$kegg_org, id_source = m$ensembl$id_source %||% "ensembl",
    expand = expand_go, refresh = refresh)
}

if (enabled(m$uniprot)) {
  log_msg("== module: uniprot ==")
  prepare_uniprot(
    taxon_id = m$uniprot$taxon_id, out_dir = out_dir, cache_dir = cache_dir,
    reviewed = if (is.null(m$uniprot$reviewed)) TRUE else isTRUE(m$uniprot$reviewed),
    kegg_org = m$uniprot$kegg_org, expand = expand_go, refresh = refresh)
}

log_msg("run_all complete -> ", out_dir)
