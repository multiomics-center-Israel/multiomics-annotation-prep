#!/usr/bin/env Rscript
# Non-model GO: build GO enrichment files from a per-gene GO table (Trinotate),
# auto-downloading term names and expanding the hierarchy.
# Usage:
#   Rscript scripts/run_go.R --go-table examples/trinotate_go.txt \
#           --out results --cache data [--no-expand] [--refresh]

suppressWarnings(suppressMessages(library(optparse)))

.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(.this), "..", "R", "bootstrap.R"))

opt_list <- list(
  make_option("--go-table", type = "character", dest = "go_table",
              help = "per-gene GO table (gene<TAB>GO ids)"),
  make_option("--out",      type = "character", default = "results", help = "output dir [%default]"),
  make_option("--cache",    type = "character", default = "data",    help = "cache dir [%default]"),
  make_option("--no-expand", action = "store_true", default = FALSE, dest = "no_expand",
              help = "do NOT expand GO to parental terms"),
  make_option("--refresh",   action = "store_true", default = FALSE,
              help = "force re-download of GO term files")
)
opt <- parse_args(OptionParser(option_list = opt_list))
if (is.null(opt$go_table)) stop("--go-table is required")

prepare_go(
  go_table  = opt$go_table,
  out_dir   = opt$out,
  cache_dir = opt$cache,
  refresh   = opt$refresh,
  expand    = !opt$no_expand
)
log_msg("GO done -> ", opt$out)
