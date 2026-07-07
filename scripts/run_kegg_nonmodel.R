#!/usr/bin/env Rscript
# Non-model KEGG: auto-download KEGG REST tables + process a KAAS result.
# Usage:
#   Rscript scripts/run_kegg_nonmodel.R --kaas examples/query.ko.txt \
#           --out results --cache data [--refresh] [--no-strip-isoform]

suppressWarnings(suppressMessages(library(optparse)))

# --- locate & load the repo modules ----------------------------------------
.this <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
source(file.path(dirname(.this), "..", "R", "bootstrap.R"))

# --- options ----------------------------------------------------------------
opt_list <- list(
  make_option("--kaas",   type = "character", help = "KAAS query.ko.txt (gene<TAB>KO)"),
  make_option("--out",    type = "character", default = "results", help = "output dir [%default]"),
  make_option("--cache",  type = "character", default = "data",    help = "download cache dir [%default]"),
  make_option("--genes",  type = "character", default = NULL,
              help = "optional gene list file for the descriptive annotation (one per line)"),
  make_option("--refresh", action = "store_true", default = FALSE,
              help = "force re-download of KEGG REST files"),
  make_option("--no-strip-isoform", action = "store_true", default = FALSE,
              dest = "no_strip", help = "keep _iN isoform suffixes (do not collapse to gene)")
)
opt <- parse_args(OptionParser(option_list = opt_list))
if (is.null(opt$kaas)) stop("--kaas is required")

genes <- if (!is.null(opt$genes)) readLines(opt$genes, warn = FALSE) else NULL

prepare_kegg_nonmodel(
  kaas_file     = opt$kaas,
  out_dir       = opt$out,
  cache_dir     = opt$cache,
  refresh       = opt$refresh,
  strip_isoform = !opt$no_strip,
  gene_list     = genes
)
log_msg("KEGG (non-model) done -> ", opt$out)
