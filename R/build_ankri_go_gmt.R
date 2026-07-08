#!/usr/bin/env Rscript
#' Build a GO gene-set (GMT) for the Serge Ankri E. histolytica proteomics run
#'
#' The NCBI GO annotation for this assembly (GCF_000208925.1) is a GAF that keys
#' every annotation on an NCBI GeneID, while the DIA-NN pg_matrix (and therefore
#' the DE FeatureID) keys proteins on RefSeq XP_ accessions. The two never match
#' directly, so this builder bridges them through the genome GTF, whose CDS rows
#' carry both a GeneID db_xref and the protein_id (XP_) accession:
#'
#'   GAF  (GeneID -> GO)  x  GTF (GeneID -> XP_)  ->  (XP_ -> GO)
#'
#' It then re-keys members onto the pg_matrix Protein.Group strings (some are
#' semicolon-joined accession groups) so the GMT members equal the DE FeatureID,
#' and writes two files the proteomics pipeline can use for GO enrichment:
#'
#'   1. protein_to_GO.tsv    long mapping, one (Protein.Group, GO_id) pair per row
#'   2. GO_ehistolytica.gmt  GMT: GO_id <tab> "name [aspect]" <tab> members...
#'
#' Single host genome here, so (unlike the Elad combined-genome builder) members
#' carry no species suffix -- they are the bare Protein.Group strings.
#'
#' Usage (from the repo root, inside the pipeline env):
#'   conda activate mofa2_env
#'   Rscript utils/build_ankri_go_gmt.R

suppressWarnings(suppressMessages({
    library(GO.db)
    library(AnnotationDbi)
}))

`%||%` <- function(x, y) if (is.null(x)) y else x

# ------------------------------------------------------------------------------
# Source the vendored GMT helpers (read_gmt/write_gmt/validate_gmt/
# filter_gmt_by_size) from gmt_utils.R next to this script.
# ------------------------------------------------------------------------------
script_dir <- tryCatch(
    dirname(sys.frame(1)$ofile),
    error = function(e) {
        args <- commandArgs(trailingOnly = FALSE)
        file_arg <- grep("^--file=", args, value = TRUE)
        if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg))) else "."
    }
)
source(file.path(script_dir, "gmt_utils.R"))

# ------------------------------------------------------------------------------
# Inputs / outputs
# ------------------------------------------------------------------------------
# The project's GAF/GTF/pg_matrix are not stored in this repo. Point
# ANNPREP_DATA_ROOT at the analysis project's data/ dir (e.g. multiomics-core's
# data/), or drop the files under this repo's own data/ and leave it unset.
data_root <- Sys.getenv("ANNPREP_DATA_ROOT", file.path(dirname(script_dir), "data"))
data_dir  <- file.path(data_root, "Serge_Ankri_June2026", "proteomics")
gaf_file  <- file.path(data_dir, "GCF_000208925.1_JCVI_ESG2_1.0_gene_ontology.gaf")
gtf_file  <- file.path(data_dir, "GCF_000208925.1_JCVI_ESG2_1.0_genomic.gtf")
pg_matrix <- file.path(data_dir, "98858-81_hystolytica_MBR.pg_matrix.tsv")
out_map   <- file.path(data_dir, "protein_to_GO.tsv")
out_gmt   <- file.path(data_dir, "GO_ehistolytica.gmt")

# Runtime enrichment size window, mirrored here only for the validation report
# and for trimming the written GMT. Keep in sync with config pathway.min/max_size.
MIN_SIZE <- 10
MAX_SIZE <- 500

# ------------------------------------------------------------------------------
# 1. GeneID -> XP_ accession map, from the genome GTF CDS rows
# ------------------------------------------------------------------------------

#' Map NCBI GeneIDs to RefSeq protein accessions from a genome GTF.
#'
#' Reads the CDS rows and pairs each `db_xref "GeneID:<id>"` with the row's
#' `protein_id "<XP_...>"`. A GeneID may yield several protein isoforms.
#'
#' @param path Path to the NCBI RefSeq genomic GTF.
#' @return A data.frame with columns GeneID (character) and accession (XP_...),
#'   one row per distinct pair.
parse_gtf_gene_to_protein <- function(path) {
    raw <- readLines(path, warn = FALSE)
    raw <- raw[grepl('protein_id "', raw, fixed = TRUE) &
               grepl('GeneID:', raw, fixed = TRUE)]
    if (!length(raw)) stop(sprintf("No CDS rows with GeneID+protein_id in %s", path))
    gene <- sub('.*GeneID:([0-9]+).*', '\\1', raw)
    prot <- sub('.*protein_id "([^"]+)".*', '\\1', raw)
    unique(data.frame(GeneID = gene, accession = prot, stringsAsFactors = FALSE))
}

# ------------------------------------------------------------------------------
# 2. GeneID -> GO_id pairs, from the GAF
# ------------------------------------------------------------------------------

#' Extract (GeneID, GO_id) pairs from an NCBI GAF 2.x file.
#'
#' Uses the DB Object ID (column 2, an NCBI GeneID here) and GO ID (column 5);
#' skips "!" comment lines and NOT-qualified (negated) annotations (column 4).
#'
#' @param path Path to the GAF file.
#' @return A long data.frame (GeneID, GO_id), or NULL if empty.
parse_gaf_gene_go <- function(path) {
    raw <- readLines(path, warn = FALSE)
    raw <- raw[!startsWith(raw, "!") & nzchar(raw)]
    if (!length(raw)) return(NULL)
    parts <- strsplit(raw, "\t", fixed = TRUE)
    gene <- vapply(parts, function(p) if (length(p) >= 2) p[2] else NA_character_, "")
    qual <- vapply(parts, function(p) if (length(p) >= 4) p[4] else "",           "")
    go   <- vapply(parts, function(p) if (length(p) >= 5) p[5] else NA_character_, "")
    keep <- !grepl("NOT", qual, fixed = TRUE) & !is.na(gene) & nzchar(gene) &
            !is.na(go) & nzchar(go)
    if (!any(keep)) return(NULL)
    unique(data.frame(GeneID = gene[keep], GO_id = go[keep], stringsAsFactors = FALSE))
}

# ------------------------------------------------------------------------------
# 3. Accession -> Protein.Group map, from the measured pg_matrix
# ------------------------------------------------------------------------------

#' Map each protein accession to the pg_matrix Protein.Group that contains it.
#'
#' DIA-NN reports razor/shared proteins as semicolon-joined groups (e.g.
#' "XP_a;XP_b"); the DE FeatureID equals the full group string, so members must
#' be re-keyed onto it. Restricting to these accessions also restricts the GMT
#' to the measured universe.
#'
#' @param path Path to the DIA-NN pg_matrix.tsv (first column = Protein.Group).
#' @return A data.frame with columns accession and group (the Protein.Group).
parse_pg_accession_to_group <- function(path) {
    groups <- read.delim(path, header = TRUE, quote = "", comment.char = "",
                         check.names = FALSE)[[1]]
    groups <- unique(as.character(groups))
    acc_list <- strsplit(groups, ";", fixed = TRUE)
    data.frame(
        accession = trimws(unlist(acc_list, use.names = FALSE)),
        group     = rep(groups, lengths(acc_list)),
        stringsAsFactors = FALSE
    )
}

# ------------------------------------------------------------------------------
# Build the long Protein.Group -> GO mapping
# ------------------------------------------------------------------------------
message("== Building GeneID -> protein and GeneID -> GO tables ==")
g2p <- parse_gtf_gene_to_protein(gtf_file)
message(sprintf("  GTF: %d GeneID<->accession pairs (%d GeneIDs, %d accessions)",
                nrow(g2p), length(unique(g2p$GeneID)), length(unique(g2p$accession))))

g2go <- parse_gaf_gene_go(gaf_file)
if (is.null(g2go)) stop("No usable annotations parsed from the GAF.")
message(sprintf("  GAF: %d GeneID<->GO pairs (%d GeneIDs with GO)",
                nrow(g2go), length(unique(g2go$GeneID))))

# GeneID -> GO joined to GeneID -> accession gives accession -> GO.
acc_go <- merge(g2go, g2p, by = "GeneID")[, c("accession", "GO_id")]
message(sprintf("== Bridged to accessions: %d accession<->GO pairs (%d accessions) ==",
                nrow(unique(acc_go)), length(unique(acc_go$accession))))

# Re-key onto the measured Protein.Group strings (also restricts to the universe).
acc2grp <- parse_pg_accession_to_group(pg_matrix)
feature_ids <- unique(acc2grp$group)
long <- merge(unique(acc_go), acc2grp, by = "accession")
long <- unique(data.frame(protein_id = long$group, GO_id = long$GO_id,
                          stringsAsFactors = FALSE))
message(sprintf("  Re-keyed to Protein.Group: %d pairs (%d of %d groups matched)",
                nrow(long), length(unique(long$protein_id)), length(feature_ids)))

long <- long[order(long$protein_id, long$GO_id), ]
write.table(long, out_map, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("  Wrote %s (%d unique pairs)", basename(out_map), nrow(long)))

# ------------------------------------------------------------------------------
# Invert to gene sets and name the terms via GO.db
# ------------------------------------------------------------------------------
message("== Building GO gene sets ==")
gene_sets <- lapply(split(long$protein_id, long$GO_id), unique)
go_ids <- names(gene_sets)

go_keys <- keys(GO.db, "GOID")
known <- go_ids %in% go_keys
ann <- suppressMessages(AnnotationDbi::select(
    GO.db, keys = go_ids[known], columns = c("TERM", "ONTOLOGY"), keytype = "GOID"))
term_map <- setNames(ann$TERM, ann$GOID)
ont_map  <- setNames(ann$ONTOLOGY, ann$GOID)

descriptions <- ifelse(
    go_ids %in% names(term_map) & !is.na(term_map[go_ids]),
    sprintf("%s [%s]", term_map[go_ids], ont_map[go_ids]),
    go_ids  # fall back to the bare GO id for obsolete / secondary terms
)
names(descriptions) <- go_ids
attr(gene_sets, "descriptions") <- descriptions

n_unnamed <- sum(!(go_ids %in% names(term_map)) | is.na(term_map[go_ids]))
message(sprintf("  %d GO terms (%d without a current GO.db name -> kept as GO id)",
                length(gene_sets), n_unnamed))
sizes <- lengths(gene_sets)
message(sprintf("  Set sizes: min %d, median %d, max %d; %d within the %d-%d window",
                min(sizes), as.integer(median(sizes)), max(sizes),
                sum(sizes >= MIN_SIZE & sizes <= MAX_SIZE), MIN_SIZE, MAX_SIZE))

# ------------------------------------------------------------------------------
# Trim to the enrichment size window and write the GMT
# ------------------------------------------------------------------------------
# Members are already restricted to the measured universe, so these sizes are
# final: sets outside [MIN_SIZE, MAX_SIZE] are exactly the ones the pipeline would
# drop at runtime (fgsea min/maxSize; ORA gene-set size filter). NOTE: if the
# config's pathway.min_size/max_size change, re-run this builder with matching
# bounds.
kept <- filter_gmt_by_size(gene_sets, min_size = MIN_SIZE, max_size = MAX_SIZE)
message(sprintf("  %d sets -> %d after size %d-%d (dropped %d); writing %s",
                length(gene_sets), length(kept), MIN_SIZE, MAX_SIZE,
                length(gene_sets) - length(kept), basename(out_gmt)))
write_gmt(kept, out_gmt, verbose = TRUE)

# ------------------------------------------------------------------------------
# Coverage report against the actual DIA-NN features + round-trip check
# ------------------------------------------------------------------------------
message("== Validating against pg_matrix features ==")
invisible(validate_gmt(gene_sets, feature_ids,
                       min_coverage = 0.05,
                       min_pathway_size = MIN_SIZE,
                       max_pathway_size = MAX_SIZE,
                       verbose = TRUE))

rt <- read_gmt(out_gmt)
message(sprintf("Round-trip: read_gmt() loaded %d sets from %s.",
                length(rt), basename(out_gmt)))

message("Done.")
