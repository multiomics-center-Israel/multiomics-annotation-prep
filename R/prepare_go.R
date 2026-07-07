# ---------------------------------------------------------------------------
# prepare_go.R
# ---------------------------------------------------------------------------
# Build GO enrichment files from a per-gene GO table (e.g. Trinotate output).
# Optionally expands direct GO annotations to parental terms so enrichment can
# use every level of the hierarchy (as in Neat_Annotation's Expand_GO.Rmd).
#
# Input
#   go_table : two columns  <gene id>\t<GO ids>
#              GO ids can be separated by comma / semicolon / backtick and may
#              carry Trinotate suffixes (GO:0005634^cellular_component^nucleus);
#              only the GO:####### tokens are extracted.
# Outputs (in out_dir)
#   GO2gene_BP.tab / _MF.tab / _CC.tab   (GO, Gene)   -> enrichment
#   GO2name_BP.tab / _MF.tab / _CC.tab   (GO, Term)   -> enrichment
# ---------------------------------------------------------------------------

prepare_go <- function(go_table,
                       out_dir   = "results",
                       cache_dir = "data",
                       refresh   = FALSE,
                       expand    = TRUE) {

  stopifnot(file.exists(go_table))
  ensure_dir(out_dir)

  ## 1. parse gene -> GO ids -------------------------------------------------
  log_msg("reading GO table: ", go_table)
  raw <- utils::read.delim(go_table, header = FALSE, sep = "\t",
                           stringsAsFactors = FALSE, quote = "",
                           col.names = c("gene", "go"), fill = TRUE)
  gene2go <- do.call(rbind, lapply(seq_len(nrow(raw)), function(i) {
    gos <- unique(regmatches(raw$go[i],
                             gregexpr("GO:[0-9]{7}", raw$go[i]))[[1]])
    if (!length(gos)) return(NULL)
    data.frame(Gene = raw$gene[i], GO = gos, stringsAsFactors = FALSE)
  }))
  gene2go <- unique(gene2go)
  log_msg("direct gene-GO pairs: ", nrow(gene2go))

  ## 2. optional expansion to ancestral terms --------------------------------
  if (expand) gene2go <- expand_go(gene2go)

  ## 3. term names + namespaces ---------------------------------------------
  terms <- go_term_table(cache_dir, refresh = refresh)
  gene2go <- merge(gene2go, terms, by.x = "GO", by.y = "go", all.x = TRUE)
  gene2go$namespace[is.na(gene2go$namespace)] <- ""

  ## 4. write one pair of files per namespace -------------------------------
  for (ns in c("BP", "MF", "CC")) {
    sub <- gene2go[gene2go$namespace == ns, , drop = FALSE]
    if (!nrow(sub)) { log_msg("no ", ns, " terms, skipping"); next }
    write_go2gene(list(sub$GO, sub$Gene),
                  file.path(out_dir, paste0("GO2gene_", ns, ".tab")))
    nm <- unique(sub[, c("GO", "name")])
    write_go2name(list(nm$GO, nm$name),
                  file.path(out_dir, paste0("GO2name_", ns, ".tab")))
  }
  invisible(out_dir)
}

# Expand gene->GO to include all ancestor terms.
# Uses GO.db ancestor maps (deterministic, offline) by default; falls back to
# clusterProfiler::buildGOmap if GO.db is unavailable.
expand_go <- function(gene2go, method = c("go.db", "clusterprofiler")) {
  method <- match.arg(method)
  if (method == "go.db" &&
      requireNamespace("GO.db", quietly = TRUE) &&
      requireNamespace("AnnotationDbi", quietly = TRUE)) {
    return(expand_go_godb(gene2go))
  }
  if (requireNamespace("clusterProfiler", quietly = TRUE)) {
    log_msg("expanding GO with clusterProfiler::buildGOmap ...")
    gomap <- data.frame(GO = gene2go$GO, Gene = gene2go$Gene,
                        stringsAsFactors = FALSE)
    expanded <- tryCatch(clusterProfiler::buildGOmap(gomap),
                         error = function(e) {
                           log_msg("buildGOmap failed (", conditionMessage(e),
                                   "), falling back to GO.db")
                           NULL
                         })
    if (!is.null(expanded)) {
      # be defensive about the returned column names/positions
      cn  <- names(expanded)
      gcol <- if ("GO" %in% cn) "GO" else grep("GO", cn, value = TRUE)[1]
      ecol <- setdiff(cn, gcol)[1]
      out  <- unique(data.frame(Gene = expanded[[ecol]], GO = expanded[[gcol]],
                                stringsAsFactors = FALSE))
      log_msg("expanded gene-GO pairs: ", nrow(out))
      return(out)
    }
  }
  expand_go_godb(gene2go)
}

expand_go_godb <- function(gene2go) {
  if (!requireNamespace("GO.db", quietly = TRUE) ||
      !requireNamespace("AnnotationDbi", quietly = TRUE)) {
    log_msg("GO.db unavailable; skipping expansion (using direct terms only)")
    return(gene2go)
  }
  log_msg("expanding GO with GO.db ancestor maps ...")
  anc <- c(as.list(GO.db::GOBPANCESTOR),
           as.list(GO.db::GOMFANCESTOR),
           as.list(GO.db::GOCCANCESTOR))
  rows <- lapply(seq_len(nrow(gene2go)), function(i) {
    g  <- gene2go$Gene[i]; go <- gene2go$GO[i]
    a  <- anc[[go]]; a <- a[!is.na(a) & a != "all"]
    data.frame(Gene = g, GO = c(go, a), stringsAsFactors = FALSE)
  })
  out <- unique(do.call(rbind, rows))
  log_msg("expanded gene-GO pairs: ", nrow(out))
  out
}
