# ---------------------------------------------------------------------------
# prepare_kegg_nonmodel.R
# ---------------------------------------------------------------------------
# Non-model organism KEGG annotation, driven by a KAAS result file.
# R re-implementation of Retrieve_KEGG_annot_nonmodel_org.pl, with the KEGG
# reference tables downloaded automatically (see download_kegg.R).
#
# Inputs
#   kaas_file : KAAS output, 2 columns  <transcript/gene id>\t<K number>
#               (rows without a K number are ignored)
# Outputs (in out_dir)
#   KEGG_pathway2gene.tab   (v1, index)     -> enrichment
#   KEGG_pathway2name.tab   (pathway, info) -> enrichment
#   KEGG_annot_genes.txt                    -> descriptive annotation
# ---------------------------------------------------------------------------

prepare_kegg_nonmodel <- function(kaas_file,
                                   out_dir       = "results",
                                   cache_dir     = "data",
                                   refresh       = FALSE,
                                   strip_isoform = TRUE,
                                   gene_list     = NULL) {

  stopifnot(file.exists(kaas_file))
  ensure_dir(out_dir)

  ## 1. gene -> {KO} from the KAAS result -----------------------------------
  log_msg("reading KAAS result: ", kaas_file)
  kaas <- utils::read.delim(kaas_file, header = FALSE, sep = "\t",
                            stringsAsFactors = FALSE, quote = "",
                            col.names = c("contig", "ko"), fill = TRUE)
  kaas <- kaas[!is.na(kaas$ko) & nzchar(kaas$ko), , drop = FALSE]
  if (strip_isoform)
    kaas$contig <- sub("_i\\d+$", "", kaas$contig)   # collapse isoforms to gene
  kaas$ko <- sub("^ko:", "", kaas$ko)
  kaas <- unique(kaas)

  gene2ko <- split(kaas$ko, kaas$contig)             # named list gene -> KOs
  log_msg("genes with >=1 KO: ", length(gene2ko))

  ## 2. KEGG reference tables (auto-download + parse) ------------------------
  kegg  <- download_kegg_rest(cache_dir, refresh = refresh)
  ko2info <- parse_ko_to_name(kegg$ko_to_name)
  ko2path <- parse_ko_to_path(kegg$ko_to_path)
  path2nm <- parse_pathway_names(kegg$pathway_names)

  ## 3. enrichment file: pathway -> name ------------------------------------
  write_pathway2name(list(names(path2nm), unname(path2nm)),
                     file.path(out_dir, "KEGG_pathway2name.tab"))

  ## 4. enrichment file: pathway -> gene ------------------------------------
  p2g_path <- character(0)
  p2g_gene <- character(0)
  for (gene in names(gene2ko)) {
    paths <- unique(unlist(ko2path[gene2ko[[gene]]], use.names = FALSE))
    paths <- paths[!is.na(paths)]
    if (length(paths)) {
      p2g_path <- c(p2g_path, paths)
      p2g_gene <- c(p2g_gene, rep(gene, length(paths)))
    }
  }
  write_pathway2gene(list(p2g_path, p2g_gene),
                     file.path(out_dir, "KEGG_pathway2gene.tab"))

  ## 5. descriptive annotation ----------------------------------------------
  genes <- if (is.null(gene_list)) names(gene2ko) else gene_list
  annot <- build_kegg_annot(genes, gene2ko, ko2info, ko2path, path2nm)
  annot_file <- file.path(out_dir, "KEGG_annot_genes.txt")
  write_tab(annot, annot_file)

  invisible(list(
    pathway2gene = file.path(out_dir, "KEGG_pathway2gene.tab"),
    pathway2name = file.path(out_dir, "KEGG_pathway2name.tab"),
    annotation   = annot_file
  ))
}

# per-gene descriptive table (values joined with " | " as in the Perl version)
build_kegg_annot <- function(genes, gene2ko, ko2info, ko2path, path2nm) {
  j <- function(x) paste(x, collapse = " | ")
  rows <- lapply(genes, function(gene) {
    kos <- gene2ko[[gene]]
    if (is.null(kos)) kos <- character(0)
    names  <- vapply(kos, function(k) if (!is.null(ko2info[[k]])) ko2info[[k]]$names else "", "")
    titles <- vapply(kos, function(k) if (!is.null(ko2info[[k]])) ko2info[[k]]$title else "", "")
    ecs    <- vapply(kos, function(k) if (!is.null(ko2info[[k]])) ko2info[[k]]$ec    else "", "")
    paths  <- unique(unlist(ko2path[kos], use.names = FALSE))
    paths  <- paths[!is.na(paths)]
    pnames <- unname(path2nm[paths]); pnames[is.na(pnames)] <- ""
    data.frame(
      Gene              = gene,
      KEGG_ID           = j(kos),
      KEGG_names        = j(names),
      KEGG_description  = j(titles),
      EC_number         = j(ecs),
      Pathway_IDs       = j(paths),
      Pathway_names     = j(pnames),
      stringsAsFactors  = FALSE
    )
  })
  do.call(rbind, rows)
}
