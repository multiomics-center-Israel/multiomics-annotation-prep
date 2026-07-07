# ---------------------------------------------------------------------------
# download_kegg_org.R  -  KEGG for a specific organism code (model organisms)
# ---------------------------------------------------------------------------
# For organisms that HAVE a KEGG organism code (e.g. hsa, mmu, dme) we can pull
# gene -> pathway directly, instead of going through KAAS/KO. Optionally convert
# KEGG gene IDs to NCBI GeneID or UniProt so they match pipeline gene IDs.
# ---------------------------------------------------------------------------

# Produce KEGG_pathway2gene.tab + KEGG_pathway2name.tab for a KEGG organism.
#   kegg_org  : KEGG organism code, e.g. "mmu"
#   id_source : "kegg" (default), "ncbi-geneid", or "uniprot"
# Returns a named list of output paths.
prepare_kegg_by_org <- function(kegg_org,
                                out_dir   = "results",
                                cache_dir = "data",
                                refresh   = FALSE,
                                id_source = "kegg") {
  ensure_dir(out_dir)

  ## gene -> pathway  (org-specific pathway ids, e.g. mmu00010) --------------
  link_url  <- sprintf("https://rest.kegg.jp/link/pathway/%s", kegg_org)
  link_file <- cached_download(link_url,
                               sprintf("kegg_%s_gene2path.txt", kegg_org),
                               cache_dir, refresh)
  gp <- utils::read.delim(link_file, header = FALSE, sep = "\t",
                          stringsAsFactors = FALSE, quote = "",
                          col.names = c("gene", "path"))
  gp$path <- sub("^path:", "", gp$path)

  ## optional id conversion --------------------------------------------------
  if (id_source != "kegg") {
    conv_url  <- sprintf("https://rest.kegg.jp/conv/%s/%s", id_source, kegg_org)
    conv_file <- cached_download(conv_url,
                                 sprintf("kegg_%s_to_%s.txt", kegg_org, id_source),
                                 cache_dir, refresh)
    cv <- utils::read.delim(conv_file, header = FALSE, sep = "\t",
                            stringsAsFactors = FALSE, quote = "",
                            col.names = c("kegg", "ext"))
    cv$ext <- sub("^[^:]+:", "", cv$ext)                # strip prefix
    map <- stats::setNames(cv$ext, cv$kegg)
    gp$gene <- unname(map[gp$gene])
    gp <- gp[!is.na(gp$gene), , drop = FALSE]
  } else {
    gp$gene <- sub(sprintf("^%s:", kegg_org), "", gp$gene)
  }

  write_pathway2gene(list(gp$path, gp$gene),
                     file.path(out_dir, "KEGG_pathway2gene.tab"))

  ## pathway -> name ---------------------------------------------------------
  pn_url  <- sprintf("https://rest.kegg.jp/list/pathway/%s", kegg_org)
  pn_file <- cached_download(pn_url,
                             sprintf("kegg_%s_pathway_names.txt", kegg_org),
                             cache_dir, refresh)
  pn <- utils::read.delim(pn_file, header = FALSE, sep = "\t",
                          stringsAsFactors = FALSE, quote = "",
                          col.names = c("pathway", "info"))
  pn$pathway <- sub("^path:", "", pn$pathway)
  write_pathway2name(list(pn$pathway, pn$info),
                     file.path(out_dir, "KEGG_pathway2name.tab"))

  invisible(list(
    pathway2gene = file.path(out_dir, "KEGG_pathway2gene.tab"),
    pathway2name = file.path(out_dir, "KEGG_pathway2name.tab")
  ))
}
