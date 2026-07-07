# ---------------------------------------------------------------------------
# prepare_ensembl.R  -  model organisms with a genome in Ensembl (via BioMart)
# ---------------------------------------------------------------------------
# Fetches gene descriptions + GO annotations from Ensembl/BioMart, and (if a
# KEGG organism code is supplied) KEGG pathway files via download_kegg_org.R.
# Output formats are identical to the non-model modules.
# ---------------------------------------------------------------------------

# dataset  : biomaRt dataset, e.g. "mmusculus_gene_ensembl"
# mart      : "genes" (Ensembl) or an Ensembl Genomes mart name
# host      : optional host override (e.g. an archive: "https://may2024.archive.ensembl.org")
# kegg_org  : optional KEGG organism code (e.g. "mmu") to also build KEGG files
prepare_ensembl <- function(dataset,
                            out_dir   = "results",
                            cache_dir = "data",
                            mart      = "ensembl",
                            host      = NULL,
                            kegg_org  = NULL,
                            id_source = "ensembl",   # for the KEGG id conversion
                            expand    = TRUE,
                            refresh   = FALSE) {

  if (!requireNamespace("biomaRt", quietly = TRUE))
    stop("Package 'biomaRt' is required for the Ensembl module. See R/install_deps.R")
  ensure_dir(out_dir)

  ## connect -----------------------------------------------------------------
  log_msg("connecting to BioMart: ", mart, " / ", dataset)
  ensembl <- if (is.null(host))
    biomaRt::useEnsembl(biomart = mart, dataset = dataset)
  else
    biomaRt::useMart(biomart = mart, dataset = dataset, host = host)

  ## 1. descriptive annotation ----------------------------------------------
  desc <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "external_gene_name",
                   "gene_biotype", "description"),
    mart = ensembl)
  names(desc) <- c("Gene", "gene_name", "gene_biotype", "description")
  write_tab(unique(desc), file.path(out_dir, "Annotation.tab"))

  ## 2. GO enrichment files --------------------------------------------------
  go <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "go_id", "namespace_1003"),
    mart = ensembl)
  go <- go[nzchar(go$go_id), , drop = FALSE]
  gene2go <- data.frame(Gene = go$ensembl_gene_id, GO = go$go_id,
                        stringsAsFactors = FALSE)
  gene2go <- unique(gene2go)
  log_msg("Ensembl direct gene-GO pairs: ", nrow(gene2go))

  if (expand) gene2go <- expand_go(gene2go)
  terms <- go_term_table(cache_dir, refresh = refresh)
  gene2go <- merge(gene2go, terms, by.x = "GO", by.y = "go", all.x = TRUE)
  gene2go$namespace[is.na(gene2go$namespace)] <- ""
  for (ns in c("BP", "MF", "CC")) {
    sub <- gene2go[gene2go$namespace == ns, , drop = FALSE]
    if (!nrow(sub)) next
    write_go2gene(list(sub$GO, sub$Gene),
                  file.path(out_dir, paste0("GO2gene_", ns, ".tab")))
    nm <- unique(sub[, c("GO", "name")])
    write_go2name(list(nm$GO, nm$name),
                  file.path(out_dir, paste0("GO2name_", ns, ".tab")))
  }

  ## 3. KEGG (optional, if organism has a KEGG code) -------------------------
  if (!is.null(kegg_org))
    prepare_kegg_by_org(kegg_org, out_dir, cache_dir,
                        refresh = refresh, id_source = id_source)

  invisible(out_dir)
}
