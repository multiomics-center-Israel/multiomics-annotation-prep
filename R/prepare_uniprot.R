# ---------------------------------------------------------------------------
# prepare_uniprot.R  -  proteomes available in UniProt (proteomics projects)
# ---------------------------------------------------------------------------
# Downloads a per-protein annotation table from the UniProt REST API and builds
# descriptive + GO enrichment files keyed by UniProt accession. Optionally adds
# KEGG pathway files (mapped to UniProt accessions) via download_kegg_org.R.
# ---------------------------------------------------------------------------

# taxon_id  : NCBI taxonomy id, e.g. 10090 (mouse)
# reviewed  : TRUE = Swiss-Prot only (recommended); FALSE = include TrEMBL
# kegg_org  : optional KEGG organism code (e.g. "mmu") to also build KEGG files
prepare_uniprot <- function(taxon_id,
                            out_dir   = "results",
                            cache_dir = "data",
                            reviewed  = TRUE,
                            kegg_org  = NULL,
                            expand    = TRUE,
                            refresh   = FALSE) {

  ensure_dir(out_dir)

  ## 1. download the UniProt table ------------------------------------------
  query <- sprintf("organism_id:%s%s", taxon_id,
                   if (reviewed) "+AND+reviewed:true" else "")
  fields <- "accession,gene_names,protein_name,go_id"
  url <- sprintf(
    "https://rest.uniprot.org/uniprotkb/stream?query=%s&format=tsv&fields=%s",
    query, fields)
  up_file <- cached_download(url,
                             sprintf("uniprot_%s.tsv", taxon_id),
                             cache_dir, refresh)
  up <- utils::read.delim(up_file, header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE, quote = "",
                          check.names = FALSE)
  # normalise column names (UniProt returns friendly headers)
  names(up) <- c("Gene", "gene_names", "protein_name", "go_ids")[seq_len(ncol(up))]
  log_msg("UniProt entries: ", nrow(up))

  ## 2. descriptive annotation ----------------------------------------------
  desc <- data.frame(
    Gene         = up$Gene,
    gene_name    = up$gene_names,
    protein_name = up$protein_name,
    stringsAsFactors = FALSE)
  write_tab(unique(desc), file.path(out_dir, "Annotation.tab"))

  ## 3. GO enrichment files --------------------------------------------------
  gene2go <- do.call(rbind, lapply(seq_len(nrow(up)), function(i) {
    gos <- unique(regmatches(up$go_ids[i],
                             gregexpr("GO:[0-9]{7}", up$go_ids[i]))[[1]])
    if (!length(gos)) return(NULL)
    data.frame(Gene = up$Gene[i], GO = gos, stringsAsFactors = FALSE)
  }))
  gene2go <- unique(gene2go)
  log_msg("UniProt direct gene-GO pairs: ", nrow(gene2go))

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

  ## 4. KEGG (optional), mapped to UniProt accessions ------------------------
  if (!is.null(kegg_org))
    prepare_kegg_by_org(kegg_org, out_dir, cache_dir,
                        refresh = refresh, id_source = "uniprot")

  invisible(out_dir)
}
