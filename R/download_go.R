# ---------------------------------------------------------------------------
# download_go.R  -  obtain GO term names and namespaces
# ---------------------------------------------------------------------------
# Preferred source is the Bioconductor GO.db package (versioned, offline).
# Fallback is the official go-basic.obo, downloaded and cached.
# ---------------------------------------------------------------------------

GO_OBO_URL <- "http://purl.obolibrary.org/obo/go/go-basic.obo"

NS_MAP <- c(biological_process = "BP",
            molecular_function = "MF",
            cellular_component = "CC")

# Returns data.frame(go, name, namespace) where namespace in {BP, MF, CC}.
go_term_table <- function(cache_dir, refresh = FALSE, prefer_godb = TRUE) {
  if (prefer_godb && requireNamespace("GO.db", quietly = TRUE) &&
      requireNamespace("AnnotationDbi", quietly = TRUE)) {
    log_msg("GO terms from GO.db (version ", as.character(utils::packageVersion("GO.db")), ")")
    ids  <- AnnotationDbi::keys(GO.db::GO.db)
    tab  <- AnnotationDbi::select(GO.db::GO.db, keys = ids,
                                  columns = c("TERM", "ONTOLOGY"), keytype = "GOID")
    data.frame(go = tab$GOID, name = tab$TERM,
               namespace = tab$ONTOLOGY, stringsAsFactors = FALSE)
  } else {
    obo <- cached_download(GO_OBO_URL, "go-basic.obo", cache_dir, refresh)
    parse_obo_terms(obo)
  }
}

# Minimal OBO parser: id / name / namespace per [Term] stanza.
parse_obo_terms <- function(path) {
  lines <- readLines(path, warn = FALSE)
  go <- name <- ns <- character(0)
  cur_id <- cur_name <- cur_ns <- NA_character_
  in_term <- FALSE
  flush <- function() {
    if (in_term && !is.na(cur_id)) {
      go   <<- c(go, cur_id)
      name <<- c(name, if (is.na(cur_name)) "" else cur_name)
      ns   <<- c(ns, if (is.na(cur_ns)) "" else cur_ns)
    }
  }
  for (ln in lines) {
    if (ln == "[Term]") {
      flush(); in_term <- TRUE
      cur_id <- cur_name <- cur_ns <- NA_character_
    } else if (grepl("^\\[", ln)) {
      flush(); in_term <- FALSE
    } else if (in_term) {
      if (startsWith(ln, "id: "))            cur_id   <- sub("^id: ", "", ln)
      else if (startsWith(ln, "name: "))     cur_name <- sub("^name: ", "", ln)
      else if (startsWith(ln, "namespace: "))cur_ns   <- sub("^namespace: ", "", ln)
    }
  }
  flush()
  ns2 <- unname(NS_MAP[ns]); ns2[is.na(ns2)] <- ""
  data.frame(go = go, name = name, namespace = ns2, stringsAsFactors = FALSE)
}
