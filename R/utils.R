# ---------------------------------------------------------------------------
# utils.R  -  shared helpers: logging, caching, downloads, and output writers
# ---------------------------------------------------------------------------

## ---- misc ------------------------------------------------------------------

# null-coalescing operator (used by the config-driven runner)
`%||%` <- function(a, b) if (is.null(a)) b else a

## ---- logging --------------------------------------------------------------

log_msg <- function(...) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[%s] %s", ts, paste0(..., collapse = "")))
}

## ---- filesystem -----------------------------------------------------------

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

## ---- cached download -------------------------------------------------------
# Downloads `url` to `<cache_dir>/<dest>` unless it already exists (and refresh
# is FALSE). Returns the local path. Used for all remote resource files so a
# re-run is offline and reproducible.

cached_download <- function(url, dest, cache_dir, refresh = FALSE) {
  ensure_dir(cache_dir)
  local <- file.path(cache_dir, dest)
  if (file.exists(local) && !refresh) {
    log_msg("cache hit: ", dest)
    return(local)
  }
  log_msg("downloading: ", url)
  ok <- tryCatch({
    utils::download.file(url, local, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) {
    log_msg("download.file failed (", conditionMessage(e), "), trying httr ...")
    FALSE
  })
  if (!ok || !file.exists(local) || file.info(local)$size == 0) {
    if (!requireNamespace("httr", quietly = TRUE))
      stop("Download failed and 'httr' is not installed: ", url)
    resp <- httr::GET(url, httr::write_disk(local, overwrite = TRUE),
                      httr::timeout(300))
    httr::stop_for_status(resp)
  }
  if (!file.exists(local) || file.info(local)$size == 0)
    stop("Downloaded file is empty: ", url)
  local
}

## ---- output writers -------------------------------------------------------
# Isolated here so the exact output format can be changed in ONE place if the
# downstream pipeline (multiomic-core) expects a different layout.

# tab-delimited, no quotes, no row names — the clusterProfiler .tab convention
write_tab <- function(df, path) {
  ensure_dir(dirname(path))
  utils::write.table(df, path, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = TRUE)
  log_msg("wrote: ", path, "  (", nrow(df), " rows)")
  invisible(path)
}

# KEGG_pathway2gene.tab  (columns literally named v1 / index)
write_pathway2gene <- function(path2gene, out_file) {
  df <- data.frame(v1 = path2gene[[1]], index = path2gene[[2]],
                   stringsAsFactors = FALSE)
  df <- df[order(df$v1, df$index), , drop = FALSE]
  write_tab(df, out_file)
}

# KEGG_pathway2name.tab  (columns pathway / info)
write_pathway2name <- function(path2name, out_file) {
  df <- data.frame(pathway = path2name[[1]], info = path2name[[2]],
                   stringsAsFactors = FALSE)
  write_tab(df, out_file)
}

# GO2gene_<NS>.tab  (columns GO / Gene)
write_go2gene <- function(go2gene, out_file) {
  df <- data.frame(GO = go2gene[[1]], Gene = go2gene[[2]],
                   stringsAsFactors = FALSE)
  df <- unique(df[order(df$GO, df$Gene), , drop = FALSE])
  write_tab(df, out_file)
}

# GO2name_<NS>.tab  (columns GO / Term)
write_go2name <- function(go2name, out_file) {
  df <- data.frame(GO = go2name[[1]], Term = go2name[[2]],
                   stringsAsFactors = FALSE)
  df <- unique(df)
  write_tab(df, out_file)
}
