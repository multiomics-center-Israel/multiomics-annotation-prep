# ---------------------------------------------------------------------------
# download_kegg.R  -  fetch the reference tables from the KEGG REST API
# ---------------------------------------------------------------------------
# These are the same endpoints used manually in Neat_Annotation, here fetched
# automatically and cached:
#   ko  -> name    : https://rest.kegg.jp/list/ko
#   ko  -> pathway : https://rest.kegg.jp/link/pathway/ko
#   path-> name    : https://rest.kegg.jp/list/pathway
# ---------------------------------------------------------------------------

KEGG_REST <- list(
  ko_to_name    = "https://rest.kegg.jp/list/ko",
  ko_to_path    = "https://rest.kegg.jp/link/pathway/ko",
  pathway_names = "https://rest.kegg.jp/list/pathway"
)

# Returns a named list of local file paths (downloaded + cached).
download_kegg_rest <- function(cache_dir, refresh = FALSE) {
  list(
    ko_to_name    = cached_download(KEGG_REST$ko_to_name,
                                    "kegg_ko_to_name.txt",    cache_dir, refresh),
    ko_to_path    = cached_download(KEGG_REST$ko_to_path,
                                    "kegg_ko_to_path.txt",    cache_dir, refresh),
    pathway_names = cached_download(KEGG_REST$pathway_names,
                                    "kegg_pathway_names.txt", cache_dir, refresh)
  )
}

# ---- parsers ---------------------------------------------------------------

# ko -> list(names, title, ec).  Input line:  "K00844\tHK; hexokinase [EC:2.7.1.1]"
parse_ko_to_name <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  out <- vector("list", length(lines))
  ids  <- character(length(lines))
  for (i in seq_along(lines)) {
    parts <- strsplit(lines[i], "\t", fixed = TRUE)[[1]]
    ko   <- sub("^ko:", "", parts[1])
    info <- if (length(parts) >= 2) parts[2] else ""
    if (grepl(";", info, fixed = TRUE)) {
      sp    <- strsplit(info, "; ", fixed = TRUE)[[1]]
      names <- sp[1]
      title <- paste(sp[-1], collapse = "; ")
    } else {
      names <- ""
      title <- info
    }
    ec <- ""
    m  <- regmatches(title, regexpr("\\[EC:[^]]*\\]", title))
    if (length(m)) {
      ec    <- sub("^\\[(EC:[^]]*)\\]$", "\\1", m)
      title <- trimws(sub("\\s*\\[EC:[^]]*\\]", "", title))
    }
    ids[i]  <- ko
    out[[i]] <- list(names = names, title = title, ec = ec)
  }
  names(out) <- ids
  out
}

# ko -> character vector of "mapNNNNN" pathways.  Input:  "ko:K00844\tpath:map00010"
parse_ko_to_path <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  sp  <- strsplit(lines, "\t", fixed = TRUE)
  ko  <- sub("^ko:",   "", vapply(sp, `[`, character(1), 1))
  pth <- sub("^path:", "", vapply(sp, `[`, character(1), 2))
  keep <- grepl("^map", pth)                      # reference (map) pathways only
  split(pth[keep], ko[keep])
}

# pathway -> name.  Input:  "map00010\tGlycolysis / Gluconeogenesis"
parse_pathway_names <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  sp <- strsplit(lines, "\t", fixed = TRUE)
  id   <- sub("^path:", "", vapply(sp, `[`, character(1), 1))
  name <- vapply(sp, function(x) if (length(x) >= 2) x[2] else "", character(1))
  stats::setNames(name, id)
}
