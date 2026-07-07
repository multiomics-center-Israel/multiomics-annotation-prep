# bootstrap.R  -  locate the repo root and source all R modules.
# Sourced by the scripts in scripts/.

.find_repo_root <- function() {
  # 1) if run via Rscript, use the script's own path
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(file_arg)) {
    d <- normalizePath(dirname(file_arg))
    # scripts/ live one level below the repo root
    if (basename(d) == "scripts") return(dirname(d))
    if (dir.exists(file.path(d, "R"))) return(d)
    if (dir.exists(file.path(dirname(d), "R"))) return(dirname(d))
  }
  # 2) fall back to the working directory
  if (dir.exists("R")) return(normalizePath("."))
  stop("Could not locate repo root (no 'R/' directory found).")
}

REPO_ROOT <- .find_repo_root()

.source_modules <- function(root = REPO_ROOT) {
  files <- c("utils.R",
             "download_kegg.R", "download_kegg_org.R", "prepare_kegg_nonmodel.R",
             "download_go.R", "prepare_go.R",
             "prepare_ensembl.R", "prepare_uniprot.R")
  for (f in files) source(file.path(root, "R", f))
}

.source_modules()
