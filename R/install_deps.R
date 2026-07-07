# Install all dependencies for multiomic-annotation-prep.
# Run once:  Rscript -e 'source("R/install_deps.R")'

cran_pkgs <- c("httr", "yaml", "optparse")
bioc_pkgs <- c("KEGGREST", "clusterProfiler", "GO.db", "AnnotationDbi", "biomaRt")

install_cran <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
}

install_bioc <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) BiocManager::install(missing, update = FALSE, ask = FALSE)
}

message("Installing CRAN packages ...")
install_cran(cran_pkgs)
message("Installing Bioconductor packages ...")
install_bioc(bioc_pkgs)
message("Done. All dependencies installed.")
