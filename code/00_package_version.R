############################################################
## 00_session_info.R
## Load all packages used in the pipeline and print versions
############################################################

## List of packages used across the full pipeline
pkgs <- c(
  # Core data handling & plotting
  "dplyr",
  "tidyr",
  "purrr",
  "tibble",
  "stringr",
  "ggplot2",
  "ggrepel",
  "patchwork",
  
  # Differential expression & stats
  "DESeq2",
  "metap",
  
  # Enrichment & pathway analysis
  "clusterProfiler",
  "enrichplot",
  "fgsea",
  "msigdbr",
  
  # Annotation
  "org.Mm.eg.db",
  "biomaRt",
  
  # Multivariate analysis
  "mixOmics",
  "FactoMineR",
  
  # Heatmaps / utilities
  "pheatmap",
  
  # I/O
  "openxlsx"
)

## ---------------------------------------------------------
## Check availability
## ---------------------------------------------------------

installed <- rownames(installed.packages())
missing   <- setdiff(pkgs, installed)

if (length(missing) > 0) {
  message("The following packages are NOT installed:")
  message(paste(" -", missing, collapse = "\n"))
} else {
  message("All required packages are installed.")
}

## ---------------------------------------------------------
## Load packages
## ---------------------------------------------------------

suppressPackageStartupMessages({
  invisible(lapply(pkgs[pkgs %in% installed], library, character.only = TRUE))
})

## ---------------------------------------------------------
## Print versions
## ---------------------------------------------------------

pkg_versions <- data.frame(
  package = pkgs[pkgs %in% installed],
  version = sapply(pkgs[pkgs %in% installed], function(p) {
    as.character(packageVersion(p))
  }),
  stringsAsFactors = FALSE
)

cat("\n===== Package versions used in this pipeline =====\n")
print(pkg_versions, row.names = FALSE)

cat("\n===== R session info =====\n")
print(sessionInfo())


## ---------------------------------------------------------
## Save package versions table
## ---------------------------------------------------------
pkg_versions <- data.frame(
  package = pkgs,
  version = sapply(pkgs, function(p) as.character(packageVersion(p))),
  stringsAsFactors = FALSE
)

write.table(
  pkg_versions,
  file = "package_versions.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

## ---------------------------------------------------------
## Save full session info
## ---------------------------------------------------------
session_file <- "session_info.txt"
sink(session_file)
cat("===== R session info =====\n\n")
print(sessionInfo())
sink()

cat("Session info saved to:", session_file, "\n")

## ---------------------------------------------------------
##  Report
## ---------------------------------------------------------

cat("\n===== Recorded package versions =====\n")
print(
  data.frame(
    package = pkgs,
    version = sapply(pkgs, function(p) as.character(packageVersion(p))),
    row.names = NULL
  )
)
