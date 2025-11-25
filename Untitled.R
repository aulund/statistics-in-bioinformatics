## Install CRAN & Bioconductor packages needed for the course

cran_pkgs <- c(
  "dplyr", "tibble", "umap", "igraph", "Seurat", "pheatmap", "ggplot2",
  "ggrepel", "Rtsne", "enrichR", "openxlsx", "RColorBrewer", "gridExtra",
  "stringr", "knitr", "prettydoc", "caret", "pROC", "C50", "statmod",
  "e1071", "reshape", "MASS", "neuralnet", "pvclust", "randomForest",
  "ncvreg", "devtools", "rmarkdown", "prettydoc", "knitr", "stats"
)

bioc_pkgs <- c(
  "DESeq2", "edgeR", "limma", "bluster", "pathview",
  "clusterProfiler", "org.Hs.eg.db", "enrichplot"
)

## Install CRAN packages
to_install_cran <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(to_install_cran)) {
  install.packages(to_install_cran)
}

## Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
to_install_bioc <- setdiff(bioc_pkgs, rownames(installed.packages()))
if (length(to_install_bioc)) {
  BiocManager::install(to_install_bioc)
}

## install github package BioStudies:

library(devtools)
devtools::install_gitlab(repo = "wolftower/biostudies")
