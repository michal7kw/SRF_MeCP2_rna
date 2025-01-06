#!/usr/bin/env Rscript

# List of required packages
packages <- c(
    "DESeq2", "tidyverse", "pheatmap", "EnhancedVolcano",
    "clusterProfiler", "org.Hs.eg.db", "DOSE", "enrichplot",
    "ComplexHeatmap", "BiocManager"
)

# Function to install packages if not already installed
install_if_missing <- function(package) {
    if (!requireNamespace(package, quietly = TRUE)) {
        if (package %in% c("DESeq2", "clusterProfiler", "org.Hs.eg.db", 
                          "DOSE", "enrichplot", "ComplexHeatmap")) {
            if (!requireNamespace("BiocManager", quietly = TRUE))
                install.packages("BiocManager")
            BiocManager::install(package)
        } else {
            install.packages(package)
        }
    }
}

# Install packages
sapply(packages, install_if_missing) 