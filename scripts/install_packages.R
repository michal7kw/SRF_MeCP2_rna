#!/usr/bin/env Rscript

# List of required packages
packages <- c("DESeq2", "tidyverse")

# Function to install packages if not already installed
install_if_missing <- function(package) {
    if (!requireNamespace(package, quietly = TRUE)) {
        if (package == "DESeq2") {
            if (!requireNamespace("BiocManager", quietly = TRUE))
                install.packages("BiocManager")
            BiocManager::install("DESeq2")
        } else {
            install.packages(package)
        }
    }
}

# Install packages
sapply(packages, install_if_missing) 