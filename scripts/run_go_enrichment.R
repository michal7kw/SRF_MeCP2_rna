#!/usr/bin/env Rscript

# --- Argument Parsing ---
# Use commandArgs for basic argument handling
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript run_go_enrichment.R <input_csv> <cell_type> <output_dir>", call. = FALSE)
}

input_file <- args[1]
cell_type <- args[2]
output_dir <- args[3]

message("--- Independent GO Enrichment Analysis ---")
message("Input File: ", input_file)
message("Cell Type: ", cell_type)
message("Output Directory: ", output_dir)

# --- Load Libraries ---
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db) # Mouse annotation
  library(DOSE)
  library(enrichplot)
  library(dplyr)
  library(AnnotationDbi)
  library(stringr)
  library(ggplot2)
})

# --- Create Output Directory ---
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Read Input Data ---
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file, call. = FALSE)
}
message("Reading input file: ", input_file)
deg_results_df <- read.csv(input_file)

# --- Prepare Gene List ---
message("Filtering significant genes (padj < 0.05)...")
genes_df <- deg_results_df %>%
  filter(!is.na(padj) & padj < 0.05) %>%
  arrange(desc(log2FoldChange)) # Keep arrangement for consistency, though not strictly needed for enrichGO

# --- Check for Significant Genes ---
if (nrow(genes_df) == 0) {
  warning(paste("No significant genes found for", cell_type, "at padj < 0.05. Skipping GO analysis."), immediate. = TRUE)
  # Create empty files to avoid errors downstream if needed
  file.create(file.path(output_dir, paste0(cell_type, "_GO_dotplot.png")))
  file.create(file.path(output_dir, paste0(cell_type, "_GO_enrichment.csv")))
  quit(save = "no", status = 0) # Exit cleanly
}
message(nrow(genes_df), " significant genes found.")

# --- Convert GENE SYMBOLS to ENTREZ IDs ---
message("Mapping SYMBOL to ENTREZ IDs for ", cell_type, "...")
genes_df$ENTREZ <- tryCatch({
  mapIds(org.Mm.eg.db,
         keys = genes_df$gene,
         keytype = "SYMBOL",
         column = "ENTREZID",
         multiVals = "first")
}, error = function(e) {
  warning("Error mapping SYMBOL to ENTREZ IDs: ", e$message, immediate. = TRUE)
  return(rep(NA, nrow(genes_df)))
})
message("Finished mapping IDs for ", cell_type, ". Mapped ", sum(!is.na(genes_df$ENTREZ)), " genes.")

# --- Filter out Unmapped Genes ---
genes_df_mapped <- genes_df[!is.na(genes_df$ENTREZ), ]

if (nrow(genes_df_mapped) == 0) {
  warning(paste("No genes could be mapped to ENTREZ IDs for", cell_type, ". Skipping GO analysis."), immediate. = TRUE)
  # Create empty files
  file.create(file.path(output_dir, paste0(cell_type, "_GO_dotplot.png")))
  file.create(file.path(output_dir, paste0(cell_type, "_GO_enrichment.csv")))
  quit(save = "no", status = 0) # Exit cleanly
}
message(nrow(genes_df_mapped), " genes successfully mapped to ENTREZ IDs.")

# --- GO Enrichment Analysis ---
message("Starting GO Enrichment Analysis (enrichGO) for ", cell_type, "...")
go_results <- tryCatch({
  enrichGO(gene = genes_df_mapped$ENTREZ,
           OrgDb = org.Mm.eg.db,
           ont = "BP", # Biological Process
           pAdjustMethod = "BH",
           pvalueCutoff = 0.05, # Use 0.05 for enrichment analysis itself
           readable = TRUE) # Map back to gene symbols in results
}, error = function(e) {
  warning("Error during GO enrichment analysis: ", e$message, immediate. = TRUE)
  return(NULL)
})
message("Finished GO Enrichment Analysis (enrichGO) for ", cell_type, ".")

# --- Process and Plot GO Results ---
plot_path <- file.path(output_dir, paste0(cell_type, "_GO_dotplot.png"))
plot_saved_successfully <- FALSE

if (!is.null(go_results) && nrow(go_results@result) > 0) {
  message("Simplifying GO results...")
  go_results_simple <- simplify(go_results, cutoff = 0.7, by = "p.adjust", select_fun = min)

  # Filter out specific GO terms *after* simplification
  if (nrow(go_results_simple@result) > 0) {
      message("Filtering out GO:0007059 and GO:0006913...")
      go_results_simple@result <- go_results_simple@result %>%
          dplyr::filter(!ID %in% c("GO:0007059", "GO:0006913"))
  }

  if (nrow(go_results_simple@result) > 0) {
    message("Generating GO dotplot...")
    plot_device_opened <- FALSE
    tryCatch({
      png(plot_path, width = 10, height = 8, units = "in", res = 300)
      plot_device_opened <- TRUE
    }, error = function(e) {
      warning("Could not open PNG device: ", plot_path, " Error: ", e$message, immediate. = TRUE)
    })

    if (plot_device_opened) {
      # Wrap labels before plotting
      go_results_simple@result$Description <- stringr::str_wrap(go_results_simple@result$Description, width = 40)

      # Generate the plot using the filtered results, showing top 10
      p_dot <- dotplot(go_results_simple, showCategory = 10, x = "GeneRatio", color = "p.adjust", size = "Count") +
               ggplot2::theme(axis.text.y = ggplot2::element_text(size = 16)) # Increased font size yet again

      plot_printed <- FALSE
      tryCatch({
        print(p_dot)
        plot_printed <- TRUE
      }, error = function(e) {
        warning("Error printing dotplot: ", e$message, immediate. = TRUE)
        tryCatch({ # Attempt placeholder
          plot.new()
          text(0.5, 0.5, "Error generating dotplot.")
        }, error = function(e2) {
          warning("Failed to plot placeholder text: ", e2$message, immediate. = TRUE)
        })
      })

      dev.off() # Close PNG device

      if (plot_printed) {
        plot_saved_successfully <- TRUE
        message("GO dotplot saved to: ", plot_path)
      } else {
        warning("Plot file may be empty or contain an error message.", immediate. = TRUE)
      }
    } # End if plot_device_opened
  } else {
    warning("No significant GO terms remaining after simplification and filtering. Plot not generated.", immediate. = TRUE)
  }
} else {
  warning("No significant GO enrichment results found. Plot not generated.", immediate. = TRUE)
}

# --- Save GO Results Table ---
if (!is.null(go_results)) {
  go_df <- as.data.frame(go_results)
  if (nrow(go_df) > 0) {
      csv_path <- file.path(output_dir, paste0(cell_type, "_GO_enrichment.csv"))
      write.csv(go_df, csv_path, row.names = FALSE)
      message("GO enrichment results table saved to: ", csv_path)
  } else {
      message("GO enrichment result table is empty, not saving CSV.")
       # Create empty file if needed for consistency
      file.create(file.path(output_dir, paste0(cell_type, "_GO_enrichment.csv")))
  }

} else {
    message("No GO results object generated, skipping CSV save.")
     # Create empty file if needed for consistency
    file.create(file.path(output_dir, paste0(cell_type, "_GO_enrichment.csv")))
}


message("--- Independent GO Enrichment Analysis Finished for ", cell_type, " ---")