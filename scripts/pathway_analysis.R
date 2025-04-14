#!/usr/bin/env Rscript

# Load libraries
library(clusterProfiler)
# library(org.Hs.eg.db) # Switched to Mouse
library(org.Mm.eg.db) # Added Mouse annotation
library(DOSE)
library(enrichplot)
library(dplyr) # Ensure dplyr is loaded for the pipe operator and filter/arrange
library(AnnotationDbi) # Needed for mapIds
library(stringr) # For wrapping text labels

# --- Function Definition ---
perform_pathway_analysis <- function(deg_results_df, cell_type, output_dir) {
    # Create output directory
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    # Prepare gene list (using the input dataframe directly)
    genes_df <- deg_results_df %>%
        filter(!is.na(padj) & padj < 0.05) %>%
        arrange(desc(log2FoldChange))

    # Check if there are any significant genes
    if (nrow(genes_df) == 0) {
        warning(paste("No significant genes found for", cell_type, "at padj < 0.05. Skipping pathway analysis."))
        # Create empty files to avoid errors in downstream processes if needed
        file.create(file.path(output_dir, paste0(cell_type, "_pathway_analysis.pdf")))
        file.create(file.path(output_dir, paste0(cell_type, "_GO_enrichment.csv")))
        file.create(file.path(output_dir, paste0(cell_type, "_GSEA_results.csv")))
        return(NULL)
    }

    # Convert GENE SYMBOLS to ENTREZ IDs
    # Handle potential errors during mapping (using gene symbols now)
    message("Mapping SYMBOL to ENTREZ IDs for ", cell_type, "...")
    genes_df$ENTREZ <- tryCatch({
        mapIds(org.Mm.eg.db, # Use Mouse DB
               keys=genes_df$gene, # Use the 'gene' column directly
               keytype="SYMBOL",   # Keytype is SYMBOL
               column="ENTREZID", # Target column
               multiVals="first") # Take the first mapping if multiple exist
    }, error = function(e) {
        warning("Error mapping SYMBOL to ENTREZ IDs: ", e$message)
        return(rep(NA, nrow(genes_df)))
    })
    message("Finished mapping IDs for ", cell_type, ". Mapped ", sum(!is.na(genes_df$ENTREZ)), " genes.")

    # Filter out genes where mapping failed
    genes_df_mapped <- genes_df[!is.na(genes_df$ENTREZ), ]

    if (nrow(genes_df_mapped) == 0) {
        warning(paste("No genes could be mapped to ENTREZ IDs for", cell_type, ". Skipping pathway analysis."))
        # Create empty files
        file.create(file.path(output_dir, paste0(cell_type, "_pathway_analysis.pdf")))
        file.create(file.path(output_dir, paste0(cell_type, "_GO_enrichment.csv")))
        file.create(file.path(output_dir, paste0(cell_type, "_GSEA_results.csv")))
        return(NULL)
    }

    # Create ranked gene list for GSEA
    ranked_genes <- genes_df_mapped$log2FoldChange
    names(ranked_genes) <- genes_df_mapped$ENTREZ
    ranked_genes <- na.omit(ranked_genes) # Should already be handled, but good practice
    ranked_genes <- sort(ranked_genes, decreasing = TRUE) # GSEA requires sorted list

    # GO Enrichment Analysis
    message("Starting GO Enrichment Analysis (enrichGO) for ", cell_type, "...")
    go_results <- tryCatch({
        enrichGO(gene = genes_df_mapped$ENTREZ,
                 OrgDb = org.Mm.eg.db, # Use Mouse DB
                 ont = "BP", # Biological Process
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 readable = TRUE) # Automatically maps ENTREZID to gene symbols
    }, error = function(e) {
        warning("Error during GO enrichment analysis: ", e$message)
        return(NULL)
    })
    message("Finished GO Enrichment Analysis (enrichGO) for ", cell_type, ".")
    # GSEA
    message("Starting GSEA Analysis (gseGO) for ", cell_type, "...")
    gsea_results <- tryCatch({
        gseGO(geneList = ranked_genes,
              OrgDb = org.Mm.eg.db, # Use Mouse DB
              ont = "BP",
              pvalueCutoff = 0.05,
              verbose = FALSE,
              pAdjustMethod = "BH")
    }, error = function(e) {
        warning("Error during GSEA analysis: ", e$message)
        return(NULL)
    })
    message("Finished GSEA Analysis (gseGO) for ", cell_type, ".")

    # Prepare plot path (PNG)
    plot_path <- file.path(output_dir, paste0(cell_type, "_pathway_analysis.png"))
    plot_saved_successfully <- FALSE # Flag to track if plot was saved

    # Plot GO enrichment if results exist
    if (!is.null(go_results) && nrow(go_results@result) > 0) {
        # Simplify results first
        go_results_simple <- simplify(go_results, cutoff=0.7, by="p.adjust", select_fun=min)

        # Filter out specific GO terms before selecting top N for visualization
        if (nrow(go_results_simple@result) > 0) {
             go_results_simple@result <- go_results_simple@result %>%
                 dplyr::filter(!ID %in% c("GO:0007059", "GO:0006913"))
        }

        # Check if any terms remain *after* simplification
        if (nrow(go_results_simple@result) > 0) {
            # Attempt to open PNG device *only if* there are simplified results to plot
            plot_device_opened <- FALSE
            tryCatch({
                png(plot_path, width = 10, height = 8, units = "in", res = 300) # Reverted height
                plot_device_opened <- TRUE
            }, error = function(e) {
                warning("Could not open PNG device: ", plot_path, " Error: ", e$message)
            })

            if (plot_device_opened) {
                # Wrap labels before plotting
                go_results_simple@result$Description <- stringr::str_wrap(go_results_simple@result$Description, width = 40) # Adjust width as needed

                # Generate the plot using the modified results
                # Generate the plot using the filtered and modified results, showing top 10
                p_dot <- dotplot(go_results_simple, showCategory=10, x = "GeneRatio", color = "p.adjust", size = "Count") +
                         ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10)) # Increased font size

                # Try to print the plot
                plot_printed <- FALSE
                tryCatch({
                    print(p_dot)
                    plot_printed <- TRUE
                }, error = function(e) {
                    warning("Error printing dotplot: ", e$message)
                    # Attempt to plot placeholder text if printing fails
                    tryCatch({
                        plot.new()
                        text(0.5, 0.5, "Error generating dotplot.")
                    }, error = function(e2) {
                        warning("Failed to plot placeholder text: ", e2$message)
                    })
                })

                # Close the device
                dev.off()

                # Set success flag only if the plot was actually printed without error
                if (plot_printed) {
                    plot_saved_successfully <- TRUE
                } else {
                     warning("Plot file may be empty or contain an error message due to plotting issues.")
                     # Optionally delete the empty/error file:
                     # tryCatch(file.remove(plot_path), warning = function(w) warning("Could not remove potentially empty plot file: ", plot_path))
                }
            } # End if plot_device_opened
        } else {
            warning("No significant GO terms remaining after simplification. Plot not generated.")
        }
    } else {
        warning("No significant GO enrichment results found. Plot not generated.")
    }


    # Save tables
    if (!is.null(go_results)) {
        go_df <- as.data.frame(go_results)
        write.csv(go_df,
                  file.path(output_dir, paste0(cell_type, "_GO_enrichment.csv")),
                  row.names = FALSE)
        message("GO enrichment results saved.")
    }
    if (!is.null(gsea_results)) {
        gsea_df <- as.data.frame(gsea_results)
        # Add leading edge gene symbols if readable=TRUE wasn't used or failed for GSEA
         if (!"leadingGeneSymbols" %in% colnames(gsea_df) && "core_enrichment" %in% colnames(gsea_df)) {
             gsea_df$leadingGeneSymbols <- sapply(gsea_df$core_enrichment, function(x) {
                 if(is.na(x) || x == "") return(NA)
                 genes <- strsplit(x, "/")[[1]]
                 symbols <- tryCatch({
                     mapIds(org.Mm.eg.db, # Use Mouse DB
                            keys = genes,
                            keytype = "ENTREZID",
                            column = "SYMBOL",
                            multiVals = "list")
                 }, error = function(e) { list() }) # Return empty list on error
                 paste(unlist(symbols), collapse = "/")
             })
         }
        write.csv(gsea_df,
                  file.path(output_dir, paste0(cell_type, "_GSEA_results.csv")),
                  row.names = FALSE)
        message("GSEA results saved.")
    }
    
    return(plot_saved_successfully) # Return the status
}

# --- Main Execution Block ---
if (sys.nframe() == 0) { # Check if the script is run directly
    # Define file paths and parameters directly
    neu_file <- "DEA_NEU.csv" # Updated input file
    nsc_file <- "DEA_NSC.csv" # Updated input file
    output_dir <- "results/07_standalone_pathway_analysis" # Define a common output directory

    # --- Analysis for NEU ---
    message("--- Starting Pathway Analysis for NEU ---")
    if (file.exists(neu_file)) {
        message("Reading input file: ", neu_file)
        neu_deg_results <- read.csv(neu_file)
        message("Performing pathway analysis for cell type: NEU")
        neu_success <- perform_pathway_analysis(neu_deg_results, "NEU", output_dir) # Updated cell type name
        if (neu_success) {
            message("NEU pathway analysis plot successfully saved.")
        } else {
            message("NEU pathway analysis plot generation failed or was skipped.")
        }
    } else {
        warning("NEU input file not found: ", neu_file)
    }
    message("--- Finished Pathway Analysis for NEU ---")
    message("\n") # Add a newline for separation

    # --- Analysis for NSC ---
    message("--- Starting Pathway Analysis for NSC ---")
    if (file.exists(nsc_file)) {
        message("Reading input file: ", nsc_file)
        nsc_deg_results <- read.csv(nsc_file)
        message("Performing pathway analysis for cell type: NSC")
        nsc_success <- perform_pathway_analysis(nsc_deg_results, "NSC", output_dir) # Updated cell type name
        if (nsc_success) {
            message("NSC pathway analysis plot successfully saved.")
        } else {
            message("NSC pathway analysis plot generation failed or was skipped.")
        }
    } else {
        warning("NSC input file not found: ", nsc_file)
    }
    message("--- Finished Pathway Analysis for NSC ---")

    message("\nPathway analysis script finished for both cell types.")
}