# Script to generate volcano plots from DESeq2 results
# with specific styling and handling for the Mecp2 gene.

library(tidyverse)
library(ggrepel)
library(org.Hs.eg.db) # Ensure this is installed if needed for gene symbols

setwd("/home/michal/Github/SRF_MeCP2_rna")


# --- Configuration ---
RESULTS_DIR <- "results/05_deseq2"
OUTPUT_DIR <- "results/06_results"
CELL_TYPES <- c("NPC", "Neuron")
P_THRESHOLD <- 0.05
FC_THRESHOLD <- 0.5 # Log2 fold change threshold
GENE_TO_LABEL <- "MECP2" # Use the gene symbol as it appears in your results

# --- Function Definition ---

#' Create a Volcano Plot with Specific Styling
#'
#' Generates a volcano plot similar to the provided example image,
#' handling a specific gene (e.g., Mecp2) separately for labeling
#' and axis scaling.
#'
#' @param results_df Data frame containing differential expression results.
#'        Must include columns: 'log2FoldChange', 'padj', 'gene_symbol'.
#' @param cell_type Character string for the cell type (used in title).
#' @param output_file Path where the PDF plot will be saved.
#' @param p_threshold Adjusted p-value significance threshold.
#' @param fc_threshold Absolute log2 fold change significance threshold.
#' @param gene_to_label The specific gene symbol to label.
#'
#' @return The ggplot object.
create_styled_volcano <- function(results_df, cell_type, output_file,
                                   p_threshold, fc_threshold, gene_to_label) {

    # Ensure gene symbols exist
    if (!"gene_symbol" %in% colnames(results_df)) {
        stop("Results data frame must contain a 'gene_symbol' column.")
    }
    
    # --- Filter out the specific gene BEFORE any processing ---
    results_df <- results_df %>% filter(gene_symbol != gene_to_label)
    message(paste("Filtered out", gene_to_label, "for", cell_type, "plot."))

    # Remove rows with NA p-values which cannot be plotted
    results_df <- results_df %>% filter(!is.na(padj))

    # Define significance categories based on example image
    results_df <- results_df %>%
        mutate(
            significance = case_when(
                padj < p_threshold & abs(log2FoldChange) >= fc_threshold ~ paste("p-value and log2 FC"),
                padj < p_threshold & abs(log2FoldChange) < fc_threshold ~ "p-value",
                padj >= p_threshold & abs(log2FoldChange) >= fc_threshold ~ "Log2 FC",
                TRUE ~ "NS" # Not Significant
            ),
            # Ensure significance is a factor with the desired order for the legend
            significance = factor(significance, levels = c("NS", "Log2 FC", "p-value", "p-value and log2 FC"))
        )

    # Identify the specific gene to label - REMOVED as gene is filtered out
    # label_data <- results_df %>% filter(gene_symbol == gene_to_label)
    # 
    # if (nrow(label_data) == 0) {
    #     warning(paste("Gene to label", gene_to_label, "not found in the results for", cell_type))
    # }

    # Define colors similar to the example image
    color_map <- c(
        "NS" = "grey80", # Lighter grey for non-significant
        "Log2 FC" = "#90EE90", # Light green
        "p-value" = "#ADD8E6", # Light blue
        "p-value and log2 FC" = "#FF6347"  # Tomato red
    )
    
    # Calculate y-axis limit based on non-labeled significant genes to avoid distortion
    # Use the 99th percentile of -log10(padj) among significant genes
    # (gene_to_label is already removed from results_df)
    y_limit_data <- results_df %>%
      filter(significance != "NS", padj > 0) %>% # Exclude p=0 for quantile calculation
      pull(padj)
      
    # Handle case where there are few/no other significant genes
    if (length(y_limit_data) < 10) {
        # Fallback: use max -log10(p) of all (remaining) genes, or 10 if none exist
         max_y_base <- max(c(10, -log10(results_df$padj[results_df$padj > 0])), na.rm = TRUE)
    } else {
         # Use 99th percentile of significant genes
         max_y_base <- -log10(quantile(y_limit_data, 0.01, na.rm = TRUE)) # Changed from 0.001 to 0.01
    }

    # Add buffer
    max_y <- max_y_base * 1.1 

    # Handle Inf cases if p-values are extremely small
    if (is.infinite(max_y)) {
         # Find max finite -log10(p) and add a buffer, or set arbitrary high limit if all are Inf
         finite_y_values <- -log10(results_df$padj[results_df$padj > 0])
         max_y <- if (length(finite_y_values) > 0) max(finite_y_values, na.rm=TRUE) * 1.2 else 350
    }
    # Ensure max_y is at least moderately high if data range is very small
    max_y <- max(max_y, 10) 

    # --- Add an upper cap to max_y to prevent extreme outliers from dominating ---
    max_y <- min(max_y, 75) # Cap the y-axis limit at 75


    # Create the plot
    p <- ggplot(results_df, aes(x = log2FoldChange, y = -log10(padj))) +
        # Plot points, setting alpha and size
        geom_point(aes(color = significance), alpha = 0.6, size = 1.5) +
        
        # Add the specific gene label using ggrepel - REMOVED
        # geom_text_repel(
        #     data = label_data,
        #     aes(label = gene_symbol),
        #     size = 4,
        #     fontface = "bold",
        #     nudge_y = (max_y * 0.02), # Nudge slightly relative to y-axis scale
        #     nudge_x = 0.5,
        #     box.padding = 0.5,
        #     point.padding = 0.5,
        #     segment.color = 'grey50',
        #     max.overlaps = Inf # Ensure this label is shown
        # ) +

        # Add significance threshold lines
        geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "grey50") +
        geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "grey50") +

        # Set colors
        scale_color_manual(values = color_map) +

        # Set theme similar to example
        theme_bw(base_size = 14) +
        theme(
            legend.position = "top", # Position legend at the top
            legend.title = element_blank(), # Remove legend title
            panel.grid.major = element_line(colour = "grey90"),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold") # Center title
        ) +

        # Add labels and title
        labs(
            x = expression(Log[2]~fold~change),
            y = expression(-Log[10]~italic(P)),
            title = sprintf("%s: Differential Expression (Mecp2 vs Control)", cell_type)
        ) +
        
        # Use coord_cartesian to zoom without removing data points outside limits
        # Set xlim based on data range, ylim based on calculated max_y
        coord_cartesian(
            ylim = c(0, max_y),
            xlim = c(min(results_df$log2FoldChange, na.rm=TRUE)-0.5, max(results_df$log2FoldChange, na.rm=TRUE)+0.5) # Add slight padding to xlim
            )

    # Save the plot
    ggsave(output_file, p, width = 8, height = 7, device = "pdf")
    
    print(paste("Saved volcano plot to:", output_file))

    return(p)
}

# --- Main Script Logic ---

# Create output directory if it doesn't exist
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Loop through each cell type and generate the plot
for (ct in CELL_TYPES) {
    results_file <- file.path(RESULTS_DIR, sprintf("%s_differential_expression.csv", ct))
    output_plot_file <- file.path(OUTPUT_DIR, sprintf("%s_styled_volcano_plot_noMecp2.pdf", ct)) # Changed filename

    if (file.exists(results_file)) {
        message(paste("Processing:", results_file))
        
        # Read the results data
        res_data <- read.csv(results_file)

        # Create the volcano plot
        create_styled_volcano(
            results_df = res_data,
            cell_type = ct,
            output_file = output_plot_file,
            p_threshold = P_THRESHOLD,
            fc_threshold = FC_THRESHOLD,
            gene_to_label = GENE_TO_LABEL
        )
    } else {
        warning(paste("Results file not found, skipping:", results_file))
    }
}

message("Volcano plot generation complete.") 