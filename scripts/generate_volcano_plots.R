# Script to generate volcano plots from DESeq2 results
# with specific styling and handling for the Mecp2 gene.

library(tidyverse)
library(ggrepel)
library(org.Hs.eg.db) # Ensure this is installed if needed for gene symbols

setwd("D:\\Github\\SRF_MeCP2_rna\\")


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
#' @param exclude_target_gene Logical, whether to filter out the gene_to_label before plotting.
#'
#' @return The ggplot object.
create_styled_volcano <- function(results_df, cell_type, output_file,
                                   p_threshold, fc_threshold, gene_to_label,
                                   exclude_target_gene = TRUE) { # Default to excluding

    # Ensure gene symbols exist
    if (!"gene_symbol" %in% colnames(results_df)) {
        stop("Results data frame must contain a 'gene_symbol' column.")
    }

    # --- Store original data before potential filtering ---
    original_results_df <- results_df

    # --- Conditionally filter out the specific gene ---
    if (exclude_target_gene) {
        results_df <- results_df %>% filter(gene_symbol != gene_to_label)
        message(paste("Filtered out", gene_to_label, "for", cell_type, "plot."))
    } else {
        # If including the gene, prepare its data for labeling later
        label_data <- original_results_df %>% filter(gene_symbol == gene_to_label)
        if (nrow(label_data) == 0) {
            warning(paste("Gene to label", gene_to_label, "not found in the results for", cell_type))
        }
    }

    # Remove rows with NA p-values which cannot be plotted (from the potentially filtered df)
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

    # Define colors similar to the example image
    color_map <- c(
        "NS" = "#404040",      # Very dark grey for non-significant
        "Log2 FC" = "green",         # Intense Green
        "p-value" = "blue",          # Intense Blue
        "p-value and log2 FC" = "red" # Intense Red
    )
    
    # --- Calculate y-axis limit --- 
    if (exclude_target_gene) {
        # Case 1: Target gene EXCLUDED - Focus on other significant genes
        y_limit_data_df <- results_df # Use the filtered data (target gene already removed)
        
        y_limit_data <- y_limit_data_df %>%
          filter(significance != "NS", padj > 0) %>% # Exclude p=0 for quantile calculation
          pull(padj)
          
        # Handle case where there are few/no other significant genes
        if (length(y_limit_data) < 10) {
            # Fallback: use max -log10(p) of remaining genes, or 10 if none exist
             max_y_base <- max(c(10, -log10(results_df$padj[results_df$padj > 0])), na.rm = TRUE)
        } else {
             # Use 99th percentile of significant genes (excluding target)
             max_y_base <- -log10(quantile(y_limit_data, 0.01, na.rm = TRUE)) 
        }
        
        # Add buffer
        max_y <- max_y_base * 1.1 
        
        # Handle Inf cases 
        if (is.infinite(max_y)) {
             finite_y_values <- -log10(results_df$padj[results_df$padj > 0])
             max_y <- if (length(finite_y_values) > 0) max(finite_y_values, na.rm=TRUE) * 1.2 else 350
        }
        # Ensure minimum limit and apply cap
        max_y <- max(max_y, 10)
        max_y <- min(max_y, 75) # Apply cap only when target gene is excluded
        
    } else {
        # Case 2: Target gene INCLUDED - Ensure it's visible
        # Use the *original* data frame before any filtering for y-limit calculation
        y_vals_all <- -log10(original_results_df$padj)
        y_vals_all <- y_vals_all[!is.na(y_vals_all) & is.finite(y_vals_all)] # Remove NA/Inf
        
        if (length(y_vals_all) == 0) {
            max_y_base <- 10 # Default if no valid p-values
        } else {
            max_y_base <- max(y_vals_all, na.rm = TRUE)
        }
        
        # Add buffer (slightly larger buffer might be good)
        max_y <- max_y_base * 1.15 
        
        # Ensure minimum limit (no upper cap needed here)
        max_y <- max(max_y, 10) 
    }
    
    # Handle potential remaining Inf case for max_y (if all pvals were 0 initially)
    if (is.infinite(max_y)) {
        max_y <- 350 # Set arbitrary high limit if calculation resulted in Inf
    }

    # Create the plot
    p <- ggplot(results_df, aes(x = log2FoldChange, y = -log10(padj)))
    
    # Add layers sequentially
    p <- p +
        geom_point(aes(color = significance), alpha = 0.6, size = 1.5) # Plot points
    
    p <- p +
        geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "grey50") # Significance lines
    
    p <- p +
        geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "grey50")
    
    p <- p +
        scale_color_manual(values = color_map) # Set colors
        
    p <- p +
        theme_bw(base_size = 14) # Set theme
        
    p <- p +
        theme(
            legend.position = "top", 
            legend.title = element_blank(),
            panel.grid.major = element_line(colour = "grey90"),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold")
        )
        
    p <- p +
        labs( # Add labels and title
            x = expression(Log[2]~fold~change),
            y = expression(-Log[10]~italic(P)),
            title = sprintf("%s: Differential Expression (Mecp2 vs Control)", cell_type)
        )

    # Use coord_cartesian to zoom without removing data points outside limits
    # Set xlim based on data range (use original_results_df if target gene included, otherwise results_df)
    # Set ylim based on calculated max_y
    plot_data_for_lims <- if (exclude_target_gene) results_df else original_results_df
    x_limits <- c(min(plot_data_for_lims$log2FoldChange, na.rm=TRUE)-0.5, max(plot_data_for_lims$log2FoldChange, na.rm=TRUE)+0.5) # Add slight padding to xlim

    p <- p + coord_cartesian(
        ylim = c(0, max_y),
        xlim = x_limits
    )

    # --- Add the specific gene label using ggrepel ONLY if included ---
    if (!exclude_target_gene && nrow(label_data) > 0) {
        p <- p + geom_text_repel(
            data = label_data,
            aes(label = gene_symbol),
            size = 4,
            fontface = "bold",
            color = "black", # Ensure label is visible
            nudge_y = (max_y * 0.02), # Nudge slightly relative to y-axis scale
            nudge_x = (x_limits[2] - x_limits[1]) * 0.05, # Nudge based on x-axis range
            box.padding = 0.5,
            point.padding = 0.5,
            segment.color = 'grey50',
            max.overlaps = Inf # Ensure this label is shown
        )
    }

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
    output_plot_file_excluded <- file.path(OUTPUT_DIR, sprintf("%s_styled_volcano_plot_noMecp2.pdf", ct))
    output_plot_file_included <- file.path(OUTPUT_DIR, sprintf("%s_styled_volcano_plot_withMecp2.pdf", ct)) # New filename

    if (file.exists(results_file)) {
        message(paste("Processing:", results_file))
        
        # Read the results data
        res_data <- read.csv(results_file)

        # --- Create the volcano plot EXCLUDING Mecp2 ---
        message(paste("Generating plot excluding", GENE_TO_LABEL))
        create_styled_volcano(
            results_df = res_data,
            cell_type = ct,
            output_file = output_plot_file_excluded,
            p_threshold = P_THRESHOLD,
            fc_threshold = FC_THRESHOLD,
            gene_to_label = GENE_TO_LABEL,
            exclude_target_gene = TRUE # Explicitly exclude
        )

        # --- Create the volcano plot INCLUDING Mecp2 ---
        message(paste("Generating plot including", GENE_TO_LABEL))
        create_styled_volcano(
            results_df = res_data,
            cell_type = ct,
            output_file = output_plot_file_included,
            p_threshold = P_THRESHOLD,
            fc_threshold = FC_THRESHOLD,
            gene_to_label = GENE_TO_LABEL,
            exclude_target_gene = FALSE # Explicitly include
        )

    } else {
        warning(paste("Results file not found, skipping:", results_file))
    }
}

message("Volcano plot generation complete.") 