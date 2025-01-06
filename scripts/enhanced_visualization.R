library(EnhancedVolcano)
library(ComplexHeatmap)

create_enhanced_plots <- function(deg_results, normalized_counts, 
                                sample_info, cell_type, output_dir) {
    # Create output directory
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Volcano plot
    pdf(file.path(output_dir, paste0(cell_type, "_volcano_plot.pdf")))
    EnhancedVolcano(deg_results,
                    lab = deg_results$gene_id,
                    x = 'log2FoldChange',
                    y = 'padj',
                    title = paste(cell_type, "Differential Expression"),
                    pCutoff = 0.05,
                    FCcutoff = 1,
                    pointSize = 3.0,
                    labSize = 6.0)
    dev.off()
    
    # Top DEGs heatmap
    top_genes <- deg_results %>%
        filter(padj < 0.05) %>%
        top_n(50, abs(log2FoldChange)) %>%
        pull(gene_id)
    
    counts_matrix <- normalized_counts[top_genes, ]
    
    pdf(file.path(output_dir, paste0(cell_type, "_top_DEGs_heatmap.pdf")))
    Heatmap(scale(t(scale(t(counts_matrix)))),
            column_annotation = HeatmapAnnotation(
                CellType = sample_info$cell_type,
                Condition = sample_info$condition
            ),
            name = "Z-score",
            show_row_names = TRUE,
            cluster_columns = TRUE)
    dev.off()
} 