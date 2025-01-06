library(tidyverse)
library(pheatmap)
library(DESeq2)

# Function to create QC visualizations
create_qc_plots <- function(dds, output_dir) {
    # Create output directory
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # 1. Sample correlation heatmap
    vst_counts <- vst(dds, blind=FALSE)
    sample_cors <- cor(assay(vst_counts))
    
    pdf(file.path(output_dir, "sample_correlation_heatmap.pdf"))
    pheatmap(sample_cors,
             main="Sample Correlation Heatmap",
             annotation_col=data.frame(
                 CellType=dds$cell_type,
                 Condition=dds$condition,
                 row.names=colnames(dds)
             ))
    dev.off()
    
    # 2. Count distribution plots
    counts_df <- as.data.frame(counts(dds, normalized=TRUE)) %>%
        rownames_to_column("gene") %>%
        gather(key="sample", value="counts", -gene)
    
    ggplot(counts_df, aes(x=log2(counts + 1))) +
        geom_density() +
        facet_wrap(~sample) +
        theme_bw() +
        ggtitle("Distribution of Normalized Counts")
    ggsave(file.path(output_dir, "count_distributions.pdf"))
    
    # 3. Number of detected genes per sample
    detected_genes <- colSums(counts(dds) > 0)
    pdf(file.path(output_dir, "detected_genes.pdf"))
    barplot(detected_genes, 
            main="Number of Detected Genes per Sample",
            las=2)
    dev.off()
} 