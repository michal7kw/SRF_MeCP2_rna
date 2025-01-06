library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(tidyverse)
library(DESeq2)

# Read DESeq2 results and counts
read_deg_results <- function(file_path) {
    read.csv(file_path) %>%
        mutate(
            significant = padj < 0.05 & abs(log2FoldChange) > 1,
            regulation = case_when(
                significant & log2FoldChange > 1 ~ "Up",
                significant & log2FoldChange < -1 ~ "Down",
                TRUE ~ "NS"
            )
        )
}

# Read results
npc_results <- read_deg_results("90-1115332406/05_deseq2/NPC_differential_expression.csv")
neuron_results <- read_deg_results("90-1115332406/05_deseq2/Neuron_differential_expression.csv")

# Read counts and sample info
counts_data <- read.table("90-1115332406/04_counts/all_samples_counts.txt", 
                         header=TRUE, row.names=1, skip=1)
count_cols <- grep("X90\\.1115332406\\.03_alignment\\.h[NP][GM][123]\\.Aligned\\.sortedByCoord\\.out\\.bam$", 
                  colnames(counts_data), value=TRUE)
counts_data <- counts_data[, count_cols]
colnames(counts_data) <- sub("X90\\.1115332406\\.03_alignment\\.(h[NP][GM][123])\\..*", "\\1", 
                            colnames(counts_data))

# Create sample information
sample_info <- data.frame(
    sample = colnames(counts_data)
) %>%
    mutate(
        condition = factor(case_when(
            grepl("G", substr(sample, 3, 4)) ~ "G",
            grepl("M", substr(sample, 3, 4)) ~ "M"
        ), levels=c("G", "M")),
        cell_type = factor(ifelse(substr(sample, 2, 2) == "P", "NPC", "Neuron"))
    )

# Create DESeq2 object for normalized counts
dds <- DESeqDataSetFromMatrix(
    countData = counts_data,
    colData = sample_info,
    design = ~ cell_type + condition
)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

# Create output directories
pathway_dir <- "90-1115332406/06_results/pathway_analysis"
viz_dir <- "90-1115332406/06_results/enhanced_visualization"
dir.create(pathway_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(viz_dir, recursive = TRUE, showWarnings = FALSE)

# Function for pathway analysis
perform_pathway_analysis <- function(deg_results, cell_type) {
    # Convert ENSEMBL IDs to ENTREZ IDs
    genes_df <- deg_results %>%
        filter(!is.na(padj) & padj < 0.05) %>%
        arrange(desc(log2FoldChange))
    
    genes_df$ENTREZ <- mapIds(org.Hs.eg.db,
                             keys=sub("\\..*", "", genes_df$gene_id),
                             keytype="ENSEMBL",
                             column="ENTREZID")
    
    # Create ranked gene list for GSEA
    ranked_genes <- genes_df$log2FoldChange
    names(ranked_genes) <- genes_df$ENTREZ
    ranked_genes <- na.omit(ranked_genes)
    
    # GO Enrichment Analysis
    go_results <- enrichGO(gene = genes_df$ENTREZ[!is.na(genes_df$ENTREZ)],
                          OrgDb = org.Hs.eg.db,
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05)
    
    # GSEA
    gsea_results <- gseGO(geneList = ranked_genes,
                         OrgDb = org.Hs.eg.db,
                         ont = "BP",
                         pvalueCutoff = 0.05)
    
    # Save results
    pdf(file.path(pathway_dir, paste0(cell_type, "_pathway_analysis.pdf")))
    print(dotplot(go_results, showCategory=20, 
                 title=paste0(cell_type, " GO Enrichment Analysis")))
    print(gseaplot2(gsea_results, geneSetID=1:3, 
                    title=paste0(cell_type, " Top 3 GSEA Pathways")))
    dev.off()
    
    # Save tables
    write.csv(as.data.frame(go_results), 
              file.path(pathway_dir, paste0(cell_type, "_GO_enrichment.csv")))
    write.csv(as.data.frame(gsea_results), 
              file.path(pathway_dir, paste0(cell_type, "_GSEA_results.csv")))
}

# Function for enhanced visualization
create_enhanced_plots <- function(deg_results, cell_type) {
    # Volcano plot
    pdf(file.path(viz_dir, paste0(cell_type, "_volcano_plot.pdf")))
    print(EnhancedVolcano(deg_results,
                         lab = deg_results$gene_id,
                         x = 'log2FoldChange',
                         y = 'padj',
                         title = paste(cell_type, "Differential Expression"),
                         pCutoff = 0.05,
                         FCcutoff = 1,
                         pointSize = 3.0,
                         labSize = 6.0))
    dev.off()
    
    # Top DEGs heatmap
    top_genes <- deg_results %>%
        filter(padj < 0.05) %>%
        top_n(50, abs(log2FoldChange)) %>%
        pull(gene_id)
    
    # Get normalized counts for top genes
    cell_type_samples <- sample_info$cell_type == cell_type
    counts_subset <- normalized_counts[top_genes, cell_type_samples]
    
    # Scale the data
    scaled_counts <- t(scale(t(counts_subset)))
    
    # Create column annotation
    column_ha <- HeatmapAnnotation(
        Condition = sample_info$condition[cell_type_samples]
    )
    
    pdf(file.path(viz_dir, paste0(cell_type, "_top_DEGs_heatmap.pdf")))
    print(Heatmap(scaled_counts,
                  name = "Z-score",
                  show_row_names = TRUE,
                  cluster_columns = TRUE,
                  top_annotation = column_ha))
    dev.off()
}

# Run analyses for each cell type
for (cell_type in c("NPC", "Neuron")) {
    results_df <- if(cell_type == "NPC") npc_results else neuron_results
    
    # Perform pathway analysis
    perform_pathway_analysis(results_df, cell_type)
    
    # Create enhanced visualizations
    create_enhanced_plots(results_df, cell_type)
} 