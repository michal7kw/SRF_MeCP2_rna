library(DESeq2)
library(tidyverse)
library(pheatmap)
library(org.Hs.eg.db)
library(ggrepel)  # For non-overlapping text labels
library(dplyr)

# Add at the beginning of the script, after library imports
log_file <- "90-1115332406/06_results/deseq2_analysis.log"
dir.create(dirname(log_file), recursive=TRUE, showWarnings=FALSE)

# Function to log messages to both console and file
log_message <- function(...) {
    msg <- paste(...)
    # Print to console
    print(msg)
    # Append to log file
    cat(paste0(msg, "\n"), file=log_file, append=TRUE)
}

# Read counts data
counts_data <- read.table("90-1115332406/04_counts/all_samples_counts.txt", 
                         header=TRUE, row.names=1, skip=1)

# Print column names for debugging
log_message("Original column names:")
log_message(paste(colnames(counts_data), collapse=", "))

# Keep only the count columns
count_cols <- grep("X90\\.1115332406\\.03_alignment\\.h[NP][GM][123]\\.Aligned\\.sortedByCoord\\.out\\.bam$", 
                  colnames(counts_data), value=TRUE)

counts_data <- counts_data[, count_cols]

# Clean up sample names
colnames(counts_data) <- sub("X90\\.1115332406\\.03_alignment\\.(h[NP][GM][123])\\..*", "\\1", 
                            colnames(counts_data))

# Create sample information
sample_info <- data.frame(
    sample = colnames(counts_data)
) %>%
    mutate(
        cell_type = factor(case_when(
            grepl("^hP", sample) ~ "NPC",
            grepl("^hN", sample) ~ "Neuron",
            TRUE ~ NA_character_
        )),
        condition = factor(case_when(
            grepl("G\\d$", sample) ~ "G",
            grepl("M\\d$", sample) ~ "M",
            TRUE ~ NA_character_
        ), levels=c("G", "M"))
    )

# Add verification step with more detailed output
verify_samples <- function(sample_info) {
    log_message("\nVerifying sample grouping:")
    print(sample_info)
    
    log_message("\nSample counts by group:")
    print(table(sample_info$cell_type, sample_info$condition))
    
    # Verify each combination has exactly 3 replicates
    for(ct in c("NPC", "Neuron")) {
        for(cond in c("G", "M")) {
            n_samples <- sum(sample_info$cell_type == ct & sample_info$condition == cond)
            if(n_samples != 3) {
                stop(sprintf("Error: %s-%s has %d samples (expected 3)", ct, cond, n_samples))
            }
        }
    }
    log_message("All sample groups have correct number of replicates")
}

verify_samples(sample_info)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
    countData = counts_data,
    colData = sample_info,
    design = ~ cell_type + condition + cell_type:condition
)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Function to get results with more detailed output
get_de_results <- function(dds, cell_type_value, output_file) {
    log_message(paste("\n=== Analyzing", cell_type_value, "==="))
    
    # Subset data for cell type
    cell_type_samples <- dds$cell_type == cell_type_value
    dds_subset <- dds[, cell_type_samples]
    
    # Reset the design for the subset to just look at condition effect
    design(dds_subset) <- ~ condition
    
    # Rerun DESeq on subset
    dds_subset <- DESeq(dds_subset)
    
    # Get results comparing Mecp2 vs GFP
    res <- results(dds_subset, contrast=c("condition", "M", "G"), alpha=0.05)
    
    # Print summary of results
    print(paste("\nSummary of differential expression for", cell_type_value))
    print(summary(res))
    
    # Convert to dataframe and add gene symbols
    res_df <- as.data.frame(res) %>%
        rownames_to_column("gene_id") %>%
        mutate(
            gene_id_clean = sub("\\..*", "", gene_id),
            gene_symbol = mapIds(org.Hs.eg.db,
                               keys = gene_id_clean,
                               keytype = "ENSEMBL",
                               column = "SYMBOL",
                               multiVals = "first")
        ) %>%
        # Add regulation status
        mutate(
            regulation = case_when(
                padj < 0.05 & log2FoldChange >= 0.5 ~ "Upregulated",
                padj < 0.05 & log2FoldChange <= -0.5 ~ "Downregulated",
                TRUE ~ "Not Significant"
            )
        ) %>%
        arrange(padj)
    
    # Print counts of differential expression
    log_message("\nDifferential Expression Counts:")
    log_message(capture.output(table(res_df$regulation)))
    
    # Save full results
    write.csv(res_df, output_file, row.names=FALSE)
    
    # Create and save lists of up/down regulated genes
    up_genes <- res_df %>%
        filter(regulation == "Upregulated") %>%
        dplyr::select(gene_symbol, gene_id, log2FoldChange, padj) %>%
        arrange(desc(log2FoldChange))
    
    down_genes <- res_df %>%
        filter(regulation == "Downregulated") %>%
        dplyr::select(gene_symbol, gene_id, log2FoldChange, padj) %>%
        arrange(log2FoldChange)
    
    # Save separate files for up/down regulated genes
    write.csv(up_genes, 
              sub("\\.csv$", "_upregulated.csv", output_file), 
              row.names=FALSE)
    write.csv(down_genes, 
              sub("\\.csv$", "_downregulated.csv", output_file), 
              row.names=FALSE)
    
    return(res_df)
}

# Enhanced volcano plot function
create_volcano_plot <- function(res_df, cell_type, output_file) {
    # Create custom color scheme
    colors <- c(
        "Upregulated" = "#FF4B4B",  # Red for upregulated
        "Downregulated" = "#4B4BFF", # Blue for downregulated
        "Not Significant" = "grey"
    )
    
    # Create volcano plot
    p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
        geom_point(aes(color = regulation), alpha = 0.6, size = 1) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
        scale_color_manual(values = colors) +
        theme_bw() +
        theme(
            legend.position = "right",
            plot.title = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 12),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 10)
        ) +
        labs(
            title = sprintf("%s: Mecp2 vs Control", cell_type),
            x = "log2 Fold Change",
            y = "-log10 Adjusted P-value",
            color = "Regulation Status"
        )
    
    # Add labels for top genes
    top_genes <- res_df %>%
        filter(regulation != "Not Significant") %>%
        group_by(regulation) %>%
        dplyr::top_n(10, abs(log2FoldChange)) %>%
        ungroup()
    
    p <- p + geom_text_repel(
        data = top_genes,
        aes(label = gene_symbol),
        size = 3,
        max.overlaps = 20,
        box.padding = 0.5
    )
    
    # Save plot
    ggsave(output_file, p, width = 10, height = 8)
    
    return(p)
}

# Create output directories
dir.create("90-1115332406/05_deseq2", showWarnings=FALSE, recursive=TRUE)
dir.create("90-1115332406/06_results", showWarnings=FALSE, recursive=TRUE)

# Get results for each comparison
npc_results <- get_de_results(dds, "NPC", "90-1115332406/05_deseq2/NPC_differential_expression.csv")
neuron_results <- get_de_results(dds, "Neuron", "90-1115332406/05_deseq2/Neuron_differential_expression.csv")

# Create QC plots
qc_dir <- "90-1115332406/06_results/qc"
dir.create(qc_dir, recursive=TRUE, showWarnings=FALSE)

# Sample correlation heatmap
vst_counts <- vst(dds, blind=FALSE)
sample_cors <- cor(assay(vst_counts))

# Get gene symbols for all genes - add error checking
gene_symbols <- tryCatch({
    mapIds(org.Hs.eg.db,
           keys = sub("\\..*", "", rownames(vst_counts)),
           keytype = "ENSEMBL",
           column = "SYMBOL",
           multiVals = "first")
}, error = function(e) {
    warning("Error mapping gene symbols: ", e$message)
    return(rownames(vst_counts))
})

# Replace NAs with original IDs
gene_symbols[is.na(gene_symbols)] <- rownames(vst_counts)[is.na(gene_symbols)]

# Update row names with gene symbols
rownames(vst_counts) <- gene_symbols

pdf(file.path(qc_dir, "sample_correlation_heatmap.pdf"))
pheatmap(sample_cors,
         main="Sample Correlation Heatmap",
         annotation_col=data.frame(
             CellType=dds$cell_type,
             Condition=dds$condition,
             row.names=colnames(dds)
         ))
dev.off()

# Count distribution plots - use gene symbols
counts_df <- as.data.frame(counts(dds, normalized=TRUE)) %>%
    rownames_to_column("gene_id") %>%
    mutate(
        gene_symbol = tryCatch({
            mapIds(org.Hs.eg.db,
                   keys = sub("\\..*", "", gene_id),
                   keytype = "ENSEMBL",
                   column = "SYMBOL",
                   multiVals = "first")
        }, error = function(e) {
            warning("Error mapping gene symbols: ", e$message)
            return(gene_id)
        }),
        gene_name = coalesce(gene_symbol, gene_id)
    ) %>%
    gather(key="sample", value="counts", -gene_id, -gene_symbol, -gene_name)

pdf(file.path(qc_dir, "count_distributions.pdf"))
ggplot(counts_df, aes(x=log2(counts + 1))) +
    geom_density() +
    facet_wrap(~sample) +
    theme_bw() +
    ggtitle("Distribution of Normalized Counts")
dev.off()

# Detected genes plot
detected_genes <- colSums(counts(dds) > 0)
pdf(file.path(qc_dir, "detected_genes.pdf"))
barplot(detected_genes, main="Number of Detected Genes per Sample", las=2)
dev.off()

# MA plots with gene symbols
pdf("90-1115332406/06_results/MA_plots.pdf")
par(mfrow=c(2,1))

# Function to add gene symbols to MA plot
plot_ma_with_symbols <- function(res, main) {
    # Plot basic MA plot first
    plotMA(res, main=main)
    
    # Get significant genes
    sig_idx <- which(res$padj < 0.05 & abs(res$log2FoldChange) > 1)
    
    # Only try to add labels if there are significant genes
    if (length(sig_idx) > 0) {
        gene_ids <- rownames(res)[sig_idx]
        
        # Check if we have valid gene IDs
        if (length(gene_ids) > 0) {
            gene_symbols <- mapIds(org.Hs.eg.db,
                                 keys = sub("\\..*", "", gene_ids),
                                 keytype = "ENSEMBL",
                                 column = "SYMBOL",
                                 multiVals = "first")
            
            # Replace NA symbols with original IDs
            gene_symbols[is.na(gene_symbols)] <- gene_ids[is.na(gene_symbols)]
            
            # Add labels
            with(res[sig_idx, ], 
                 text(log2FoldChange, log2(baseMean), 
                      labels=gene_symbols, 
                      cex=0.8, pos=4))
        }
    }
}

# NPC comparison
res_npc <- results(dds[, dds$cell_type == "NPC"], 
                  contrast=c("condition", "M", "G"))
plot_ma_with_symbols(res_npc, "NPC: MeCP2 vs Control")

# Neuron comparison
res_neuron <- results(dds[, dds$cell_type == "Neuron"], 
                     contrast=c("condition", "M", "G"))
plot_ma_with_symbols(res_neuron, "Neuron: MeCP2 vs Control")

dev.off()

# PCA plot with gene symbols
pdf("90-1115332406/06_results/PCA_plot.pdf", width=10, height=8)
pcaData <- plotPCA(vst_counts, intgroup=c("cell_type", "condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x=PC1, y=PC2, color=cell_type, shape=condition)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle("PCA Plot of Samples") +
    theme_bw()
dev.off()

# Save session info
writeLines(capture.output(sessionInfo()), 
           "90-1115332406/06_results/session_info.txt") 

# Add after getting results
for (cell_type in c("NPC", "Neuron")) {
    results_file <- sprintf("90-1115332406/05_deseq2/%s_differential_expression.csv", cell_type)
    volcano_file <- sprintf("90-1115332406/06_results/%s_volcano_plot.pdf", cell_type)
    
    # Get results
    res_df <- get_de_results(dds, cell_type, results_file)
    
    # Create volcano plot
    create_volcano_plot(res_df, cell_type, volcano_file)
} 