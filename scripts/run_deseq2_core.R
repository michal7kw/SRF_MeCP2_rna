library(DESeq2)
library(tidyverse)
library(pheatmap)

# Read counts data
counts_data <- read.table("90-1115332406/04_counts/all_samples_counts.txt", 
                         header=TRUE, row.names=1, skip=1)

# Print column names for debugging
print("Original column names:")
print(colnames(counts_data))

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
        group = substr(sample, 2, 4),
        condition = factor(case_when(
            grepl("G", substr(sample, 3, 4)) ~ "G",
            grepl("M", substr(sample, 3, 4)) ~ "M"
        ), levels=c("G", "M")),
        cell_type = factor(ifelse(substr(sample, 2, 2) == "P", "NPC", "Neuron"))
    )

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
    countData = counts_data,
    colData = sample_info,
    design = ~ cell_type + condition
)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Function to get results
get_de_results <- function(dds, cell_type_value, output_file) {
    dir.create(dirname(output_file), showWarnings=FALSE, recursive=TRUE)
    
    cell_type_samples <- dds$cell_type == cell_type_value
    dds_subset <- dds[, cell_type_samples]
    
    res <- results(dds_subset, contrast=c("condition", "M", "G"), alpha=0.05)
    
    res_df <- as.data.frame(res) %>%
        rownames_to_column("gene_id") %>%
        arrange(padj)
    
    write.csv(res_df, output_file, row.names=FALSE)
    return(res_df)
}

# Create output directories
dir.create("90-1115332406/05_deseq2", showWarnings=FALSE, recursive=TRUE)
dir.create("90-1115332406/06_results", showWarnings=FALSE, recursive=TRUE)

# Get results for each comparison
npc_results <- get_de_results(dds, "NPC", "90-1115332406/05_deseq2/NPC_differential_expression.csv")
neuron_results <- get_de_results(dds, "Neuron", "90-1115332406/05_deseq2/Neuron_differential_expression.csv")

# Create MA plots
pdf("90-1115332406/06_results/MA_plots.pdf")
par(mfrow=c(2,1))
plotMA(results(dds[, dds$cell_type == "NPC"], contrast=c("condition", "M", "G")), 
       main="NPC: MeCP2 vs Control")
plotMA(results(dds[, dds$cell_type == "Neuron"], contrast=c("condition", "M", "G")), 
       main="Neuron: MeCP2 vs Control")
dev.off()

# Create QC plots
qc_dir <- "90-1115332406/06_results/qc"
dir.create(qc_dir, recursive=TRUE, showWarnings=FALSE)

# Sample correlation heatmap
vst_counts <- vst(dds, blind=FALSE)
sample_cors <- cor(assay(vst_counts))
pdf(file.path(qc_dir, "sample_correlation_heatmap.pdf"))
pheatmap(sample_cors,
         main="Sample Correlation Heatmap",
         annotation_col=data.frame(
             CellType=dds$cell_type,
             Condition=dds$condition,
             row.names=colnames(dds)
         ))
dev.off()

# Count distribution plots
counts_df <- as.data.frame(counts(dds, normalized=TRUE)) %>%
    rownames_to_column("gene") %>%
    gather(key="sample", value="counts", -gene)

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

# PCA plot
pdf("90-1115332406/06_results/PCA_plot.pdf")
plotPCA(vst_counts, intgroup=c("cell_type", "condition"))
dev.off()

# Save session info
writeLines(capture.output(sessionInfo()), 
           "90-1115332406/06_results/session_info.txt") 