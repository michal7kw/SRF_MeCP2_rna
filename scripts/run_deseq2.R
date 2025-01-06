library(DESeq2)
library(tidyverse)

# Read counts data
counts_data <- read.table("90-1115332406/04_counts/all_samples_counts.txt", 
                         header=TRUE, row.names=1, skip=1)

# Print column names for debugging
print("Original column names:")
print(colnames(counts_data))

# Keep only the count columns which match our expected pattern
count_cols <- grep("X90\\.1115332406\\.03_alignment\\.h[NP][GM][123]\\.Aligned\\.sortedByCoord\\.out\\.bam$", 
                  colnames(counts_data), value=TRUE)

# Verify we found all expected columns
print("Found count columns:")
print(count_cols)
print(paste("Number of count columns:", length(count_cols)))

counts_data <- counts_data[, count_cols]

# Clean up sample names to match our expected format
colnames(counts_data) <- sub("X90\\.1115332406\\.03_alignment\\.(h[NP][GM][123])\\..*", "\\1", colnames(counts_data))

# Print cleaned column names
print("Column names after cleanup:")
print(colnames(counts_data))

# Verify sample names match our expected list
expected_samples <- c("hPG1", "hPG2", "hPG3", "hPM1", "hPM2", "hPM3",
                     "hNG1", "hNG2", "hNG3", "hNM1", "hNM2", "hNM3")

# Print comparison of actual vs expected samples
print("Expected samples:")
print(expected_samples)
print("Actual samples:")
print(sort(colnames(counts_data)))

# Verify all expected samples are present
missing_samples <- setdiff(expected_samples, colnames(counts_data))
if(length(missing_samples) > 0) {
    print("Counts data column names:")
    print(colnames(counts_data))
    print("Missing samples:")
    print(missing_samples)
    stop("Missing samples: ", paste(missing_samples, collapse=", "))
}

# Create sample information with correct order matching the counts data
sample_info <- data.frame(
    sample = colnames(counts_data)
) %>%
    mutate(
        group = substr(sample, 2, 4),  # Extract PG1, PM2, etc.
        condition = factor(case_when(
            grepl("G", substr(sample, 3, 4)) ~ "G",
            grepl("M", substr(sample, 3, 4)) ~ "M"
        ), levels=c("G", "M")),  # G or M
        cell_type = factor(ifelse(substr(sample, 2, 2) == "P", "NPC", "Neuron"))
    )

# Print sample information for debugging
print("Sample information:")
print(sample_info)

# Verify condition assignment
print("\nCondition table:")
print(table(sample_info$condition, sample_info$cell_type))

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
    countData = counts_data,
    colData = sample_info,
    design = ~ cell_type + condition
)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Function to get results for a specific comparison
get_de_results <- function(dds, cell_type_value, output_file) {
    # Create output directory if it doesn't exist
    dir.create(dirname(output_file), showWarnings=FALSE, recursive=TRUE)
    
    # Subset data for cell type
    cell_type_samples <- dds$cell_type == cell_type_value
    dds_subset <- dds[, cell_type_samples]
    
    # Get results (M vs G)
    res <- results(dds_subset, 
                  contrast=c("condition", "M", "G"),
                  alpha=0.05)
    
    # Convert to data frame and add gene names
    res_df <- as.data.frame(res) %>%
        rownames_to_column("gene_id") %>%
        arrange(padj)
    
    # Write results
    write.csv(res_df, output_file, row.names=FALSE)
    
    return(res_df)
}

# Create output directories
dir.create("90-1115332406/05_deseq2", showWarnings=FALSE, recursive=TRUE)
dir.create("90-1115332406/06_results", showWarnings=FALSE, recursive=TRUE)

# Get results for each comparison
npc_results <- get_de_results(dds, "NPC", 
                            "90-1115332406/05_deseq2/NPC_differential_expression.csv")
neuron_results <- get_de_results(dds, "Neuron", 
                                "90-1115332406/05_deseq2/Neuron_differential_expression.csv")

# Create MA plots
pdf("90-1115332406/06_results/MA_plots.pdf", width=10, height=12)
par(mfrow=c(2,1))

# NPC comparison
res_npc <- results(dds[, dds$cell_type == "NPC"], 
                  contrast=c("condition", "M", "G"))
plotMA(res_npc, main="NPC: MeCP2 vs Control")

# Neuron comparison
res_neuron <- results(dds[, dds$cell_type == "Neuron"], 
                     contrast=c("condition", "M", "G"))
plotMA(res_neuron, main="Neuron: MeCP2 vs Control")

dev.off()

# Print summary of results
print("Summary of NPC results:")
summary(res_npc)
print("\nSummary of Neuron results:")
summary(res_neuron) 