#!/bin/bash
#SBATCH --job-name=snakemake
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --partition=workq
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/snakemake.err"
#SBATCH --output="logs/snakemake.out"

# Create logs directory if it doesn't exist
mkdir -p logs

# Set working directory
WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_rna"
cd $WORKING_DIR || exit 1

# # Install required R packages
# Rscript -e '
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install(c(
#     "DESeq2",
#     "clusterProfiler",
#     "org.Hs.eg.db",
#     "DOSE",
#     "enrichplot",
#     "GO.db",
#     "EnhancedVolcano",
#     "ComplexHeatmap"
# ))
# install.packages(c("tidyverse", "pheatmap"))
# '

# Download reference files if needed
# bash scripts/get_reference.sh

snakemake --unlock

# Run snakemake
snakemake \
    --cores 16 \
    --rerun-incomplete \
    --keep-going \
    --latency-wait 120 \
    --restart-times 3 \
    --executor slurm \
    --default-resources \
        slurm_account="kubacki.michal" \
        "slurm_partition='workq'" \
        mem_mb=32000 \
        runtime=720 \
    --jobs 12
