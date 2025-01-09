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

# snakemake --unlock

# Run snakemake
snakemake \
    --cores 16 \
    --rerun-incomplete \
    --keep-going \
    --latency-wait 120 \
    --restart-times 3 \
    --executor slurm \
    --verbose \
    --printshellcmds \
    --default-resources \
        slurm_account="kubacki.michal" \
        "slurm_partition='workq'" \
        mem_mb=32000 \
        runtime=720 \
    --jobs 12 \
    2>&1 | tee logs/snakemake_full.log
