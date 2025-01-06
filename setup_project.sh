#!/bin/bash

# Create project directory structure
mkdir -p 90-1115332406/{01_qc,02_trimmed,03_alignment,04_counts,05_deseq2,06_results}
mkdir -p scripts

# Make sure the R script is executable
chmod +x scripts/run_deseq2.R 