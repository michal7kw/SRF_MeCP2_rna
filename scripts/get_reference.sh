#!/bin/bash

# Create reference directory
mkdir -p reference

cd reference

# Download genome
if [ ! -f "GRCh38.primary_assembly.genome.fa" ]; then
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz
    gunzip GRCh38.primary_assembly.genome.fa.gz
fi

# Download annotation
if [ ! -f "gencode.v38.annotation.gtf" ]; then
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
    gunzip gencode.v38.annotation.gtf.gz
fi 