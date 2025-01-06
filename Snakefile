import os
from pathlib import Path

# Configuration
SAMPLES = ["hPG1", "hPG2", "hPG3", "hPM1", "hPM2", "hPM3",
          "hNG1", "hNG2", "hNG3", "hNM1", "hNM2", "hNM3"]

# Paths
FASTQ_DIR = "90-1115332406/00_fastq"
QC_DIR = "90-1115332406/01_qc"
TRIM_DIR = "90-1115332406/02_trimmed"
ALIGN_DIR = "90-1115332406/03_alignment"
COUNTS_DIR = "90-1115332406/04_counts"
DESEQ_DIR = "90-1115332406/05_deseq2"
RESULTS_DIR = "90-1115332406/06_results"

# New output directories
PATHWAY_DIR = f"{RESULTS_DIR}/pathway_analysis"
ENHANCED_VIZ_DIR = f"{RESULTS_DIR}/enhanced_visualization"
QC_PLOTS_DIR = f"{RESULTS_DIR}/qc"

# Reference genome info
GENOME_DIR = "reference"
GENOME_FASTA = "reference/GRCh38.primary_assembly.genome.fa"
GENOME_GTF = "reference/gencode.v38.annotation.gtf"

# Load config
configfile: "config/config.yaml"

rule all:
    input:
        # Core analysis outputs
        f"{GENOME_DIR}/star_index",
        expand(f"{QC_DIR}/{{sample}}_R{{read}}_001_fastqc.html", sample=SAMPLES, read=[1,2]),
        expand(f"{QC_DIR}/{{sample}}_R{{read}}_001_fastqc.zip", sample=SAMPLES, read=[1,2]),
        f"{COUNTS_DIR}/all_samples_counts.txt",
        # DESeq2 core outputs
        f"{DESEQ_DIR}/NPC_differential_expression.csv",
        f"{DESEQ_DIR}/Neuron_differential_expression.csv",
        f"{RESULTS_DIR}/MA_plots.pdf",
        f"{RESULTS_DIR}/PCA_plot.pdf",
        f"{QC_PLOTS_DIR}/sample_correlation_heatmap.pdf",
        f"{QC_PLOTS_DIR}/count_distributions.pdf",
        f"{QC_PLOTS_DIR}/detected_genes.pdf",
        f"{RESULTS_DIR}/session_info.txt",
        # Advanced analysis outputs
        expand(f"{PATHWAY_DIR}/{{celltype}}_pathway_analysis.pdf", celltype=["NPC", "Neuron"]),
        expand(f"{PATHWAY_DIR}/{{celltype}}_GO_enrichment.csv", celltype=["NPC", "Neuron"]),
        expand(f"{PATHWAY_DIR}/{{celltype}}_GSEA_results.csv", celltype=["NPC", "Neuron"]),
        expand(f"{ENHANCED_VIZ_DIR}/{{celltype}}_volcano_plot.pdf", celltype=["NPC", "Neuron"]),
        expand(f"{ENHANCED_VIZ_DIR}/{{celltype}}_top_DEGs_heatmap.pdf", celltype=["NPC", "Neuron"])

rule fastqc:
    input:
        r1 = f"{FASTQ_DIR}/{{sample}}_R1_001.fastq.gz",
        r2 = f"{FASTQ_DIR}/{{sample}}_R2_001.fastq.gz"
    output:
        html_r1 = f"{QC_DIR}/{{sample}}_R1_001_fastqc.html",
        html_r2 = f"{QC_DIR}/{{sample}}_R2_001_fastqc.html",
        zip_r1 = f"{QC_DIR}/{{sample}}_R1_001_fastqc.zip",
        zip_r2 = f"{QC_DIR}/{{sample}}_R2_001_fastqc.zip"
    resources:
        mem_mb=16000,
        runtime=120
    shell:
        """
        # Clean up any existing files
        rm -f {QC_DIR}/{wildcards.sample}_R*_fastqc.* || true
        
        # Create output directory
        mkdir -p {QC_DIR}
        
        # Run FastQC with increased memory
        fastqc -t 2 --memory {resources.mem_mb}M --outdir {QC_DIR} {input.r1} {input.r2}
        
        # Verify outputs exist
        for f in {output}; do
            if [ ! -f "$f" ]; then
                echo "Error: Expected output file $f was not created"
                exit 1
            fi
        done
        """

rule trim_reads:
    input:
        r1 = f"{FASTQ_DIR}/{{sample}}_R1_001.fastq.gz",
        r2 = f"{FASTQ_DIR}/{{sample}}_R2_001.fastq.gz"
    output:
        r1 = f"{TRIM_DIR}/{{sample}}_R1_trimmed.fastq.gz",
        r2 = f"{TRIM_DIR}/{{sample}}_R2_trimmed.fastq.gz",
        r1_unpaired = f"{TRIM_DIR}/{{sample}}_R1_unpaired.fastq.gz",
        r2_unpaired = f"{TRIM_DIR}/{{sample}}_R2_unpaired.fastq.gz"
    resources:
        mem_mb=32000,
        runtime=240
    shell:
        """
        # Create output directory
        mkdir -p {TRIM_DIR}
        
        # Get adapter path
        ADAPTER_PATH=$(dirname $(which trimmomatic))/adapters/TruSeq3-PE.fa
        
        trimmomatic PE \
            {input.r1} {input.r2} \
            {output.r1} {output.r1_unpaired} \
            {output.r2} {output.r2_unpaired} \
            ILLUMINACLIP:$ADAPTER_PATH:2:30:10 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """

rule star_index:
    input:
        fasta = GENOME_FASTA,
        gtf = GENOME_GTF
    output:
        directory(f"{GENOME_DIR}/star_index")
    threads: 8
    resources:
        mem_mb=64000,
        runtime=240
    shell:
        """
        mkdir -p {output}
        STAR --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang 100
        """

rule star_align:
    input:
        r1 = f"{TRIM_DIR}/{{sample}}_R1_trimmed.fastq.gz",
        r2 = f"{TRIM_DIR}/{{sample}}_R2_trimmed.fastq.gz",
        index = f"{GENOME_DIR}/star_index"
    output:
        bam = f"{ALIGN_DIR}/{{sample}}/Aligned.sortedByCoord.out.bam"
    threads: 8
    resources:
        mem_mb=64000,
        runtime=240
    shell:
        """
        mkdir -p {ALIGN_DIR}/{wildcards.sample}

        STAR --runThreadN {threads} \
            --genomeDir {input.index} \
            --readFilesIn {input.r1} {input.r2} \
            --readFilesCommand zcat \
            --outFileNamePrefix {ALIGN_DIR}/{wildcards.sample}/ \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts \
            --limitBAMsortRAM {resources.mem_mb}000000
        """

rule featureCounts:
    input:
        bams = expand(f"{ALIGN_DIR}/{{sample}}/Aligned.sortedByCoord.out.bam", sample=SAMPLES)
    output:
        counts = f"{COUNTS_DIR}/all_samples_counts.txt"
    threads: 8
    resources:
        mem_mb=32000,
        runtime=60
    shell:
        """
        featureCounts -T {threads} \
            -p -t exon -g gene_id \
            -a {GENOME_GTF} \
            -o {output.counts} \
            {input.bams}
        """

rule deseq2_core:
    input:
        counts = f"{COUNTS_DIR}/all_samples_counts.txt"
    output:
        # Core outputs
        npc_de = f"{DESEQ_DIR}/NPC_differential_expression.csv",
        neuron_de = f"{DESEQ_DIR}/Neuron_differential_expression.csv",
        ma_plots = f"{RESULTS_DIR}/MA_plots.pdf",
        pca_plot = f"{RESULTS_DIR}/PCA_plot.pdf",
        # QC outputs
        sample_corr = f"{QC_PLOTS_DIR}/sample_correlation_heatmap.pdf",
        count_dist = f"{QC_PLOTS_DIR}/count_distributions.pdf",
        detected_genes = f"{QC_PLOTS_DIR}/detected_genes.pdf",
        # Session info
        session_info = f"{RESULTS_DIR}/session_info.txt"
    resources:
        mem_mb=32000,
        runtime=60
    script:
        "scripts/run_deseq2_core.R"

rule advanced_analysis:
    input:
        npc_de = f"{DESEQ_DIR}/NPC_differential_expression.csv",
        neuron_de = f"{DESEQ_DIR}/Neuron_differential_expression.csv",
        counts = f"{COUNTS_DIR}/all_samples_counts.txt"
    output:
        # Pathway analysis outputs
        pathway_plots = expand(f"{PATHWAY_DIR}/{{celltype}}_pathway_analysis.pdf", celltype=["NPC", "Neuron"]),
        go_results = expand(f"{PATHWAY_DIR}/{{celltype}}_GO_enrichment.csv", celltype=["NPC", "Neuron"]),
        gsea_results = expand(f"{PATHWAY_DIR}/{{celltype}}_GSEA_results.csv", celltype=["NPC", "Neuron"]),
        # Enhanced visualization outputs
        volcano_plots = expand(f"{ENHANCED_VIZ_DIR}/{{celltype}}_volcano_plot.pdf", celltype=["NPC", "Neuron"]),
        deg_heatmaps = expand(f"{ENHANCED_VIZ_DIR}/{{celltype}}_top_DEGs_heatmap.pdf", celltype=["NPC", "Neuron"])
    resources:
        mem_mb=32000,
        runtime=60
    script:
        "scripts/run_advanced_analysis.R" 