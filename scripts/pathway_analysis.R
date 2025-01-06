library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)

perform_pathway_analysis <- function(deg_results, cell_type, output_dir) {
    # Create output directory
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Prepare gene list
    genes_df <- deg_results %>%
        filter(!is.na(padj) & padj < 0.05) %>%
        arrange(desc(log2FoldChange))
    
    # Convert ENSEMBL IDs to ENTREZ IDs
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
    pdf(file.path(output_dir, paste0(cell_type, "_pathway_analysis.pdf")))
    
    # Plot GO enrichment
    print(dotplot(go_results, showCategory=20, title="GO Enrichment Analysis"))
    
    # Plot GSEA results
    print(gseaplot2(gsea_results, geneSetID=1:3, title="Top 3 GSEA Pathways"))
    
    dev.off()
    
    # Save tables
    write.csv(as.data.frame(go_results), 
              file.path(output_dir, paste0(cell_type, "_GO_enrichment.csv")))
    write.csv(as.data.frame(gsea_results), 
              file.path(output_dir, paste0(cell_type, "_GSEA_results.csv")))
} 