# Create summary tables

library(KEGGREST)

source("config.R")

# Load DESeq2 results
load(file.path(output_dir, "deseq2_results.RData"))

# Summary table
summary_df <- data.frame(
  Comparison = comparison_labels[names(results_list)],
  Total_DEGs = sapply(results_list, function(x) sum(x$padj < threshold_padj & abs(x$log2FoldChange) > threshold_log2fc, na.rm = TRUE)),
  Upregulated = sapply(results_list, function(x) sum(x$log2FoldChange > threshold_log2fc & x$padj < threshold_padj, na.rm = TRUE)),
  Downregulated = sapply(results_list, function(x) sum(x$log2FoldChange < -threshold_log2fc & x$padj < threshold_padj, na.rm = TRUE))
)

print(summary_df)
write.csv(summary_df, file.path(output_dir, "tables", "DEG_summary.csv"), row.names = FALSE)

# Create supplementary table with all DEGs including pathway info
supp_table <- data.frame()

for (comp_name in names(results_list)) {
  res <- results_list[[comp_name]]
  sig_genes <- rownames(res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1, ])
  
  if (length(sig_genes) > 0) {
    res_sub <- res[sig_genes, ]
    
    # Get pathway info for each gene
    cat("Getting pathway info for", comp_name, "-", length(sig_genes), "genes\n")
    pathways <- sapply(sig_genes, function(gene) {
      tryCatch({
        Sys.sleep(0.2)
        pw <- keggGet(paste0(kegg_organism, ":", gene))[[1]]$PATHWAY
        if (is.null(pw)) return(NA)
        paste(names(pw), collapse = "; ")
      }, error = function(e) { NA })
    })
    
    temp_df <- data.frame(
      Comparison = comparison_labels[comp_name],
      Gene = sig_genes,
      Gene_Name = names_dict$Name[match(sig_genes, names_dict$Rv)],
      log2FoldChange = res_sub$log2FoldChange,
      padj = res_sub$padj,
      baseMean = res_sub$baseMean,
      Pathway = pathways,
      stringsAsFactors = FALSE
    )
    supp_table <- rbind(supp_table, temp_df)
  }
}

supp_table <- supp_table[order(supp_table$Comparison, supp_table$padj), ]
write.csv(supp_table, file.path(output_dir, "tables", "Supplementary_Table_All_DEGs.csv"), row.names = FALSE)

cat("Summary tables complete!\n")
