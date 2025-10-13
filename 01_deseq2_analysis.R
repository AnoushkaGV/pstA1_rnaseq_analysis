# DESeq2 differential expression analysis

library(DESeq2)

source("config.R")

# Load data
cc_matrix <- read.csv(file.path(data_dir, "cc_matrix.csv"), header=T)
rownames(cc_matrix) <- cc_matrix$X
cc_matrix$X <- NULL

# Load gene names dictionary
names_dict <- read.table(file.path(data_dir, "rv_dict.csv"), sep=",", header=T)

# Sample metadata
sample_info <- data.frame(
  row.names = colnames(cc_matrix),
  condition = factor(c(rep("pstA1_ctrl", 3), rep("pstA1_RIF_30min", 3), rep("pstA1_RIF_60min", 3),
                       rep("WT_ctrl", 3), rep("WT_RIF_30min", 3), rep("WT_RIF_60min", 3)))
)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = cc_matrix, colData = sample_info, design = ~ condition)
dds <- DESeq(dds)

# Define all compariP0+r\P0+r\P0+r\P0+r\P0+r\P0+r\P0+r\sons
comparisons <- list(
  "WT_30min_vs_ctrl" = c("condition", "WT_RIF_30min", "WT_ctrl"),
  "WT_60min_vs_ctrl" = c("condition", "WT_RIF_60min", "WT_ctrl"),
  "pstA1_30min_vs_ctrl" = c("condition", "pstA1_RIF_30min", "pstA1_ctrl"),
  "pstA1_60min_vs_ctrl" = c("condition", "pstA1_RIF_60min", "pstA1_ctrl"),
  "pstA1_vs_WT_ctrl" = c("condition", "pstA1_ctrl", "WT_ctrl"),
  "pstA1_vs_WT_30min" = c("condition", "pstA1_RIF_30min", "WT_RIF_30min"),
  "pstA1_vs_WT_60min" = c("condition", "pstA1_RIF_60min", "WT_RIF_60min")
)

# Function to perform DESeq2 analysis and get results
get_deseq_results <- function(comparison_name, contrast) {
  res <- results(dds, contrast = contrast)
  res_df <- as.data.frame(res)
  res_df$significant <- with(res_df, padj < threshold_padj & abs(log2FoldChange) > threshold_log2fc)
  res_df$comparison <- comparison_name
  return(res_df)
}

# Run all comparisons
results_list <- list()
for (comp_name in names(comparisons)) {
  cat("Running DESeq2 for:", comp_name, "\n")
  results_list[[comp_name]] <- get_deseq_results(comp_name, comparisons[[comp_name]])
  
  # Print summary
  res <- results_list[[comp_name]]
  cat("Total:", sum(res$padj < threshold_padj & abs(res$log2FoldChange) > threshold_log2fc, na.rm = TRUE), 
      "Up:", sum(res$log2FoldChange > threshold_log2fc & res$padj < threshold_padj, na.rm = TRUE), 
      "Down:", sum(res$log2FoldChange < -threshold_log2fc & res$padj < threshold_padj, na.rm = TRUE), "\n\n")
}

# Save objects for other scripts
save(dds, results_list, names_dict, comparisons, file = file.path(output_dir, "deseq2_results.RData"))

cat("DESeq2 analysis complete!\n")
