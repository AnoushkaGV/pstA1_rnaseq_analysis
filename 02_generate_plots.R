# Generate plots: MA plots, heatmap, and Venn diagrams

library(DESeq2)
library(pheatmap)
library(ggplot2)
library(ggvenn)

source("config.R")

# Load DESeq2 results
load(file.path(output_dir, "deseq2_results.RData"))

# 1. MA PLOTS
pdf(file.path(output_dir, "figures", "MA_plots.pdf"), width = 12, height = 10)
par(mfrow = c(3, 3), mar = c(5, 5, 4, 2), cex.main = 1.2, cex.lab = 1.3, cex.axis = 1.1)
for (i in 1:length(comparisons)) {
  comp_name <- names(comparisons)[i]
  res <- results(dds, contrast = comparisons[[comp_name]])
  plotMA(res, ylim = c(-6, 6), 
         main = comparison_labels[comp_name], 
         xlab = "Mean of Normalized Counts",
         ylab = "Log2 Fold Change",
         cex.main = 1.2,
         colSig = "red",
         colNonSig = "gray60",
         alpha = 0.05)
}
dev.off()

# 2. HEATMAP WITH AVERAGED REPLICATES
# Get all significant genes across all comparisons
all_significant <- unique(unlist(lapply(results_list, function(x) {
  rownames(x[!is.na(x$padj) & x$significant, ])
})))

cat("Total unique significant genes across all comparisons:", length(all_significant), "\n")

# Get normalized counts and subset
normalized_counts <- counts(dds, normalized = TRUE)
significant_counts <- normalized_counts[all_significant, , drop = FALSE]

# Log transform
significant_counts_log <- log2(significant_counts + 1)

# Average replicates and organize columns
pstA1_ctrl_cols <- grep("pstA1_ctrl", colnames(significant_counts_log))
pstA1_30min_cols <- grep("pstA1_RIF_30min", colnames(significant_counts_log))
pstA1_60min_cols <- grep("pstA1_RIF_60min", colnames(significant_counts_log))
wt_ctrl_cols <- grep("WT_ctrl", colnames(significant_counts_log))
wt_30min_cols <- grep("WT_RIF_30min", colnames(significant_counts_log))
wt_60min_cols <- grep("WT_RIF_60min", colnames(significant_counts_log))

# Create averaged data
averaged_data <- data.frame(
  "pstA1_Control" = rowMeans(significant_counts_log[, pstA1_ctrl_cols]),
  "pstA1_RIF_30min" = rowMeans(significant_counts_log[, pstA1_30min_cols]),
  "pstA1_RIF_60min" = rowMeans(significant_counts_log[, pstA1_60min_cols]),
  "WT_Control" = rowMeans(significant_counts_log[, wt_ctrl_cols]),
  "WT_RIF_30min" = rowMeans(significant_counts_log[, wt_30min_cols]),
  "WT_RIF_60min" = rowMeans(significant_counts_log[, wt_[2;2R60min_cols])
)

# Mean center the data
mean_centered_log2 <- averaged_data - rowMeans(averaged_data)

# Add gene names
gene_names_mapped <- names_dict$Name[match(rownames(mean_centered_log2), names_dict$Rv)]
gene_names_final <- ifelse(is.na(gene_names_mapped) | gene_names_mapped == "", 
                           rownames(mean_centered_log2), gene_names_mapped)
rownames(mean_centered_log2) <- gene_names_final

col_annotation_full <- data.frame(
  Strain = c(rep("pstA1", 3), rep("WT", 3)),
  Treatment = rep(c("Control", "RIF 30min", "RIF 60min"), 2)
)
rownames(col_annotation_full) <- colnames(mean_centered_log2)

annotation_colors <- list(Strain = strain_colors, Treatment = treatment_colors)

pdf(file.path(output_dir, "figures", "heatmap_all_DEGs.pdf"), width[3;1R = 10, height = 12)
pheatmap(mean_centered_log2,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = ifelse(nrow(mean_centered_log2) < 80, TRUE, FALSE),
         annotation_col = col_annotation_full,
         annotation_colors = annotation_colors,
         color = colorRampPalette(c("navy", "blue", "lightblue", "white", "pink", "red", "darkred"))(100),
         fontsize_row = 12,
         fontsize_col = 12,
         fontsize = 12,
         angle_col = 45,
         cellwidth = 40,
         cellheight = 10,
         breaks = seq(-6, 6, length.out = 101))
dev.off()

# 3. Venn diagrams
get_degs <- function(comparison) {
  res <- results_list[[comparison]]
  sig <- res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]
  up <- rownames(sig[sig$log2FoldChange > 1, ])
  down <- rownames(sig[sig$log2FoldChange < -1, ])
  return(list(up = up, down = down))
}

# pstA1 vs WT across timepoints
ctrl <- get_degs("pstA1_vs_WT_ctrl")
min30 <- get_degs("pstA1_vs_WT_30min")
min60 <- get_degs("pstA1_vs_WT_60min")

# RIF response at 60min
pstA1_rif <- get_degs("pstA1_60min_vs_ctrl")
wt_rif <- get_degs("WT_60min_vs_ctrl")

# Venn lists
up_venn_3way <- list(Ctrl = ctrl$up, `30min` = min30$up, `60min` = min60$up)
down_venn_3way <- list(Ctrl = ctrl$down, `30min` = min30$down, `60min` = min60$down)
up_venn_rif <- list(WT = wt_rif$up, pstA1 = pstA1_rif$up)
down_venn_rif <- list(WT = wt_rif$down, pstA1 = pstA1_rif$down)

# Plot all 4 Venns
pdf(file.path(output_dir, "figures", "venn_diagr[>41;2500;0c]10;rgb:dcaa/dcab/dcaa\]11;rgb:158e/193a/1e75\ams.pdf"), width = 10, height = 10)

print(ggvenn(up_venn_3way, fill_color = rep(NA, 3), stroke_size = 2, 
       set_name_size = 6, text_size = 6, show_percentage = FALSE) + 
  ggtitle("Upregulated DEGs: pstA1 vs WT") +
  theme(plot.title = element_text(size=16, face="bold", hjust=0.5)))

print(ggvenn(down_venn_3way, fill_color = rep(NA, 3), stroke_size = 2,
       set_name_size = 6, text_size = 6, show_percentage = FALSE) + 
  ggtitle("Downregulated DEGs: pstA1 vs WT") +
  theme(plot.title = element_text(size=16, face="bold", hjust=0.5)))

print(ggvenn(up_venn_rif, fill_color = rep(NA, 2), stroke_size = 2,
       set_name_size = 6, text_size = 6, show_percentage = FALSE) + 
  ggtitle("Upregulated DEGs: WT vs pstA1 (60min)") +
  theme(plot.title = element_text(size=16, face="bold", hjust=0.5)))

print(ggvenn(down_venn_rif, fill_color = rep(NA, 2), stroke_size = 2,
       set_name_size = 6, text_size = 6, show_percentage = FALSE) + 
  ggtitle("Downregulated DEGs: WT vs pstA1 (60min)") +
  theme(plot.title = element_text(size=16, face="bold", hjust=0.5)))

dev.off()

cat("All plots generated!\n")
