# Configuration file for RNA-seq analysis

# Paths
data_dir <- "data"
output_dir <- "results"

# Thresholds
threshold_log2fc <- 1
threshold_padj <- 0.05

# KEGG settings
kegg_organism <- "mtu"
kegg_delay <- 0.5
top_pathways <- 5

# Colors
strain_colors <- c("pstA1" = "#1F78B4", "WT" = "#E31A1C")
treatment_colors <- c("Control" = "#33A02C", "RIF 30min" = "#FF7F00", "RIF 60min" = "#6A3D9A")
plot_colors <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F")

# Comparison labels
comparison_labels <- c(
  "WT_30min_vs_ctrl" = "WT 30min vs Control",
  "WT_60min_vs_ctrl" = "WT 60min vs Control", 
  "pstA1_30min_vs_ctrl" = "pstA1 30min vs Control",
  "pstA1_60min_vs_ctrl" = "pstA1 60min vs Control",
  "pstA1_vs_WT_ctrl" = "pstA1 vs WT (Control)",
  "pstA1_vs_WT_30min" = "pstA1 vs WT (30min)",
  "pstA1_vs_WT_60min" = "pstA1 vs WT (60min)"
)

# Create output directories
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
