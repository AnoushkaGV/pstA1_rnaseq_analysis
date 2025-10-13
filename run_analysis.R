# Master script to run complete RNA-seq analysis pipeline

cat("Starting pstA1 RNA-seq analysis pipeline\n")
cat("==========================================\n\n")

start_time <- Sys.time()

# 1. DESeq2 analysis
cat("Step 1: Running DESeq2 differential expression analysis...\n")
source("01_deseq2_analysis.R")
cat("\n")

# 2. Generate plots
cat("Step 2: Generating visualizations...\n")
source("02_generate_plots.R")
cat("\n")

# 3. Pathway analysis
cat("Step 3: Running KEGG pathway enrichment analysis...\n")
cat("Note: This step may take 10-30 minutes due to KEGG API rate limiting\n")
source("03_pathway_analysis.R")
cat("\n")

# 4. Summary tables
cat("Step 4: Creating summary tables...\n")
source("04_summary_tables.R")
cat("\n")

# 5. Save session info
cat("Step 5: Saving session information...\n")
sink(file.path("results", "session_info.txt"))
cat("pstA1 RNA-seq Analysis\n")
cat("Analysis completed:", as.character(Sys.time()), "\n\n")
cat("===========================================\n")
cat("Session Information\n")
cat("===========================================\n\n")
print(sessionInfo())
sink()

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "mins")

cat("\n==========================================\n")
cat("Pipeline complete!\n")
cat("Total time:", round(elapsed, 2), "minutes\n")
cat("Results saved to: results/\n")
cat("Session info saved to: results/session_info.txt\n")
cat("==========================================\n")
