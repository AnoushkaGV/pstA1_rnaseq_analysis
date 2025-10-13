# KEGG pathway enrichment analysis

library(KEGGREST)
library(ggplot2)
library(tidyr)
library(dplyr)

source("config.R")

# Load DESeq2 results
load(file.path(output_dir, "deseq2_results.RData"))

# Pathway analysis function
perform_pathway_analysis <- function(significant_genes, comparison_name) {
  if (length(significant_genes) == 0) {
    cat("No significant genes for", comparison_name, "\n")
    return(NULL)
  }
  
  genes <- paste(kegg_organism, ":", significant_genes, sep = "")
  
  # Get pathways for each gene
  get_pathways <- function(gene) {
    tryCatch({
      Sys.sleep(kegg_delay)  # Rate limiting
      pathways <- keggGet(gene)[[1]]$PATHWAY
      return(pathways)
    }, error = function(e) {
      return(NULL)
    })
  }
  
  pathways_list <- lapply(genes, get_pathways)
  pathways_df <- data.frame(Gene = genes, Pathways = I(pathways_list))
  
  # Count pathways
  pathways_long <- pathways_df %>%
    unnest(Pathways) %>%
    filter(!is.null(Pathways)) %>%
    group_by(Pathways) %>%
    summarise(Count = n(), .groups = 'drop') %>%
    arrange(desc(Count)) %>%
    head(top_pathways)
  
  return(pathways_long)
}

# Run pathway analysis for all comparisons
pathway_results <- list()

for (comp_name in names(results_list)) {
  sig_genes <- rownames(results_list[[comp_name]][!is.na(results_list[[comp_name]]$padj) & 
                                                    results_list[[comp_name]]$significant, ])
  
  cat("Running pathway analysis for", comp_name, "with", length(sig_genes), "genes\n")
  
  pathways_result <- perform_pathway_analysis(sig_genes, comp_name)
  pathway_results[[comp_name]] <- pathways_result
}

# Pathway analysis plots
pathway_plots <- list()

for (i in 1:length(pathway_results)) {
  comp_name <- names(pathway_results)[i]
  pathways_result <- pathway_results[[comp_name]]
  
  if (!is.null(pathways_result) && nrow(pathways_result) > 0) {
    # Clean pathway names for display
    pathways_result$Pathways_clean <- gsub(paste0(kegg_organism, "\\d+"), "", pathways_result$Pathways)
    pathways_result$Pathways_clean <- gsub("^\\s+|\\s+$", "", pathways_result$Pathways_clean)
    pathways_result$Pathways_clean <- gsub("Biosynthesis of ", "", pathways_result$Pathways_clean)
    pathways_result$Pathways_clean <- gsub("Two-component system", "Two component system", pathways_result$Pathways_clean)
    pathways_result$Pathways_clean <- paste0(toupper(substr(pathways_result$Pathways_clean, 1, 1)), 
                                             substr(pathways_result$Pathways_clean, 2, nchar(pathways_result$Pathways_clean)))
    
    # Plot
    p <- ggplot(pathways_result, aes(x = reorder(Pathways_clean, Count), y = Count)) +
      geom_col(fill = plot_colors[i], width = 0.7) +
      coord_flip() +
      labs(title = comparison_labels[comp_name],
           x = "", 
           y = "Number of DEGs") +
      theme_classic() +
      theme(
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        axis.text = element_text(size = 14, color = "black", face = "bold"),
        axis.title.x = element_text(size = 14, color = "black", face = "bold"),
        axis.title.y = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1)
      ) +
      geom_text(aes(label = Count), hjust = -0.2, size = 4.5, color = "black", face = "bold") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks = function(x) seq(0, max(x), by = 1))
    
    pathway_plots[[comp_name]] <- p
    
    ggsave(file.path(output_dir, "figures", paste0(comp_name, "_pathway.pdf")), 
           plot = p, width = 8, height = 6)
  }
}

cat("Pathway analysis complete!\n")
