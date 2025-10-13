# pstA1 RNA-seq Analysis

RNA-seq differential expression analysis for *Mycobacterium tuberculosis* pstA1 mutant response to rifampicin treatment.

## Citation

## Requirements

**R version:** >= 4.2.0

**R packages:**
```r
# From Bioconductor
BiocManager::install(c("DESeq2", "KEGGREST"))

# From CRAN
install.packages(c("pheatmap", "ggplot2", "tidyr", "dplyr", "ggvenn"))
```

## Installation

```bash
git clone https://github.com/AnoushkaGV/pstA1.git
cd pstA1
```
## Data Availability

Raw sequencing data and processed count matrices are available at NCBI Gene Expression Omnibus (GEO):
- **GEO Accession:**
- **Link:**

To run this analysis pipeline:
1. Download count matrix and gene annotation files from GEO
2. Place files in the `data/` directory as described below

## Data Structure

Place your input files in the `data/` directory:

```
data/
├── cc_matrix.csv    # Raw count matrix (genes × samples)
└── rv_dict.csv      # Gene ID to gene name mapping
```

**cc_matrix.csv format:**
- First column: Gene IDs (Rv numbers)
- Remaining columns: Sample counts
- Column names: `pstA1_ctrl_1`, `pstA1_ctrl_2`, `pstA1_ctrl_3`, `pstA1_RIF_30min_1`, etc.

**rv_dict.csv format:**
```
Rv,Name
Rv0001,dnaA
Rv0002,dnaN
```

## Usage

Run the complete analysis:

```bash
Rscript run_analysis.R
```

Or run scripts individually:

```r
source("01_deseq2_analysis.R")    # DESeq2 analysis
source("02_generate_plots.R")     # MA plots, heatmap, Venn diagrams
source("03_pathway_analysis.R")   # KEGG pathway enrichment
source("04_summary_tables.R")     # Summary and supplementary tables
```

## Output

Results are saved in the `results/` directory:

```
results/
├── deseq2_results.RData          # R objects for downstream analysis
├── figures/
│   ├── MA_plots.pdf
│   ├── heatmap_all_DEGs.pdf
│   ├── venn_diagrams.pdf
│   └── [comparison]_pathway.pdf
├── tables/
│   ├── DEG_summary.csv
│   └── Supplementary_Table_All_DEGs.csv
└── session_info.txt
```

## Analysis Parameters

Default thresholds (modify in `config.R`):
- Log2 fold change: |LFC| > 1
- Adjusted p-value: padj < 0.05

## Experimental Design

- **Strains:** Wild-type (WT) and pstA1 mutant
- **Treatment:** Rifampicin (RIF)
- **Time points:** 0 (control), 30 min, 60 min
- **Replicates:** 3 biological replicates per condition

## Comparisons

1. Time course within strain (WT and pstA1: 30min vs ctrl, 60min vs ctrl)
2. Strain differences at each time point (pstA1 vs WT at ctrl, 30min, 60min)

## License

MIT License

## Contact

G V Anoushka Chinmayi 
Johns Hopkins University

## References

Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*. 2014;15(12):550.
Kanehisa M, Goto S. KEGG: Kyoto Encyclopedia of Genes and Genomes. *Nucleic Acids Research*. 2000;28(1):27-30.
