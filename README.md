# pstA1 RNA-seq Analysis

RNA-seq differential expression analysis for *Mycobacterium tuberculosis* pstA1 mutant response to rifampicin treatment.

## Citation

If you use this code or data, please cite:
Biorxiv link
## Authors
G V Anoushka Chinmayi

**Data:**
```
Gene Expression Omnibus accession GSE[ACCESSION_NUMBER]
```

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
# 1. Clone the repository
git clone https://github.com/AnoushkaGV/pstA1.git
cd pstA1

# 2. Create data directory
mkdir data

# 3. Download data files from GEO (see Data Availability section above)
# Place count_matrix.csv and rv_dict.csv in the data/ directory
```

## Data Availability

**All sequencing data and count matrices are available on GEO (Gene Expression Omnibus):**

- **GEO Accession:** GSE[ACCESSION_NUMBER]
- **Link:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE[ACCESSION_NUMBER]

### Required Files from GEO:

Download the following files from GEO and place them in a `data/` directory:

1. **count_matrix.csv** - Raw count matrix (genes × samples)
   - Format: First column = Gene IDs (Rv numbers), remaining columns = sample counts
   - Column names: `pstA1_ctrl_1`, `pstA1_ctrl_2`, `pstA1_ctrl_3`, `pstA1_RIF_30min_1`, etc.

2. **rv_dict.csv** - Gene ID to gene name mapping
   - Format: Two columns (Rv, Name)
   - Example:
   ```
   Rv,Name
   Rv0001,dnaA
   Rv0002,dnaN
   ```

### Data Structure:

After downloading from GEO, organize files as:
```
pstA1/
├── data/
│   ├── count_matrix.csv
│   └── rv_dict.csv
└── ...
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

