# Huntington's Disease RNA-Seq Explorer

An R Shiny application for exploring bulk RNA-Seq data from the Huntington's Disease dataset GSE64810 (post-mortem BA9 prefrontal cortex: 20 HD patients vs 49 neurologically normal controls).

## Dataset

- GEO accession: [GSE64810](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810)
- Reference: Labadorf et al. (2015), PLoS One, 10(12):e0143563. PMID: 26636579
- Tissue: Human post-mortem dorsolateral prefrontal cortex (Brodmann Area 9)
- Samples: 69 total (20 HD, 49 controls)
- Genes: 28,087

## Repository contents

```
.
├── app.R              # Shiny application (UI + server)
├── preprocessing.R    # Script that converts raw GEO files into the CSVs the app loads
├── sample_info.csv    # Sample metadata (69 rows x 14 columns)
├── de_results.csv     # DESeq2 differential expression results (28,087 rows x 10 columns)
└── README.md
```

Note: `norm_counts.csv` (29 MB) is not included due to file size. To generate it, download `GSE64810_mlhd_DESeq2_norm_counts_adjust.txt.gz` from the [GEO supplementary files](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810), place it in the same folder as `preprocessing.R`, and run `preprocessing.R`.

## Running the app

1. Install required R packages:

```r
install.packages(c("shiny", "tidyverse", "DT", "pheatmap", "patchwork", "colourpicker"))
```

2. Open `app.R` in RStudio and click Run App, or from the Console:

```r
shiny::runApp("app.R")
```

3. In each tab, upload the corresponding CSV file.

## App tabs

### 1. Samples
Upload `sample_info.csv`. Three sub-tabs: a type/summary table, a sortable data table, and histograms of continuous variables (user-selectable column).

### 2. Counts
Upload `norm_counts.csv`. Two sliders filter genes by variance percentile and minimum non-zero samples. Four sub-tabs: filter summary, diagnostic scatter plots (median vs variance; median vs zeros), clustered heatmap with optional log2 transform, and PCA (scatter of two user-chosen PCs, or beeswarm of top N PCs).

### 3. DE (Differential Expression)
Upload `de_results.csv`. User picks x/y axis columns, base/highlight colors, and a significance threshold (10^X). Two sub-tabs: volcano plot and a searchable, sortable results table.

### 4. Gene Viz (Individual Gene Visualization)
Upload both `norm_counts.csv` and `sample_info.csv`. User picks a grouping variable, a gene (searchable dropdown), and a plot type (Bar, Boxplot, Violin, or Beeswarm). The Plot button triggers rendering.

## Input validation

All file inputs are validated: the extension must be `.csv` or `.tsv`, and required columns must be present. If validation fails, a clear error message is shown in place of the plot/table.

## Author

Shafeer Faizan — Boston University, Spring 2026
