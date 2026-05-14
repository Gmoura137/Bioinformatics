# RNA-Seq Differential Expression & Heatmap Analysis

A complete RNA-Seq analysis workflow in **R** for identifying differentially expressed genes, generating clustered heatmaps, performing PCA analysis, and visualizing gene expression patterns.

This repository contains scripts for:

* Differential expression analysis using `DESeq2`
* Heatmap generation with hierarchical clustering
* PCA visualization
* Z-score normalization
* RNA-Seq expression analysis in treated vs control samples

---

# Repository Structure

```bash
.
├── DEseq_2024.R
├── RNAseq.R
├── hela_cells_expression.R
├── airway_scaledcounts.csv
├── airway_metadata.csv
├── RNAseq_mat_top20.csv
├── de_result.all.csv
├── de_result.filtered.csv
├── normalized_counts.csv
├── heatmap_de.jpg
├── rna_seq.jpg
└── README.md
```

---

# Overview

This project demonstrates a complete RNA-Seq data analysis pipeline using publicly available airway dataset examples.

The workflow includes:

1. Loading RNA-Seq count data
2. Metadata preprocessing
3. Differential expression analysis
4. Filtering significant genes
5. Sample clustering
6. Heatmap visualization
7. PCA dimensionality reduction
8. Exporting publication-quality figures

---

# Workflow Overview

```text
Raw Count Matrix
        ↓
Metadata Processing
        ↓
DESeq2 Dataset Construction
        ↓
Normalization & Dispersion Estimation
        ↓
Differential Expression Analysis
        ↓
Gene Filtering
        ↓
PCA & Hierarchical Clustering
        ↓
Heatmap Generation
        ↓
Export Results & Figures
```

---

# Scripts Description

## 1. `DEseq_2024.R`

Performs differential gene expression analysis using the `DESeq2` package.

### Features

* Imports count matrix and metadata
* Creates DESeq2 dataset
* Filters low-count genes
* Performs normalization
* Computes differential expression statistics
* Generates:

  * PCA plots
  * Dispersion plots
  * Heatmaps
  * Z-score normalized heatmaps
* Exports filtered differentially expressed genes

### Output Files

| File                     | Description             |
| ------------------------ | ----------------------- |
| `de_result.all.csv`      | Complete DESeq2 results |
| `de_result.filtered.csv` | Significant DE genes    |
| `normalized_counts.csv`  | Normalized count matrix |

---

## 2. `RNAseq.R`

Creates clustered heatmaps from RNA-Seq expression matrices.

### Features

* Row scaling (Z-score normalization)
* Hierarchical clustering
* Sample annotations
* Custom color palettes
* Publication-quality heatmaps
* JPEG export

### Output

| File          | Description             |
| ------------- | ----------------------- |
| `rna_seq.jpg` | Final clustered heatmap |

---

## 3. `hela_cells_expression.R`

Performs manual preprocessing and clustering analysis for gene expression data.

### Features

* Removes NA values
* Removes zero-count genes
* Calculates treated/control averages
* Computes log2 fold changes
* Selects top differentially expressed genes
* Generates clustered heatmaps

### Output

| File             | Description             |
| ---------------- | ----------------------- |
| `heatmap_de.jpg` | Heatmap of top DE genes |

---

# Installation

## Clone Repository

```bash
git clone https://github.com/yourusername/rnaseq-analysis.git
cd rnaseq-analysis
```

---

# Required R Packages

Install all required packages before running the analysis scripts.

```r
install.packages(c(
  "pheatmap",
  "tidyverse",
  "ggplotify",
  "heatmaply",
  "ggrepel",
  "RColorBrewer"
))

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
```

---

# Input Data

## `airway_scaledcounts.csv`

RNA-Seq raw count matrix.

* Rows = genes
* Columns = samples

---

## `airway_metadata.csv`

Sample metadata containing treatment information.

### Example

| Sample     | dex     |
| ---------- | ------- |
| SRR1039508 | control |
| SRR1039509 | treated |

---

## `RNAseq_mat_top20.csv`

Expression matrix containing the top 20 differentially expressed genes.

---

# Usage

## Step 1 — Set Working Directory

Update the working directory inside each R script.

```r
setwd("path/to/project")
```

---

## Step 2 — Run Differential Expression Analysis

```r
source("DEseq_2024.R")
```

---

## Step 3 — Generate Heatmaps

```r
source("RNAseq.R")
```

---

## Step 4 — Run HeLa Cell Expression Analysis

```r
source("hela_cells_expression.R")
```

---

# Differential Expression Criteria

Genes are considered significantly differentially expressed if:

```text
Adjusted p-value (padj) < 0.05
AND
|log2FoldChange| > 1
```

---

# Visualization

## Heatmaps

The workflow generates clustered heatmaps using normalized gene expression values.

### Included Features

* Hierarchical clustering
* Sample annotations
* Z-score normalization
* Gene expression scaling
* Publication-quality formatting

### Color Scheme

| Color | Meaning                 |
| ----- | ----------------------- |
| Blue  | Low expression          |
| White | Intermediate expression |
| Red   | High expression         |

---

## PCA Analysis

Principal Component Analysis (PCA) is used to visualize sample clustering and assess variation between treated and control groups.

---

# Example Outputs

## Differential Expression Results

```text
Gene ID
log2FoldChange
pvalue
padj
```

---

## Heatmap Outputs

* Clustered genes
* Clustered samples
* Treatment annotations
* Expression intensity visualization

---

# Key Packages Used

| Package        | Purpose                          |
| -------------- | -------------------------------- |
| `DESeq2`       | Differential expression analysis |
| `pheatmap`     | Heatmap visualization            |
| `tidyverse`    | Data wrangling                   |
| `ggplot2`      | Plotting and visualization       |
| `RColorBrewer` | Color palette generation         |
| `heatmaply`    | Interactive heatmaps             |

---

# Applications

This workflow can be applied to:

* RNA-Seq analysis
* Cancer genomics
* Drug response studies
* Biomarker discovery
* Gene expression profiling
* Transcriptomics research
* Machine learning preprocessing for genomics

---

# Notes

* Ensure metadata sample names match count matrix column names.
* `DESeq2` requires raw integer count data.
* Heatmaps perform best after normalization and scaling.
* Removing low-count genes improves statistical power.
* PCA helps identify batch effects and sample outliers.

---

# Future Improvements

Potential future additions include:

* Volcano plot generation
* Interactive heatmaps
* GO enrichment analysis
* KEGG pathway analysis
* Automated HTML/PDF reports
* Shiny dashboard integration
* Batch effect correction
* Multi-factor experimental design support

---

# License

This project is licensed under the MIT License.

---

# Author

RNA-Seq Bioinformatics Workflow in R using:

* `DESeq2`
* `pheatmap`
* `tidyverse`
* `ggplot2`
* `heatmaply`

Developed for RNA-Seq differential expression and visualization analysis.
