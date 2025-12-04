
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scAnalysisToolkit

<!-- badges: start -->

<!-- badges: end -->

## Overview

**For internal use by the DeCarvalho and Constantinides labs at Scripps
Research.**

scAnalysisToolkit provides a comprehensive suite of helper functions for
single-cell RNA-seq and ATAC-seq analysis workflows. The package
includes helper tools and wrappers for:

- **Clustering and Annotation**: Seurat-based clustering workflows and
  cell type annotation with SingleR
- **Differential Expression**: Pseudobulk, Wilcoxon, and bulk RNA-seq
  differential expression analysis
- **Differential Accessibility**: ATAC-seq peak accessibility analysis
- **Pathway Enrichment**: GO term and Metascape enrichment analysis
- **Repertoire Analysis**: T-cell receptor repertoire analysis and
  visualization
- **Cell Communication**: Cell-cell communication analysis with CellChat
- **Gene Programs**: cNMF gene expression program analysis
- **Motif Analysis**: Motif enrichment in accessible chromatin regions
- **Visualizations**: Volcano plots, scatterplots, heatmaps, and circos
  diagrams

The package is designed to streamline end-to-end single-cell analysis
with consistent interfaces and automated output generation, enabling
iterative frictionless analysis of datasets. It builds on functionality
provided by packages such as Seurat, Signac, SingleR, DESeq2, and
CellChat.

**Note**: This package operates with Seurat objects in R only (no
support for scran or scanpy).

## Installation

You can install the development version of scAnalysisToolkit from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("EduardAnsaldo/scAnalysisToolkit")
```

## Usage

Section in construction…

The package provides high-level wrapper functions for common single-cell
analysis workflows:

``` r
library(scAnalysisToolkit)

# Clustering workflow
clustering_results <- perform_seurat_clustering(
  seurat_object,
  dimensions = 30,
  resolutions = c(0.1, 0.25, 0.5)
)

# Cell type annotation
annotated_seurat <- annotate_seurat_with_SingleR_Eduard(
  seurat_object,
  local_path = "./results",
  database = "ImmGen"
)

# Differential expression analysis
de_results <- pseudobulk(
  seurat_object,
  comparison = "condition",
  group1 = "WT",
  group2 = "KO",
  cluster = "Tcells"
)
```
