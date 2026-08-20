
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scAnalysisToolkit

<!-- badges: start -->

<!-- badges: end -->

## Overview

scAnalysisToolkit provides a comprehensive suite of helper functions for
single-cell RNA-seq and ATAC-seq analysis workflows. The package
includes helper tools and wrappers for:

- **Clustering and Annotation**: Seurat-based clustering workflows and
  cell type annotation with SingleR
- **Differential Expression**: Pseudobulk, Wilcoxon, and bulk RNA-seq
  differential expression analysis
- **Differential Accessibility**: ATAC-seq peak accessibility analysis
- **Pathway Enrichment**: MSigDB, GO term and Metascape enrichment
  analysis, with optional collapse of redundant terms into gene-overlap
  clusters
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

# ...with pathway enrichment, collapsing redundant terms into gene-overlap
# clusters and adding a pair of plots for each topic of interest
de_results <- pseudobulk(
  seurat_object,
  comparison = "condition",
  group1 = "WT",
  group2 = "KO",
  cluster = "Tcells",
  run_pathway_enrichment = "ClusterProfiler",
  collapse_redundancy = TRUE,
  topics = c("lipid", "immune")
)
```

An enrichment result of a few hundred significant terms is mostly
near-duplicates: the same genes reported under `GOBP_`, `REACTOME_`,
`WP_` and `HALLMARK_` names. `collapse_redundancy` groups those terms by
Jaccard similarity of their enriched genes followed by Louvain
clustering, which works across mixed collections where
`clusterProfiler::simplify()` does not. Each enrichment directory then
also gets a cluster bar plot, a term-overlap network and both cluster
tables. The pieces are usable on their own for interactive work:

``` r
clusters <- collapse_pathway_redundancy(de_results$enrichment$results_up)
clusters
#> 137 terms -> 6 clusters (jaccard >= 0.25, louvain, seed 1)
#>   1. Cholesterol metabolism with bloch and kandutsch... (39 terms, p.adj 5e-24)
#>   2. Defense response to virus (67 terms, p.adj 1.7e-11)

plot_pathway_clusters(clusters)     # one bar per cluster
plot_pathway_cluster_map(clusters)  # terms linked by gene overlap
plot_pathway_overlap(clusters)      # the redundancy itself, as a heatmap


# Is a topic one redundant block or several distinct signals?
plot_pathway_topic(clusters, "lipid")
plot_pathway_topic_genes(clusters, "lipid")
```

With `collapse_redundancy = TRUE` the enrichment call prints nothing at
all — not even the over-representation bar plots. It writes its figures
and tables to disk and hands every plot back, so you decide which ones
appear and where:

``` r
enrichment <- de_results$enrichment

enrichment$plot_up                          # over-representation bar plot
enrichment$redundancy_plots_up$clusters     # one bar per collapsed cluster
enrichment$redundancy_plots_up$network      # terms linked by gene overlap
enrichment$redundancy_plots_up$topics$lipid$terms
enrichment$redundancy_plots_up$topics$lipid$genes
```

`plot_pathway_cluster_map()` can also return an interactive plotly
network, where every node carries its full name, collection, adjusted
p-value, gene count and cluster on hover. That is worth asking for,
because the static version can only label one term per cluster. It has
to be the value of a top-level expression for knitr to embed it:

``` r
plot_pathway_cluster_map(de_results$enrichment$redundancy_up, interactive = TRUE)
```
