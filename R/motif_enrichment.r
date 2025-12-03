## Motif Enrichment Analysis Functions
## Functions for finding and analyzing enriched motifs in accessible regions

#' Find Enriched Transcription Factor Binding Motifs
#'
#' Identifies transcription factor binding motifs enriched in differentially
#' accessible chromatin regions for a specific cluster. Uses Signac's FindMotifs
#' to perform motif enrichment analysis and generates visualization of top motifs.
#' Requires ATAC-seq data with motif annotations.
#'
#' @param seurat_object Seurat object with ATAC assay and motif annotations
#' @param da_results Data frame; differential accessibility results from FindAllMarkers
#'   containing columns: peak, cluster, pct.1, avg_log2FC, p_val_adj
#' @param cluster_id Character; cluster identifier to analyze
#' @param pct_threshold Numeric; minimum percentage of cells with accessible peak.
#'   Default 0.1
#' @param log2fc_threshold Numeric; minimum log2 fold change for peak filtering.
#'   Default 1
#' @param n_motifs_to_plot Integer; number of top motifs to visualize. Default 10
#'
#' @return List with elements or NULL if no peaks pass filtering:
#'   \item{enriched_motifs}{Data frame; motif enrichment results with p-values}
#'   \item{plot}{ggplot; visualization of top enriched motif logos}
#'   \item{n_peaks}{Integer; number of peaks used for enrichment analysis}
#'
#' @examples
#' motifs <- find_enriched_motifs(seurat_atac, da_peaks, cluster_id = "0",
#'                               pct_threshold = 0.1, log2fc_threshold = 1,
#'                               n_motifs_to_plot = 15)
#' @export
find_enriched_motifs <- function(seurat_object,
                  da_results,
                  cluster_id,
                  pct_threshold = 0.1,
                  log2fc_threshold = 1,
                  n_motifs_to_plot = 10) {

  # Check for required packages
  if (!requireNamespace("Signac", quietly = TRUE)) {
    stop("Package 'Signac' is required but not installed. Please install it with:\n",
         "  install.packages('Signac')",
         call. = FALSE)
  }

  # Filter top accessible peaks for the specified cluster
  top_accessible_peaks <- da_results |>
  filter(cluster == cluster_id) |>
  filter(pct.1 > pct_threshold & avg_log2FC > log2fc_threshold) |>
  arrange((p_val_adj)) |>
  pull(peak) |>
  unique()

  # Print number of peaks
  message(sprintf("Number of top accessible peaks for cluster %s: %d",
          cluster_id, length(top_accessible_peaks)))

  if (length(top_accessible_peaks)  ==  0) {
    message("No peaks passed the filtering criteria.")
    return(NULL)

  }
  # Find enriched motifs
  enriched_motifs <- FindMotifs(
  object = seurat_object,
  features = top_accessible_peaks
  )

  # Create motif plot
  motif_plot <- MotifPlot(
  object = seurat_object,
  motif = rownames(enriched_motifs)[1:min(n_motifs_to_plot, nrow(enriched_motifs))]
  ) + labs(title = paste0('Enriched Motifs in Cluster ', cluster_id),
         subtitle = paste0('Top ', min(n_motifs_to_plot, nrow(enriched_motifs)), ' motifs'))

  print(motif_plot)

  # Return both results and plot
  return(list(
  enriched_motifs = enriched_motifs,
  plot = motif_plot,
  n_peaks = length(top_accessible_peaks)
  ))
}
