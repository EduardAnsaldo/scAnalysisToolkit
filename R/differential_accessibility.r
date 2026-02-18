## Differential Accessibility Analysis Functions
## Functions for ATAC-seq analysis

#' Find Top Differentially Accessible Peaks Per Cluster
#'
#' Identifies differentially accessible chromatin regions for each cluster using
#' Seurat's FindAllMarkers on ATAC-seq data. Maps peaks to nearest genes and
#' optionally performs pathway enrichment analysis. Saves ranked peak lists
#' (top 10, 50, 100) to TSV files.
#'
#' @param seurat Seurat object with ATAC assay and cluster assignments
#' @param n_genes_to_plot Integer; number of top peaks per cluster for visualization.
#'   Default 3
#' @param grouping_var Character; metadata column name for cluster identity.
#'   Default 'seurat_clusters_atac'
#' @param object_annotations Character; string to append to output file names. Default ''
#' @param tables_path Character; path to save output tables. Default 'results/tables/'
#' @param figures_path Character; path to save figures. Default 'results/figures/'
#' @param results_path Character; path for pathway enrichment results. Default 'results/'
#' @param run_pathway_enrichment Character vector; enrichment methods to run
#'   ('Metascape', 'ClusterProfiler', or NULL). Default NULL
#' @param ... Additional arguments passed to ClusterProfiler analysis
#'
#' @return List with elements:
#'   \item{ClusterProfiler_results}{ClusterProfiler enrichment results or NULL}
#'   \item{metascape_results}{Metascape enrichment results or NULL}
#'   \item{topn}{Data frame; top n peaks per cluster}
#'   \item{da_results}{Data frame; all differential accessibility results}
#'
#' @export
top_peaks_per_cluster <- function (seurat, n_genes_to_plot = 3, grouping_var = 'seurat_clusters_atac', object_annotations = '', tables_path = 'results/tables/', figures_path = 'results/figures/', results_path = 'results/', run_pathway_enrichment = NULL, ...) {

    # Set the identity class for clustering
    Idents(seurat) <- grouping_var
    # Set the default assay to 'ATAC'
    DefaultAssay(seurat) <- 'ATAC'

    da_peaks <- FindAllMarkers(seurat, test.use = 'LR', only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, latent.vars = 'nCount_ATAC')

    genes <- Signac::ClosestFeature(seurat, rownames(da_peaks))
    da_peaks <- left_join(da_peaks, genes, by = c('gene' = 'query_region'))
    da_peaks <- da_peaks |>
        select(-tx_id, -gene_id) |>
        rename(peak = gene, gene = gene_name)


    #Add gene annotations:
    da_peaks <- da_peaks |>
                    left_join(y= unique(annotations[,c('gene_name', 'description')]),
                        by = c('gene' = 'gene_name'))   
    
    if (all(is.numeric(unique(da_peaks$cluster)))) {
        da_peaks <- da_peaks |>
            mutate(cluster = fct_inseq(cluster))
    }

    #Top10 markers
    da_peaks |>
        group_by(cluster) |>
        arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) |>
        slice_head(n = 10) -> top10

    #Top25 markers
    da_peaks |>
        group_by(cluster) |>
        arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) |>
        slice_head(n = 50) -> top50

    #Top100 markers
    da_peaks |>
        group_by(cluster) |>
        arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) |>
        slice_head(n = 100) -> top100

    #Topn markers
    da_peaks |>
        group_by(cluster) |>
        arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) |>
        slice_head(n = n_genes_to_plot) -> topn

    # Save the top markers to files
    write.table(top100,file=here::here(tables_path, paste0('top100_peaks', '_',object_annotations, ".tsv")), sep="\t",row.names = FALSE)
    # write.table(top25,file=here::where(path,'top25',object_annotations, ".tsv"), sep="\t",row.names = FALSE)
    write.table(top10,file=here::here(tables_path, paste0('top10_peaks', '_',object_annotations, ".tsv")), sep="\t",row.names = FALSE)

    top100_genes_per_cluster <- top100 |>
        group_by(cluster) |>
        summarise(genes = str_flatten_comma(gene))

    write.table(top100_genes_per_cluster,
                file = here::here(tables_path, paste0('top100_peak_gene_names_per_cluster_', object_annotations, ".tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)


    metascape_results <- NULL
    ClusterProfiler_results <- NULL
        # Run pathway enrichment analysis
    if (!is.null(run_pathway_enrichment)) {
        if ('Metascape' %in% run_pathway_enrichment) {
            metascape_results <- Metascape_functional_analysis_cluster_identification(seurat, top100, identities = grouping_var, path=results_path, object_annotations = object_annotations)
        }
        if ('ClusterProfiler' %in% run_pathway_enrichment) {
            ClusterProfiler_results <- GO_functional_analysis_cluster_identification(seurat, da_peaks, path=results_path, object_annotations = object_annotations, identities = grouping_var, top_gene_number = 100, ...)
        }

    }

    return(list(ClusterProfiler_results = ClusterProfiler_results, metascape_results = metascape_results, topn = topn, da_results =  da_peaks ) )

    }


#' Core Differential Accessibility Analysis
#'
#' Performs differential accessibility testing on ATAC-seq data using Seurat's
#' FindMarkers with Logistic regression tests with nCount_ATAC as latent variable. Maps differentially accessible peaks
#' to nearest genes and saves results to CSV files. Requires sufficient cells
#' per comparison group.
#'
#' @param seurat Seurat object with ATAC assay and chromatin accessibility data
#' @param comparison Character; metadata column name for comparison groups
#' @param group1 Character; pattern to match for first group (control)
#' @param group2 Character; pattern to match for second group (treatment)
#' @param cluster Character; cluster or cell type identifier. Default 'all_clusters'
#' @param path Character; output directory path. Default './'
#' @param FC_threshold Numeric; log2 fold change threshold for filtering. Default 0.3
#' @param p_value_threshold Numeric; adjusted p-value threshold. Default 0.05
#' @param min_fraction Numeric; minimum fraction of cells with accessible peak. Default 0.05
#' @param minimum_cell_number Integer; minimum cells required per group. Default 30
#' @param ... Additional arguments (unused)
#'
#' @return List with elements:
#'   \item{all_count}{Integer; number of significant differentially accessible peaks}
#'   \item{UP_count}{Integer; number of more accessible peaks}
#'   \item{DOWN_count}{Integer; number of less accessible peaks}
#'   \item{results}{Data frame; complete DA results with nearest genes and annotations}
#'
#' @export
run_differential_accessibility_FindMarkers <- function(seurat, comparison, group1, group2,
                                         cluster = 'all_clusters', path = './', FC_threshold = 0.3,
                                         p_value_threshold = 0.05, min_fraction = 0.05,
                                         minimum_cell_number = 30, ...) {

    # Set Paths
    gene_lists_path <- here::here(path, 'gene_lists')
    dir.create(gene_lists_path, showWarnings = FALSE, recursive = TRUE)

    message('Cluster: ', cluster)

    group1 <- fixed(group1)
    group2 <- fixed(group2)

    Idents(seurat) <- comparison

    message('Number of cells in group 1: ', seurat@meta.data |> filter(str_detect(!!as.name(comparison), group1)) |> nrow())
    message('Number of cells in group 2: ', seurat@meta.data |> filter(str_detect(!!as.name(comparison), group2)) |> nrow())

    # Check there are enough cells
    n_group1 <- seurat@meta.data |> filter(str_detect(!!as.name(comparison), group1)) |> nrow()
    n_group2 <- seurat@meta.data |> filter(str_detect(!!as.name(comparison), group2)) |> nrow()

    if (n_group1 < minimum_cell_number | n_group2 < minimum_cell_number) {
        warning('Not enough cells per group (minimum ', minimum_cell_number, ' required)')
        return(list(
            all_count = 'Not enough cells',
            UP_count = 'Not enough cells',
            DOWN_count = 'Not enough cells',
            results = NULL
        ))
    }

    DefaultAssay(seurat) <- 'ATAC'
    # Run Wilcox test
    results <- FindMarkers(object = seurat, ident.1 = group1, ident.2 = group2,
                          assay = 'ATAC', slot = 'data', test.use = 'LR', latent.vars = 'nCount_ATAC', min.pct = min_fraction, only.pos = F, logfc.threshold = 0.5)


    closest_genes <- Signac::ClosestFeature(seurat, rownames(results))


    results <- results  |>
        rownames_to_column('gene') |>
        left_join(closest_genes, by = c('gene' = 'query_region'))
    results <- results |>
        select(-tx_id, -gene_id) |>
        rename(peak = gene, genes = gene_name)

    #Add gene annotations:
    results <- results |>
                left_join(y= unique(annotations[,c('gene_name', 'description')]),
                    by = c('genes' = 'gene_name')) |>
                rename(
                    log2FoldChange = avg_log2FC,
                    padj = p_val_adj,
                    pvalue = p_val
                ) |>
                mutate(log2FoldChange = log2FoldChange * -1)

    # Filter results
    results_filtered <- filter(
        results,
        padj < p_value_threshold &
            (log2FoldChange >= FC_threshold | log2FoldChange <= -1 * FC_threshold)
    ) |>
        arrange(padj)

    results_filtered_UP <- filter(results_filtered, log2FoldChange >= FC_threshold)
    results_filtered_DOWN <- filter(results_filtered, log2FoldChange <= -1 * FC_threshold)

    # Write results to CSV files
    write.csv(results |> arrange(padj),
             file = here::here(gene_lists_path, paste('ALL_PEAKS_DA_Analysis', cluster, 'Wilcox', group2, 'vs', group1, '.csv', sep = '_')))
    write.csv(results_filtered_UP |> arrange(padj, desc(log2FoldChange)),
             file = here::here(gene_lists_path, paste('PEAKS_UP_in', group2, cluster, 'Wilcox', group2, 'vs', group1, '.csv', sep = '_')))
    write.csv(results_filtered_DOWN |> arrange(padj, log2FoldChange),
             file = here::here(gene_lists_path, paste('PEAKS_UP_in_DOWN_in', group2, cluster, 'Wilcox', group2, 'vs', group1, '.csv', sep = '_')))

    # Return counts and results
    list(
        all_count = nrow(results_filtered),
        UP_count = nrow(results_filtered_UP),
        DOWN_count = nrow(results_filtered_DOWN),
        results = results
    )
}


#' Complete Differential Accessibility Analysis Pipeline
#'
#' Orchestrates a full ATAC-seq differential accessibility analysis workflow including
#' statistical testing, volcano plot generation for both peak regions and nearest genes,
#' and optional pathway enrichment analysis. Creates two volcano plots: one labeled
#' with peak coordinates and one with nearest gene names.
#'
#' @param seurat Seurat object with ATAC assay and chromatin accessibility data
#' @param comparison Character; metadata column name for comparison groups
#' @param group1 Character; pattern to match for first group (control)
#' @param group2 Character; pattern to match for second group (treatment)
#' @param is_integrated_subset Logical; currently unused, maintained for consistency.
#'   Default FALSE
#' @param cluster Character; cluster or cell type identifier. Default 'all_clusters'
#' @param path Character; output directory path. Default './'
#' @param FC_threshold Numeric; log2 fold change threshold for filtering. Default 0.3
#' @param p_value_threshold Numeric; adjusted p-value threshold. Default 0.05
#' @param min_fraction Numeric; minimum fraction of cells with accessible peak. Default 0.05
#' @param max_overlaps Integer; maximum label overlaps in plots. Default 15
#' @param label_size Numeric; size of peak/gene labels in plots. Default 5
#' @param pathways_of_interest Character vector; specific pathways to analyze. Default NULL
#' @param label_threshold Numeric; accessibility threshold for labeling peaks. Default 100000
#' @param gene_lists_to_plot Named list; custom gene lists to visualize. Default NULL
#' @param colors Character vector of length 2; colors for DOWN and UP peaks.
#'   Default c('green4', 'darkorchid4')
#' @param minimum_cell_number Integer; minimum cells required per group. Default 30
#' @param run_pathway_enrichment Character; method for pathway enrichment or NULL. Default NULL
#' @param ... Additional arguments passed to core analysis and visualization functions
#'
#' @return List with elements:
#'   \item{all_count}{Integer; number of significant differentially accessible peaks}
#'   \item{UP_count}{Integer; number of more accessible peaks}
#'   \item{DOWN_count}{Integer; number of less accessible peaks}
#'   \item{results}{Data frame; complete DA results}
#'
#' @export
run_differential_accessibility <- function(seurat, comparison, group1, group2, is_integrated_subset = FALSE,
                                      cluster = 'all_clusters', path = './', FC_threshold = 0.3,
                                      p_value_threshold = 0.05, min_fraction = 0.05, max_overlaps = 15, label_size = 5,
                                      pathways_of_interest = NULL, label_threshold = 100000,
                                      gene_lists_to_plot = NULL,
                                      colors = c('green4', 'darkorchid4'),
                                      minimum_cell_number = 30, run_pathway_enrichment = NULL, ...) {

    local_figures_path <- here::here(path, 'figures')
    dir.create(local_figures_path, showWarnings = FALSE, recursive = TRUE)

    group1 <- fixed(group1)
    group2 <- fixed(group2)

    # Set colors for the plot
    my_colors <- c(colors, "gray")
    names(my_colors) <- c("DOWN", "UP", "NO")

    # Run core DE analysis
    da_results <- run_differential_accessibility_FindMarkers(
        seurat = seurat,
        comparison = comparison,
        group1 = group1,
        group2 = group2,
        cluster = cluster,
        path = path,
        FC_threshold = FC_threshold,
        p_value_threshold = p_value_threshold,
        min_fraction = min_fraction,
        minimum_cell_number = minimum_cell_number,
        ...
    )

    # Check if analysis was successful
    if (is.null(da_results$results)) {
        return(list(
            all_count = da_results$all_count,
            UP_count = da_results$UP_count,
            DOWN_count = da_results$DOWN_count
        ))
    }

    results_regions <- da_results$results |>
        select(-genes) |>
        rename(genes = peak)

    # Without gene names
    volcanoplot_output <- volcano_plot(
        results = results_regions,
        group1 = group1,
        group2 = group2,
        cluster = cluster,
        my_colors = my_colors,
        local_figures_path = local_figures_path,
        FC_threshold = FC_threshold,
        p_value_threshold = p_value_threshold,
        max_overlaps = 1,
        label_size = label_size,
        label_threshold = label_threshold,
        test_type = 'Wilcox_ATAC',
        ...
    )

    # With gene names
    volcanoplot_output <- volcano_plot(
        results = da_results$results,
        group1 = group1,
        group2 = group2,
        cluster = cluster,
        my_colors = my_colors,
        local_figures_path = local_figures_path,
        FC_threshold = FC_threshold,
        p_value_threshold = p_value_threshold,
        max_overlaps = max_overlaps,
        label_size = label_size,
        label_threshold = label_threshold,
        test_type = 'Wilcox_ATAC_closest_genes',
        ...
    )

    # Overrepresentation analysis

    run_DEG_functional_analysis(
        results = da_results$results,
        method = run_pathway_enrichment,
        grouping_var = cluster,
        path = path,
        FC_threshold = FC_threshold,
        p_value_threshold = p_value_threshold,
        group1 = group1,
        group2 = group2,
        pathways_of_interest = pathways_of_interest,
        comparison = comparison,
        ...
    )

    return(list(
        all_count = da_results$all_count,
        UP_count = da_results$UP_count,
        DOWN_count = da_results$DOWN_count,
        results = da_results$results
        # volcanoplot = volcanoplot_output
    ))
}


#' Compare Motif Activity Between Two Clusters
#'
#' This function performs pairwise comparison of transcription factor motif activities
#' between two specified clusters using ChromVAR scores. It identifies differentially
#' active motifs and generates visualization plots including volcano plots, motif
#' sequence logos, UMAP feature plots, and violin plots.
#'
#' @param seurat A Seurat object containing chromVAR and ATAC assays
#' @param cluster_id Character or numeric. The primary cluster identifier to compare
#' @param other_cluster_id Character or numeric. The comparison cluster identifier
#' @param identities Character. The metadata column to use for cell identities.
#'   Default is 'seurat_clusters_atac'
#' @param reduction Character. The dimensional reduction to use for feature plots.
#'   Default is 'umap.wnn'
#' @param topn Numeric. Number of top motifs to display for up and down-regulated
#'   activities. Default is 9
#' @param JASPARConnect A JASPAR database connection object for retrieving motif
#'   information and sequence logos
#' @param sequential_palette Color palette for feature plots
#' @param ... Additional arguments passed to volcano_plot function
#'
#' @return A data frame of differential motif activities with columns:
#'   \itemize{
#'     \item avg_diff: Average difference in chromVAR scores
#'     \item p_val: P-value from likelihood ratio test
#'     \item p_val_adj: Adjusted p-value
#'     \item pct.1: Percentage of cells in cluster_id
#'     \item pct.2: Percentage of cells in other_cluster_id
#'     \item motif_name: Human-readable motif name from JASPAR
#'   }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Runs FindMarkers using likelihood ratio test on chromVAR scores,
#'         controlling for ATAC fragment counts
#'   \item Identifies top up-regulated motifs (avg_diff > 0, padj < 0.05)
#'   \item Identifies top down-regulated motifs (avg_diff < 0, padj < 0.05)
#'   \item Generates a volcano plot of all differential motif activities
#'   \item For both up and down-regulated sets:
#'     \itemize{
#'       \item Displays motif sequence logos
#'       \item Creates UMAP feature plots showing motif activity
#'       \item Generates violin plots comparing activity between clusters
#'     }
#'   \item Annotates results with motif names from JASPAR database
#' }
#'
#' @note This function requires the following packages: Seurat, Signac, TFBSTools,
#'   patchwork, dplyr, purrr, scCustomize
#'
#' @examples
#' \dontrun{
#' # Compare motif activities between cluster 1 and cluster 2
#' diff_motifs <- compare_motive_activity_pairwise(
#'   seurat = my_seurat,
#'   cluster_id = 1,
#'   other_cluster_id = 2,
#'   JASPARConnect = jaspar_db,
#'   sequential_palette = viridis::viridis(100)
#' )
#' }
#'
#' @export

# Custom function, to be incorporated into the package later
compare_motive_activity_pairwise <- function(seurat, cluster_id, other_cluster_id, figures_path, identities='seurat_clusters_atac', reduction = 'umap.wnn', topn=9, JASPARConnect, sequential_palette , ...) {
  
  # Set the default assay and identities
  DefaultAssay(seurat) <- 'chromvar'
  Idents(seurat) <- identities  
  
  differential_activity <- FindMarkers(
    object = seurat,
    ident.1 = cluster_id,
    ident.2 = other_cluster_id,
    assay = 'chromvar',
    slot = 'data',
    only.pos = FALSE,
    mean.fxn = rowMeans,
    fc.name = 'avg_diff',
    # latent.vars = 'nCount_ATAC'
  )

  

  top_motif_activities_up <- differential_activity |> 
    filter(avg_diff > 0) |> 
    filter(p_val_adj < 0.05) |>
    arrange(p_val_adj) |>
    slice_head(n = topn) |>
    rownames_to_column(var = 'gene') |>
    pull(gene) |> 
    unique()

  top_motif_activities_down <- differential_activity |> 
    filter(avg_diff < 0) |> 
    filter(p_val_adj < 0.05) |>
    arrange(p_val_adj) |>
    slice_head(n = topn) |>
    rownames_to_column(var = 'gene') |>
    pull(gene) |> 
    unique()
  
  differential_activity_cluster <- differential_activity |> 
    rownames_to_column(var = 'gene') |>
    rename(log2FoldChange = avg_diff, pvalue = p_val, padj = p_val_adj) |>
    mutate(genes = map_chr(gene, ID_to_symbol)) 


  plot <- volcano_plot(differential_activity_cluster, group1 = other_cluster_id, group2 = cluster_id, cluster = cluster_id, local_figures_path = figures_path, FC_threshold = 0.25, p_value_threshold = 0.05, test_type = 'Wilcox_ChromVar_motif', ...)

    # Create a list of motif sets to iterate over
    motif_sets <- list(
      up = list(motifs = top_motif_activities_up, label = "UP"),
      down = list(motifs = top_motif_activities_down, label = "DOWN")
    )
    
    # Iterate over both up and down motif sets
    walk(motif_sets, function(motif_set) {
      motifs <- motif_set$motifs
      label <- motif_set$label
      
      # Create motif plot
      motif_plot <- MotifPlot(
        object = seurat,
        motifs = motifs,
        assay = 'ATAC'
      ) + labs(title = paste0('Top Motif Activities ', label, ' in Cluster ', cluster_id))
      
    print(motif_plot)
    ggsave(
      filename = here::here(figures_path, paste0('motif_sequences_', label, '_cluster_', cluster_id, '_vs_', other_cluster_id, '.pdf')),
      plot = motif_plot,
      width = 12,
      height = 8
    )
      
      # Get motif names
      top_motif_names <- map(motifs, function(motif_id) {
        motif_matrix <- getMatrixByID(
          x = JASPARConnect, 
          ID = motif_id   
        )
        return(motif_matrix@name)
      }) 
      
      seurat_local <- subset(seurat, idents = c(cluster_id, other_cluster_id))
      
      # Create feature plots
      plots <- map2(motifs, top_motif_names, function(TF_ID, TF_name) {
        p <- FeaturePlot_scCustom(seurat, feature = TF_ID, reduction = reduction, colors_use = sequential_palette) + 
          labs(title = paste0(TF_name)) + 
          NoAxes() + 
          NoLegend()
        return(p)
      })
    combined_plot <- wrap_plots(plots, axes = 'collect', axis_titles = 'collect', guides = 'auto') + 
        plot_annotation(title = paste0('Top Motif Activities ', label, ' in Cluster ', cluster_id))
    print(combined_plot)
    ggsave(
      filename = here::here(figures_path, paste0('motif_activities_', label, '_cluster_', cluster_id, '_vs_', other_cluster_id, '.pdf')),
      plot = combined_plot,
      width = 12,
      height = 10
    )
      
      # Create violin plots
      violin_plots <- map2(motifs, top_motif_names, function(TF_ID, TF_name) {
        p <- VlnPlot(seurat_local, assay = 'chromvar', features = TF_ID, group.by = identities, pt.size = 0) + 
          labs(title = paste0(TF_name)) + 
          NoLegend() + 
          theme(axis.title.x = element_blank(), axis.title.y = element_blank())
        return(p)
      })
    combined_violin_plot <- wrap_plots(violin_plots, axes = 'collect', axis_titles = 'collect', guides = 'auto') + 
        plot_annotation(title = paste0('Top Motif Activities ', label, ' in Cluster ', cluster_id))
    print(combined_violin_plot)
    ggsave(
      filename = here::here(figures_path, paste0('motif_activities_violin_', label, '_cluster_', cluster_id, '_vs_', other_cluster_id, '.pdf')),
      plot = combined_violin_plot,
      width = 12,
      height = 10
    )
    })
  
    # Update differential_activity rownames with motif names
    differential_activity <- differential_activity |>
      rownames_to_column(var = 'motif_id') |>
      mutate(
        motif_name = map_chr(motif_id, function(motif_id) {
            motif_matrix <- getMatrixByID(
              x = JASPARConnect,
              ID = motif_id
            )
            return(motif_matrix@name)
          })) |> 
      column_to_rownames(var = 'motif_id')      
  return(differential_activity)
}
