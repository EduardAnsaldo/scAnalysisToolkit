## Seurat Clustering Functions
## Functions for performing Seurat clustering and dimensionality reduction

#' Perform Complete Seurat Clustering Workflow
#'
#' Executes a comprehensive Seurat clustering pipeline including SCTransform
#' normalization, PCA, UMAP dimensionality reduction, clustering at multiple
#' resolutions, marker gene identification, and optional SingleR cell type annotation.
#' Generates extensive visualizations and optionally performs pathway enrichment.
#'
#' @param seurat Seurat object to cluster
#' @param npcs Integer; number of principal components to compute. Default 100
#' @param dimensions Integer; number of dimensions to use for UMAP and clustering. Default 30
#' @param resolutions Numeric vector; clustering resolutions to test. Default c(0.05, 0.1, 0.15, 0.25, 0.375, 0.5)
#' @param filter_variable Character or FALSE; metadata column to filter by. Default FALSE
#' @param keep_values Character vector; values to keep when filtering. Default NULL
#' @param plot_elbow Logical; whether to plot elbow plot. Default TRUE
#' @param seed Integer; random seed for reproducibility. Default 3514L
#' @param verbose Logical; whether to print detailed messages. Default FALSE
#' @param selected_resolution Numeric; resolution for detailed analysis. If NULL, only
#'   comparison plots are generated. Default NULL
#' @param object_annotations Character; string to append to output file names. Default 'full_object'
#' @param figures_path Character; path to save figures. Default NULL
#' @param tables_path Character; path to save tables. Default NULL
#' @param results_path Character; path for enrichment results. Default NULL
#' @param run_pathway_enrichment Character vector or FALSE; enrichment methods to run.
#'   Default FALSE
#' @param n_genes_to_plot Integer; number of marker genes per cluster to plot. Default 3
#' @param run_singler Logical; whether to run SingleR cell type annotation. Default FALSE
#' @param singler_database Character; SingleR reference database. Default 'ImmGen'
#'
#' @return List with elements:
#'   \item{seurat}{Processed Seurat object with clusters}
#'   \item{resolution_plot}{Grid of UMAP plots for all resolutions}
#'   \item{detailed_plots}{List of detailed plots for selected resolution or NULL}
#'   \item{top_genes_plot}{Marker gene analysis results or NULL}
#'   \item{dimensions}{Number of dimensions used}
#'   \item{resolutions}{Resolutions tested}
#'   \item{selected_resolution}{Resolution selected for detailed analysis}
#'
#' @examples
#' results <- perform_seurat_clustering(seurat, npcs = 50, dimensions = 30,
#'                                     selected_resolution = 0.25,
#'                                     run_pathway_enrichment = c("Metascape"))
#' @export
perform_seurat_clustering <- function(
    seurat,
    npcs = 100,
    dimensions = 30,
    resolutions = c(0.05, 0.1, 0.15, 0.25, 0.375, 0.5),
    filter_variable = FALSE,
    keep_values = NULL,
    plot_elbow = TRUE,
    seed = 3514L,
    verbose = FALSE,
    selected_resolution = NULL,
    object_annotations = 'full_object',
    figures_path = NULL,
    tables_path = NULL,
    results_path = NULL,
    run_pathway_enrichment = FALSE,
    n_genes_to_plot = 3,
    run_singler = FALSE,
    singler_database = 'ImmGen'
) {
    # Filter object
    if (!is.null(filter_variable)) {
        if (is.null(keep_values)) {
            stop("keep_values must be provided when filter_variable is specified")
        }
        if (!filter_variable %in% colnames(seurat@meta.data)) {
            stop(paste0("Column '", filter_variable, "' not found in seurat metadata"))
        }
        seurat <- subset(seurat, subset = !!sym(filter_variable) %in% keep_values)
        p1 <- DimPlot(seurat, label = T) + ggtitle(paste0('Filtered by ', filter_variable, ' values: ', paste(keep_values, collapse = ', '))) + theme(legend.position = "none")
        print(p1)
    }

    # Normalize RNA data with SCTransform
    seurat <- SCTransform(seurat, verbose = verbose)

    # Dimensionality reduction
    seurat <- RunPCA(seurat, npcs = npcs, verbose = verbose)

    # Plot elbow plot if requested
    if (plot_elbow) {
        print(ElbowPlot(seurat, ndims = npcs))
    }

    # Run UMAP
    seurat <- RunUMAP(seurat, dims = 1:dimensions, verbose = verbose, seed.use = seed)

    # Find neighbors and clusters
    seurat <- FindNeighbors(seurat, dims = 1:dimensions, verbose = verbose)
    seurat <- FindClusters(seurat, resolution = resolutions, verbose = verbose, algorithm = 4, random.seed = seed)

    # Create UMAP plots for each resolution
    p <- resolutions |>
        map(\(resolution) {
            Idents(seurat) <- paste0('SCT_snn_res.', resolution)
            DimPlot(seurat, reduction = "umap", label = TRUE) +
                ggtitle(paste0('R ', resolution)) +
                theme(legend.position = "none")
        })

    # Arrange plots in grid
    resolution_plot <- gridExtra::grid.arrange(grobs = p)

    # If a specific resolution is selected, generate detailed plots
    detailed_plots <- NULL
    top_genes_plot <- NULL
    if (!is.null(selected_resolution)) {
        Idents(seurat) <- paste0('SCT_snn_res.', selected_resolution)
        seurat[['seurat_clusters']] <- Idents(seurat)

        # Generate detailed plots
        p1 <- DimPlot(seurat, reduction = "umap", label = TRUE) +
            ggtitle(paste0('R ', selected_resolution)) +
            theme(legend.position = "none")

        if (!is.null(figures_path)) {
            ggsave(
                paste0('UMAP_clusters_R_', selected_resolution, '_', object_annotations, '.pdf'),
                plot = p1,
                path = figures_path
            )
        }

        p2 <- DimPlot(seurat, reduction = "umap", label = TRUE, split.by = 'Groups', ncol = 2) +
            ggtitle(paste0('R ', selected_resolution)) +
            theme(legend.position = "none")

        if (!is.null(figures_path)) {
            ggsave(
                paste0('UMAP_clusters_by_group_R_', selected_resolution, '_', object_annotations, '.pdf'),
                plot = p2,
                path = figures_path,
                width = 5,
                height = 5
            )
        }

        p3 <- DimPlot(seurat, reduction = "umap", label = FALSE, group.by = 'Groups') +
            ggtitle(paste0('R ', selected_resolution)) +
            theme(legend.position = "right")

        detailed_plots <- list(cluster_plot = p1, group_split_plot = p2, groups_plot = p3)
        walk(detailed_plots, print)

        # Top genes per cluster
        top_genes_plot <- top_genes_per_cluster(
            seurat,
            n_genes_to_plot = n_genes_to_plot,
            object_annotations = object_annotations,
            tables_path = tables_path,
            figures_path = figures_path,
            results_path = results_path,
            run_pathway_enrichment = run_pathway_enrichment
        )

        if (!is.null(figures_path)) {
            ggsave(
                paste0(figures_path, 'DotPlot_Top', n_genes_to_plot, '_per_cluster_', object_annotations, '.pdf'),
                plot = top_genes_plot$plot,
                width = 8,
                height = 14
            )
        }

        # SingleR annotations if requested
        if (run_singler) {
            local_path <- paste0(figures_path, object_annotations, '_cell_type_annotations')
            unlink(local_path, recursive = TRUE)
            dir.create(local_path)

            # Normalize and scale RNA data
            seurat <- NormalizeData(seurat, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)
            seurat <- ScaleData(seurat, assay = 'RNA')

            # Run SingleR annotations
            seurat <- annotate_seurat_with_SingleR_Eduard(
                seurat,
                local_path,
                database = singler_database,
                annotation_basis = 'cluster_coarse',
                split_by_groups = FALSE
            )

            annotate_seurat_with_SingleR_Eduard(
                seurat,
                local_path,
                database = singler_database,
                annotation_basis = 'cell_coarse',
                split_by_groups = FALSE
            )
        }
    }

    # Return list with processed seurat object and plots
    list(
        seurat = seurat,
        resolution_plot = resolution_plot,
        detailed_plots = detailed_plots,
        top_genes_plot = top_genes_plot,
        dimensions = dimensions,
        resolutions = resolutions,
        selected_resolution = selected_resolution
    )
}

## Cell Counting and Statistics Functions
## Functions for extracting and visualizing cell count statistics

#' Extract and Visualize Cell Count Statistics
#'
#' Calculates cell counts and frequencies for each cluster/cell type across samples
#' and experimental groups. Generates frequency tables and bar plots with error bars
#' showing mean frequency and variability across replicates.
#'
#' @param seurat Seurat object with cluster assignments and metadata
#' @param grouping_var Unquoted column name; metadata column for grouping cells
#'   (e.g., seurat_clusters, cell_type)
#' @param figures_path Character; path to save figures
#' @param tables_path Character; path to save frequency and count tables
#' @param object_annotations Character; string to append to output file names. Default ''
#'
#' @return Invisibly returns NULL. Creates output files:
#'   \itemize{
#'     \item Frequency table CSV: Percentage of each cluster per sample
#'     \item Count table CSV: Raw cell counts per cluster per sample
#'     \item Bar plot PDF: Mean frequency with error bars and individual points
#'   }
#'
#' @examples
#' extract_cell_counts(seurat, seurat_clusters, "./figures", "./tables",
#'                    object_annotations = "experiment1")
#' @export
extract_cell_counts <- function(seurat, grouping_var, figures_path, tables_path, object_annotations='') {
    cell_counts <- FetchData(seurat, vars = c(englue("{{grouping_var}}"), "Samples", "Groups"))
    cell_counts <- arrange(cell_counts, Samples)

    counts <- cell_counts %>% add_count(Samples, name='total_cell_count_by_sample')

    counts <- counts %>%
        dplyr::count(  {{grouping_var}},  , Samples, Groups,  total_cell_count_by_sample,name='cluster_count')  |>
            mutate(frequency_within_sample=cluster_count*100/total_cell_count_by_sample)  |>
            mutate(Samples = as.character(Samples)) |>
            arrange(Samples, desc(Samples))

    frequency_table <- counts |>
        arrange(Samples) |>
        pivot_wider(id_cols = {{grouping_var}},  names_from = 'Samples', values_from = frequency_within_sample)
    write.csv(frequency_table,file=here::here(tables_path, paste0(englue("frequency per {{grouping_var}} per condition "), object_annotations, ".csv")), row.names=F)

    count_table <- counts |>
        arrange(Samples) |>
        pivot_wider(id_cols = {{grouping_var}}, names_from = 'Samples', values_from = cluster_count)
    write.csv(count_table,file = here::here(tables_path, paste0(englue("counts per {{grouping_var}} per condition "), object_annotations, ".csv")),row.names=F)

    # # Barplot of proportion of cells in each cluster by sample
    # plot1 <- ggplot(seurat@meta.data) +
    #     geom_bar(aes(x=Groups, fill={{grouping_var}}), position=position_fill()) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    # ggsave(plot = plot1, filename = paste0(figures_path, englue("frequency per {{grouping_var}} per group"), object_annotations, ".pdf"))
    # print(plot1)

    # counts <- cell_counts %>% add_count(Groups, name='total_cell_count_by_sample')
    # counts <- counts %>%
    #     dplyr::count({{grouping_var}}, Groups, total_cell_count_by_sample,name='frequency_within_sample')  |>
    #         mutate(frequency_within_sample=frequency_within_sample*100/total_cell_count_by_sample)

    # Barplot of proportion of cells in each cluster by sample
    plot2 <- ggplot(counts, aes(x={{grouping_var}} |> fct_reorder(frequency_within_sample) |> fct_rev(), y = frequency_within_sample, fill=Groups)) +
        geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), show.legend = TRUE) +
        stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.9), aes(fill = Groups, alpha = 0.5)) +
        stat_summary(fun.data = mean_se, fun.args = list(mult = 1),
            geom = "errorbar", width = 0.2, position = position_dodge(width = 0.9) ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        theme_classic() +
        labs(x = 'Cell type', y = 'Frequency (%)', title = englue('Frequency per {{grouping_var}}'))+
        guides(alpha = 'none')+
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    ggsave(plot = plot2, filename =  paste0(englue('frequency_per_{{grouping_var}}_per_sample '), object_annotations, '.pdf'), width = 12, height = 6, path = figures_path)
    print(plot2)

#     # Tidyplot of proportion of cells in each cluster by sample
# #     counts <- mutate(
# #         {{grouping_var}} := fct_reorder({{grouping_var}}frequency_within_sample),
# #         {{grouping_var}} := fct_rev({{grouping_var}}frequency_within_sample)
# #   )
#     plot2 <- tidyplot(counts, x = {{grouping_var}}, y = frequency_within_sample, color = Groups)  |>
#       add_data_points_beeswarm() |>
#       add_mean_bar(alpha = 0.4) |>
#       add_sem_errorbar() |>
#       add_test_asterisks(test = 't_test', p.adjust.method = 'BH')  |>
#       save_plot(filename = here::here(figures_path, paste0(englue('frequency_per_{{grouping_var}}_per_sample_tidyplot '), object_annotations, '.pdf')))

}

