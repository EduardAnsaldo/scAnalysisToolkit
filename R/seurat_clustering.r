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
#' @param filter_ig Logical; whether to filter out immunoglobulin genes from marker
#'   gene analysis. Default FALSE
#' @param filter_tcr Logical; whether to filter out T cell receptor genes from marker
#'   gene analysis. Default FALSE
#' @param quiet_vdj Logical; whether to remove VDJ genes (Igh, Igk, Igl, Tra, Trb, Trd, Trg)
#'   from variable features after SCTransform. Default FALSE
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
    singler_database = 'ImmGen',
    filter_ig = FALSE,
    filter_tcr = FALSE,
    quiet_vdj = FALSE
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

    # Optionally quiet VDJ genes from variable features
    if (quiet_vdj) {
        if (verbose) {
            message("Variable features before VDJ filtering: ", length(VariableFeatures(seurat)))
        }
        VariableFeatures(seurat) <- VariableFeatures(seurat)[!stringr::str_detect(VariableFeatures(seurat), "Igh|Igl|Igk|Tra|Trb|Trd|Trg")]
        if (verbose) {
            message("Variable features after VDJ filtering: ", length(VariableFeatures(seurat)))
        }
    }

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
            run_pathway_enrichment = run_pathway_enrichment,
            filter_ig = filter_ig,
            filter_tcr = filter_tcr
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
#' @param plot_grouping Character; type of plot to generate. "Samples" creates a grouped
#'   bar plot with individual sample points (default), "Groups" creates a stacked bar plot
#'   showing the composition of experimental groups within each cluster/cell type. Default "Samples"
#'
#' @return Invisibly returns NULL. Creates output files:
#'   \itemize{
#'     \item Frequency table CSV: Percentage of each cluster per sample
#'     \item Count table CSV: Raw cell counts per cluster per sample
#'     \item Bar plot PDF: Mean frequency with error bars and individual points (Samples mode)
#'       or stacked bar plot showing composition (Groups mode)
#'   }
#'
#' @export
extract_cell_counts <- function(seurat, grouping_var, figures_path, tables_path, object_annotations='', plot_grouping = "Samples") {
    cell_counts <- FetchData(seurat, vars = c(englue("{{grouping_var}}"), "Samples", "Groups"))
    cell_counts <- arrange(cell_counts, Samples)

    counts <- cell_counts |> add_count(Samples, name='total_cell_count_by_sample')

    counts <- counts |>
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

    # Barplot of proportion of cells in each cluster
    if (plot_grouping == "Samples") {
        # Grouped bar plot with individual sample points
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
    } else if (plot_grouping == "Groups") {
        # Stacked bar plot showing composition of Groups within each grouping_var
        plot2 <- ggplot(seurat@meta.data, aes(x={{grouping_var}}, fill=Groups)) +
            geom_bar(position=position_fill()) +
            scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0, 0.01))) +
            theme_classic() +
            labs(x = englue('{{grouping_var}}'), y = 'Frequency (%)', title = englue('Composition of Groups per {{grouping_var}}'), fill = 'Groups') +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
        ggsave(plot = plot2, filename =  paste0(englue('composition_of_groups_per_{{grouping_var}} '), object_annotations, '.pdf'), width = 10, height = 6, path = figures_path)
        print(plot2)
    } else {
        stop("plot_grouping must be either 'Samples' or 'Groups'")
    }
}
#' Process Multimodal Seurat Object with RNA and ATAC Data
#'
#' This function performs comprehensive processing of multimodal single-cell data
#' including RNA and ATAC-seq modalities. It handles normalization, dimensionality
#' reduction, clustering, and weighted nearest neighbor (WNN) analysis. The function
#' can operate in interactive mode, stopping at key decision points to allow
#' examination of diagnostic plots.
#'
#' @param seurat A Seurat object containing both RNA and ATAC assays
#' @param figures_path Character string specifying the directory path where
#'   generated plots should be saved. Required if `generate_plots = TRUE`.
#' @param object_annotations Character string used as prefix for saved plot filenames.
#'   Required if `generate_plots = TRUE`.
#' @param dimensions_RNA Integer specifying the number of principal components to
#'   use for RNA data. If NULL, function stops after PCA to allow examination of
#'   elbow plot.
#' @param dimensions_ATAC Integer specifying the upper dimension for ATAC LSI
#'   components (dimension 1 is excluded). If NULL, function stops after LSI to
#'   allow examination of elbow plot.
#' @param resolution_RNA Numeric value specifying the clustering resolution for
#'   RNA data. If NULL, function stops after showing resolution exploration plots.
#' @param resolution_ATAC Numeric value specifying the clustering resolution for
#'   ATAC data. If NULL, function stops after showing resolution exploration plots.
#' @param resolution_WNN Numeric value specifying the clustering resolution for
#'   weighted nearest neighbor analysis. If NULL, function stops after showing
#'   resolution exploration plots.
#' @param resolutions_RNA Numeric vector of resolution values to explore for RNA
#'   clustering. Default: c(0.375, 0.4, 0.475, 0.5, 0.55, 0.6, 0.625, 0.75, 1)
#' @param resolutions_ATAC Numeric vector of resolution values to explore for ATAC
#'   clustering. Default: c(0.1, 0.25, 0.3, 0.325, 0.375, 0.5, 0.625, 0.75, 1)
#' @param resolutions_WNN Numeric vector of resolution values to explore for WNN
#'   clustering. Default: c(0.1, 0.25, 0.375, 0.5, 0.625, 0.75, 1, 1.25, 1.5)
#' @param npcs Integer specifying the number of principal components to compute
#'   for RNA PCA. Default: 50
#' @param ndims_lsi Integer specifying the number of LSI dimensions to compute
#'   for ATAC data. Default: 40
#' @param random_seed Integer seed for reproducibility in clustering algorithms.
#'   Default: 3514
#' @param remove_ig Logical indicating whether to filter out immunoglobulin-related
#'   features from both RNA and ATAC data. Default: FALSE
#' @param generate_plots Logical indicating whether to generate and save final
#'   plots. Requires `figures_path` and `object_annotations` to be specified.
#'   Default: TRUE
#' @param skip_preprocessing Logical indicating whether to skip RNA/ATAC normalization
#'   and dimensionality reduction steps. Use TRUE when rerunning after examining
#'   elbow plots. Default: FALSE
#' @param skip_umap_and_neighbors Logical indicating whether to skip UMAP generation
#'   and neighbor finding steps. Use TRUE when rerunning after selecting dimensions
#'   to only redo clustering. Default: FALSE
#'
#' @return A list containing:
#'   \item{seurat}{The processed Seurat object with added reductions and clustering}
#'   \item{stage}{Character string indicating completion stage: "dimensions",
#'     "resolutions", or "complete"}
#'   \item{missing}{Character vector of missing parameters (only if stage != "complete")}
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item RNA Processing: SCTransform normalization, removal of immunoglobulin
#'     genes from variable features, and PCA
#'   \item ATAC Processing: Feature selection, removal of immunoglobulin loci,
#'     TF-IDF normalization, and LSI dimensionality reduction
#'   \item Interactive dimension selection based on elbow plots
#'   \item UMAP generation for both modalities
#'   \item Clustering exploration across multiple resolutions
#'   \item Weighted Nearest Neighbor (WNN) analysis integrating both modalities
#'   \item Optional generation of comprehensive visualization plots
#' }
#'
#' When `remove_ig = TRUE`, the function automatically filters out immunoglobulin-related features:
#' \itemize{
#'   \item RNA: Genes matching patterns "Igh", "Igl", or "Igk"
#'   \item ATAC: Peaks overlapping genomic regions chr12:113222388-115937574 (Igh),
#'     chr6:67532620-70703738 (Igk), and chr16:18845608-19079594 (Igl)
#' }
#'
#' @note This function uses the Leiden algorithm (algorithm = 4) for clustering.
#'   The first LSI component is excluded from ATAC analysis as it typically
#'   correlates with sequencing depth.
#'
#' @examples
#' \dontrun{
#' # Interactive mode - examine plots at each stage
#' result <- process_multimodal_seurat(seurat_obj)
#' 
#' # After examining elbow plots, specify dimensions and skip preprocessing
#' result <- process_multimodal_seurat(
#'   seurat_obj,
#'   dimensions_RNA = 30,
#'   dimensions_ATAC = 25,
#'   skip_preprocessing = TRUE
#' )
#' 
#' # After examining clustering, specify resolutions and skip earlier steps
#' result <- process_multimodal_seurat(
#'   seurat_obj,
#'   dimensions_RNA = 30,
#'   dimensions_ATAC = 25,
#'   resolution_RNA = 0.5,
#'   resolution_ATAC = 0.3,
#'   resolution_WNN = 0.5,
#'   skip_preprocessing = TRUE,
#'   skip_umap_and_neighbors = TRUE
#' )
#' }
#'
#' @export
process_multimodal_seurat <- function(
  seurat,
  figures_path = NULL,
  object_annotations = NULL,
  # Dimension parameters (NULL = stop and show plots)
  dimensions_RNA = NULL,
  dimensions_ATAC = NULL,
  # Resolution parameters (NULL = stop and show plots)
  resolution_RNA = NULL,
  resolution_ATAC = NULL,
  resolution_WNN = NULL,
  # Resolution exploration ranges
  resolutions_RNA = c(0.375, 0.4, 0.475, 0.5, 0.55, 0.6, 0.625, 0.75, 1),
  resolutions_ATAC = c(0.1, 0.25, 0.3, 0.325, 0.375, 0.5, 0.625, 0.75, 1),
  resolutions_WNN = c(0.1, 0.25, 0.375, 0.5, 0.625, 0.75, 1, 1.25, 1.5),
  # Other parameters
  npcs = 50,
  ndims_lsi = 40,
  random_seed = 3514,
  remove_ig = FALSE,
  genome = "GRCm39",
  generate_plots = TRUE,
  # Skip parameters for iterative workflows
  skip_preprocessing = FALSE,
  skip_umap_and_neighbors = FALSE
) {
  
  # Define Ig regions based on genome
  ig_regions <- list(
    GRCm38 = list(
      Igh = 'chr12-113258768-116009954',
      Igk = 'chr6-67555636-70726754',
      Igl = 'chr16-19026858-19260844'
    ),
    GRCm39 = list(
      Igh = 'chr12-113222388-115937574',
      Igk = 'chr6-67532620-70703738',
      Igl = 'chr16-18845608-19079594'
    )
  )
  
  # Validate genome parameter
  if (!genome %in% names(ig_regions)) {
    stop("Genome must be one of: ", paste(names(ig_regions), collapse = ", "))
  }
  
  # ========== RNA Processing ==========
  if (!skip_preprocessing) {
    message("Processing RNA data...")
    DefaultAssay(seurat) <- 'RNA'
    seurat <- SCTransform(seurat)
    
    if (remove_ig) {
      message("Variable features before filtering: ", length(VariableFeatures(seurat)))
      VariableFeatures(seurat) <- VariableFeatures(seurat)[!str_detect(VariableFeatures(seurat), "Igh|Igl|Igk")]
      message("Variable features after filtering: ", length(VariableFeatures(seurat)))
    } else {
      message("Variable features (Ig filtering disabled): ", length(VariableFeatures(seurat)))
    }
    
    seurat <- RunPCA(seurat, npcs = npcs)
    
    # Always show elbow plot
    p <- ElbowPlot(seurat, reduction = 'pca', ndims = npcs) + 
      labs(title = "RNA PCA Elbow Plot")
    print(p)
  } else {
    message("Skipping RNA preprocessing (already completed)")
  }
  
  # ========== ATAC Processing ==========
  if (!skip_preprocessing) {
    message("\nProcessing ATAC data...")
    DefaultAssay(seurat) <- 'ATAC'
    seurat <- FindTopFeatures(seurat, min.cutoff = 5)
    
    if (remove_ig) {
      message("Variable features before filtering: ", length(VariableFeatures(seurat)))
      message("Using ", genome, " genome coordinates for Ig region filtering")
      
      # Get Ig regions for the specified genome
      selected_regions <- ig_regions[[genome]]
      Ig_peaks <- IntersectMatrix(
        seurat@assays$ATAC$counts, 
        regions = c(selected_regions$Igh, selected_regions$Igk, selected_regions$Igl)
      ) |> rownames()
      
      VariableFeatures(seurat) <- VariableFeatures(seurat)[!(VariableFeatures(seurat) %in% Ig_peaks)]
      message("Variable features after filtering: ", length(VariableFeatures(seurat)))
    } else {
      message("Variable features (Ig filtering disabled): ", length(VariableFeatures(seurat)))
    }
    
    seurat <- RunTFIDF(seurat)
    seurat <- RunSVD(seurat)
    
    # Always show elbow plot and depth correlation
    p1 <- ElbowPlot(seurat, reduction = 'lsi', ndims = ndims_lsi) + 
      labs(title = "ATAC LSI Elbow Plot") + 
      scale_y_log10()
    print(p1)
    print(DepthCor(seurat))
  } else {
    message("Skipping ATAC preprocessing (already completed)")
  }
  
  # Stop if dimensions not provided
  if (is.null(dimensions_RNA) || is.null(dimensions_ATAC)) {
    message("\n=== STOPPING: Please examine the elbow plots and provide missing dimensions ===")
    if (is.null(dimensions_RNA)) {
      message("Missing: dimensions_RNA")
    }
    if (is.null(dimensions_ATAC)) {
      message("Missing: dimensions_ATAC (typically exclude the first dimension)")
    }
    message("\nTo continue, rerun with skip_preprocessing = TRUE and provide the dimension parameters")
    return(list(seurat = seurat, stage = "dimensions"))
  }
  
  # ========== RNA UMAP ==========
  if (!skip_umap_and_neighbors) {
    message("\nUsing ", dimensions_RNA, " dimensions for RNA")
    DefaultAssay(seurat) <- 'RNA'
    seurat <- RunUMAP(seurat, reduction = 'pca', dims = 1:dimensions_RNA, 
              assay = 'SCT', reduction.name = 'umap.rna', 
              reduction.key = 'rnaUMAP_')
  } else {
    message("\nSkipping RNA UMAP generation (already completed)")
  }
  
  # ========== ATAC UMAP ==========
  if (!skip_umap_and_neighbors) {
    message("Using dimensions 2:", dimensions_ATAC, " for ATAC (excluding first dimension)")
    DefaultAssay(seurat) <- 'ATAC'
    seurat <- RunUMAP(seurat, reduction = 'lsi', dims = 2:dimensions_ATAC, 
              assay = 'ATAC', reduction.name = 'umap.atac', 
              reduction.key = 'atacUMAP_')
  } else {
    message("Skipping ATAC UMAP generation (already completed)")
  }
  
  # Track whether we need to stop
  missing_resolutions <- c()
  
  # ========== RNA Clustering ==========
  if (!skip_umap_and_neighbors) {
    message("\nClustering RNA data...")
    DefaultAssay(seurat) <- 'RNA'
    Idents(seurat) <- 'Groups'
    seurat <- FindNeighbors(seurat, reduction = 'pca', dims = 1:dimensions_RNA)
    seurat <- FindClusters(seurat, graph.name = 'SCT_snn', resolution = resolutions_RNA, 
                verbose = FALSE, algorithm = 4, random.seed = random_seed)
  } else {
    message("\nSkipping RNA neighbor finding (already completed)")
    message("Re-clustering RNA data...")
    DefaultAssay(seurat) <- 'RNA'
    seurat <- FindClusters(seurat, graph.name = 'SCT_snn', resolution = resolutions_RNA, 
                verbose = FALSE, algorithm = 4, random.seed = random_seed)
  }
  
  # Always show RNA resolution plots
  p <- list()
  for (i in seq_along(resolutions_RNA)) {
    Idents(seurat) <- paste0('SCT_snn_res.', resolutions_RNA[i])
    p[[i]] <- DimPlot(seurat, reduction = "umap.rna", label = TRUE) + 
      ggtitle(paste0('R ', resolutions_RNA[i])) + 
      theme(legend.position = "none")
  }
  plot <- wrap_plots(p, axes = 'collect', guides = 'collect')
  print(plot)
  
  if (is.null(resolution_RNA)) {
    missing_resolutions <- c(missing_resolutions, "resolution_RNA")
  } else {
    message("Using resolution ", resolution_RNA, " for RNA")
    Idents(seurat) <- paste0('SCT_snn_res.', resolution_RNA)
    seurat$seurat_clusters_rna <- Idents(seurat)
    seurat@meta.data <- seurat@meta.data |>
      dplyr::select(-starts_with('SCT_snn_res')) |>
      mutate(seurat_clusters_rna = fct_inseq(seurat_clusters_rna))
  }
  
  # ========== ATAC Clustering ==========
  if (!skip_umap_and_neighbors) {
    message("\nClustering ATAC data...")
    DefaultAssay(seurat) <- 'ATAC'
    seurat <- FindNeighbors(seurat, reduction = 'lsi', dims = 2:dimensions_ATAC)
    seurat <- FindClusters(seurat, graph.name = 'ATAC_snn', resolution = resolutions_ATAC, 
                verbose = FALSE, algorithm = 4, random.seed = random_seed)
  } else {
    message("\nSkipping ATAC neighbor finding (already completed)")
    message("Re-clustering ATAC data...")
    DefaultAssay(seurat) <- 'ATAC'
    seurat <- FindClusters(seurat, graph.name = 'ATAC_snn', resolution = resolutions_ATAC, 
                verbose = FALSE, algorithm = 4, random.seed = random_seed)
  }
  
  # Always show ATAC resolution plots
  p <- list()
  for (i in seq_along(resolutions_ATAC)) {
    Idents(seurat) <- paste0('ATAC_snn_res.', resolutions_ATAC[i])
    p[[i]] <- DimPlot(seurat, reduction = "umap.atac", label = TRUE) + 
      ggtitle(paste0('R ', resolutions_ATAC[i])) + 
      theme(legend.position = "none")
  }
  plot <- wrap_plots(p, axes = 'collect', guides = 'collect')
  print(plot)
  
  if (is.null(resolution_ATAC)) {
    missing_resolutions <- c(missing_resolutions, "resolution_ATAC")
  } else {
    message("Using resolution ", resolution_ATAC, " for ATAC")
    Idents(seurat) <- paste0('ATAC_snn_res.', resolution_ATAC)
    seurat$seurat_clusters_atac <- Idents(seurat)
    seurat@meta.data <- seurat@meta.data |>
      dplyr::select(-starts_with('ATAC_snn_res')) |>
      mutate(seurat_clusters_atac = fct_inseq(seurat_clusters_atac))
  }
  
  # ========== WNN Analysis ==========
  if (!skip_umap_and_neighbors) {
    message("\nPerforming Weighted Nearest Neighbor Analysis...")
    seurat <- FindMultiModalNeighbors(
      seurat,
      reduction.list = list("pca", "lsi"),
      dims.list = list(1:dimensions_RNA, 2:dimensions_ATAC),
      verbose = TRUE
    )
    
    seurat <- RunUMAP(seurat, nn.name = "weighted.nn", assay = "SCT", 
              reduction.name = "umap.wnn", reduction.key = "wnnUMAP_")
    
    seurat <- FindClusters(seurat, graph.name = 'wsnn', resolution = resolutions_WNN, 
                verbose = FALSE, algorithm = 4, random.seed = random_seed)
  } else {
    message("\nSkipping WNN neighbor finding (already completed)")
    message("Re-clustering WNN data...")
    seurat <- FindClusters(seurat, graph.name = 'wsnn', resolution = resolutions_WNN, 
                verbose = FALSE, algorithm = 4, random.seed = random_seed)
  }
  
  # Always show WNN resolution plots
  p <- list()
  for (i in seq_along(resolutions_WNN)) {
    Idents(seurat) <- paste0('wsnn_res.', resolutions_WNN[i])
    p[[i]] <- DimPlot(seurat, reduction = "umap.wnn", label = TRUE) + 
      ggtitle(paste0('R ', resolutions_WNN[i])) + 
      theme(legend.position = "none")
  }
  plot <- wrap_plots(p, axes = 'collect', guides = 'collect')
  print(plot)
  
  if (is.null(resolution_WNN)) {
    missing_resolutions <- c(missing_resolutions, "resolution_WNN")
  } else {
    message("Using resolution ", resolution_WNN, " for WNN")
    Idents(seurat) <- paste0('wsnn_res.', resolution_WNN)
    seurat$seurat_clusters_wnn <- Idents(seurat)
    seurat@meta.data <- seurat@meta.data |>
      dplyr::select(-ends_with('seurat_clusters')) |>
      dplyr::select(-starts_with('wsnn_res')) |>
      mutate(seurat_clusters_wnn = fct_inseq(seurat_clusters_wnn))
  }
  
  # Stop if any resolutions are missing
  if (length(missing_resolutions) > 0) {
    message("\n=== STOPPING: Please examine clustering results and provide missing resolutions ===")
    message("Missing: ", paste(missing_resolutions, collapse = ", "))
    message("\nTo continue, rerun with skip_preprocessing = TRUE and skip_umap_and_neighbors = TRUE")
    return(list(seurat = seurat, stage = "resolutions", missing = missing_resolutions))
  }
  
  # ========== Generate Final Plots ==========
  if (generate_plots && !is.null(figures_path) && !is.null(object_annotations)) {
    message("\nGenerating and saving final plots...")
    
    # RNA plots
    Idents(seurat) <- 'seurat_clusters_rna'
    p1 <- DimPlot(seurat, reduction = 'umap.rna', label = TRUE) + 
      labs(title = paste0('RNA UMAP R ', resolution_RNA)) + 
      theme(legend.position = "none")
    print(p1)
    ggsave(file = paste0(figures_path, object_annotations, '_RNA_UMAP_R', resolution_RNA, '.pdf'))
    
    p2 <- DimPlot(seurat, reduction = "umap.rna", group.by = 'Groups', label = TRUE, shuffle = TRUE) + 
      labs(title = paste0('RNA UMAP R ', resolution_RNA), subtitle = 'RNA Leiden Clusters')
    print(p2)
    ggsave(file = paste0(figures_path, object_annotations, '_RNA_UMAP_R', resolution_RNA, '_by_groups.pdf'))
    
    # ATAC plots
    Idents(seurat) <- 'seurat_clusters_atac'
    p3 <- DimPlot(seurat, reduction = 'umap.atac', label = TRUE) + 
      labs(title = paste0('ATAC UMAP R ', resolution_ATAC)) + 
      theme(legend.position = "none")
    print(p3)
    ggsave(file = paste0(figures_path, object_annotations, '_ATAC_UMAP_R', resolution_ATAC, '.pdf'))
    
    p4 <- DimPlot(seurat, reduction = "umap.atac", group.by = 'Groups', label = TRUE, shuffle = TRUE) + 
      labs(title = paste0('ATAC UMAP R ', resolution_ATAC), subtitle = 'ATAC Leiden Clusters')
    print(p4)
    ggsave(file = paste0(figures_path, object_annotations, '_ATAC_UMAP_R', resolution_ATAC, '_by_groups.pdf'))
    
    # WNN plots
    Idents(seurat) <- 'seurat_clusters_wnn'
    p5 <- DimPlot(seurat, reduction = 'umap.wnn', label = TRUE) + 
      labs(title = paste0('WNN (joint) UMAP R ', resolution_WNN)) + 
      theme(legend.position = "none")
    print(p5)
    ggsave(file = paste0(figures_path, object_annotations, '_WNN_UMAP_R', resolution_WNN, '.pdf'))
    
    p6 <- DimPlot(seurat, reduction = "umap.wnn", group.by = 'Groups', label = TRUE, shuffle = TRUE) + 
      labs(title = paste0('WNN (joint) UMAP R ', resolution_WNN), subtitle = 'WNN Leiden Clusters')
    print(p6)
    ggsave(file = paste0(figures_path, object_annotations, '_WNN_UMAP_R', resolution_WNN, '_by_groups.pdf'))
    
    # Split plots
    Idents(seurat) <- 'seurat_clusters_rna'
    p7 <- DimPlot(seurat, reduction = 'umap.rna', label = TRUE, split.by = 'Groups') + 
      labs(title = 'RNA UMAP', subtitle = 'RNA Leiden Clusters') + 
      theme(legend.position = "none")
    print(p7)
    ggsave(file = paste0(figures_path, object_annotations, '_RNA_UMAP_split_by_groups.pdf'), width = 8, height = 5)
    
    Idents(seurat) <- 'seurat_clusters_atac'
    p8 <- DimPlot(seurat, reduction = 'umap.atac', label = TRUE, split.by = 'Groups') + 
      labs(title = 'ATAC UMAP', subtitle = 'ATAC Leiden Clusters') + 
      theme(legend.position = "none")
    print(p8)
    ggsave(file = paste0(figures_path, object_annotations, '_ATAC_UMAP_split_by_groups.pdf'), width = 8, height = 5)
    
    Idents(seurat) <- 'seurat_clusters_wnn'
    p9 <- DimPlot(seurat, reduction = 'umap.wnn', label = TRUE, split.by = 'Groups') + 
      labs(title = 'WNN (joint) UMAP', subtitle = 'WNN Leiden Clusters') + 
      theme(legend.position = "none")
    print(p9)
    ggsave(file = paste0(figures_path, object_annotations, '_WNN_UMAP_split_by_groups.pdf'), width = 8, height = 5)
  }
  
  message("\n=== COMPLETE: All processing finished ===")
  return(list(seurat = seurat, stage = "complete"))
}
