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

    sequential_palette_dotplot <- grDevices::hcl.colors(n = 20,'YlGn',rev = T)

    # Set the identity class for clustering
    Idents(seurat) <- grouping_var
    # Set the default assay to 'ATAC'
    DefaultAssay(seurat) <- 'ATAC'

    da_peaks <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)

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
        arrange(desc(avg_log2FC), .by_group = TRUE) |>
        slice_head(n = 10) -> top10

    #Top25 markers
    da_peaks |>
        group_by(cluster) |>
        arrange(desc(avg_log2FC), .by_group = TRUE) |>
        slice_head(n = 50) -> top50

    #Top100 markers
    da_peaks |>
        group_by(cluster) |>
        arrange(desc(avg_log2FC), .by_group = TRUE) |>
        slice_head(n = 100) -> top100

    #Topn markers
    da_peaks |>
        group_by(cluster) |>
        arrange(desc(avg_log2FC), .by_group = TRUE) |>
        slice_head(n = n_genes_to_plot) -> topn

    # Save the top markers to files
    write.table(top100,file=here::here(tables_path, paste0('top100_peaks', '_',object_annotations, ".tsv")), sep="\t",row.names = FALSE)
    # write.table(top25,file=here::here(path,'top25',object_annotations, ".tsv"), sep="\t",row.names = FALSE)
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
#' FindMarkers with Wilcoxon rank-sum test. Maps differentially accessible peaks
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

    print(paste('Cluster', cluster))

    group1 <- fixed(group1)
    group2 <- fixed(group2)

    Idents(seurat) <- comparison

    print('number of cells in group 1')
    print(seurat@meta.data |> filter(str_detect(!!as.name(comparison), group1)) |> nrow())
    print('number of cells in group 2')
    print(seurat@meta.data |> filter(str_detect(!!as.name(comparison), group2)) |> nrow())

    # Check there are enough cells
    n_group1 <- seurat@meta.data |> filter(str_detect(!!as.name(comparison), group1)) |> nrow()
    n_group2 <- seurat@meta.data |> filter(str_detect(!!as.name(comparison), group2)) |> nrow()

    if (n_group1 < minimum_cell_number | n_group2 < minimum_cell_number) {
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
                          assay = 'ATAC', slot = 'data', test.use = 'wilcox', min.pct = min_fraction, only.pos = F, logfc.threshold = 0.5)
    print(head(results |> filter(p_val < 0.05 & (avg_log2FC)  < -0.5)))


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
    write.csv(results_filtered_UP |> arrange(desc(log2FoldChange)),
             file = here::here(gene_lists_path, paste('PEAKS_UP_in', group2, cluster, 'Wilcox', group2, 'vs', group1, '.csv', sep = '_')))
    write.csv(results_filtered_DOWN |> arrange(log2FoldChange),
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
