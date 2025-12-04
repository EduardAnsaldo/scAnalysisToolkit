## Differential Expression Wrapper Functions
## Functions that orchestrate DE analysis with visualization and functional analysis

#' Complete Pseudobulk DE Analysis Pipeline
#'
#' Orchestrates a full pseudobulk differential expression analysis workflow including
#' statistical testing with DESeq2, scatter and volcano plot generation, and optional
#' pathway enrichment analysis. Returns both statistical results and visualizations.
#'
#' @param scRNAseq Seurat object containing single-cell RNA-seq data
#' @param comparison Character; metadata column name for comparison groups
#' @param group1 Character; pattern to match for first group (control)
#' @param group2 Character; pattern to match for second group (treatment)
#' @param cluster Character; cluster or cell type identifier. Default 'all_clusters'
#' @param path Character; output directory path. Default './'
#' @param FC_threshold Numeric; log2 fold change threshold for filtering. Default 0.3
#' @param p_value_threshold Numeric; adjusted p-value threshold. Default 0.05
#' @param expression_threshold_for_gene_list Numeric; minimum average normalized counts
#'   for filtering gene lists. Default 20
#' @param minimum_cell_number Integer; minimum cells required per group. Default 10
#' @param run_pathway_enrichment Character; method for pathway enrichment ('Metascape',
#'   'ClusterProfiler', or NULL). Default NULL
#' @param genes_to_exclude Character vector; genes to exclude from analysis. Default c()
#' @param ... Additional arguments passed to visualization and enrichment functions
#'
#' @return List with elements:
#'   \item{all_count}{Integer; number of significant DEGs}
#'   \item{UP_count}{Integer; number of upregulated genes}
#'   \item{DOWN_count}{Integer; number of downregulated genes}
#'   \item{results}{Data frame; complete DE results}
#'   \item{scatterplot}{ggplot; scatterplot object}
#'   \item{volcanoplot}{ggplot; volcano plot object}
#'
#' @export
pseudobulk <- function(scRNAseq, comparison, group1, group2, cluster = 'all_clusters',
                       path = './', FC_threshold = 0.3, p_value_threshold = 0.05,
                       expression_threshold_for_gene_list = 20, minimum_cell_number = 10,
                       run_pathway_enrichment = NULL, genes_to_exclude = c(), ...) {

    local_figures_path <- here::here(path, 'figures')
    dir.create(local_figures_path, showWarnings = FALSE, recursive = TRUE)

  group1 <- fixed(group1)
  group2 <- fixed(group2)

    # Run core DE analysis
    de_results <- pseudobulk_de(
        scRNAseq = scRNAseq,
        comparison = comparison,
        group1 = group1,
        group2 = group2,
        cluster = cluster,
        path = path,
        FC_threshold = FC_threshold,
        p_value_threshold = p_value_threshold,
        expression_threshold_for_gene_list = expression_threshold_for_gene_list,
        minimum_cell_number = minimum_cell_number,
        genes_to_exclude = genes_to_exclude,
        ...
    )

    # Check if analysis was successful
    if (is.null(de_results$results)) {
        return(list(
            all_count = de_results$all_count,
            UP_count = de_results$UP_count,
            DOWN_count = de_results$DOWN_count
        ))
    }

    # Generate scatter and volcano plots
    scatterplot_output <- scatterplot(
        results = de_results$results,
        group1 = group1,
        group2 = group2,
        local_figures_path = local_figures_path,
        FC_threshold = FC_threshold,
        p_value_threshold = p_value_threshold,
        cluster = cluster,
        test_type = 'Pseudobulk',
        ...
    )

    volcanoplot_output <- volcano_plot(
        results = de_results$results,
        group1 = group1,
        group2 = group2,
        cluster = cluster,
        local_figures_path = local_figures_path,
        FC_threshold = FC_threshold,
        p_value_threshold = p_value_threshold,
        test_type = 'Pseudobulk',
        ...
    )

    # Overrepresentation analysis
    run_DEG_functional_analysis(
        results = de_results$results,
        method = run_pathway_enrichment,
        grouping_var = cluster,
        path = path,
        FC_threshold = FC_threshold,
        p_value_threshold = p_value_threshold,
        group1 = group1,
        group2 = group2,
        ...
    )

    return(list(
        all_count = de_results$all_count,
        UP_count = de_results$UP_count,
        DOWN_count = de_results$DOWN_count,
        results = de_results$results,
        scatterplot = scatterplot_output,
        volcanoplot = volcanoplot_output
    ))
}


#' Complete Wilcoxon DE Analysis Pipeline for RNA Assay
#'
#' Orchestrates a full Wilcoxon differential expression analysis workflow for the
#' RNA assay including statistical testing, scatter and volcano plot generation,
#' pathway enrichment analysis, and optional custom gene list visualization.
#'
#' @param scRNAseq Seurat object containing single-cell RNA-seq data
#' @param comparison Character; metadata column name for comparison groups
#' @param group1 Character; pattern to match for first group (control)
#' @param group2 Character; pattern to match for second group (treatment)
#' @param cluster Character; cluster or cell type identifier. Default 'all_clusters'
#' @param path Character; output directory path. Default './'
#' @param FC_threshold Numeric; log2 fold change threshold for filtering. Default 0.3
#' @param p_value_threshold Numeric; adjusted p-value threshold. Default 0.05
#' @param max_overlaps Integer; maximum label overlaps in plots. Default 15
#' @param label_size Numeric; size of gene labels in plots. Default 5
#' @param pathways_of_interest Character vector; specific pathways to analyze. Default NULL
#' @param label_threshold Numeric; expression threshold for labeling genes. Default 100000
#' @param distance_from_diagonal_threshold Numeric; minimum distance from diagonal
#'   for labeling in scatterplot. Default 0.7
#' @param gene_lists_to_plot Named list; custom gene lists to visualize. Default NULL
#' @param expression_threshold_for_gene_list Numeric; minimum expression for gene list
#'   filtering. Default 20
#' @param colors Character vector of length 2; colors for DOWN and UP genes.
#'   Default c('green4', 'darkorchid4')
#' @param minimum_cell_number Integer; minimum cells required per group. Default 30
#' @param run_pathway_enrichment Logical; whether to run pathway enrichment. Default TRUE
#'
#' @return Named vector with elements:
#'   \item{all_count}{Integer; number of significant DEGs}
#'   \item{UP_count}{Integer; number of upregulated genes}
#'   \item{DOWN_count}{Integer; number of downregulated genes}
#'
#' @export
DEG_FindMarkers_RNA_assay <- function(scRNAseq, comparison, group1, group2, cluster = 'all_clusters',
                                      path = './', FC_threshold = 0.3, p_value_threshold = 0.05,
                                      max_overlaps = 15, label_size = 5, pathways_of_interest = NULL,
                                      label_threshold = 100000, distance_from_diagonal_threshold = 0.7,
                                      gene_lists_to_plot = NULL, expression_threshold_for_gene_list = 20,
                                      colors = c('green4', 'darkorchid4'), minimum_cell_number = 30,
                                      run_pathway_enrichment = TRUE) {

    local_figures_path <- here::here(path, 'figures')
    dir.create(local_figures_path, showWarnings = FALSE, recursive = TRUE)

    group1 <- fixed(group1)
    group2 <- fixed(group2)

    # Set colors for the plot
    my_colors <- create_deg_colors(colors)

    # Run core DE analysis
    de_results <- DEG_FindMarkers_RNA_assay_de(
        scRNAseq = scRNAseq,
        comparison = comparison,
        group1 = group1,
        group2 = group2,
        cluster = cluster,
        path = path,
        FC_threshold = FC_threshold,
        p_value_threshold = p_value_threshold,
        minimum_cell_number = minimum_cell_number
    )

    # Check if analysis was successful
    if (is.null(de_results$results)) {
        return(c(
            all_count = de_results$all_count,
            UP_count = de_results$UP_count,
            DOWN_count = de_results$DOWN_count
        ))
    }

    # Generate scatter and volcano plots
    scatterplot(results = de_results$results, group1 = group1, group2 = group2,
               cluster = cluster, my_colors = my_colors,
               local_figures_path = local_figures_path,
               FC_threshold = FC_threshold, p_value_threshold = p_value_threshold,
               max_overlaps = max_overlaps, label_size = label_size,
               label_threshold = label_threshold,
               distance_from_diagonal_threshold = distance_from_diagonal_threshold,
               test_type = 'Wilcox')

    volcano_plot(results = de_results$results, group1 = group1, group2 = group2,
                cluster = cluster, my_colors = my_colors,
                local_figures_path = local_figures_path,
                FC_threshold = FC_threshold, p_value_threshold = p_value_threshold,
                max_overlaps = max_overlaps, label_size = label_size,
                label_threshold = label_threshold, test_type = 'Wilcox')

    # Pathway enrichment analysis
    run_standard_enrichment(
        results = de_results$results,
        run_pathway_enrichment = run_pathway_enrichment,
        pathways_of_interest = pathways_of_interest,
        cluster = cluster,
        path = path,
        group1 = group1,
        group2 = group2,
        comparison = comparison,
        FC_threshold = FC_threshold
    )

    # Plot specific gene lists
    plot_gene_lists(
        results = de_results$results,
        gene_lists_to_plot = gene_lists_to_plot,
        group1 = group1,
        group2 = group2,
        cluster = cluster,
        my_colors = my_colors,
        local_figures_path = local_figures_path,
        FC_threshold = FC_threshold,
        p_value_threshold = p_value_threshold,
        test_type = 'Wilcox',
        max_overlaps = max_overlaps,
        label_size = label_size,
        label_threshold = label_threshold,
        distance_from_diagonal_threshold = distance_from_diagonal_threshold
    )

    return(c(
        all_count = de_results$all_count,
        UP_count = de_results$UP_count,
        DOWN_count = de_results$DOWN_count
    ))
}


#' Complete Wilcoxon DE Analysis Pipeline for SCT Assay
#'
#' Orchestrates a full Wilcoxon differential expression analysis workflow for
#' SCTransform-normalized data including statistical testing, scatter and volcano
#' plot generation, pathway enrichment analysis, and optional custom gene list
#' visualization. Handles integrated datasets appropriately.
#'
#' @param scRNAseq Seurat object containing SCT-normalized single-cell RNA-seq data
#' @param comparison Character; metadata column name for comparison groups
#' @param group1 Character; pattern to match for first group (control)
#' @param group2 Character; pattern to match for second group (treatment)
#' @param is_integrated_subset Logical; if TRUE, skips UMI recorrection. Default FALSE
#' @param cluster Character; cluster or cell type identifier. Default 'all_clusters'
#' @param path Character; output directory path. Default './'
#' @param FC_threshold Numeric; log2 fold change threshold for filtering. Default 0.3
#' @param p_value_threshold Numeric; adjusted p-value threshold. Default 0.05
#' @param min_fraction Numeric; minimum fraction of cells expressing a gene. Default 0.01
#' @param max_overlaps Integer; maximum label overlaps in plots. Default 15
#' @param label_size Numeric; size of gene labels in plots. Default 5
#' @param pathways_of_interest Character vector; specific pathways to analyze. Default NULL
#' @param label_threshold Numeric; expression threshold for labeling genes. Default 100000
#' @param distance_from_diagonal_threshold Numeric; minimum distance from diagonal
#'   for labeling in scatterplot. Default 0.7
#' @param gene_lists_to_plot Named list; custom gene lists to visualize. Default NULL
#' @param expression_threshold_for_gene_list Numeric; minimum expression for gene list
#'   filtering. Default 20
#' @param colors Character vector of length 2; colors for DOWN and UP genes.
#'   Default c('green4', 'darkorchid4')
#' @param minimum_cell_number Integer; minimum cells required per group. Default 30
#' @param run_pathway_enrichment Character; method for pathway enrichment or NULL. Default NULL
#' @param ... Additional arguments passed to core analysis and visualization functions
#'
#' @return List with elements:
#'   \item{all_count}{Integer; number of significant DEGs}
#'   \item{UP_count}{Integer; number of upregulated genes}
#'   \item{DOWN_count}{Integer; number of downregulated genes}
#'   \item{results}{Data frame; complete DE results}
#'   \item{scatterplot}{ggplot; scatterplot object}
#'   \item{volcanoplot}{ggplot; volcano plot object}
#'
#' @export
DEG_FindMarkers_SCT_assay <- function(scRNAseq, comparison, group1, group2, is_integrated_subset = FALSE,
                                      cluster = 'all_clusters', path = './', FC_threshold = 0.3,
                                      p_value_threshold = 0.05, min_fraction = 0.01, max_overlaps = 15, label_size = 5,
                                      pathways_of_interest = NULL, label_threshold = 100000,
                                      distance_from_diagonal_threshold = 0.7, gene_lists_to_plot = NULL,
                                      expression_threshold_for_gene_list = 20, colors = c('green4', 'darkorchid4'),
                                      minimum_cell_number = 30, run_pathway_enrichment = NULL, ...) {

    local_figures_path <- here::here(path, 'figures')
    dir.create(local_figures_path, showWarnings = FALSE, recursive = TRUE)

    group1 <- fixed(group1)
    group2 <- fixed(group2)

    # Set colors for the plot
    my_colors <- create_deg_colors(colors)

    # Run core DE analysis
    de_results <- DEG_FindMarkers_SCT_assay_de(
        scRNAseq = scRNAseq,
        comparison = comparison,
        group1 = group1,
        group2 = group2,
        is_integrated_subset = is_integrated_subset,
        cluster = cluster,
        path = path,
        FC_threshold = FC_threshold,
        p_value_threshold = p_value_threshold,
        min_fraction = min_fraction,
        expression_threshold_for_gene_list = expression_threshold_for_gene_list,
        minimum_cell_number = minimum_cell_number,
        ...
    )

    # Check if analysis was successful
    if (is.null(de_results$results)) {
        return(list(
            all_count = de_results$all_count,
            UP_count = de_results$UP_count,
            DOWN_count = de_results$DOWN_count
        ))
    }

    # Generate scatter and volcano plots
    scatterplot_output <- scatterplot(
        results = de_results$results,
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
        distance_from_diagonal_threshold = distance_from_diagonal_threshold,
        test_type = 'Wilcox',
        ...
    )

    volcanoplot_output <- volcano_plot(
        results = de_results$results,
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
        test_type = 'Wilcox',
        ...
    )

    # Overrepresentation analysis
    run_DEG_functional_analysis(
        results = de_results$results,
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

    # Plotting individual genes of interest
    if (!is.null(gene_lists_to_plot)) {
        for (gene_list in names(gene_lists_to_plot)) {
            genes_to_plot <- gene_lists_to_plot[[gene_list]]
            print(genes_to_plot)

            # Generate scatter and volcano plots for gene lists
            results_scatter <- scatterplot(
                genes_to_plot = genes_to_plot,
                gene_list_name = gene_list,
                results = de_results$results,
                group1 = group1,
                group2 = group2,
                cluster = cluster,
                my_colors = my_colors,
                local_figures_path = local_figures_path,
                FC_threshold = FC_threshold,
                p_value_threshold = p_value_threshold,
                max_overlaps = 15,
                label_size = label_size,
                label_threshold = label_threshold,
                distance_from_diagonal_threshold = distance_from_diagonal_threshold,
                test_type = 'Wilcox'
            )

            volcano_plot(
                genes_to_plot = genes_to_plot,
                gene_list_name = gene_list,
                results_scatter = results_scatter,
                group1 = group1,
                group2 = group2,
                cluster = cluster,
                my_colors = my_colors,
                local_figures_path = local_figures_path,
                FC_threshold = FC_threshold,
                p_value_threshold = p_value_threshold,
                max_overlaps = 15,
                label_size = label_size,
                label_threshold = label_threshold,
                test_type = 'Wilcox'
            )
        }
    }

    return(list(
        all_count = de_results$all_count,
        UP_count = de_results$UP_count,
        DOWN_count = de_results$DOWN_count,
        results = de_results$results,
        scatterplot = scatterplot_output,
        volcanoplot = volcanoplot_output
    ))
}


#' Complete Bulk RNA-seq DE Analysis Pipeline
#'
#' Orchestrates a full bulk RNA-seq differential expression analysis workflow using
#' DESeq2 including statistical testing, scatter and volcano plot generation, pathway
#' enrichment analysis, and optional custom gene list visualization.
#'
#' @param counts_table Data frame; count matrix with genes as rows and samples as columns,
#'   must include a 'genes' column
#' @param comparison Character; column name in metadata for comparison groups. Default 'Groups'
#' @param group1 Character; pattern to match for first group (control)
#' @param group2 Character; pattern to match for second group (treatment)
#' @param cluster Character; identifier for the dataset. Default ''
#' @param path Character; output directory path. Default './'
#' @param FC_threshold Numeric; log2 fold change threshold for filtering. Default 0.3
#' @param p_value_threshold Numeric; adjusted p-value threshold. Default 0.05
#' @param max_overlaps Integer; maximum label overlaps in plots. Default 15
#' @param label_size Numeric; size of gene labels in plots. Default 5
#' @param pathways_of_interest Character vector; specific pathways to analyze. Default NULL
#' @param label_threshold Numeric; expression threshold for labeling genes. Default 100000
#' @param distance_from_diagonal_threshold Numeric; minimum distance from diagonal
#'   for labeling in scatterplot. Default 0.4
#' @param gene_lists_to_plot Named list; custom gene lists to visualize. Default NULL
#' @param expression_threshold_for_gene_list Numeric; minimum average normalized counts
#'   for filtering gene lists. Default 20
#' @param minimum_cell_number Integer; minimum replicates required (not used for bulk). Default 10
#' @param run_pathway_enrichment Logical; whether to run pathway enrichment. Default FALSE
#' @param ... Additional arguments passed to core analysis and visualization functions
#'
#' @return Named vector with elements:
#'   \item{all_count}{Integer; number of significant DEGs}
#'   \item{UP_count}{Integer; number of upregulated genes}
#'   \item{DOWN_count}{Integer; number of downregulated genes}
#'
#' @export
bulk_analysis <- function(counts_table, comparison = 'Groups', group1, group2, cluster = '',
                         path = './', FC_threshold = 0.3, p_value_threshold = 0.05,
                         max_overlaps = 15, label_size = 5, pathways_of_interest = NULL,
                         label_threshold = 100000, distance_from_diagonal_threshold = 0.4,
                         gene_lists_to_plot = NULL, expression_threshold_for_gene_list = 20,
                         minimum_cell_number = 10, run_pathway_enrichment = FALSE, ...) {

    local_figures_path <- here::here(path, 'figures')
    dir.create(local_figures_path, showWarnings = FALSE, recursive = TRUE)

    group1 <- fixed(group1)
    group2 <- fixed(group2)

    # Run core DE analysis
    de_results <- bulk_analysis_de(
        counts_table = counts_table,
        comparison = comparison,
        group1 = group1,
        group2 = group2,
        cluster = cluster,
        path = path,
        FC_threshold = FC_threshold,
        p_value_threshold = p_value_threshold,
        expression_threshold_for_gene_list = expression_threshold_for_gene_list,
        ...
    )

    # Check if analysis was successful
    if (is.null(de_results$results)) {
        return(c(
            all_count = de_results$all_count,
            UP_count = de_results$UP_count,
            DOWN_count = de_results$DOWN_count
        ))
    }

    # Generate scatter and volcano plots
    scatterplot(results = de_results$results, group1 = group1, group2 = group2,
               cluster = cluster, local_figures_path = local_figures_path,
               FC_threshold = FC_threshold, p_value_threshold = p_value_threshold,
               max_overlaps = max_overlaps, label_size = label_size,
               label_threshold = label_threshold,
               distance_from_diagonal_threshold = distance_from_diagonal_threshold,
               test_type = 'Bulk', ...)

    volcano_plot(results = de_results$results, group1 = group1, group2 = group2,
                cluster = cluster, local_figures_path = local_figures_path,
                FC_threshold = FC_threshold, p_value_threshold = p_value_threshold,
                max_overlaps = max_overlaps, label_size = label_size,
                label_threshold = label_threshold, test_type = 'Bulk', ...)

    # Pathway enrichment analysis
    run_standard_enrichment(
        results = de_results$results,
        run_pathway_enrichment = run_pathway_enrichment,
        pathways_of_interest = pathways_of_interest,
        cluster = cluster,
        path = path,
        group1 = group1,
        group2 = group2,
        comparison = comparison,
        FC_threshold = FC_threshold
    )

    # Plot specific gene lists
    plot_gene_lists(
        results = de_results$results,
        gene_lists_to_plot = gene_lists_to_plot,
        group1 = group1,
        group2 = group2,
        cluster = cluster,
        my_colors = create_deg_colors(),
        local_figures_path = local_figures_path,
        FC_threshold = FC_threshold,
        p_value_threshold = p_value_threshold,
        test_type = 'Bulk',
        max_overlaps = max_overlaps,
        label_size = label_size,
        label_threshold = label_threshold,
        distance_from_diagonal_threshold = distance_from_diagonal_threshold,
        ...
    )

    return(c(
        all_count = de_results$all_count,
        UP_count = de_results$UP_count,
        DOWN_count = de_results$DOWN_count
    ))
}
