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
#' @param shrink_LFC Logical; if TRUE, performs log fold change shrinkage for more accurate estimates, especially for low count genes. Default FALSE
#' @param shrink_type Character; shrinkage estimator passed to \code{DESeq2::lfcShrink}, one of
#'   'apeglm', 'ashr' or 'normal'. Default NULL, which selects 'apeglm' when no \code{contrast} is
#'   given and 'ashr' when one is, because \code{lfcShrink} rejects a contrast with 'apeglm'.
#'   Requires the corresponding package to be installed.
#' @param design_formula Formula; custom design formula for DESeq2 analysis. Default NULL (uses ~condition)
#' @param parameter_2 Character; name of the second parameter in the design formula, if applicable. Default NULL
#' @param group1_parameter_2 Character; value of the second parameter for the first group. Default NULL
#' @param group2_parameter_2 Character; value of the second parameter for the second group. Default NULL
#' @param contrast Character vector or list; specifies the contrast for DESeq2 results. Default NULL (uses the default contrast based on the default ~condition design) See DESeq2 documentation and ?DESeq2::results for details on specifying contrasts.
#' @param FC_threshold Numeric; log2 fold change threshold for filtering. Default 0.3  
#' @param p_value_threshold Numeric; adjusted p-value threshold. Default 0.05
#' @param max_overlaps Integer; maximum label overlaps in plots. Default 15
#' @param expression_threshold_for_gene_list Numeric; minimum average normalized counts
#'   for filtering gene lists. Default 20
#' @param minimum_cell_number Integer; minimum cells required per group. Default 10
#' @param run_pathway_enrichment Character; method for pathway enrichment ('Metascape',
#'   'ClusterProfiler', or NULL). Default NULL
#' @param genes_to_exclude Character vector of regular expressions; genes matching any
#'   pattern are excluded from analysis. Plain gene names work as literal patterns.
#'   Default c()
#' @param exclude_vdj Logical; if TRUE, also exclude VDJ (immunoglobulin and T-cell
#'   receptor) genes. Default FALSE
#' @param ... Additional arguments passed to visualization and enrichment functions. Note that the
#'   plot label arguments (\code{plot_subtitle}, \code{x_axis_label}, \code{up_label},
#'   \code{down_label}, \code{x_column}, \code{y_column}, \code{x_group_label},
#'   \code{y_group_label}, \code{scatter_applicable}, \code{file_suffix}) are set from the contrast
#'   and must not be passed here.
#'
#' @section Plot labels:
#' Plot subtitles, the volcano x-axis title and its direction arrows, and the scatterplot axes
#' are derived from \code{design_formula} and \code{contrast}, so they describe the comparison
#' DESeq2 actually performed rather than always reading as \code{group2} vs \code{group1}. Under
#' an interaction design the scatterplot is replaced by a placeholder, since an
#' expression-vs-expression plot cannot represent a difference of fold changes. Contrasts that
#' cannot be classified fall back to neutral, direction-agnostic wording and emit a warning rather
#' than guessing. Output filenames gain a contrast-specific suffix so several contrasts on one
#' cluster do not overwrite each other.
#'
#' @return List with elements:
#'   \item{all_count}{Integer; number of significant DEGs}
#'   \item{UP_count}{Integer; number of upregulated genes}
#'   \item{DOWN_count}{Integer; number of downregulated genes}
#'   \item{results}{Data frame; complete DE results}
#'   \item{dds}{The fitted DESeqDataSet}
#'   \item{contrast_labels}{List; the plot labels describing the comparison the contrast encodes}
#'   \item{scatterplot}{ggplot; scatterplot object}
#'   \item{volcanoplot}{ggplot; volcano plot object}
#'   \item{enrichment}{Enrichment analysis results (plot_up and plot_down), or NULL}
#'
#' @export
pseudobulk <- function(scRNAseq, comparison, group1, group2, cluster = 'all_clusters',
                       path = './', shrink_LFC = FALSE, design_formula = NULL, parameter_2 = NULL, group1_parameter_2 = NULL, group2_parameter_2 = NULL,
                       contrast = NULL, shrink_type = NULL, FC_threshold = 0.3, p_value_threshold = 0.05,
                       max_overlaps = 15,
                       expression_threshold_for_gene_list = 20, minimum_cell_number = 10,
                       run_pathway_enrichment = NULL, genes_to_exclude = c(),
                       exclude_vdj = FALSE, ...) {

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
        shrink_LFC = shrink_LFC,
        design_formula = design_formula,
        parameter_2 = parameter_2,
        group1_parameter_2 = group1_parameter_2,
        group2_parameter_2 = group2_parameter_2,
        contrast = contrast,
        shrink_type = shrink_type,
        FC_threshold = FC_threshold,
        p_value_threshold = p_value_threshold,
        expression_threshold_for_gene_list = expression_threshold_for_gene_list,
        minimum_cell_number = minimum_cell_number,
        genes_to_exclude = genes_to_exclude,
        exclude_vdj = exclude_vdj,
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

    # Labels describing the comparison the DESeq2 contrast actually encodes
    contrast_labels <- de_results$contrast_labels

    # Generate scatter and volcano plots. The label arguments are derived from the contrast, not
    # caller-supplied: x_column / y_column / scatter_applicable in particular select which Avg_
    # columns exist to be plotted, so overriding them would name columns absent from results.
    scatterplot_output <- scatterplot(
        results = de_results$results,
        group1 = group1,
        group2 = group2,
        local_figures_path = local_figures_path,
        FC_threshold = FC_threshold,
        p_value_threshold = p_value_threshold,
        cluster = cluster,
        max_overlaps = max_overlaps,
        test_type = 'Pseudobulk',
        x_column = contrast_labels$x_column,
        y_column = contrast_labels$y_column,
        x_group_label = contrast_labels$x_group_label,
        y_group_label = contrast_labels$y_group_label,
        plot_subtitle = contrast_labels$subtitle,
        scatter_applicable = contrast_labels$scatter_applicable,
        file_suffix = contrast_labels$file_tag,
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
        max_overlaps = max_overlaps,
        test_type = 'Pseudobulk',
        plot_subtitle = contrast_labels$subtitle,
        x_axis_label = contrast_labels$volcano_x,
        up_label = contrast_labels$up_label,
        down_label = contrast_labels$down_label,
        file_suffix = contrast_labels$file_tag,
        ...
    )

    # Overrepresentation analysis. Use the contrast's own groups so enrichment output is
    # labelled with the same comparison as the plots.
    enrichment_output <- run_DEG_functional_analysis(
        results = de_results$results,
        method = run_pathway_enrichment,
        grouping_var = cluster,
        path = path,
        FC_threshold = FC_threshold,
        p_value_threshold = p_value_threshold,
        group1 = if (is.na(contrast_labels$denominator)) group1 else contrast_labels$denominator,
        group2 = if (is.na(contrast_labels$numerator)) group2 else contrast_labels$numerator,
        # The enrichment output directories are named from one group plus the cluster only, so two
        # contrasts sharing those would collide -- and create_analysis_directory() unlinks first,
        # deleting the earlier result rather than merging.
        file_suffix = contrast_labels$file_tag,
        ...
    )

    return(list(
        all_count = de_results$all_count,
        UP_count = de_results$UP_count,
        DOWN_count = de_results$DOWN_count,
        results = de_results$results,
        dds = de_results$dds,
        contrast_labels = contrast_labels,
        scatterplot = scatterplot_output,
        volcanoplot = volcanoplot_output,
        enrichment = enrichment_output
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
#' @param minimum_cell_number Integer; minimum cells required per group. Default 30
#' @param run_pathway_enrichment Character; method for pathway enrichment or NULL. Default NULL
#' @param exclude_sex_chr_genes Logical; if TRUE, removes mouse sex chromosome genes
#'   from results to avoid false positives from sex imbalance. Default FALSE
#' @param genes_to_exclude Character vector of regular expressions; genes matching any
#'   pattern are excluded from results. Plain gene names work as literal patterns.
#'   Default c()
#' @param exclude_vdj Logical; if TRUE, also exclude VDJ (immunoglobulin and T-cell
#'   receptor) genes. Default FALSE
#' @param ... Additional arguments passed to the underlying DE function
#'   (\code{DEG_FindMarkers_RNA_assay_de}) and the scatter/volcano plotting functions.
#'
#' @return List with elements:
#'   \item{all_count}{Integer; number of significant DEGs}
#'   \item{UP_count}{Integer; number of upregulated genes}
#'   \item{DOWN_count}{Integer; number of downregulated genes}
#'   \item{results}{Data frame; complete DE results}
#'   \item{enrichment}{Enrichment analysis results (plot_up and plot_down), or NULL}
#'
#' @export
DEG_FindMarkers_RNA_assay <- function(scRNAseq, comparison, group1, group2, cluster = 'all_clusters',
                                      path = './', FC_threshold = 0.3, p_value_threshold = 0.05,
                                      min_fraction = 0.01,
                                      max_overlaps = 15, label_size = 5, pathways_of_interest = NULL,
                                      label_threshold = 100000, distance_from_diagonal_threshold = 0.7,
                                      gene_lists_to_plot = NULL, expression_threshold_for_gene_list = 20,
                                      minimum_cell_number = 30,
                                      run_pathway_enrichment = NULL,
                                      exclude_sex_chr_genes = FALSE,
                                      genes_to_exclude = c(), exclude_vdj = FALSE, ...) {

    local_figures_path <- here::here(path, 'figures')
    dir.create(local_figures_path, showWarnings = FALSE, recursive = TRUE)

    group1 <- fixed(group1)
    group2 <- fixed(group2)

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
        min_fraction = min_fraction,
        minimum_cell_number = minimum_cell_number,
        exclude_sex_chr_genes = exclude_sex_chr_genes,
        genes_to_exclude = genes_to_exclude,
        exclude_vdj = exclude_vdj,
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
               cluster = cluster,
               local_figures_path = local_figures_path,
               FC_threshold = FC_threshold, p_value_threshold = p_value_threshold,
               max_overlaps = max_overlaps, label_size = label_size,
               label_threshold = label_threshold,
               distance_from_diagonal_threshold = distance_from_diagonal_threshold,
               test_type = 'Wilcox', ...)

    volcano_plot(results = de_results$results, group1 = group1, group2 = group2,
                cluster = cluster,
                local_figures_path = local_figures_path,
                FC_threshold = FC_threshold, p_value_threshold = p_value_threshold,
                max_overlaps = max_overlaps, label_size = label_size,
                label_threshold = label_threshold, test_type = 'Wilcox', ...)

    # Overrepresentation analysis
    enrichment_output <- run_DEG_functional_analysis(
        results = de_results$results,
        method = run_pathway_enrichment,
        grouping_var = cluster,
        path = path,
        FC_threshold = FC_threshold,
        p_value_threshold = p_value_threshold,
        group1 = group1,
        group2 = group2
    )

    # Plot specific gene lists
    plot_gene_lists(
        results = de_results$results,
        gene_lists_to_plot = gene_lists_to_plot,
        group1 = group1,
        group2 = group2,
        cluster = cluster,
        local_figures_path = local_figures_path,
        FC_threshold = FC_threshold,
        p_value_threshold = p_value_threshold,
        test_type = 'Wilcox',
        max_overlaps = max_overlaps,
        label_size = label_size,
        label_threshold = label_threshold,
        distance_from_diagonal_threshold = distance_from_diagonal_threshold
    )

    return(list(
        all_count = de_results$all_count,
        UP_count = de_results$UP_count,
        DOWN_count = de_results$DOWN_count,
        results = de_results$results,
        enrichment = enrichment_output
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
#' @param minimum_cell_number Integer; minimum cells required per group. Default 30
#' @param run_pathway_enrichment Character; method for pathway enrichment or NULL. Default NULL
#' @param exclude_sex_chr_genes Logical; if TRUE, removes mouse sex chromosome genes
#'   from results to avoid false positives from sex imbalance. Default FALSE
#' @param genes_to_exclude Character vector of regular expressions; genes matching any
#'   pattern are excluded from results. Plain gene names work as literal patterns.
#'   Default c()
#' @param exclude_vdj Logical; if TRUE, also exclude VDJ (immunoglobulin and T-cell
#'   receptor) genes. Default FALSE
#' @param ... Additional arguments passed to core analysis and visualization functions
#'
#' @return List with elements:
#'   \item{all_count}{Integer; number of significant DEGs}
#'   \item{UP_count}{Integer; number of upregulated genes}
#'   \item{DOWN_count}{Integer; number of downregulated genes}
#'   \item{results}{Data frame; complete DE results}
#'   \item{scatterplot}{ggplot; scatterplot object}
#'   \item{volcanoplot}{ggplot; volcano plot object}
#'   \item{enrichment}{Enrichment analysis results (plot_up and plot_down), or NULL}
#'
#' @export
DEG_FindMarkers_SCT_assay <- function(scRNAseq, comparison, group1, group2, is_integrated_subset = FALSE,
                                      cluster = 'all_clusters', path = './', FC_threshold = 0.3,
                                      p_value_threshold = 0.05, min_fraction = 0.01, max_overlaps = 15, label_size = 5,
                                      pathways_of_interest = NULL, label_threshold = 100000,
                                      distance_from_diagonal_threshold = 0.7, gene_lists_to_plot = NULL,
                                      expression_threshold_for_gene_list = 20,
                                      minimum_cell_number = 30, run_pathway_enrichment = NULL,
                                      exclude_sex_chr_genes = FALSE,
                                      genes_to_exclude = c(), exclude_vdj = FALSE, ...) {

    local_figures_path <- here::here(path, 'figures')
    dir.create(local_figures_path, showWarnings = FALSE, recursive = TRUE)

    group1 <- fixed(group1)
    group2 <- fixed(group2)

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
        exclude_sex_chr_genes = exclude_sex_chr_genes,
        genes_to_exclude = genes_to_exclude,
        exclude_vdj = exclude_vdj,
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

    # Generate scatter and volcano plots. Guarded with tryCatch so a plotting
    # failure cannot prevent the UP/DOWN counts from being returned.
    scatterplot_output <- tryCatch(
        scatterplot(
            results = de_results$results,
            group1 = group1,
            group2 = group2,
            cluster = cluster,
            local_figures_path = local_figures_path,
            FC_threshold = FC_threshold,
            p_value_threshold = p_value_threshold,
            max_overlaps = max_overlaps,
            label_size = label_size,
            label_threshold = label_threshold,
            distance_from_diagonal_threshold = distance_from_diagonal_threshold,
            test_type = 'Wilcox',
            ...
        ),
        error = function(e) { warning('scatterplot failed for ', cluster, ': ', conditionMessage(e)); NULL }
    )

    volcanoplot_output <- tryCatch(
        volcano_plot(
            results = de_results$results,
            group1 = group1,
            group2 = group2,
            cluster = cluster,
            local_figures_path = local_figures_path,
            FC_threshold = FC_threshold,
            p_value_threshold = p_value_threshold,
            max_overlaps = max_overlaps,
            label_size = label_size,
            label_threshold = label_threshold,
            test_type = 'Wilcox',
            ...
        ),
        error = function(e) { warning('volcano_plot failed for ', cluster, ': ', conditionMessage(e)); NULL }
    )

    # Overrepresentation analysis
    enrichment_output <- run_DEG_functional_analysis(
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

    # Plot specific gene lists
    plot_gene_lists(
        results = de_results$results,
        gene_lists_to_plot = gene_lists_to_plot,
        group1 = group1,
        group2 = group2,
        cluster = cluster,
        local_figures_path = local_figures_path,
        FC_threshold = FC_threshold,
        p_value_threshold = p_value_threshold,
        test_type = 'Wilcox',
        max_overlaps = max_overlaps,
        label_size = label_size,
        label_threshold = label_threshold,
        distance_from_diagonal_threshold = distance_from_diagonal_threshold,
        ...
    )

    return(list(
        all_count = de_results$all_count,
        UP_count = de_results$UP_count,
        DOWN_count = de_results$DOWN_count,
        results = de_results$results,
        scatterplot = scatterplot_output,
        volcanoplot = volcanoplot_output,
        enrichment = enrichment_output
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
#' @param genes_to_exclude Character vector of regular expressions; genes matching any
#'   pattern are excluded from analysis. Plain gene names work as literal patterns.
#'   Default c()
#' @param exclude_vdj Logical; if TRUE, also exclude VDJ (immunoglobulin and T-cell
#'   receptor) genes. Default FALSE
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
                         minimum_cell_number = 10, run_pathway_enrichment = FALSE,
                         genes_to_exclude = c(), exclude_vdj = FALSE, ...) {

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
        genes_to_exclude = genes_to_exclude,
        exclude_vdj = exclude_vdj,
        ...
    )

    # Check if analysis was successful
    if (is.null(de_results$results)) {
        return(c(
            all_count = de_results$all_count,
            UP_count = de_results$UP_count,
            DOWN_count = de_results$DOWN_count,
            results = de_results$results
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

    # Overrepresentation analysis
    enrichment_output <- run_DEG_functional_analysis(
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

    # Plot specific gene lists
    plot_gene_lists(
        results = de_results$results,
        gene_lists_to_plot = gene_lists_to_plot,
        group1 = group1,
        group2 = group2,
        cluster = cluster,
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

    return(list(
        all_count = de_results$all_count,
        UP_count = de_results$UP_count,
        DOWN_count = de_results$DOWN_count,
        results = de_results$results,
        enrichment = enrichment_output
    ))
}
