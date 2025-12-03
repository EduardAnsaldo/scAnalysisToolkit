## Shared Utility Functions
## Functions used across multiple analysis scripts

#' Create or Recreate Analysis Directory
#'
#' Removes an existing directory and creates a fresh one. Useful for ensuring
#' clean output directories for analysis results.
#'
#' @param path Character string specifying the directory path to create
#' @param recursive Logical; should directories be created recursively? Default TRUE
#' @param showWarnings Logical; should warnings be shown? Default FALSE
#'
#' @return Character string of the created path
#'
#' @examples
#' create_analysis_directory("results/figures")
create_analysis_directory <- function(path, recursive = TRUE, showWarnings = FALSE) {
    unlink(path, recursive = TRUE)
    dir.create(path, recursive = recursive, showWarnings = showWarnings)
    return(path)
}

#' Create "No Data" Placeholder Plot
#'
#' Generates a blank plot with a text label, used when there are insufficient
#' results to visualize (e.g., no significant genes found).
#'
#' @param label Character string to display on the plot. Default 'N/A'
#' @param text_size Numeric; size of the text label. Default 24
#'
#' @return A ggplot object
#'
#' @examples
#' create_no_data_plot(label = "No significant results")
create_no_data_plot <- function(label = 'N/A', text_size = 24) {
    ggplot() +
        theme_void() +
        geom_text(aes(0, 0, label = label)) +
        theme(text = element_text(size = text_size)) +
        xlab(NULL)
}

#' Set Color Scale for Enrichment Plots
#'
#' Sets the global color scale option for enrichplot package visualizations
#' using a viridis color palette.
#'
#' @return NULL (invisibly); sets global option
#'
#' @examples
#' set_enrichment_color_scale()
set_enrichment_color_scale <- function() {
    color_scale <- viridis::viridis(n = 4, direction = -1)
    options(enrichplot.colours = color_scale)
}

#' Create Color Palette for Differential Expression Plots
#'
#' Generates a named color vector for DE visualizations with colors for
#' downregulated, upregulated, and non-significant genes.
#'
#' @param colors Character vector of length 2 with colors for DOWN and UP.
#'   Default c('green4', 'darkorchid4')
#'
#' @return Named character vector with colors for "DOWN", "UP", and "NO"
#'
#' @examples
#' my_colors <- create_deg_colors(c("blue", "red"))
create_deg_colors <- function(colors = c('green4', 'darkorchid4')) {
    my_colors <- c(colors, "gray")
    names(my_colors) <- c("DOWN", "UP", "NO")
    return(my_colors)
}

#' Add Gene Labels Layer for DE Plots
#'
#' Creates a geom_text_repel layer for adding gene labels to scatterplots
#' and volcano plots. Automatically adjusts nudge direction based on whether
#' genes are up or downregulated.
#'
#' @param label_column Character; name of the column containing gene labels
#' @param diffexpressed Character vector; differential expression status per gene
#' @param label_size Numeric; size of text labels. Default 5
#' @param max_overlaps Integer; maximum overlaps allowed for labels. Default 15
#' @param nudge_x Numeric; horizontal nudge for labels. Default NULL (auto-calculated)
#' @param nudge_y Numeric; vertical nudge for labels. Default NULL (auto-calculated)
#'
#' @return A geom_text_repel layer
#'
#' @examples
#' add_gene_labels_layer("genes_to_label", de_status, label_size = 4)
add_gene_labels_layer <- function(label_column, diffexpressed, label_size = 5,
                                   max_overlaps = 15, nudge_x = NULL, nudge_y = NULL) {
    # Default nudge based on direction if not specified
    if (is.null(nudge_x)) {
        nudge_x <- ifelse(diffexpressed == 'UP', -0.4, 1.25)
    }
    if (is.null(nudge_y)) {
        nudge_y <- ifelse(diffexpressed == 'UP', 1.25, -0.4)
    }

    geom_text_repel(
        size = label_size,
        box.padding = 0.35,
        show.legend = FALSE,
        max.overlaps = max_overlaps,
        max.time = 10,
        max.iter = 10000000,
        nudge_x = nudge_x,
        nudge_y = nudge_y,
        fontface = 'italic',
        aes(label = !!sym(label_column), segment.size = 0.3, segment.alpha = 0.4, segment.curvature = 0)
    )
}

#' Validate Cell Counts for Differential Expression
#'
#' Checks if each comparison group has sufficient cells for differential
#' expression analysis. Returns validation status and cell counts.
#'
#' @param seurat Seurat object
#' @param comparison Character; metadata column name for comparison
#' @param group1 Character; pattern to match for group 1
#' @param group2 Character; pattern to match for group 2
#' @param minimum_cell_number Integer; minimum cells required per group
#'
#' @return List with elements:
#'   \item{valid}{Logical; TRUE if both groups have enough cells}
#'   \item{n_group1}{Integer; number of cells in group 1}
#'   \item{n_group2}{Integer; number of cells in group 2}
#'   \item{message}{Character; error message if validation fails}
#'
#' @examples
#' validation <- validate_cell_counts(seurat, "condition", "WT", "KO", 30)
validate_cell_counts <- function(seurat, comparison, group1, group2, minimum_cell_number) {
    n_group1 <- seurat@meta.data |>
        filter(str_detect(!!as.name(comparison), group1)) |>
        nrow()
    n_group2 <- seurat@meta.data |>
        filter(str_detect(!!as.name(comparison), group2)) |>
        nrow()

    if (n_group1 < minimum_cell_number | n_group2 < minimum_cell_number) {
        return(list(
            valid = FALSE,
            n_group1 = n_group1,
            n_group2 = n_group2,
            message = 'Not enough cells'
        ))
    }

    return(list(
        valid = TRUE,
        n_group1 = n_group1,
        n_group2 = n_group2
    ))
}

#' Run Standard Pathway Enrichment Analysis
#'
#' Executes Metascape functional analysis and optional pathway-specific
#' analysis on differential expression results.
#'
#' @param results Data frame; differential expression results
#' @param run_pathway_enrichment Logical; whether to run enrichment
#' @param pathways_of_interest Character vector; specific pathways to analyze
#' @param cluster Character; cluster identifier
#' @param path Character; output directory path
#' @param group1 Character; first comparison group
#' @param group2 Character; second comparison group
#' @param comparison Character; comparison metadata column name
#' @param FC_threshold Numeric; fold change threshold
#'
#' @return NULL (invisibly); creates output files
#'
#' @examples
#' run_standard_enrichment(de_results, TRUE, NULL, "cluster_1",
#'                        "./results", "WT", "KO", "condition", 0.5)
run_standard_enrichment <- function(results, run_pathway_enrichment,
                                    pathways_of_interest = NULL, cluster, path,
                                    group1, group2, comparison = NULL, FC_threshold) {
    if (!run_pathway_enrichment) return(invisible(NULL))

    Metascape_functional_analysis(results, grouping_var = cluster,
                                  group1 = group1, group2 = group2,
                                  path = path, FC_threshold = FC_threshold)

    if (!is.null(pathways_of_interest) && !is.null(comparison)) {
        pathways_of_interest_analysis(results = results,
                                     pathways_of_interest = pathways_of_interest,
                                     cluster = cluster, path = path,
                                     group1 = group1, group2 = group2,
                                     comparison = comparison)
    }
}

#' Plot Scatter and Volcano Plots for Gene Lists
#'
#' Iterates through named lists of genes and generates corresponding
#' scatterplots and volcano plots for each gene list.
#'
#' @param results Data frame; differential expression results
#' @param gene_lists_to_plot Named list; each element is a vector of gene names
#' @param group1 Character; first comparison group
#' @param group2 Character; second comparison group
#' @param cluster Character; cluster identifier
#' @param my_colors Named character vector; colors for plot
#' @param local_figures_path Character; path to save figures
#' @param FC_threshold Numeric; fold change threshold
#' @param p_value_threshold Numeric; p-value threshold
#' @param test_type Character; type of test performed
#' @param max_overlaps Integer; max label overlaps. Default 15
#' @param label_size Numeric; label text size. Default 5
#' @param label_threshold Numeric; expression threshold for labeling. Default 100000
#' @param distance_from_diagonal_threshold Numeric; distance threshold. Default 0.7
#' @param ... Additional arguments passed to plot functions
#'
#' @return NULL (invisibly); creates plot files
#'
#' @examples
#' gene_lists <- list(immune = c("Il2", "Ifng"), metabolic = c("Hk2", "Pkm"))
#' plot_gene_lists(results, gene_lists, "WT", "KO", "cluster_1",
#'                my_colors, "./figures", 0.5, 0.05, "Wilcox")
plot_gene_lists <- function(results, gene_lists_to_plot, group1, group2,
                           cluster, my_colors, local_figures_path,
                           FC_threshold, p_value_threshold, test_type,
                           max_overlaps = 15, label_size = 5,
                           label_threshold = 100000,
                           distance_from_diagonal_threshold = 0.7, ...) {
    if (is.null(gene_lists_to_plot)) return(invisible(NULL))

    for (gene_list in names(gene_lists_to_plot)) {
        genes_to_plot <- gene_lists_to_plot[[gene_list]]
        print(genes_to_plot)

        scatterplot(genes_to_plot = genes_to_plot, gene_list_name = gene_list,
                   results = results, group1 = group1, group2 = group2,
                   cluster = cluster, my_colors = my_colors,
                   local_figures_path = local_figures_path,
                   FC_threshold = FC_threshold, p_value_threshold = p_value_threshold,
                   max_overlaps = max_overlaps, label_size = label_size,
                   label_threshold = label_threshold,
                   distance_from_diagonal_threshold = distance_from_diagonal_threshold,
                   test_type = test_type, ...)

        volcano_plot(genes_to_plot = genes_to_plot, gene_list_name = gene_list,
                    results = results, group1 = group1, group2 = group2,
                    cluster = cluster, my_colors = my_colors,
                    local_figures_path = local_figures_path,
                    FC_threshold = FC_threshold, p_value_threshold = p_value_threshold,
                    max_overlaps = max_overlaps, label_size = label_size,
                    label_threshold = label_threshold,
                    test_type = test_type, ...)
    }
}
