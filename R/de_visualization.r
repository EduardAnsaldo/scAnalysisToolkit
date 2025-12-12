## Differential Expression Visualization Functions
## Scatterplot and volcano plot for DE results

#' Create Scatterplot for Differential Expression Results
#'
#' Generates a scatterplot comparing average expression between two groups,
#' with genes colored by differential expression status. Labels are added
#' intelligently based on distance from diagonal and expression levels.
#'
#' @param results Data frame containing DE results with columns: genes, log2FoldChange,
#'   padj, pvalue, Avg_group1, Avg_group2
#' @param group1 Character; name of first comparison group
#' @param group2 Character; name of second comparison group
#' @param local_figures_path Character; directory path to save the plot
#' @param FC_threshold Numeric; log2 fold change threshold for significance
#' @param p_value_threshold Numeric; adjusted p-value threshold for significance
#' @param cluster Character; cluster or cell type identifier. Default 'all_clusters'
#' @param my_colors Character vector of length 3; colors for DOWN, UP, NO.
#'   Default c('green4', 'darkorchid4', 'gray')
#' @param max_overlaps Integer; maximum label overlaps allowed. Default 15
#' @param label_size Numeric; size of gene labels. Default 5
#' @param label_threshold Numeric; expression threshold for secondary labels. Default 10000
#' @param distance_from_diagonal_threshold Numeric; minimum distance from diagonal
#'   for labeling. Default 0.5
#' @param test_type Character; type of test ('Wilcox', 'Pseudobulk', or 'Bulk')
#' @param genes_to_plot Character vector; specific genes to highlight. Default NULL
#' @param pt_size Numeric; point size. Default 1.3
#' @param ... Additional arguments passed to plotting functions
#'
#' @return A ggplot object
#'
#' @export
scatterplot <- function (results, group1, group2, local_figures_path, FC_threshold, p_value_threshold, cluster = 'all_clusters', my_colors = c('green4', 'darkorchid4', 'gray'), max_overlaps = 15, label_size = 5, label_threshold = 10000, distance_from_diagonal_threshold = 0.5, test_type = c('Wilcox', 'Pseudobulk', 'Bulk'), genes_to_plot = NULL, pt_size = 1.3, ...) {

    # Set colors for the plot
    names(my_colors) <- c("DOWN", "UP", "NO")

    #Determine test type
    if (test_type == 'Pseudobulk') {
        axis_test <- 'Average Normalized Counts'
    } else if (test_type == 'Bulk') {
        axis_test <- 'Average Normalized Counts'
    } else if (test_type == 'Wilcox') {
        axis_test <- 'Average CPMs'
    }
    # If genes_to_plot is provided, plot only those genes as described in the first paragraph

    ## prepare for visualization
    results_scatter <- results |>
        drop_na(pvalue) |>
        mutate(
            log10_pval = log10(padj+10^-300)*-1,
            distance_from_diagonal =  (abs((log10(!!sym(paste0('Avg_', group2))+1)) - (log10(!!sym(paste0('Avg_', group1))+1)))/sqrt(2))) |>
        mutate(
            genes_to_label_first = ifelse(
                (log2FoldChange >= FC_threshold | log2FoldChange <= -1 * FC_threshold) &
                (padj < p_value_threshold) &
                (distance_from_diagonal > distance_from_diagonal_threshold) &
                ((!!sym(paste0('Avg_', group2)) > 100) | (!!sym(paste0('Avg_', group1)) > 100)),
                genes, NA
            ),
            genes_to_label_second = ifelse(
                (log2FoldChange >= FC_threshold | log2FoldChange <= -1 * FC_threshold) &
                (padj < p_value_threshold) &
                ((!!sym(paste0('Avg_', group2)) > label_threshold) | (!!sym(paste0('Avg_', group1)) > label_threshold)) &
                is.na(genes_to_label_first),
                genes, NA
            ),
            genes_to_label = ifelse(
                (log2FoldChange >= FC_threshold | log2FoldChange <= -1 * FC_threshold) &
                (padj < p_value_threshold) &
                is.na(genes_to_label_first) &
                is.na(genes_to_label_second),
                genes, NA
            ),
            diffexpressed = case_when(
                log2FoldChange >= FC_threshold & padj < p_value_threshold ~ 'UP',
                log2FoldChange <= -1 * FC_threshold & padj < p_value_threshold ~ "DOWN",
                TRUE ~ 'NO'
            )
        ) |>
        mutate(
            diffexpressed = factor(diffexpressed, levels = c('NO', 'DOWN', 'UP'))
        ) |>
        arrange(diffexpressed)

      # Replace genes to label with provided gene list if applicable
      if (!is.null(genes_to_plot)) {
        results_scatter <- results_scatter |>
            mutate(
                genes_to_label_first = ifelse(
                    genes %in% genes_to_plot,
                    genes, NA
                ),
                genes_to_label_second = NA,
                genes_to_label = NA)
              }

    # Scatterplot
    limx <- results_scatter |> pull(paste0('Avg_', group1)) |> max()
    limy <- results_scatter |> pull(paste0('Avg_', group2)) |> max()
    mylims <- max(limx, limy)*6

    scatter_plot <- results_scatter |>
        ggplot(aes(x = !!sym(paste0('Avg_', group1)), y = !!sym(paste0('Avg_', group2)), col = diffexpressed))+
            geom_point(size=pt_size, stroke = 0) +
            geom_abline(slope = 1, intercept = 0)+
            add_gene_labels_layer('genes_to_label_first', results_scatter$diffexpressed,
                                 label_size, max_overlaps) +
            add_gene_labels_layer('genes_to_label_second', results_scatter$diffexpressed,
                                 label_size, max_overlaps = 10) +
            add_gene_labels_layer('genes_to_label', results_scatter$diffexpressed,
                                 label_size, max_overlaps = 10) +
        scale_colour_manual(values=my_colors)+
        theme(text=element_text(size=20), legend.position="none")+
        labs(title=paste0(test_type, ' DEGs in ', str_replace(cluster,pattern = '_',replace = ' ') ),
                    x=paste0(axis_test, ' in ', group1),
                    y=paste0(axis_test, ' in ',  group2))+
        theme_classic(base_size = 28, base_line_size=1) +
        theme(legend.position="none",
            title = element_text(size=15),
            axis.text= element_text(size=10),
            axis.title= element_text(size=13))+
       scale_x_log10(limits =  c(0.5, mylims), expand = expansion(mult = c(0.01, 0.1)))+
       scale_y_log10(limits =  c(0.5, mylims), expand = expansion(mult = c(0.01, 0.1)))

    ggsave(plot = scatter_plot, filename = paste0(test_type,'_scatter_DEG_in_', cluster, '.pdf'), path = local_figures_path)
    print(scatter_plot)
    return(scatter_plot)
}

#' Create Volcano Plot for Differential Expression Results
#'
#' Generates a volcano plot displaying log2 fold change vs. -log10 adjusted p-value,
#' with genes colored by differential expression status. Supports both RNA-seq and
#' ATAC-seq data. Automatically removes extreme non-significant genes to improve
#' visualization.
#'
#' @param results Data frame containing DE results with columns: genes, log2FoldChange,
#'   padj, pvalue
#' @param group1 Character; name of first comparison group
#' @param group2 Character; name of second comparison group
#' @param cluster Character; cluster or cell type identifier
#' @param local_figures_path Character; directory path to save the plot
#' @param FC_threshold Numeric; log2 fold change threshold for significance
#' @param p_value_threshold Numeric; adjusted p-value threshold for significance
#' @param max_overlaps Integer; maximum label overlaps allowed. Default 15
#' @param label_size Numeric; size of gene labels. Default 5
#' @param my_colors Character vector of length 3; colors for DOWN, UP, NO.
#'   Default c('green4', 'darkorchid4', 'gray')
#' @param test_type Character; type of test. Options: 'Wilcox', 'Pseudobulk', 'Bulk',
#'   'Wilcox_ATAC', 'Wilcox_ATAC_closest_genes'
#' @param genes_to_plot Character vector; specific genes to highlight. Default NULL
#' @param pt_size Numeric; point size. Default 1.5
#' @param ... Additional arguments passed to plotting functions
#'
#' @return A ggplot object
#'
#' @export
volcano_plot <- function (results, group1, group2, cluster, local_figures_path, FC_threshold, p_value_threshold, max_overlaps = 15, label_size = 5, my_colors = c('green4', 'darkorchid4', 'gray'), test_type = c('Wilcox', 'Pseudobulk', 'Bulk', 'Wilcox_ATAC', 'Wilcox_ATAC_closest_genes', 'Wilcox_ChromVar_motif'), genes_to_plot = NULL, pt_size = 1.5, ...) {

    # determine plot title
    test_type <- match.arg(test_type)

    if (test_type == 'Wilcox_ATAC') {
        plot_title <- paste0('Differentially Accessible Regions in ', str_replace(cluster,pattern = '_',replace = ' ') )
    } else if (test_type == 'Wilcox_ATAC_closest_genes') {
        plot_title <- paste0('Genes Closest to Differentially Accessible Regions in ', str_replace(cluster,pattern = '_',replace = ' ') )
    } else if (test_type == 'Wilcox_ChromVar_motif') {
        plot_title <- paste0('Differentially Active Motifs in ', str_replace(cluster,pattern = '_',replace = ' ') )
    } else {
        plot_title <- paste0(test_type, ' DEGs in ', str_replace(cluster,pattern = '_',replace = ' ') )
    }

    # Set colors for the plot
    names(my_colors) <- c("DOWN", "UP", "NO")

    nudge_x <- 0
    nudge_y <- 0

    ## prepare for visualization
    results_volcano <- results |>
        drop_na(pvalue) |>
        mutate(
            log10_pval = log10(padj+10^-600)*-1) |>
        mutate(
            genes_to_label_UP = ifelse(
                (log2FoldChange >= FC_threshold) &
                (padj < p_value_threshold),
                genes, NA
            ),
            genes_to_label_DOWN = ifelse(
                (log2FoldChange <= -1 * FC_threshold) &
                (padj < p_value_threshold),
                genes, NA
            ),
            diffexpressed = case_when(
                log2FoldChange >= FC_threshold & padj < p_value_threshold ~ 'UP',
                log2FoldChange <= -1*FC_threshold & padj < p_value_threshold ~ "DOWN",
                TRUE ~ 'NO'
            )
        ) |>
        mutate(
            diffexpressed = factor(diffexpressed, levels = c('NO', 'DOWN', 'UP'))

        ) |>
        arrange(diffexpressed)

      # Replace genes to label with provided gene list if applicable
      if (!is.null(genes_to_plot)) {
        results_volcano <- results_volcano |>
            mutate(
                genes_to_label_UP = ifelse((
                    genes %in% genes_to_plot) & (log2FoldChange >= FC_threshold),
                    genes, NA
                ),
                genes_to_label_DOWN = ifelse((
                    genes %in% genes_to_plot) & (log2FoldChange <= -1 * FC_threshold),
                    genes, NA
                )
            )
            nudge_x <- 3
            nudge_y <- 3
      }

    # Remove non-significant genes that would bias the plot visualization
    initial_number_of_genes <- nrow(results_volcano)
    max_FC_up_significant <- results_volcano |> filter(diffexpressed != 'NO') |> dplyr::select(log2FoldChange) |> max(na.rm = T)
    min_FC_up_significant <- results_volcano |> filter(diffexpressed != 'NO') |> dplyr::select(log2FoldChange) |> min(na.rm = T)
    if (min_FC_up_significant > -3 | is.na(min_FC_up_significant)) {
        min_FC_up_significant <- -3
    }
    if (max_FC_up_significant  < 3 | is.na(max_FC_up_significant)) {
        max_FC_up_significant <- 3
    }
    results_volcano <- results_volcano |> filter(!(diffexpressed == 'NO' & (log2FoldChange < min_FC_up_significant | log2FoldChange > max_FC_up_significant))) |>
        arrange(diffexpressed)
    final_number_of_genes <- nrow(results_volcano)
    print(paste('Removed', initial_number_of_genes-final_number_of_genes, 'non-significant genes that would bias the plot visualization'))

    volcano_plot <- results_volcano |>
        arrange(desc(padj)) |>
        ggplot(aes(x=log2FoldChange, y=log10_pval,  col=diffexpressed)) +
        geom_point(size=pt_size, stroke = 0) +
        geom_text_repel(
            size=label_size,
            box.padding = 0.35,
            show.legend = FALSE,
            max.overlaps = max_overlaps,
            max.time = 10,
            max.iter = 10000000,
            nudge_x = nudge_x,
            nudge_y = nudge_y,
            aes(label = genes_to_label_UP,segment.size=0.5, segment.alpha=0.8, segment.curvature=0),
            fontface = 'italic') +
        geom_text_repel(
            size=label_size,
            box.padding = 0.35,
            show.legend = FALSE,
            max.overlaps = max_overlaps,
            max.time = 10,
            max.iter = 10000000,
            nudge_x = -1*nudge_x,
            nudge_y = nudge_y,
            fontface = 'italic',
            aes(label = genes_to_label_DOWN, segment.size=0.5, segment.alpha=0.8, segment.curvature=0)) +
        scale_colour_manual(values=my_colors)+
        geom_vline(xintercept= FC_threshold, col="lavenderblush2", linetype=2, size=0.5) +
        geom_vline(xintercept=-FC_threshold, col="lavenderblush2", linetype=2, size=0.5) +
        geom_hline(yintercept=-1*log10(p_value_threshold), col="lavenderblush2", linetype=2, size=0.5)+
        theme(text=element_text(size=20), legend.position="none")+
        labs(title=plot_title,
                    x=paste('Average log2 FC (', group2, '/', group1, ')', sep=''),
                    y= '-Log10 Adj. p-value')+
        theme_classic(base_size = 28, base_line_size=1) +
        theme(legend.position="none",
            title = element_text(size=15),
            axis.text= element_text(size=10),
            axis.title= element_text(size=13),
            )         +
            scale_y_continuous(n.breaks = 8, expand = expansion(mult = c(0.01, 0.1))) +
            scale_x_continuous(n.breaks = 8)
    ggsave(plot = volcano_plot, filename = paste0(test_type, '_volcano_in_', cluster, '.pdf'), path = local_figures_path)
    print(volcano_plot)
    return(volcano_plot)
}
