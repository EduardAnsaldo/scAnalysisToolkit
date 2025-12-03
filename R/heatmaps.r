## Heatmap Visualization Functions
## Functions for creating heatmaps of pathway gene expression

#' Create Pathway Gene Expression ComplexHeatmap::Heatmap (Genes on X-axis)
#'
#' Generates a heatmap showing z-scored gene expression across samples or groups.
#' Aggregates expression by grouping variable, scales values per gene, and displays
#' genes on the x-axis with samples/groups on the y-axis. Uses a diverging color
#' palette centered at zero.
#'
#' @param genes_to_plot Character vector; gene names to include in heatmap
#' @param seurat Seurat object with gene expression data
#' @param pathway_name Character; title for the heatmap
#' @param color_palette Character vector; diverging color palette for heatmap.
#'   Default diverging_palette_2
#' @param grouping_var Character; metadata column name for aggregating expression
#'   (e.g., 'Samples', 'cell_type'). Default 'Samples'
#'
#' @return ggplot object with heatmap
#'
#' @examples
#' genes <- c("Il2", "Ifng", "Tnf", "Il10")
#' plot_pathways_heatmap(genes, seurat, "Cytokine Signaling",
#'                      grouping_var = "Samples")
#' @export
plot_pathways_heatmap <- function(genes_to_plot, seurat, pathway_name, color_palette = diverging_palette_2, grouping_var = 'Samples') {

    Aggregated_expression <- AggregateExpression(seurat, group.by = grouping_var, return.seurat = T)
    data_to_plot <- Aggregated_expression[['RNA']]$data |>
        as.data.frame() |>
        rownames_to_column('gene') |>
        as_tibble() |>
        filter(gene %in% genes_to_plot) |>
        mutate(gene = factor(gene, levels = genes_to_plot)) |>
        pivot_longer(cols = -gene, names_to = 'Sample', values_to = 'expression') |>
        mutate(Sample = factor(Sample, levels = colnames(Aggregated_expression[['RNA']]))) |>
        group_by(gene) |>
        mutate(scaled_expression = scale(expression)[,1]) |>
        ungroup() |>
        mutate(Sample = fct_rev(factor(Sample)))

    # Calculate symmetric limits for the color scale
    max_abs <- max(abs(data_to_plot$scaled_expression), na.rm = TRUE)

    plot <- ggplot(data_to_plot, aes(x = gene, y = Sample, fill = scaled_expression)) +
        geom_tile(color = "grey60", linewidth = 0.3) +
        scale_fill_gradientn(
            colors = color_palette,
            name = "z-score<br>expression",
            limits = c(-max_abs, max_abs)
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
            legend.title = element_text(),
            plot.title = element_text(hjust = 0.5, size = 14)
        ) +
        labs(x = NULL, y = NULL, title = pathway_name)

    return(plot)
}


#' Create Pathway Gene Expression ComplexHeatmap::Heatmap (Genes on Y-axis)
#'
#' Generates a heatmap showing z-scored gene expression across samples or groups.
#' Aggregates expression by grouping variable, scales values per gene, and displays
#' genes on the y-axis with samples/groups on the x-axis. Uses a diverging color
#' palette centered at zero. Alternative orientation to plot_pathways_heatmap.
#'
#' @param genes_to_plot Character vector; gene names to include in heatmap
#' @param seurat Seurat object with gene expression data
#' @param pathway_name Character; title for the heatmap
#' @param color_palette Character vector; diverging color palette for heatmap.
#'   Default diverging_palette_2
#' @param grouping_var Character; metadata column name for aggregating expression
#'   (e.g., 'Samples', 'cell_type'). Default 'Samples'
#'
#' @return ggplot object with heatmap
#'
#' @examples
#' genes <- c("Il2", "Ifng", "Tnf", "Il10")
#' plot_pathways_heatmap2(genes, seurat, "Cytokine Signaling",
#'                       grouping_var = "cell_type")
#' @export
plot_pathways_heatmap2 <- function(genes_to_plot, seurat, pathway_name, color_palette = diverging_palette_2, grouping_var = 'Samples') {

    Aggregated_expression <- AggregateExpression(seurat, group.by = grouping_var, return.seurat = T)
    data_to_plot <- Aggregated_expression[['RNA']]$data |>
        as.data.frame() |>
        rownames_to_column('gene') |>
        as_tibble() |>
        filter(gene %in% genes_to_plot) |>
        mutate(gene = factor(gene, levels = genes_to_plot)) |>
        pivot_longer(cols = -gene, names_to = 'Sample', values_to = 'expression') |>
        mutate(Sample = factor(Sample, levels = colnames(Aggregated_expression[['RNA']]))) |>
        group_by(gene) |>
        mutate(scaled_expression = scale(expression)[,1]) |>
        ungroup() |>
        mutate(Sample = fct_rev(factor(Sample)))


    # Calculate symmetric limits for the color scale
    max_abs <- max(abs(data_to_plot$scaled_expression), na.rm = TRUE)

    plot <- ggplot(data_to_plot, aes(x = Sample, y = gene, fill = scaled_expression)) +
        geom_tile(color = "grey60", linewidth = 0.3) +
        scale_fill_gradientn(
            colors = color_palette,
            name = "z-scored\nexpression",
            limits = c(-max_abs, max_abs)
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5, size = 14)
        ) +
        labs(x = NULL, y = NULL, title = pathway_name)

    return(plot)
}
