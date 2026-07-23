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
#' @param title Character; title for the heatmap
#' @param color_palette Character vector; diverging color palette for heatmap.
#'   Default grDevices::hcl.colors(n = 20,'RdBu',rev = T)
#' @param grouping_var Character; metadata column name for aggregating expression
#'   (e.g., 'Samples', 'cell_type'). Default 'Samples'
#'
#' @return ggplot object with heatmap
#'
#' @export
plot_pathways_heatmap <- function(genes_to_plot, seurat, title, color_palette = grDevices::hcl.colors(n = 20,'RdBu',rev = T), grouping_var = 'Samples') {

    Aggregated_expression <- AggregateExpression(seurat, assays = 'RNA', group.by = grouping_var, return.seurat = T)
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
        labs(x = NULL, y = NULL, title = title)

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
#' @param title Character; title for the heatmap
#' @param color_palette Character vector; diverging color palette for heatmap.
#'   Default grDevices::hcl.colors(n = 20,'RdBu',rev = T)
#' @param grouping_var Character; metadata column name for aggregating expression
#'   (e.g., 'Samples', 'cell_type'). Default 'Samples'
#'
#' @return ggplot object with heatmap
#'
#' @export
plot_pathways_heatmap2 <- function(genes_to_plot, seurat, title, color_palette = grDevices::hcl.colors(n = 20,'RdBu',rev = T), grouping_var = 'Samples') {

    Aggregated_expression <- AggregateExpression(seurat, assays = 'RNA', group.by = grouping_var, return.seurat = T)
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
        labs(x = NULL, y = NULL, title = title)

    return(plot)
}

#' Plot chromVAR motif activity heatmap by group
#'
#' Creates a heatmap visualizing the mean chromVAR motif activity (deviation
#' z-scores) across groups (e.g., samples, clusters, conditions). Unlike RNA
#' expression, chromVAR scores are already z-scores relative to a GC-matched
#' background, so they are averaged per group without further scaling. The
#' color scale is centered at 0 with symmetric limits, so white represents
#' baseline activity, red represents above-baseline, and blue represents
#' below-baseline.
#'
#' @param motifs_to_plot Character vector of motif IDs (matching rownames of
#'   the chromVAR assay) to include in the heatmap, in the order they should
#'   appear on the y-axis. If `motif_labels = TRUE`, this vector must be
#'   named, with names being the display labels (e.g., TF names) and values
#'   being the motif IDs:
#'   `c(Sox2 = "MA0143.4", "Pou5f1::Sox2" = "MA0142.1")`. Duplicated names
#'   are allowed (e.g., when the same TF binds multiple motifs).
#' @param seurat A Seurat object with a `chromvar` assay (typically created
#'   via `Signac::RunChromVAR()`).
#' @param title Character string used as the plot title.
#' @param color_palette Character vector of colors for the diverging fill
#'   gradient. Defaults to a reversed RdBu palette with 20 steps. Should be
#'   ordered low → high (negative → positive activity).
#' @param grouping_var Character string naming the metadata column in
#'   `seurat@meta.data` to group cells by (e.g., `"Samples"`,
#'   `"seurat_clusters"`). Each unique value becomes one column in the heatmap.
#' @param motif_labels Logical. If `TRUE` (default `FALSE`), uses the names of
#'   `motifs_to_plot` as y-axis labels instead of the motif IDs. Requires
#'   `motifs_to_plot` to be a named vector. Duplicated names are displayed
#'   as-is (the same label can appear multiple times on the y-axis).
#'
#' @return A `ggplot` object.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Pulls the `data` slot from the `chromvar` assay.
#'   \item Computes the mean activity per group per motif using `rowMeans()`.
#'   \item Reshapes to long format and orders motifs/samples as factors.
#'   \item Builds a tile plot with a diverging fill scale centered at 0, with
#'     symmetric limits set by the maximum absolute mean activity.
#' }
#'
#' Note: chromVAR scores can be negative, so this function uses `rowMeans()`
#' rather than `Seurat::AggregateExpression()`, which is designed for
#' log-normalized counts and would distort z-scored data.
#'
#' Internally, when `motif_labels = TRUE`, the y-axis uses the motif IDs as
#' unique factor levels (with display labels supplied via `scale_y_discrete`),
#' so duplicated TF names render correctly without being collapsed into a
#' single row. Motif IDs in `motifs_to_plot` must themselves be unique.
#'
#' @examples
#' \dontrun{
#' # With TF names as labels, including a duplicated TF name
#' plot_chromvar_heatmap(
#'   motifs_to_plot = c(Sox2 = "MA0143.4",
#'                      Sox2 = "MA0143.3",      # same TF, different motif
#'                      "Pou5f1::Sox2" = "MA0142.1",
#'                      Klf4 = "MA0039.4"),
#'   seurat = seurat_obj,
#'   title = "Pluripotency TFs",
#'   grouping_var = "seurat_clusters",
#'   motif_labels = TRUE
#' )
#' }
#'
#' @importFrom Seurat GetAssayData
#' @importFrom dplyr filter mutate bind_cols
#' @importFrom tibble as_tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom purrr map set_names
#' @importFrom forcats fct_rev
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn
#'   scale_y_discrete theme_minimal theme element_text labs
#' @export
plot_chromvar_heatmap <- function(motifs_to_plot, 
                                   seurat, 
                                   title, 
                                   color_palette = grDevices::hcl.colors(n = 20, 'RdBu', rev = TRUE), 
                                   grouping_var = 'Samples',
                                   motif_labels = FALSE) {
  
  # Validate inputs
  if (motif_labels && is.null(names(motifs_to_plot))) {
    stop("When `motif_labels = TRUE`, `motifs_to_plot` must be a named vector ",
         "(names will be used as display labels).")
  }
  if (any(duplicated(motifs_to_plot))) {
    stop("`motifs_to_plot` must contain unique motif IDs.")
  }
  
  # Pull chromVAR matrix and group vector
  chromvar_mat <- GetAssayData(seurat, assay = "chromvar", slot = "data")
  group_vec <- seurat[[grouping_var, drop = TRUE]]
  
  # Compute mean z-score per group per motif (do not re-scale: scores are
  # already z-scored relative to a GC-matched background)
  aggregated <- levels(group_vec) |> rev() |>
    set_names() |>
    map(~ rowMeans(chromvar_mat[, group_vec == .x, drop = FALSE], na.rm = TRUE)) |>
    bind_cols() |>
    as.data.frame()
  rownames(aggregated) <- rownames(chromvar_mat)
  
  # Reshape to long format. Keep motif IDs as the factor (they're unique).
  data_to_plot <- aggregated |>
    rownames_to_column('motif') |>
    as_tibble() |>
    filter(motif %in% motifs_to_plot) |>
    mutate(motif = factor(motif, levels = motifs_to_plot)) |>
    pivot_longer(cols = -motif, names_to = 'Sample', values_to = 'activity') |>
    mutate(Sample = factor(Sample, levels = colnames(aggregated))) |>
    mutate(Sample = fct_rev(Sample))
  
  # Symmetric color limits centered at 0 (baseline activity)
  max_abs <- max(abs(data_to_plot$activity), na.rm = TRUE)
  
  # Build heatmap
  plot <- ggplot(data_to_plot, aes(x = Sample, y = motif, fill = activity)) +
    geom_tile(color = "grey60", linewidth = 0.3) +
    scale_fill_gradientn(
      colors = color_palette,
      name = "Mean motif\nChromVAR\nscore (z)",
      limits = c(-max_abs, max_abs)
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, size = 14)
    ) +
    labs(x = NULL, y = NULL, title = title)
  
  # Apply display labels if requested. The factor levels remain the unique
  # motif IDs (so each motif gets its own row), but the displayed tick labels
  # are the names of motifs_to_plot — duplicates render as identical strings.
  if (motif_labels) {
    id_to_label <- set_names(names(motifs_to_plot), motifs_to_plot)
    plot <- plot + scale_y_discrete(labels = id_to_label)
  }
  
  return(plot)
}

#' Plot a heatmap of DORC accessibility across samples
#'
#' Creates a heatmap visualizing z-scored DORC (Domains of Regulatory Chromatin)
#' accessibility across samples or other grouping variables. Accessibility values
#' are aggregated within groups and then z-scored across groups for each DORC,
#' allowing relative patterns of accessibility to be compared.
#'
#' @param dorcs_to_plot Character vector of DORC names to include in the heatmap.
#'   The order of this vector determines the row order in the resulting plot.
#' @param seurat A Seurat object containing a `DORCs` assay. The assay should
#'   have a populated `data` slot (i.e. normalized accessibility values).
#' @param title Character string used as the plot title. Typically the
#'   name of the pathway or gene set the DORCs belong to.
#' @param color_palette Character vector of colors used for the fill gradient.
#'   Defaults to a reversed RdBu palette with 20 colors, suitable for
#'   diverging z-scored data.
#' @param grouping_var Character string giving the name of the metadata column
#'   in the Seurat object to aggregate accessibility by. Defaults to `'Samples'`.
#'
#' @return A `ggplot` object showing a heatmap with DORCs on the y-axis,
#'   samples (or other groups) on the x-axis, and tiles colored by z-scored
#'   accessibility. The color scale is symmetric around zero.
#'
#' @details
#' Accessibility values are aggregated within each group using
#' `Seurat::AggregateExpression()` with `return.seurat = TRUE`, then z-scored
#' across groups within each DORC using `scale()`. The color scale limits are
#' set symmetrically based on the maximum absolute z-score so that zero maps
#' to the midpoint of the palette.
#'
#' Note that depending on how the `DORCs` assay was constructed, you may need
#' to use the `counts` slot instead of `data` if no normalization has been
#' applied.
#'
#' @examples
#' \dontrun{
#' dorcs <- c("GENE1", "GENE2", "GENE3")
#' plot_dorcs_heatmap(
#'     dorcs_to_plot = dorcs,
#'     seurat = seurat_obj,
#'     title = "Inflammatory Response",
#'     grouping_var = "Samples"
#' )
#' }
#'
#' @importFrom Seurat AggregateExpression
#' @importFrom dplyr filter mutate group_by ungroup
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom tidyr pivot_longer
#' @importFrom forcats fct_rev
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn theme_minimal
#'   theme element_text labs
#'
#' @export
plot_dorcs_heatmap <- function(dorcs_to_plot, seurat, title, color_palette = grDevices::hcl.colors(n = 20, 'RdBu', rev = TRUE), grouping_var = 'Samples') {

    Aggregated_expression <- AggregateExpression(seurat, assays = 'DORCs', group.by = grouping_var, return.seurat = TRUE)
    data_to_plot <- Aggregated_expression[['DORCs']]$data |>
        as.data.frame() |>
        rownames_to_column('dorc') |>
        as_tibble() |>
        filter(dorc %in% dorcs_to_plot) |>
        mutate(dorc = factor(dorc, levels = dorcs_to_plot |> rev())) |>
        pivot_longer(cols = -dorc, names_to = 'Sample', values_to = 'accessibility') |>
        mutate(Sample = str_replace(Sample, 'g', '')) |>
        mutate(Sample = fct_inseq(Sample)) |>
        mutate(Sample = fct_rev(Sample)) |>
        group_by(dorc) |>
        mutate(scaled_accessibility = scale(accessibility)[,1]) |>
        ungroup() |>
        mutate(Sample = fct_rev(factor(Sample)))

    # Calculate symmetric limits for the color scale
    max_abs <- max(abs(data_to_plot$scaled_accessibility), na.rm = TRUE)

    plot <- ggplot(data_to_plot, aes(x = Sample, y = dorc, fill = scaled_accessibility)) +
        geom_tile(color = "grey60", linewidth = 0.3) +
        scale_fill_gradientn(
            colors = color_palette,
            name = "z-scored\naccessibility",
            limits = c(-max_abs, max_abs)
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5, size = 14)
        ) +
        labs(x = NULL, y = NULL, title = title)

    return(plot)
}


#' Plot a figR regulatory score heatmap (ggplot2 implementation)
#'
#' Visualizes TF-DORC regulatory scores as a heatmap, with DORCs on rows and
#' TFs on columns. Both axes are ordered by hierarchical clustering using
#' Pearson correlation distance by default, mirroring the default behavior of
#' the `ComplexHeatmap`-based implementation.
#'
#' Filtering behavior:
#' \itemize{
#'   \item When both `DORCs` and `TFs` are specified, `score.cut` is ignored
#'     and the heatmap shows all pairs from the user-supplied lists.
#'   \item When only `DORCs` is specified, all supplied DORCs are kept (no
#'     score filtering on rows), and TFs are filtered to those passing
#'     `score.cut` for at least one of the supplied DORCs.
#'   \item When only `TFs` is specified, all supplied TFs are kept (no score
#'     filtering on columns), and DORCs are filtered to those passing
#'     `score.cut` for at least one of the supplied TFs.
#'   \item When neither is specified, both axes are filtered by `score.cut`.
#' }
#'
#' Clustering can be disabled per-axis via `cluster_DORCs` and `cluster_TFs`.
#' When clustering is disabled and the corresponding axis was user-supplied,
#' the user's order is used but reversed (so the first element appears at the
#' top of the y-axis or the left of the x-axis, matching how the input is
#' typically read).
#'
#' DORCs may optionally be split into groups via `DORC_groups`, in which case
#' the heatmap is faceted by group along the y-axis. Grouping requires `DORCs`
#' to be specified and disables row clustering (clustering across facets is
#' not meaningful).
#'
#' Specific TFs can be visually emphasized via `highlight_TFs`, which both
#' colors/bolds the corresponding x-axis labels and (optionally) draws a
#' rectangle around the entire column for each highlighted TF. The label and
#' box colors can be set independently.
#'
#' @param figR.d A data.frame of figR results containing at least the columns
#'   `DORC`, `Motif`, and `Score`.
#' @param score.cut Numeric absolute score threshold. Applied only to the
#'   axis (or axes) that the user did not specify. Ignored when both `DORCs`
#'   and `TFs` are supplied.
#' @param DORCs Optional character vector of DORC names to restrict the plot to.
#' @param TFs Optional character vector of TF (Motif) names to restrict the
#'   plot to.
#' @param DORC_groups Optional vector the same length as `DORCs` giving a
#'   group label for each DORC. When supplied, the heatmap is faceted by group
#'   along rows. Requires `DORCs` to be non-NULL.
#' @param cluster_DORCs Logical; whether to cluster DORCs (rows). Only takes
#'   effect when `DORCs` is non-NULL. Forced to `FALSE` when `DORC_groups` is
#'   supplied. Defaults to `TRUE`.
#' @param cluster_TFs Logical; whether to cluster TFs (columns). Only takes
#'   effect when `TFs` is non-NULL. Defaults to `TRUE`.
#' @param color_palette Character vector of colors used for the fill gradient.
#'   Defaults to a reversed RdBu palette suitable for diverging scores.
#' @param limits Numeric vector of length 2 giving the fill scale limits.
#'   Defaults to `c(-2, 2)` to match the original implementation.
#' @param plot_title Optional plot title. Defaults to `NULL` (no title).
#' @param x_title Character string for the x-axis title. Defaults to `"TFs"`.
#' @param y_title Character string for the y-axis title. Defaults to `"DORCs"`.
#' @param highlight_TFs Optional character vector of TF (Motif) names to
#'   highlight. Highlighted TFs receive bold, colored x-axis labels and, if
#'   `highlight_box = TRUE`, a rectangle drawn around their column.
#' @param highlight_label_color Color used for highlighted TF axis labels.
#'   Defaults to `"red"`.
#' @param highlight_box_color Color used for the rectangles drawn around
#'   highlighted TF columns. Defaults to `"red"`.
#' @param highlight_box Logical; whether to draw a rectangle around the column
#'   of each highlighted TF. Defaults to `TRUE`. Has no effect when
#'   `highlight_TFs` is `NULL`.
#' @param highlight_box_linewidth Numeric line width for the highlight
#'   rectangles. Defaults to `1`.
#' @param group_tag Character string used as the plot tag when `DORC_groups`
#'   is supplied. Sits between the panel and the legend, rotated to read
#'   vertically. Defaults to `"DORC groups"`.
#'
#' @return A `ggplot` object.
#'
#' @importFrom dplyr filter pull mutate left_join group_by summarise n_distinct
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom tibble column_to_rownames rownames_to_column as_tibble tibble
#' @importFrom ggplot2 ggplot aes geom_tile geom_rect scale_fill_gradientn
#'   scale_x_discrete theme_minimal theme element_text element_blank
#'   element_rect labs facet_grid vars margin
#'
#' @export
plotfigRHeatmap_EA_gg <- function(figR.d,
                                  score.cut = 1,
                                  DORCs = NULL,
                                  TFs = NULL,
                                  DORC_groups = NULL,
                                  cluster_DORCs = TRUE,
                                  cluster_TFs = TRUE,
                                  color_palette = grDevices::hcl.colors(n = 20, 'RdBu', rev = TRUE),
                                  limits = c(-2, 2),
                                  plot_title = NULL,
                                  x_title = "TFs",
                                  y_title = "DORCs",
                                  highlight_TFs = NULL,
                                  highlight_label_color = "red",
                                  highlight_box_color = "red",
                                  highlight_box = TRUE,
                                  highlight_box_linewidth = 1,
                                  group_tag = "DORC groups") {

    # Validate user-supplied names up front
    if (!is.null(DORCs) && !all(DORCs %in% figR.d$DORC))
        stop("One or more DORCs specified is not a valid DORC symbol found in the data.frame")
    if (!is.null(TFs) && !all(TFs %in% figR.d$Motif))
        stop("One or more TFs specified is not a valid TF symbol found in the data.frame")
    if (!is.null(highlight_TFs) && !all(highlight_TFs %in% figR.d$Motif))
        stop("One or more highlight_TFs is not a valid TF symbol found in the data.frame")

    # Validate DORC_groups
    if (!is.null(DORC_groups)) {
        if (is.null(DORCs))
            stop("`DORC_groups` requires `DORCs` to be specified.")
        if (length(DORC_groups) != length(DORCs))
            stop("`DORC_groups` must be the same length as `DORCs`.")
        if (cluster_DORCs) {
            message("`DORC_groups` supplied; disabling DORC clustering.\n")
            cluster_DORCs <- FALSE
        }
    }

    if (!is.null(DORCs) && !is.null(TFs)) {
        message("Both DORCs and TFs specified; ignoring score.cut.\n")
        DORCsToKeep <- DORCs
        TFsToKeep   <- TFs
    } else if (!is.null(DORCs)) {
        message("DORCs specified; applying score.cut (", score.cut, ") to TFs only.\n")
        DORCsToKeep <- DORCs
        TFsToKeep   <- figR.d |>
            dplyr::filter(DORC %in% DORCsToKeep & abs(Score) >= score.cut) |>
            dplyr::pull(Motif) |> unique()
    } else if (!is.null(TFs)) {
        message("TFs specified; applying score.cut (", score.cut, ") to DORCs only.\n")
        TFsToKeep   <- TFs
        DORCsToKeep <- figR.d |>
            dplyr::filter(Motif %in% TFsToKeep & abs(Score) >= score.cut) |>
            dplyr::pull(DORC) |> unique()
    } else {
        message("Using absolute score cut-off of: ", score.cut, " ..\n")
        DORCsToKeep <- figR.d |> dplyr::filter(abs(Score) >= score.cut) |> dplyr::pull(DORC) |> unique()
        TFsToKeep   <- figR.d |> dplyr::filter(abs(Score) >= score.cut) |> dplyr::pull(Motif) |> unique()
    }

    # Build matrix (DORCs x TFs)
    net.d <- figR.d |>
        dplyr::filter(DORC %in% DORCsToKeep & Motif %in% TFsToKeep) |>
        tidyr::pivot_wider(id_cols = DORC, names_from = Motif, values_from = Score) |>
        tibble::column_to_rownames("DORC") |>
        as.matrix()

    net.d[is.na(net.d)] <- 0

    message("Plotting ", nrow(net.d), " DORCs x ", ncol(net.d), " TFs\n")

    do_cluster_rows <- if (is.null(DORCs)) TRUE else cluster_DORCs
    do_cluster_cols <- if (is.null(TFs))   TRUE else cluster_TFs

    pearson_dist <- function(x) as.dist(1 - cor(t(x), use = "pairwise.complete.obs"))

    row_order <- if (do_cluster_rows && nrow(net.d) > 1) {
        rownames(net.d)[hclust(pearson_dist(net.d))$order]
    } else if (!is.null(DORCs)) {
        rev(intersect(DORCs, rownames(net.d)))
    } else {
        rownames(net.d)
    }

    col_order <- if (do_cluster_cols && ncol(net.d) > 1) {
        colnames(net.d)[hclust(pearson_dist(t(net.d)))$order]
    } else if (!is.null(TFs)) {
        rev(intersect(TFs, colnames(net.d)))
    } else {
        colnames(net.d)
    }

    # Long-format data for ggplot
    data_to_plot <- net.d |>
        as.data.frame() |>
        tibble::rownames_to_column("DORC") |>
        tibble::as_tibble() |>
        tidyr::pivot_longer(cols = -DORC, names_to = "Motif", values_to = "Score") |>
        dplyr::mutate(
            DORC  = factor(DORC,  levels = row_order),
            Motif = factor(Motif, levels = col_order)
        )

    # Per-tick styling vectors matching col_order
    tf_label_colors <- ifelse(col_order %in% highlight_TFs, highlight_label_color, "black")
    tf_label_faces  <- ifelse(col_order %in% highlight_TFs, "bold", "plain")

    if (!is.null(DORC_groups)) {
        group_lookup <- tibble::tibble(
            DORC  = DORCs,
            Group = factor(DORC_groups, levels = unique(DORC_groups))
        )
        data_to_plot <- dplyr::left_join(data_to_plot, group_lookup, by = "DORC") |>
            dplyr::mutate(DORC = factor(DORC, levels = row_order))
    }

    p <- ggplot2::ggplot(data_to_plot, ggplot2::aes(x = Motif, y = DORC, fill = Score)) +
        ggplot2::geom_tile(color = "grey60", linewidth = 0.2) +
        ggplot2::scale_fill_gradientn(
            colors = color_palette,
            limits = limits,
            oob    = scales::squish,
            name   = "TF regulation\nscore"
        ) +
        ggplot2::labs(x = x_title, y = y_title, title = plot_title) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            panel.grid   = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
            axis.text.y  = ggplot2::element_text(face = "italic"),
            strip.text.y = ggplot2::element_text(face = "bold", angle = 0, size = 12),
            plot.title   = ggplot2::element_text(hjust = 0.5, size = 14),
            axis.text.x  = ggplot2::element_text(
                angle = 90, hjust = 1, vjust = 0.5,
                color = tf_label_colors,
                face  = tf_label_faces
            )
        ) +
        ggplot2::scale_x_discrete(labels = toupper)

    # Highlight rectangles around specified TF columns. Adjacent highlighted
    # TFs are merged into a single rectangle.
    if (!is.null(highlight_TFs) && highlight_box) {
        highlight_present <- intersect(highlight_TFs, col_order)

        if (length(highlight_present) > 0) {
            # Find positions of highlighted TFs along the x-axis, sorted
            highlight_positions <- sort(match(highlight_present, col_order))

            # Group consecutive positions into runs (gap > 1 starts a new run)
            run_id <- cumsum(c(1, diff(highlight_positions) != 1))

            # Compute xmin/xmax for each run (first and last position in run)
            runs <- data.frame(
                run  = run_id,
                pos  = highlight_positions
            ) |>
                dplyr::group_by(run) |>
                dplyr::summarise(
                    xmin = min(pos) - 0.5,
                    xmax = max(pos) + 0.5,
                    .groups = "drop"
                )

            if (!is.null(DORC_groups)) {
                # Per-facet rectangles: cross runs with each Group's height
                facet_heights <- data_to_plot |>
                    dplyr::group_by(Group) |>
                    dplyr::summarise(ymax = dplyr::n_distinct(DORC) + 0.5,
                                     .groups = "drop") |>
                    dplyr::mutate(ymin = 0.5)

                highlight_data <- tidyr::crossing(runs, facet_heights)
            } else {
                highlight_data <- runs |>
                    dplyr::mutate(
                        ymin = 0.5,
                        ymax = length(levels(data_to_plot$DORC)) + 0.5
                    )
            }

            p <- p + ggplot2::geom_rect(
                data = highlight_data,
                inherit.aes = FALSE,
                ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                color = highlight_box_color, fill = NA, linewidth = highlight_box_linewidth
            )
        }
    }

    if (!is.null(DORC_groups)) {
        p <- p + ggplot2::facet_grid(rows = ggplot2::vars(Group), scales = "free_y", space = "free_y") +
            ggplot2::labs(tag = group_tag) +
            ggplot2::theme(
                plot.tag.position = c(0.86, 0.5),
                legend.box.margin = ggplot2::margin(l = 22),
                plot.tag          = ggplot2::element_text(angle = 270)
            )
    }

    p
}