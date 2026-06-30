# =============================================================================
# Pseudotime visualization
# -----------------------------------------------------------------------------
# Provides:
#   * plot_pseudotime_heatmap() — ggplot heatmap of cells (ordered by
#                                 pseudotime) x genes, with optional rolling-mean
#                                 / LOESS smoothing and a pseudotime annotation
#                                 bar.
# =============================================================================

#' Plot a heatmap of cells x genes, ordered by pseudotime
#'
#' @param seurat_obj  A Seurat object.
#' @param genes       Character vector of genes (features) to plot.
#' @param pseudotime_col Name of the column in `seurat_obj@meta.data` used to
#'                    order cells. Default "pseudotime".
#' @param assay       Assay to pull expression from. Default DefaultAssay(seurat_obj).
#' @param slot        Slot to pull expression from ("data", "counts", "scale.data").
#'                    Default "data" (log-normalized).
#' @param normalize   How to normalize each gene across cells before plotting:
#'                    "minmax"      -> scale to \[0, 1\] using absolute min/max
#'                    "minmax_pct"  -> scale to \[0, 1\] using 1st/99th percentile
#'                                     (clipped); matches the Ma et al. /
#'                                     SHARE-seq DORC convention.
#'                    "zscore"      -> per-gene z-score (clipped at +/-3)
#'                    "none"        -> no normalization
#' @param pct_low,pct_high Percentiles used when normalize = "minmax_pct".
#' @param smooth_window If non-NULL and > 1, apply a centered rolling mean of
#'                    this window across pseudotime-ordered cells (per gene).
#' @param loess_smooth If TRUE, fit a per-gene loess on the (rolling-averaged
#'                    and normalized) values vs. the (rolling-averaged)
#'                    pseudotime, and plot the loess fit.
#' @param loess_span  span / alpha for loess. Default 0.1.
#' @param gene_order  How to order genes on the y-axis:
#'                    "input", "peak", or "hclust".
#' @param cell_subset Optional character vector of cell names to restrict to.
#' @param show_pseudotime_bar If TRUE, add a thin annotation bar on top showing
#'                    pseudotime. The bar has NO colorbar legend, only a
#'                    label naming the pseudotime column.
#' @param palette     Color palette for the heatmap (viridis/Brewer name,
#'                    vector of colors, or a ggplot Scale object).
#' @param palette_direction +1 or -1.
#' @param pseudotime_palette Palette for the pseudotime annotation bar.
#' @param legend_title Optional title for the heatmap fill legend. If NULL
#'                    (default), an automatic label is built from `normalize`
#'                    (and `loess_smooth`).
#'
#' @return A ggplot object (or a patchwork if show_pseudotime_bar = TRUE).
#'
#' @importFrom Seurat DefaultAssay GetAssayData
#' @export
plot_pseudotime_heatmap <- function(seurat_obj,
                                    genes,
                                    pseudotime_col = "pseudotime",
                                    assay = NULL,
                                    slot = "data",
                                    normalize = c("minmax", "minmax_pct",
                                                  "zscore", "none"),
                                    pct_low  = 0.01,
                                    pct_high = 0.99,
                                    smooth_window = NULL,
                                    loess_smooth = FALSE,
                                    loess_span   = 0.1,
                                    gene_order = c("input", "peak", "hclust"),
                                    cell_subset = NULL,
                                    show_pseudotime_bar = TRUE,
                                    palette = NULL,
                                    palette_direction = 1,
                                    pseudotime_palette = "magma",
                                    legend_title = NULL) {

  normalize  <- match.arg(normalize)
  gene_order <- match.arg(gene_order)
  if (is.null(assay)) assay <- DefaultAssay(seurat_obj)

  if (is.null(palette)) {
    palette <- if (normalize == "zscore") c("#2166AC", "white", "#B2182B")
               else "viridis"
  }

  # ---- 1. Check pseudotime column ------------------------------------------
  meta <- seurat_obj@meta.data
  if (!pseudotime_col %in% colnames(meta)) {
    stop("Column '", pseudotime_col, "' not found in seurat_obj@meta.data.")
  }

  # ---- 2. Select & order cells ---------------------------------------------
  cells <- if (is.null(cell_subset)) colnames(seurat_obj) else cell_subset
  pt    <- meta[cells, pseudotime_col]
  keep  <- !is.na(pt)
  cells <- cells[keep]; pt <- pt[keep]
  ord   <- order(pt)
  cells <- cells[ord];  pt <- pt[ord]

  # ---- 3. Check genes & fetch expression -----------------------------------
  missing_genes <- setdiff(genes, rownames(seurat_obj[[assay]]))
  if (length(missing_genes) > 0) {
    warning("Dropping genes not found in assay '", assay, "': ",
            paste(head(missing_genes, 5), collapse = ", "),
            if (length(missing_genes) > 5) ", ..." else "")
    genes <- setdiff(genes, missing_genes)
  }
  if (length(genes) == 0) stop("No requested genes found in assay '", assay, "'.")

  expr <- GetAssayData(seurat_obj, assay = assay, slot = slot)[genes, cells, drop = FALSE]
  expr <- as.matrix(expr)

  # ---- helper: centered rolling mean ---------------------------------------
  roll_mean <- function(x, w) {
    sm <- stats::filter(x, rep(1 / w, w), sides = 2) |> as.numeric()
    sm |> zoo::na.locf(na.rm = FALSE) |>
          zoo::na.locf(fromLast = TRUE, na.rm = FALSE)
  }

  # ---- 4. Optional smoothing along pseudotime ------------------------------
  pt_smoothed <- pt
  if (!is.null(smooth_window) && smooth_window > 1) {
    w <- as.integer(smooth_window)
    expr <- t(apply(expr, 1, roll_mean, w = w))
    rownames(expr) <- genes; colnames(expr) <- cells
    pt_smoothed <- roll_mean(pt, w)
  }

  # ---- 5. Per-gene normalization -------------------------------------------
  if (normalize == "minmax") {
    expr <- t(apply(expr, 1, function(x) {
      rng <- range(x, na.rm = TRUE)
      if (diff(rng) == 0) rep(0, length(x))
      else (x - rng[1]) / diff(rng)
    }))
    rownames(expr) <- genes; colnames(expr) <- cells

  } else if (normalize == "minmax_pct") {
    if (pct_low >= pct_high || pct_low < 0 || pct_high > 1) {
      stop("pct_low must be < pct_high and both within [0, 1].")
    }
    expr <- t(apply(expr, 1, function(x) {
      qs <- stats::quantile(x, probs = c(pct_low, pct_high),
                            na.rm = TRUE, names = FALSE)
      lo <- qs[1]; hi <- qs[2]
      if (!is.finite(lo) || !is.finite(hi) || hi - lo == 0) {
        rep(0, length(x))
      } else {
        y <- (x - lo) / (hi - lo)
        pmin(pmax(y, 0), 1)
      }
    }))
    rownames(expr) <- genes; colnames(expr) <- cells

  } else if (normalize == "zscore") {
    expr <- t(scale(t(expr)))
    expr[is.nan(expr)] <- 0
    expr[expr >  3] <-  3
    expr[expr < -3] <- -3
  }

  # ---- 6. Optional loess fit vs. (smoothed) pseudotime ---------------------
  if (loess_smooth) {
    expr <- t(apply(expr, 1, function(y) {
      ok <- is.finite(y) & is.finite(pt_smoothed)
      if (sum(ok) < 4 || length(unique(pt_smoothed[ok])) < 4) return(y)
      fit <- tryCatch(
        stats::loess(y[ok] ~ pt_smoothed[ok],
                     span = loess_span, degree = 2,
                     na.action = stats::na.exclude,
                     control = stats::loess.control(surface = "direct")),
        error = function(e) NULL
      )
      if (is.null(fit)) return(y)
      stats::predict(fit, newdata = pt_smoothed)
    }))
    rownames(expr) <- genes; colnames(expr) <- cells
    if (normalize %in% c("minmax", "minmax_pct")) {
      expr[] <- pmin(pmax(expr, 0), 1)
    }
  }

  # ---- 7. Determine gene order ---------------------------------------------
  gene_levels <- switch(
    gene_order,
    input  = genes,
    peak   = genes[order(apply(expr, 1, which.max))],
    hclust = { hc <- hclust(dist(expr)); genes[hc$order] }
  )

  # ---- 8. Long-format data for ggplot --------------------------------------
  df <- as.data.frame(expr) |>
    tibble::rownames_to_column("gene") |>
    tidyr::pivot_longer(-gene, names_to = "cell", values_to = "expression") |>
    mutate(
      cell = factor(cell, levels = cells),
      gene = factor(gene, levels = rev(gene_levels))
    )

  # ---- 9. Resolve palette --> a ggplot fill scale --------------------------
  make_fill_scale <- function(palette, name, diverging = FALSE,
                              direction = 1, guide = "colourbar") {
    if (inherits(palette, "Scale")) return(palette)

    viridis_opts <- c("viridis", "magma", "plasma", "inferno",
                      "cividis", "mako", "rocket", "turbo")

    if (is.character(palette) && length(palette) == 1) {
      if (palette %in% viridis_opts) {
        return(ggplot2::scale_fill_viridis_c(
          option = palette, direction = direction, name = name, guide = guide
        ))
      }
      if (requireNamespace("RColorBrewer", quietly = TRUE) &&
          palette %in% rownames(RColorBrewer::brewer.pal.info)) {
        info <- RColorBrewer::brewer.pal.info[palette, ]
        cols <- RColorBrewer::brewer.pal(info$maxcolors, palette)
        if (direction < 0) cols <- rev(cols)
        if (diverging) {
          return(ggplot2::scale_fill_gradientn(
            colours = cols, name = name, guide = guide,
            rescaler = function(x, to = c(0, 1), from = NULL) {
              m <- max(abs(x), na.rm = TRUE)
              scales::rescale(x, to = to, from = c(-m, m))
            }
          ))
        }
        return(ggplot2::scale_fill_gradientn(colours = cols, name = name,
                                             guide = guide))
      }
      stop("Unrecognized palette name: '", palette, "'.")
    }

    if (is.character(palette) && length(palette) >= 2) {
      cols <- if (direction < 0) rev(palette) else palette
      if (diverging) {
        return(ggplot2::scale_fill_gradientn(
          colours = cols, name = name, guide = guide,
          rescaler = function(x, to = c(0, 1), from = NULL) {
            m <- max(abs(x), na.rm = TRUE)
            scales::rescale(x, to = to, from = c(-m, m))
          }
        ))
      }
      return(ggplot2::scale_fill_gradientn(colours = cols, name = name,
                                           guide = guide))
    }

    stop("`palette` must be a Scale, a palette name, or a vector of >=2 colors.")
  }

  # ---- 10. Main heatmap ----------------------------------------------------
  fill_label <- switch(normalize,
                      minmax     = "Expression\n(min-max)",
                      minmax_pct = sprintf("Expression\n(min-max,\n%g-%g%%)",
                                           pct_low * 100, pct_high * 100),
                      zscore     = "Expression\n(z-score)",
                      none       = "Expression")
  if (loess_smooth) {
    fill_label <- paste0(fill_label, sprintf("\nloess α=%g", loess_span))
  }
  if (!is.null(legend_title)) fill_label <- legend_title

  fill_scale <- make_fill_scale(
    palette   = palette,
    name      = fill_label,
    diverging = (normalize == "zscore"),
    direction = palette_direction
  )

  p <- ggplot(df, aes(x = cell, y = gene, fill = expression)) +
    geom_raster() +
    fill_scale +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = paste0("Cells ordered by ", pseudotime_col, " →"),
         y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid   = element_blank(),
      axis.text.y  = element_text(face = "italic"),
      plot.margin  = margin(5, 10, 5, 10)
    )

  # ---- 11. Optional pseudotime annotation bar ------------------------------
  # No colorbar legend: just a label on the left naming the pseudotime column.
  if (show_pseudotime_bar) {
    bar_df <- data.frame(cell = factor(cells, levels = cells), pseudotime = pt)
    pt_scale <- make_fill_scale(pseudotime_palette,
                                name = NULL,
                                guide = "none")
    pbar <- ggplot(bar_df, aes(x = cell, y = 1, fill = pseudotime)) +
      geom_raster() +
      pt_scale +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0),
                         breaks = 1,
                         labels = pseudotime_col) +
      labs(y = NULL) +
      theme_void() +
      theme(
        legend.position = "none",
        axis.text.y     = element_text(size = 9, hjust = 1,
                                       margin = margin(r = 4)),
        plot.margin     = margin(5, 10, 0, 10)
      )

    if (!requireNamespace("patchwork", quietly = TRUE)) {
      warning("Install 'patchwork' to show the pseudotime annotation bar.")
      return(p)
    }
    return(patchwork::wrap_plots(pbar, p, ncol = 1, heights = c(0.05, 1)))
  }

  p
}
