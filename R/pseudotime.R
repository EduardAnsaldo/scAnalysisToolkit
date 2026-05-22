# =============================================================================
# Rolling-average + LOESS smoothing along pseudotime, single chosen assay.
# Input: Seurat object with assays (e.g. "RNA", "ATAC", "DORC", "chromvar")
# and a pseudotime / NN-stimulation-time metadata column.
# -----------------------------------------------------------------------------
# Provides:
#   * smooth_feature()       — single gene/feature trace + LOESS fit
#   * plot_feature()         — scatter + LOESS overlay for one feature
#   * smooth_feature_matrix()— many features, returns smoothed/normed matrix
#   * plot_feature_heatmap() — heatmap of many features, ordered by peak time
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(zoo)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ComplexHeatmap)   # install.packages("BiocManager"); BiocManager::install("ComplexHeatmap")
  library(circlize)
})

# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

#' Extract a single feature's values from a Seurat assay (internal)
#'
#' Pulls one row from `GetAssayData(seu, assay, slot)` for the requested cells
#' and returns it as a plain numeric vector. Errors if the feature is missing
#' from the assay or the assay is missing from the object.
#'
#' @param seu      A Seurat object.
#' @param feature  Character(1). Row name to extract (gene, peak, motif, ...).
#' @param assay    Character(1). Assay to pull from. Must be in `Assays(seu)`.
#' @param slot     Character(1). Slot within the assay (e.g. "data", "counts").
#' @param cells    Character vector of cell barcodes (column names of `seu`).
#'                 The returned vector is in this exact order.
#'
#' @return Numeric vector of length `length(cells)`.
#' @keywords internal
.get_feature_vec <- function(seu, feature, assay, slot, cells) {
  stopifnot(assay %in% Assays(seu))
  mat <- GetAssayData(seu, assay = assay, slot = slot)
  if (!feature %in% rownames(mat)) {
    stop(sprintf("Feature '%s' not found in assay '%s'.", feature, assay))
  }
  as.numeric(mat[feature, cells])
}

#' Rolling mean + 1st–99th percentile min–max normalization (internal)
#'
#' Applies a centered rolling mean with window size `window_size`, then
#' rescales the smoothed values so that the 1st percentile maps to 0 and
#' the 99th percentile maps to 1. Values are clipped to `[0, 1]`.
#'
#' Robust to NAs in the input (uses `na.rm = TRUE` within each window via
#' `zoo::rollapply`) and to flat features (returns all-zero vector when
#' the 1–99\% range is zero, instead of `NaN`/`Inf`).
#'
#' @param x            Numeric vector (cells ordered by pseudotime upstream).
#' @param window_size  Integer. Number of cells per sliding window.
#'
#' @return Numeric vector of length `length(x) - window_size + 1`, values
#'   in `[0, 1]`. Returns all zeros if the feature has zero dynamic range.
#' @keywords internal
.rollmean_norm <- function(x, window_size) {
  # Rolling mean that tolerates NAs within a window
  s <- zoo::rollapply(x, width = window_size, FUN = mean,
                      align = "center", na.rm = TRUE)
  # If a whole window was NA, mean(., na.rm=TRUE) returns NaN — replace with NA
  s[is.nan(s)] <- NA_real_

  q <- quantile(s, probs = c(0.01, 0.99), na.rm = TRUE)
  rng <- q[2] - q[1]

  if (!is.finite(rng) || rng <= .Machine$double.eps) {
    # Flat / all-NA feature: return a vector of zeros (no dynamics to show)
    return(rep(0, length(s)))
  }

  n <- (s - q[1]) / rng
  n <- pmin(pmax(n, 0), 1)        # clip to [0,1]
  # Carry NAs through; the caller decides whether to skip them
  n
}

#' Order cells by a pseudotime metadata column (internal)
#'
#' Subsets to cells with non-NA pseudotime and returns them sorted by it.
#' If `cells` is NULL, all cells in `seu` are considered.
#'
#' @param seu      A Seurat object.
#' @param time_col Character(1). Name of the metadata column holding the
#'                 pseudotime / stimulation-time values (must be numeric).
#' @param cells    Character vector of cell barcodes, or NULL for all cells.
#'
#' @return Named list with two elements of equal length:
#'   \describe{
#'     \item{cells}{Cell barcodes, sorted by pseudotime.}
#'     \item{pt}{The corresponding numeric pseudotime values.}
#'   }
#' @keywords internal
.ordered_cells <- function(seu, time_col, cells) {
  if (is.null(cells)) cells <- colnames(seu)
  pt <- seu@meta.data[cells, time_col]
  keep  <- !is.na(pt)
  cells <- cells[keep]
  pt    <- pt[keep]
  ord   <- order(pt)
  list(cells = cells[ord], pt = pt[ord])
}

# =============================================================================
# 1. Single feature: smoothed trace + LOESS fit
# =============================================================================

#' Smooth one feature along pseudotime and fit a LOESS curve
#'
#' For a single feature (gene, peak, or motif) in a chosen Seurat assay:
#' orders cells by pseudotime, computes a centered rolling mean over
#' `window_size` cells, min-max normalizes the result to its 1st–99th
#' percentile range, and fits a LOESS smoother with span `loess_alpha`
#' against the (also rolling-mean smoothed) pseudotime axis.
#'
#' This is the per-feature engine behind \code{\link{plot_feature}}.
#'
#' @param seu          A Seurat object. Must contain `assay` and `time_col`.
#' @param feature      Character(1). Feature to extract — must be in
#'                     `rownames(GetAssayData(seu, assay, slot))`.
#' @param time_col     Character(1). Name of the metadata column holding the
#'                     pseudotime / NN stimulation time values.
#' @param assay        Character(1). Assay to pull from. Default `"DORC"`.
#'                     Other typical choices: `"RNA"`, `"ATAC"`, `"chromvar"`.
#' @param slot         Character(1). Slot within the assay. Default `"data"`.
#' @param window_size  Integer. Cells per sliding window. Default 100.
#' @param loess_alpha  Numeric in (0, 1\]. LOESS span. Default 0.1.
#' @param cells        Optional character vector of cell barcodes to use.
#'                     NULL (default) uses all cells in `seu`.
#'
#' @return A named list:
#'   \describe{
#'     \item{feature}{The feature name (echoed back).}
#'     \item{assay}{The assay name (echoed back).}
#'     \item{data}{A data.frame with columns
#'       \code{pseudotime} (smoothed time axis),
#'       \code{value} (smoothed + 1–99\% normalized feature signal), and
#'       \code{fit} (LOESS predictions on the smoothed grid).}
#'     \item{loess}{The fitted \code{loess} object, for downstream prediction
#'       or diagnostics.}
#'   }
#'
#' @seealso \code{\link{plot_feature}}, \code{\link{smooth_feature_matrix}}.
#'
#' @examples
#' \dontrun{
#' res <- smooth_feature(seu, "IRF4", "stim_pseudotime", assay = "DORC")
#' plot_feature(res)
#' }
smooth_feature <- function(seu,
                           feature,
                           time_col,
                           assay       = "DORC",
                           slot        = "data",
                           window_size = 100,
                           loess_alpha = 0.1,
                           cells       = NULL) {

  oc   <- .ordered_cells(seu, time_col, cells)
  vec  <- .get_feature_vec(seu, feature, assay, slot, oc$cells)

  val_s <- .rollmean_norm(vec, window_size)
  pt_s  <- zoo::rollmean(oc$pt, k = window_size, align = "center")

  df <- data.frame(pseudotime = pt_s, value = val_s)

  # Drop rows where either column is non-finite (NA / NaN / Inf).
  # LOESS calls a Fortran routine that errors on any non-finite values.
  ok <- is.finite(df$pseudotime) & is.finite(df$value)
  if (sum(ok) < 10) {
    stop(sprintf(
      "Not enough finite values for '%s' in assay '%s' to fit LOESS (n = %d after smoothing). Check for NAs in the assay or in '%s'.",
      feature, assay, sum(ok), time_col
    ))
  }
  df_fit <- df[ok, , drop = FALSE]

  fit <- tryCatch(
    loess(value ~ pseudotime, data = df_fit, span = loess_alpha),
    error = function(e) {
      stop(sprintf(
        "LOESS failed for '%s' (assay '%s'): %s\nTry a larger 'loess_alpha' or a larger 'window_size'.",
        feature, assay, conditionMessage(e)
      ))
    }
  )

  # Predict back onto the full smoothed grid (NA where input was non-finite)
  df$fit <- NA_real_
  df$fit[ok] <- predict(fit)

  list(
    feature = feature,
    assay   = assay,
    data    = df,
    loess   = fit
  )
}

# =============================================================================
# 2. Single feature: plot
# =============================================================================

#' Scatter + LOESS overlay plot for one smoothed feature
#'
#' Visualizes the output of \code{\link{smooth_feature}}: smoothed
#' (rolling-mean, 1–99\% normalized) values as a point cloud with the LOESS
#' fit drawn on top, against the smoothed pseudotime axis.
#'
#' @param result      Named list returned by \code{\link{smooth_feature}}.
#' @param point_alpha Numeric in \[0, 1\]. Alpha for the point cloud.
#'                    Default 0.3.
#' @param point_size  Numeric. Point size. Default 0.7.
#' @param line_color  Character. Color for both points and LOESS line.
#'                    Default `"#1B998B"` (teal).
#'
#' @return A `ggplot` object. Print or compose with `patchwork` / `cowplot`.
#'
#' @seealso \code{\link{smooth_feature}}.
plot_feature <- function(result,
                         point_alpha = 0.3,
                         point_size  = 0.7,
                         line_color  = "#1B998B") {
  df <- result$data
  ggplot(df, aes(pseudotime, value)) +
    geom_point(alpha = point_alpha, size = point_size, color = line_color) +
    geom_line(aes(y = fit), linewidth = 1.1, color = line_color) +
    labs(
      title = paste0(result$feature, " (", result$assay, ")"),
      x = "Smoothed NN stimulation time (rolling mean)",
      y = "Normalized signal (1–99% min-max)"
    ) +
    theme_classic(base_size = 13)
}

# =============================================================================
# 3. Many features: smoothed + normalized matrix
#    Rows = features, cols = sliding windows along pseudotime.
# =============================================================================

#' Smooth many features into a matrix for heatmap plotting
#'
#' Applies the same rolling-mean + 1st–99th percentile min–max normalization
#' as \code{\link{smooth_feature}} to every feature in `features`, optionally
#' followed by a LOESS fit (with values re-clipped to \[0, 1\]). Returns a
#' matrix suitable for \code{\link{plot_feature_heatmap}}: rows are features,
#' columns are sliding-window indices along the pseudotime axis.
#'
#' Each row is normalized independently, so the heatmap shows the
#' *shape / timing* of each feature's response rather than its absolute
#' magnitude. Features requested but not present in the assay are silently
#' dropped; if none remain, the function errors.
#'
#' @param seu          A Seurat object.
#' @param features     Character vector of features to include (genes, peaks,
#'                     or motifs depending on `assay`).
#' @param time_col     Character(1). Pseudotime metadata column.
#' @param assay        Character(1). Assay to pull from. Default `"DORC"`.
#' @param slot         Character(1). Slot within the assay. Default `"data"`.
#' @param window_size  Integer. Cells per sliding window. Default 100.
#' @param loess_alpha  Numeric in (0, 1\]. LOESS span. Default 0.1.
#' @param use_loess    Logical. If TRUE (default), apply LOESS smoothing to
#'                     each row after the rolling-mean step. Set FALSE on
#'                     sparse features (e.g. raw ATAC peaks) where LOESS may
#'                     issue numerical warnings.
#' @param cells        Optional cell barcodes to restrict to. Default NULL
#'                     (all cells).
#' @param verbose      Logical. Print a status message. Default TRUE.
#'
#' @return Named list:
#'   \describe{
#'     \item{matrix}{Numeric matrix, rows = features, cols = pseudotime
#'       windows, values in \[0, 1\].}
#'     \item{pseudotime}{Numeric vector of smoothed pseudotime values,
#'       one per column of `matrix`.}
#'     \item{assay}{Echo of the assay name, used by the plotter for titles.}
#'   }
#'
#' @seealso \code{\link{plot_feature_heatmap}}.
#'
#' @examples
#' \dontrun{
#' sm <- smooth_feature_matrix(seu, head(rownames(seu[["DORC"]]), 80),
#'                             "stim_pseudotime", assay = "DORC")
#' plot_feature_heatmap(sm)
#' }
smooth_feature_matrix <- function(seu,
                                  features,
                                  time_col,
                                  assay       = "DORC",
                                  slot        = "data",
                                  window_size = 100,
                                  loess_alpha = 0.1,
                                  use_loess   = TRUE,
                                  cells       = NULL,
                                  verbose     = TRUE) {

  oc <- .ordered_cells(seu, time_col, cells)
  mat <- GetAssayData(seu, assay = assay, slot = slot)
  features <- intersect(features, rownames(mat))
  if (length(features) == 0) stop("None of the requested features are in the assay.")

  if (verbose) message("Smoothing ", length(features), " features from assay '", assay, "'")

  pt_s <- zoo::rollmean(oc$pt, k = window_size, align = "center")

  rows <- lapply(features, function(f) {
    vec   <- as.numeric(mat[f, oc$cells])
    val_s <- .rollmean_norm(vec, window_size)

    if (use_loess) {
      ok <- is.finite(pt_s) & is.finite(val_s)
      # Only attempt LOESS if there are enough finite points
      if (sum(ok) >= 10) {
        df_tmp <- data.frame(pseudotime = pt_s[ok], value = val_s[ok])
        fit <- tryCatch(
          loess(value ~ pseudotime, data = df_tmp, span = loess_alpha),
          error   = function(e) NULL,
          warning = function(w) NULL
        )
        if (!is.null(fit)) {
          out      <- val_s  # preserve length and NA positions
          out[ok]  <- predict(fit)
          out      <- pmin(pmax(out, 0), 1, na.rm = FALSE)
          return(out)
        }
      }
    }
    val_s
  })

  smoothed <- do.call(rbind, rows)
  rownames(smoothed) <- features
  colnames(smoothed) <- sprintf("w%04d", seq_len(ncol(smoothed)))

  list(matrix = smoothed, pseudotime = pt_s, assay = assay)
}

# =============================================================================
# 4. Heatmap: features × pseudotime windows, ordered by peak time
# =============================================================================

#' Heatmap of many smoothed features along pseudotime
#'
#' Renders a `ComplexHeatmap::Heatmap` from the output of
#' \code{\link{smooth_feature_matrix}}. Rows are features, columns are
#' sliding-window positions along the pseudotime axis, and cell color
#' encodes the row-normalized signal (0 = feature's 1st percentile,
#' 1 = its 99th percentile).
#'
#' Use `order_by = "peak"` for the classic "wave" presentation where
#' rows are sorted by the column index of each feature's maximum —
#' early-responding features at the top, late-responding at the bottom.
#' Use `"hclust"` to instead group features by similarity of trajectory
#' shape regardless of timing.
#'
#' All text elements expose font controls. Anything left at its default
#' (or NULL for the family arguments) inherits from `base_family` and
#' `base_size`. Note: the chosen `fontfamily` must be available to the
#' active graphics device — see `extrafont` or `cairo_pdf()` for PDFs.
#'
#' @param smoothed       List from \code{\link{smooth_feature_matrix}}.
#' @param order_by       One of `"peak"`, `"hclust"`, or `"none"`. Default
#'                       `"peak"`.
#' @param show_row_names Logical. Show feature names on the rows. Default TRUE.
#'                       Turn off for >100 features unless plotting very tall.
#' @param title          Optional column title. NULL uses an auto-generated
#'                       title from the assay name.
#' @param col_low,col_mid,col_high  Colors for the 0 / 0.5 / 1 stops of the
#'                       color ramp. Defaults are a viridis-style purple →
#'                       yellow → red.
#' @param annotate_time  Logical. Draw a line plot of pseudotime above the
#'                       heatmap to indicate how time is distributed across
#'                       columns (useful since cells often pile up at the
#'                       endpoints). Default TRUE.
#' @param base_family    Default `fontfamily` for all text. Default "Helvetica".
#' @param base_size      Default `fontsize` for all text. Default 10.
#' @param row_names_family,row_names_size,row_names_face        Row label font.
#' @param col_title_family,col_title_size,col_title_face        Column title font.
#' @param row_title_family,row_title_size,row_title_face        Row title font
#'                       (only visible if rows are split into groups).
#' @param legend_title_family,legend_title_size,legend_title_face
#'                       Color-scale legend title font.
#' @param legend_labels_family,legend_labels_size,legend_labels_face
#'                       Color-scale legend tick-label font.
#' @param anno_name_family,anno_name_size,anno_name_face
#'                       Font for the "Pseudotime" label on the top annotation.
#'
#' @details
#' Each `*_family` / `*_size` / `*_face` argument is independent; leave one
#' at its default to inherit. `fontface` accepts `"plain"`, `"italic"`,
#' `"bold"`, or `"bold.italic"`.
#'
#' @return A `ComplexHeatmap::Heatmap` object. Print to draw, or compose
#'   with `+` and pass to `draw()` for legends/annotations control.
#'
#' @seealso \code{\link{smooth_feature_matrix}}.
plot_feature_heatmap <- function(smoothed,
                                 order_by      = c("peak", "hclust", "none"),
                                 show_row_names = TRUE,
                                 title         = NULL,
                                 col_low       = "#352A87",
                                 col_mid       = "#FAFA6E",
                                 col_high      = "#A50026",
                                 annotate_time = TRUE,
                                 # ---- font controls -----------------------
                                 # Global default (used wherever a specific
                                 # *_family / *_size is left NULL):
                                 base_family   = "Helvetica",
                                 base_size     = 10,
                                 # Per-element overrides (NULL = inherit base)
                                 row_names_family   = NULL, row_names_size   = 8,  row_names_face   = "plain",
                                 col_title_family   = NULL, col_title_size   = 12, col_title_face   = "bold",
                                 row_title_family   = NULL, row_title_size   = NULL, row_title_face = "plain",
                                 legend_title_family = NULL, legend_title_size = 10, legend_title_face = "bold",
                                 legend_labels_family = NULL, legend_labels_size = 9, legend_labels_face = "plain",
                                 anno_name_family    = NULL, anno_name_size    = 9,  anno_name_face   = "plain") {

  order_by <- match.arg(order_by)
  M  <- smoothed$matrix
  pt <- smoothed$pseudotime

  # --- order rows ---
  if (order_by == "peak") {
    peaks <- apply(M, 1, which.max)
    M <- M[order(peaks), , drop = FALSE]
  } else if (order_by == "hclust") {
    M <- M[hclust(dist(M))$order, , drop = FALSE]
  }

  # --- color ---
  col_fun <- circlize::colorRamp2(c(0, 0.5, 1), c(col_low, col_mid, col_high))

  # --- helper: build a gpar that inherits from base_* when args are NULL ---
  mk_gp <- function(family, size, face) {
    grid::gpar(
      fontfamily = family %||% base_family,
      fontsize   = size   %||% base_size,
      fontface   = face   %||% "plain"
    )
  }

  # --- top annotation: pseudotime axis ---
  top_anno <- NULL
  if (annotate_time) {
    top_anno <- HeatmapAnnotation(
      Pseudotime = anno_lines(pt, gp = gpar(col = "grey30"), height = unit(1, "cm")),
      annotation_name_side = "left",
      annotation_name_gp   = mk_gp(anno_name_family, anno_name_size, anno_name_face)
    )
  }

  Heatmap(
    M,
    name              = "Norm.\nsignal",
    col               = col_fun,
    cluster_rows      = FALSE,
    cluster_columns   = FALSE,
    show_column_names = FALSE,
    show_row_names    = show_row_names,
    row_names_gp      = mk_gp(row_names_family, row_names_size, row_names_face),
    top_annotation    = top_anno,
    column_title      = title %||% paste0(smoothed$assay, " — features ordered by peak time"),
    column_title_gp   = mk_gp(col_title_family, col_title_size, col_title_face),
    row_title_gp      = mk_gp(row_title_family, row_title_size, row_title_face),
    heatmap_legend_param = list(
      title_gp  = mk_gp(legend_title_family,  legend_title_size,  legend_title_face),
      labels_gp = mk_gp(legend_labels_family, legend_labels_size, legend_labels_face)
    )
  )
}

# Compat for older R without %||%
`%||%` <- function(a, b) if (!is.null(a)) a else b

# =============================================================================
# EXAMPLE USAGE
# =============================================================================
# # --- single feature, choose the assay you want ---
# res <- smooth_feature(
#   seu         = seu,
#   feature     = "IRF4",
#   time_col    = "stim_pseudotime",
#   assay       = "DORC",         # or "RNA", "ATAC", "chromvar"
#   slot        = "data",
#   window_size = 100,
#   loess_alpha = 0.1
# )
# plot_feature(res)
#
# # --- multi-feature heatmap ---
# top_dorc_genes <- head(rownames(seu[["DORC"]]), 80)
# sm <- smooth_feature_matrix(
#   seu         = seu,
#   features    = top_dorc_genes,
#   time_col    = "stim_pseudotime",
#   assay       = "DORC",
#   window_size = 100,
#   loess_alpha = 0.1,
#   use_loess   = TRUE
# )
# plot_feature_heatmap(sm, order_by = "peak")
#
# # ChromVAR TF activity heatmap (motif names live in the chromvar assay):
# tf_motifs <- head(rownames(seu[["chromvar"]]), 50)
# sm_tf <- smooth_feature_matrix(seu, tf_motifs, "stim_pseudotime",
#                                assay = "chromvar", slot = "data")
# plot_feature_heatmap(sm_tf, order_by = "peak",
#                      title = "TF motif activity along stim time")
#
# # --- font customization example ---
# # Set everything to Arial 11 globally, but italicize row names and
# # bump the column title size:
# plot_feature_heatmap(
#   sm,
#   base_family       = "Arial",
#   base_size         = 11,
#   row_names_face    = "italic",
#   row_names_size    = 9,
#   col_title_size    = 14,
#   legend_title_size = 11
# )
# -----------------------------------------------------------------------------
