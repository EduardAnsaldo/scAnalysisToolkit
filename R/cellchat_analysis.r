## CellChat Analysis Functions
## Functions for cell-cell communication analysis

#' Run CellChat Analysis on Seurat Object
#'
#' @param seurat A Seurat object with cell type annotations
#' @param group_by Character string specifying the metadata column for cell grouping
#' @param species Character string, either "mouse" or "human"
#' @param smooth_data Logical, whether to smooth data using PPI networkw
#' @param subset_cell_groups Character vector of cell groups to subset, or NULL for all
#' @param subset_db Character string to filter database interactions, or NULL for all.
#'   Common values: "Cell-Cell Contact", "ECM-Receptor", "Secreted Signaling"
#' @param min_cells Integer, minimum number of cells required for communication
#' @param use_parallel Logical, whether to use parallel processing
#' @param n_workers Integer, number of parallel workers (only used if use_parallel = TRUE)
#' @param max_memory_mb Numeric, maximum memory for future in MB (only used if use_parallel = TRUE)
#' @param population_size Logical, whether to consider cell population size effect
#' @param nboot Integer, number of bootstrap samples for computing p-values
#' @param output_file Character string, path to save RDS file, or NULL to skip saving
#'
#' @return A CellChat object
#' @export
run_cellchat_analysis <- function(seurat,
           group_by = "major_cell_types",
           species = "mouse",
           smooth_data = FALSE,
           subset_cell_groups = NULL,
           subset_db = NULL,
           min_cells = 10,
           method = 'triMean',
           trim = NULL,
           use_parallel = TRUE,
           n_workers = 10,
           max_memory_mb = 4000,
           population_size = FALSE,
           nboot = 100,

           output_file = NULL) {

  # Check for required packages
  if (!requireNamespace("CellChat", quietly = TRUE)) {
    stop("Package 'CellChat' is required but not installed. Please install it with:\n",
         "  devtools::install_github('jinworks/CellChat')",
         call. = FALSE)
  }

  # Part I: Setup and preprocessing
  message("Starting CellChat analysis...")

  ## LogNormalize data
  message("Normalizing data...")
  seurat <- NormalizeData(seurat, assay = "RNA", normalization.method = "LogNormalize")

  ## Extract expression data and metadata from the Seurat object
  if ("Samples" %in% colnames(seurat@meta.data)) {
    seurat@meta.data <- seurat@meta.data |> mutate(samples = Samples)
  }

  message("Creating CellChat object...")
  cellchat_object <- CellChat::createCellChat(object = seurat, group.by = group_by, assay = "RNA")

  ## Set the CellChat database
  cellchat_db <- if (species == "mouse") {
    CellChat::CellChatDB.mouse
  } else if (species == "human") {
    CellChat::CellChatDB.human
  } else {
    stop("Species must be either 'mouse' or 'human'")
  }

  ## Filter database to specific interaction types if requested
  if (!is.null(subset_db)) {
    message(sprintf("Filtering database to: %s", subset_db))
    cellchat_db <- CellChat::subsetDB(cellchat_db, search = subset_db, key = 'annotation')
  }

  cellchat_object@DB <- cellchat_db

  ## Preprocessing the data for CellChat analysis
  message("Preprocessing data...")
  cellchat_object <- subsetData(cellchat_object)

  # Set up parallel processing if requested
  if (use_parallel) {
    message(sprintf("Setting up parallel processing with %d workers...", n_workers))
    future::plan("multisession", workers = n_workers)
    options(future.globals.maxSize = max_memory_mb * 1024^2)
  } else {
    message("Running in sequential mode...")
    future::plan("sequential")
  }

  cellchat_object <- CellChat::identifyOverExpressedGenes(cellchat_object)

  message("Identifying over-expressed interactions...")
  ptm <- Sys.time()
  cellchat_object <- CellChat::identifyOverExpressedInteractions(cellchat_object)
  execution_time <- Sys.time() - ptm
  message(sprintf("Execution time: %.2f mins", as.numeric(execution_time, units = "mins")))

  ## Optional: smooth data using PPI network
  if (smooth_data) {
    message("Smoothing data using PPI network...")
    ppi_network <- if (species == "mouse") CellChat::PPI.mouse else CellChat::PPI.human
    cellchat_object <- smoothData(cellchat_object, adj = ppi_network)
  }

  # Part II: Compute communication probabilities
  message(sprintf("Computing communication probabilities (nboot = %d)...", nboot))
  ptm <- Sys.time()
  cellchat_object <- CellChat::computeCommunProb(
    cellchat_object,
    type = method,
    raw.use = !smooth_data,
    population.size = population_size,
    seed.use = 3514L,
    nboot = nboot,
    trim = trim
  )
  execution_time <- Sys.time() - ptm
  message(sprintf("Execution time: %.2f mins", as.numeric(execution_time, units = "mins")))

  ## Filter out cell-cell communication
  message(sprintf("Filtering communication (min.cells = %d)...", min_cells))
  cellchat_object <- CellChat::filterCommunication(cellchat_object, min.cells = min_cells)

  ## Infer signaling pathway level communication
  message("Computing communication probabilities at pathway level...")
  ptm <- Sys.time()
  cellchat_object <- CellChat::computeCommunProbPathway(cellchat_object)
  execution_time <- Sys.time() - ptm
  message(sprintf("Execution time: %.2f mins", as.numeric(execution_time, units = "mins")))

  ## Calculate aggregated network
  message("Aggregating cell-cell communication network...")
  if (!is.null(subset_cell_groups)) {
    cellchat_object <- CellChat::aggregateNet(
      cellchat_object,
      sources.use = subset_cell_groups,
      targets.use = subset_cell_groups
    )
  } else {
    cellchat_object <- CellChat::aggregateNet(cellchat_object)
  }

  # Save if output file specified
  if (!is.null(output_file)) {
    message(sprintf("Saving CellChat object to: %s", output_file))
    saveRDS(cellchat_object, file = output_file)
  }

  future::plan("sequential")  # Reset to sequential plan
  message("CellChat analysis complete!")
  return(invisible(NULL))
}

# Example usage:
# seurat <- readRDS(here::here('data/SPF_GF_annotated.rds'))
# cellchat_object <- run_cellchat_analysis(
#   seurat,
#   group_by = "major_cell_types",
#   species = "mouse",
#   use_parallel = TRUE,
#   output_file = here::here('data/CellChat_SPF_GF_major_cell_types.rds')
# )
#
# # Example with database subsetting for cell-cell contact interactions only:
# cellchat_object <- run_cellchat_analysis(
#   seurat,
#   group_by = "major_cell_types",
#   species = "mouse",
#   subset_db = "Cell-Cell Contact",
#   use_parallel = TRUE,
#   output_file = here::here('data/CellChat_SPF_GF_cell_contact.rds')
# )

#' Bubble plot of pathway-level cell-cell communication
#'
#' Visualises the aggregated signalling pathway probabilities stored in
#' \code{object@@netP$prob}.  Each bubble represents one (sender, receiver,
#' pathway) triplet.  Bubble \strong{size} and \strong{colour} both encode
#' the communication probability; no p-value is available at the pathway level.
#'
#' This function is the pathway-level analogue of
#' \code{\link[CellChat]{netVisual_bubble}}, which operates on individual
#' ligand-receptor pairs (\code{object@@net}).  The argument interface is kept
#' as close to the original as possible so the two can be used
#' interchangeably.
#'
#' @param object A \code{\linkS4class{CellChat}} object whose
#'   \code{@@netP$prob} slot has been populated by
#'   \code{\link[CellChat]{computeCommunProbPathway}}.
#' @param sources.use Integer indices \emph{or} character names of sender cell
#'   groups to include.  \code{NULL} (default) retains all groups.
#' @param targets.use Integer indices \emph{or} character names of receiver
#'   cell groups to include.  \code{NULL} (default) retains all groups.
#' @param signaling Character vector of pathway names to display.  \code{NULL}
#'   (default) shows all pathways that survive \code{min.prob} filtering.
#'   When both \code{signaling} and \code{top.n} are supplied, \code{top.n}
#'   is applied \emph{within} the user-specified pathway set.
#' @param top.n Integer.  Retain only the \code{top.n} pathways ranked by
#'   their aggregated probability across all retained cell pairs (after
#'   \code{min.prob} and \code{signaling} filtering).  The ranking statistic
#'   is controlled by \code{top.n.by}.  Pathways are plotted top-to-bottom
#'   in descending rank order.  \code{NULL} (default) keeps all pathways.
#' @param top.n.by Character string specifying the statistic used to rank
#'   pathways when \code{top.n} is not \code{NULL}.  One of \code{"mean"}
#'   (default), \code{"max"}, or \code{"sum"}.
#' @param remove.isolate Logical.  When \code{TRUE}, cell-pair columns that
#'   carry zero probability across all retained pathways are dropped (default
#'   \code{FALSE}).
#' @param min.prob Numeric.  Bubbles whose probability is at or below this
#'   value are suppressed (default \code{0}).
#' @param sort.by.target Logical.  When \code{TRUE}, cell-pair columns are
#'   ordered by target group first, then source group.  Default \code{FALSE}
#'   orders by source first.
#' @param do.flip Logical.  Transpose the axes so that pathways appear on the
#'   x-axis and cell-group pairs on the y-axis (default \code{FALSE}).
#' @param show.lr.pairs Logical.  When \code{TRUE}, each pathway's axis tick
#'   label is replaced by the ligand-receptor pairs that make up the pathway
#'   (parsed from \code{object@@LR$LRsig$interaction_name_2}).  Ligands that
#'   share an identical receptor set are grouped and combined with \code{"/"}
#'   (e.g. \code{"TGFB1/TGFB2/TGFB3 - (TGFBR1+TGFBR2)"}); genuinely
#'   many-to-many pathways span multiple lines.  The affected axis title also
#'   changes to \code{"Ligand-receptor pairs"}.  Requires an object populated
#'   through \code{\link[CellChat]{identifyOverExpressedInteractions}}.
#'   Default \code{FALSE}.
#' @param color.use Controls bubble fill colour.  Two modes:
#'   \describe{
#'     \item{\code{NULL} (default)}{Continuous gradient controlled by
#'       \code{color.heatmap}.}
#'     \item{Named character vector}{Maps pathway names to hex colours for
#'       categorical colouring, e.g.
#'       \code{c(CXCL = "#e41a1c", VEGF = "#377eb8")}.  Pathways absent from
#'       the vector receive \code{NA} fill.}
#'   }
#' @param color.heatmap Character string naming the continuous colour palette,
#'   following \code{\link[CellChat]{netVisual_bubble}}.  Any
#'   \code{RColorBrewer} palette name (default \code{"Spectral"}); an
#'   unrecognised name falls back to a \pkg{viridis} \code{option} (e.g.
#'   \code{"magma"}).  Colours run low-to-high after CellChat's reversal
#'   (\code{direction = -1}).  Ignored when \code{color.use} is a named vector.
#' @param min.quantile,max.quantile Numeric quantile cutoffs in \code{[0, 1]}
#'   (default \code{0} / \code{1}).  Following
#'   \code{\link[CellChat]{netVisual_bubble}}, probabilities below/above these
#'   quantiles are clamped before colouring, which sets the colour-scale range
#'   and trims outliers; the two colourbar ends are labelled \code{"min"} and
#'   \code{"max"}.  Defaults leave the data unchanged.
#' @param show.legend Logical; whether to draw the legend (default \code{TRUE}).
#' @param scale.size Logical.  When \code{TRUE} (default) bubble size encodes
#'   the communication probability.  When \code{FALSE} all bubbles are drawn at
#'   a single fixed size (\code{dot.size}) and the size legend is removed;
#'   colour then remains the only probability encoding.
#' @param size.range Length-2 numeric giving the minimum and maximum bubble
#'   size passed to \code{\link[ggplot2]{scale_size_continuous}} when
#'   \code{scale.size = TRUE} (default \code{c(1, 10)}; e.g. \code{c(1, 5)} for
#'   a tighter spread).
#' @param dot.size Numeric; the single fixed bubble size used when
#'   \code{scale.size = FALSE} (default \code{4}, in the same units as
#'   \code{size.range}).  Ignored when \code{scale.size = TRUE}.
#' @param font.size Base font size in pt passed to
#'   \code{\link[ggplot2]{theme_classic}} (default \code{10}).
#' @param font.size.title Plot title font size in pt (default \code{10}).
#' @param angle.x Rotation of x-axis tick labels in degrees (default
#'   \code{90}).
#' @param vjust.x Vertical justification of x-axis labels.  \code{NULL}
#'   (default) sets \code{1} when \code{angle.x != 0}, else \code{0.5}.
#' @param hjust.x Horizontal justification of x-axis labels.  \code{NULL}
#'   (default) sets \code{1} when \code{angle.x != 0}, else \code{0.5}.
#' @param title.name Character string for the plot title.  \code{NULL}
#'   (default) produces no title.
#' @param return.data Logical.  When \code{TRUE} a named list is returned
#'   instead of a \code{ggplot} object:
#'   \describe{
#'     \item{\code{$communication}}{Tidy data frame with columns
#'       \code{source}, \code{target}, \code{pathway}, \code{prob},
#'       \code{pair}.}
#'     \item{\code{$gg.obj}}{The \code{ggplot} object.}
#'   }
#'   Default \code{FALSE}.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object, or (when
#'   \code{return.data = TRUE}) a named list; see above.
#'
#' @seealso
#'   \code{\link[CellChat]{netVisual_bubble}} for the LR-pair-level equivalent,
#'   \code{\link[CellChat]{computeCommunProbPathway}} to populate
#'   \code{@@netP$prob}.
#'
#' @examples
#' \dontrun{
#' # cellchat must have been processed with computeCommunProbPathway()
#'
#' # Default: all pathways, viridis colour scale
#' netVisual_bubble_pathways(cellchat)
#'
#' # Top 15 pathways by mean probability, magma palette
#' netVisual_bubble_pathways(
#'   cellchat,
#'   top.n         = 15,
#'   top.n.by      = "mean",
#'   color.heatmap = "magma"
#' )
#'
#' # Filter senders/receivers and pathways, top 10 by max
#' netVisual_bubble_pathways(
#'   cellchat,
#'   sources.use    = c("Macrophage", "Fibroblast"),
#'   targets.use    = c("B cell", "T cell"),
#'   top.n          = 10,
#'   top.n.by       = "max",
#'   remove.isolate = TRUE,
#'   title.name     = "Top 10 pathways"
#' )
#'
#' # Return tidy data for downstream use
#' out <- netVisual_bubble_pathways(cellchat, top.n = 20, return.data = TRUE)
#' head(out$communication)
#' out$gg.obj
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_fill_manual scale_color_manual scale_fill_gradientn scale_color_gradientn scale_size_continuous scale_x_discrete scale_y_discrete labs theme_classic theme element_text element_rect element_line element_blank unit
#' @importFrom rlang .data
#' @importFrom methods .hasSlot
#'
#' @export
netVisual_bubble_pathways <- function(
    object,
    sources.use     = NULL,
    targets.use     = NULL,
    signaling       = NULL,
    top.n           = NULL,
    top.n.by        = c("mean", "max", "sum"),
    remove.isolate  = FALSE,
    min.prob        = 0,
    sort.by.target  = FALSE,
    do.flip         = FALSE,
    show.lr.pairs   = FALSE,
    color.use       = NULL,
    color.heatmap   = "Spectral",
    min.quantile    = 0,
    max.quantile    = 1,
    show.legend     = TRUE,
    scale.size      = TRUE,
    size.range      = c(1, 10),
    dot.size        = 4,
    font.size       = 10,
    font.size.title = 10,
    angle.x         = 90,
    vjust.x         = NULL,
    hjust.x         = NULL,
    title.name      = NULL,
    return.data     = FALSE
) {
  # ---- 0. Validate scalar arguments -----------------------------------------
  top.n.by <- match.arg(top.n.by)

  if (!is.null(top.n)) {
    top.n <- as.integer(top.n)
    if (is.na(top.n) || top.n < 1L) stop("'top.n' must be a positive integer.")
  }

  if (!is.character(color.heatmap) || length(color.heatmap) != 1L) {
    stop("'color.heatmap' must be a single palette name (character string).")
  }
  if (!is.numeric(min.quantile) || !is.numeric(max.quantile) ||
      length(min.quantile) != 1L || length(max.quantile) != 1L ||
      anyNA(c(min.quantile, max.quantile)) ||
      min.quantile < 0 || max.quantile > 1 || min.quantile > max.quantile) {
    stop("'min.quantile'/'max.quantile' must be single numbers in [0, 1] with min <= max.")
  }

  if (!is.numeric(size.range) || length(size.range) != 2L || anyNA(size.range) ||
      any(size.range <= 0) || size.range[1] > size.range[2]) {
    stop("'size.range' must be a length-2 numeric of positive, non-decreasing values.")
  }
  if (!is.numeric(dot.size) || length(dot.size) != 1L || is.na(dot.size) ||
      dot.size <= 0) {
    stop("'dot.size' must be a single positive number.")
  }

  # ---- 1. Extract pathway probability array from @netP ----------------------
  if (!methods::.hasSlot(object, "netP") || is.null(object@netP$prob)) {
    stop(
      "'object@netP$prob' is empty or missing.\n",
      "Run  cellchat <- computeCommunProbPathway(cellchat)  first."
    )
  }

  prob3d        <- object@netP$prob          # [n_groups x n_groups x n_pathways]
  group_names   <- dimnames(prob3d)[[1]]
  pathway_names <- dimnames(prob3d)[[3]]

  if (is.null(group_names)) {
    stop("'dimnames(object@netP$prob)[[1]]' is NULL; cannot identify cell groups.")
  }

  # ---- 2. Resolve sources / targets to character names ----------------------
  sources_use <- .resolve_groups(sources.use, group_names, "sources.use")
  targets_use <- .resolve_groups(targets.use, group_names, "targets.use")

  # ---- 3. Optionally subset pathways via 'signaling' ------------------------
  if (!is.null(signaling)) {
    unknown <- setdiff(signaling, pathway_names)
    if (length(unknown)) {
      warning(
        "The following pathways were not found in object@netP$prob and will ",
        "be ignored: ", paste(unknown, collapse = ", ")
      )
    }
    signaling <- intersect(signaling, pathway_names)
    if (length(signaling) == 0L) {
      stop("None of the requested 'signaling' pathways exist in object@netP$prob.")
    }
    prob3d        <- prob3d[, , signaling, drop = FALSE]
    pathway_names <- signaling
  }

  # ---- 4. Reshape 3-D array to long data frame ------------------------------
  # columns: source, target, pathway, prob
  df <- reshape2::melt(
    prob3d,
    varnames   = c("source", "target", "pathway"),
    value.name = "prob",
    as.is      = TRUE
  )

  # ---- 5. Filter by source / target -----------------------------------------
  df <- df[df$source %in% sources_use & df$target %in% targets_use, , drop = FALSE]

  # ---- 6. Apply probability threshold ---------------------------------------
  df <- df[df$prob > min.prob, , drop = FALSE]

  if (nrow(df) == 0L) {
    stop(
      "No interactions survive the current filters. ",
      "Try lowering min.prob, broadening sources/targets, or check @netP$prob."
    )
  }

  # ---- 7. Retain top-N pathways by aggregate probability --------------------
  #
  # Aggregate over all *retained* cell pairs so the ranking reflects the same
  # data that will be plotted.  Pathways are ordered descending on the plot
  # (highest probability at the top of the Y-axis).
  #
  agg_fn <- switch(
    top.n.by,
    mean = function(x) mean(x, na.rm = TRUE),
    max  = function(x) max(x,  na.rm = TRUE),
    sum  = function(x) sum(x,  na.rm = TRUE)
  )

  pathway_agg   <- tapply(df$prob, df$pathway, agg_fn)
  ranked_paths  <- names(sort(pathway_agg, decreasing = TRUE))   # best first

  if (!is.null(top.n)) {
    n_keep       <- min(top.n, length(ranked_paths))
    ranked_paths <- ranked_paths[seq_len(n_keep)]
    df           <- df[df$pathway %in% ranked_paths, , drop = FALSE]
    if (nrow(df) == 0L) {
      stop("No rows remain after top.n filtering.")   # should not happen
    }
  }

  # ---- 8. Sort rows and build cell-pair label --------------------------------
  if (sort.by.target) {
    df <- df[order(df$target, df$source), , drop = FALSE]
  } else {
    df <- df[order(df$source, df$target), , drop = FALSE]
  }
  df$pair <- paste0(df$source, " \u2192 ", df$target)

  # ---- 9. Remove isolates (pairs with no surviving probability) -------------
  if (remove.isolate) {
    active_pairs <- unique(df$pair[df$prob > min.prob])
    df <- df[df$pair %in% active_pairs, , drop = FALSE]
    if (nrow(df) == 0L) {
      stop(
        "No cell pairs remain after remove.isolate = TRUE. ",
        "Lower min.prob or set remove.isolate = FALSE."
      )
    }
  }

  # ---- 10. Set factor levels ------------------------------------------------
  #
  # Pairs: preserve left-to-right display order (encounter order after sort).
  # Pathways: ranked_paths is descending (best first); rev() puts the best
  # pathway as the *last* factor level, which ggplot2 renders at the TOP of
  # the Y-axis.
  #
  df$pair    <- factor(df$pair,    levels = unique(df$pair))
  df$pathway <- factor(df$pathway, levels = rev(ranked_paths))

  # ---- 10b. Optional ligand-receptor-pair axis labels -----------------------
  # Named vector mapping each pathway (factor level) to its collapsed LR-pair
  # label; applied to the pathway axis as tick labels further below.  Computed
  # here so a missing @LR$LRsig fails fast, before any plotting.
  lr_map <- if (show.lr.pairs) .pathway_lr_labels(object, levels(df$pathway)) else NULL

  # ---- 11. Axis assignment (respects do.flip) --------------------------------
  pathway_axis_lab <- if (show.lr.pairs) "Ligand-receptor pairs" else "Signalling Pathway"
  if (do.flip) {
    x_var <- "pathway"; y_var <- "pair"
    x_lab <- pathway_axis_lab; y_lab <- "Cell-group pair"
  } else {
    x_var <- "pair"; y_var <- "pathway"
    x_lab <- "Cell-group pair"; y_lab <- pathway_axis_lab
  }

  # ---- 12. Resolve label-justification defaults -----------------------------
  if (is.null(vjust.x)) vjust.x <- if (angle.x != 0) 1 else 0.5
  if (is.null(hjust.x)) hjust.x <- if (angle.x != 0) 1 else 0.5

  # ---- 13. Build ggplot ------------------------------------------------------
  #
  # Colour handling mirrors CellChat::netVisual_bubble: clamp the plotted
  # probability to the requested quantile cutoffs (a no-op at the 0/1 defaults),
  # then span the colour scale over the clamped min-max with the two ends
  # labelled "min"/"max".  A degenerate range (all values equal) falls back to
  # the scale's default limits/breaks.
  min.cutoff <- stats::quantile(df$prob, min.quantile, na.rm = TRUE, names = FALSE)
  max.cutoff <- stats::quantile(df$prob, max.quantile, na.rm = TRUE, names = FALSE)
  df$prob[df$prob < min.cutoff] <- min.cutoff
  df$prob[df$prob > max.cutoff] <- max.cutoff

  prob_min <- min.cutoff   # == quantile(df$prob, 0) after clamping
  prob_max <- max.cutoff   # == quantile(df$prob, 1) after clamping
  if (prob_min != prob_max) {
    prob_limits <- c(prob_min, prob_max)
    prob_breaks <- c(prob_min, prob_max)
    prob_labels <- c("min", "max")
  } else {
    prob_limits <- NULL
    prob_breaks <- ggplot2::waiver()
    prob_labels <- ggplot2::waiver()
  }
  grad_colors <- .heatmap_colors(color.heatmap)

  if (!is.null(color.use) &&
      is.character(color.use) &&
      !is.null(names(color.use))) {

    # Categorical: one explicit colour per pathway --------------------------
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x    = .data[[x_var]],
        y    = .data[[y_var]],
        size = .data[["prob"]],
        fill = .data[["pathway"]],
        col  = .data[["pathway"]]
      )
    ) +
      ggplot2::geom_point(shape = 21, stroke = 0.35) +
      ggplot2::scale_fill_manual(values = color.use, name = "Pathway") +
      ggplot2::scale_color_manual(values = color.use, guide = "none")

  } else {

    # Continuous CellChat-style gradient ------------------------------------
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x    = .data[[x_var]],
        y    = .data[[y_var]],
        size = .data[["prob"]],
        fill = .data[["prob"]],
        col  = .data[["prob"]]
      )
    ) +
      ggplot2::geom_point(shape = 21, stroke = 0.35) +
      ggplot2::scale_fill_gradientn(
        colors   = grad_colors,
        na.value = "white",
        name     = "Probability",
        limits   = prob_limits,
        breaks   = prob_breaks,
        labels   = prob_labels
      ) +
      ggplot2::scale_color_gradientn(
        colors   = grad_colors,
        na.value = "white",
        guide    = "none",
        limits   = prob_limits
      )
  }

  # Size scale: encode probability across 'size.range', or collapse to a single
  # fixed 'dot.size' (legend hidden) when 'scale.size = FALSE' so colour becomes
  # the only probability encoding.
  if (scale.size) {
    p <- p + ggplot2::scale_size_continuous(range = size.range, name = "Probability")
  } else {
    p <- p + ggplot2::scale_size_continuous(range = c(dot.size, dot.size), guide = "none")
  }

  p <- p +
    ggplot2::labs(x = x_lab, y = y_lab, title = title.name) +
    ggplot2::theme_classic(base_size = font.size) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(
        size = font.size.title, face = "bold", hjust = 0.5
      ),
      axis.text.x      = ggplot2::element_text(
        angle = angle.x, vjust = vjust.x, hjust = hjust.x,
        size  = font.size * 0.9
      ),
      axis.text.y      = ggplot2::element_text(size = font.size * 0.9),
      axis.title       = ggplot2::element_text(size = font.size),
      panel.border     = ggplot2::element_rect(
        color = "grey60", fill = NA, linewidth = 0.4
      ),
      panel.grid.major = ggplot2::element_line(color = "grey92", linewidth = 0.35),
      panel.grid.minor = ggplot2::element_blank(),
      legend.key.size  = ggplot2::unit(0.45, "cm"),
      legend.text      = ggplot2::element_text(size = font.size * 0.85),
      legend.title     = ggplot2::element_text(size = font.size * 0.9, face = "bold")
    )

  if (!show.legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  # ---- 13b. Relabel pathway axis with ligand-receptor pairs -----------------
  # lr_map is named by pathway; ggplot matches names against the discrete
  # break (factor level) values, so the underlying data is left untouched.
  if (show.lr.pairs) {
    if (do.flip) {
      p <- p + ggplot2::scale_x_discrete(labels = lr_map)
    } else {
      p <- p + ggplot2::scale_y_discrete(labels = lr_map)
    }
  }

  # ---- 14. Return ------------------------------------------------------------
  if (return.data) {
    return(list(
      communication = df[, c("source", "target", "pathway", "prob", "pair")],
      gg.obj        = p
    ))
  }

  p
}


# ------------------------------------------------------------------------------
# Internal helpers
# These are not exported and generate no .Rd file.
# ------------------------------------------------------------------------------

#' Resolve a sources.use / targets.use argument to a character vector of names
#'
#' @param use        NULL, integer indices, or character names.
#' @param all_groups Character vector of all available group names.
#' @param arg_name   String used in error/warning messages.
#'
#' @return Character vector of resolved group names.
#' @noRd
.resolve_groups <- function(use, all_groups, arg_name) {
  if (is.null(use)) return(all_groups)

  if (is.numeric(use)) {
    idx <- as.integer(use)
    oob <- idx[idx < 1L | idx > length(all_groups)]
    if (length(oob)) {
      warning(
        arg_name, ": index(es) out of range and will be ignored: ",
        paste(oob, collapse = ", ")
      )
      idx <- idx[idx >= 1L & idx <= length(all_groups)]
    }
    return(all_groups[idx])
  }

  if (is.character(use)) {
    unknown <- setdiff(use, all_groups)
    if (length(unknown)) {
      warning(
        arg_name, ": group(s) not found in the object and will be ignored: ",
        paste(unknown, collapse = ", ")
      )
    }
    found <- intersect(use, all_groups)
    if (length(found) == 0L) {
      stop(arg_name, ": none of the requested groups exist in the CellChat object.")
    }
    return(found)
  }

  stop(arg_name, " must be NULL, a character vector, or an integer index vector.")
}


#' Build a CellChat-style continuous colour ramp
#'
#' Mirrors the palette logic of \code{CellChat::netVisual_bubble}: resolve a
#' base palette from an \code{RColorBrewer} name, falling back to a
#' \pkg{viridis} option, reverse it (\code{direction = -1}) and interpolate to
#' \code{n.out} colours.
#'
#' @param color.heatmap Palette name (Brewer or viridis option).
#' @param n.colors      Base number of palette colours (default 10).
#' @param direction     \code{-1} reverses the palette (default, as CellChat).
#' @param n.out         Number of interpolated colours to return (default 99).
#'
#' @return Character vector of \code{n.out} hex colours.
#' @noRd
.heatmap_colors <- function(color.heatmap, n.colors = 10L, direction = -1L, n.out = 99L) {
  cols <- tryCatch(
    suppressWarnings(RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)),
    error = function(e) scales::viridis_pal(option = color.heatmap, direction = -1)(n.colors)
  )
  if (direction == -1L) cols <- rev(cols)   # matches CellChat (incl. its viridis double-rev)
  grDevices::colorRampPalette(cols)(n.out)
}


#' Collapse a set of ligand-receptor pairs into a compact axis label
#'
#' Given the \code{interaction_name_2} strings for one pathway (each formatted
#' \code{"ligand - receptor"}, where a multi-subunit receptor looks like
#' \code{"(A+B)"}), ligands that share an identical receptor \emph{set} are
#' grouped and joined with \code{"/"}; one line is emitted per distinct
#' receptor set.  Shared cases stay on a single line
#' (\code{"TGFB1/TGFB2/TGFB3 - (TGFBR1+TGFBR2)"}); genuinely many-to-many
#' pathways span multiple newline-separated lines.
#'
#' @param inter2 Character vector of \code{interaction_name_2} strings.
#'
#' @return A single character string (possibly containing \code{"\n"}).
#' @noRd
.collapse_lr_pairs <- function(inter2) {
  inter2 <- unique(inter2)

  # Split "ligand - receptor" on the ' - ' delimiter.  \s+-\s+ tolerates the
  # stray double space seen in some mouse rows and does NOT split hyphenated
  # gene symbols such as 'H2-K1' (which carry no spaces around the hyphen).
  lig <- sub("\\s+-\\s+.*$", "", inter2)               # before first ' - '
  rec <- sub("^.*?\\s+-\\s+", "", inter2, perl = TRUE) # after  first ' - '

  ligs       <- unique(lig)
  rec_by_lig <- lapply(ligs, function(L) unique(rec[lig == L]))
  keys       <- vapply(
    rec_by_lig,
    function(r) paste(sort(r), collapse = ""),   # order-independent set key
    character(1)
  )

  entries <- vapply(unique(keys), function(k) {
    idx <- which(keys == k)
    paste(
      paste(ligs[idx],          collapse = "/"),
      paste(rec_by_lig[[idx[1]]], collapse = "/"),
      sep = " - "
    )
  }, character(1))

  paste(entries, collapse = "\n")
}


#' Map pathway names to collapsed ligand-receptor-pair labels
#'
#' Looks up each pathway's contributing LR pairs in \code{object@@LR$LRsig} and
#' collapses them via \code{\link{.collapse_lr_pairs}}.
#'
#' @param object   A CellChat object with a populated \code{@@LR$LRsig} slot.
#' @param pathways Character vector of pathway names.
#'
#' @return Named character vector (names = \code{pathways}) of labels.  A
#'   pathway with no \code{LRsig} rows falls back to its own name.
#' @noRd
.pathway_lr_labels <- function(object, pathways) {
  if (!methods::.hasSlot(object, "LR") || is.null(object@LR$LRsig)) {
    stop(
      "show.lr.pairs = TRUE requires 'object@LR$LRsig'; re-run the CellChat ",
      "pipeline through identifyOverExpressedInteractions() to populate it."
    )
  }
  lr <- object@LR$LRsig
  if (!all(c("pathway_name", "interaction_name_2") %in% colnames(lr))) {
    stop(
      "'object@LR$LRsig' lacks the required 'pathway_name' / ",
      "'interaction_name_2' columns."
    )
  }

  vapply(pathways, function(pw) {
    v <- lr$interaction_name_2[lr$pathway_name == pw]
    if (length(v) == 0L) pw else .collapse_lr_pairs(v)
  }, character(1), USE.NAMES = TRUE)
}
