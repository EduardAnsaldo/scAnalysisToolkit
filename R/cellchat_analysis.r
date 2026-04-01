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
           trim = 0.1,
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
#' @param color.use Controls bubble fill colour.  Two modes:
#'   \describe{
#'     \item{\code{NULL} (default)}{Continuous viridis-family gradient
#'       controlled by \code{viridis.option}.}
#'     \item{Named character vector}{Maps pathway names to hex colours for
#'       categorical colouring, e.g.
#'       \code{c(CXCL = "#e41a1c", VEGF = "#377eb8")}.  Pathways absent from
#'       the vector receive \code{NA} fill.}
#'   }
#' @param viridis.option Character string passed to the \code{option} argument
#'   of \code{\link[ggplot2]{scale_fill_viridis_c}}.  Accepted values:
#'   \code{"viridis"} (default), \code{"magma"}, \code{"plasma"},
#'   \code{"inferno"}, \code{"cividis"}, \code{"mako"}, \code{"rocket"},
#'   \code{"turbo"}.  Ignored when \code{color.use} is a named vector.
#' @param show.legend Logical; whether to draw the legend (default \code{TRUE}).
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
#'   top.n          = 15,
#'   top.n.by       = "mean",
#'   viridis.option = "magma"
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
#' @importFrom ggplot2 ggplot aes geom_point
#'   scale_fill_manual scale_color_manual
#'   scale_fill_viridis_c scale_color_viridis_c
#'   scale_size_continuous labs
#'   theme_classic theme element_text element_rect element_line element_blank
#'   unit
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
    color.use       = NULL,
    viridis.option  = "viridis",
    show.legend     = TRUE,
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

  viridis_opts <- c("viridis", "magma", "plasma", "inferno",
                    "cividis", "mako", "rocket", "turbo")
  if (!viridis.option %in% viridis_opts) {
    stop(
      "'viridis.option' must be one of: ",
      paste(viridis_opts, collapse = ", ")
    )
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
  df <- .array3d_to_long(prob3d)   # columns: source, target, pathway, prob

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

  # ---- 11. Axis assignment (respects do.flip) --------------------------------
  if (do.flip) {
    x_var <- "pathway"; y_var <- "pair"
    x_lab <- "Signalling Pathway"; y_lab <- "Cell-group pair"
  } else {
    x_var <- "pair"; y_var <- "pathway"
    x_lab <- "Cell-group pair"; y_lab <- "Signalling Pathway"
  }

  # ---- 12. Resolve label-justification defaults -----------------------------
  if (is.null(vjust.x)) vjust.x <- if (angle.x != 0) 1 else 0.5
  if (is.null(hjust.x)) hjust.x <- if (angle.x != 0) 1 else 0.5

  # ---- 13. Build ggplot ------------------------------------------------------
  prob_max <- max(df$prob, na.rm = TRUE)

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

    # Continuous viridis gradient -------------------------------------------
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
      ggplot2::scale_fill_viridis_c(
        option = viridis.option,
        name   = "Probability",
        limits = c(0, prob_max),
        begin  = 0.15,    # avoid very dark end on some palettes
        end    = 0.95
      ) +
      ggplot2::scale_color_viridis_c(
        option = viridis.option,
        guide  = "none",
        limits = c(0, prob_max),
        begin  = 0.15,
        end    = 0.95
      )
  }

  p <- p +
    ggplot2::scale_size_continuous(range = c(1, 10), name = "Probability") +
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


#' Melt a named 3-D array into a long data frame
#'
#' @param arr A named 3-D array with dimensions
#'   \code{[sender, receiver, pathway]}.
#'
#' @return A data frame with columns \code{source}, \code{target},
#'   \code{pathway}, and \code{prob}.
#' @noRd
.array3d_to_long <- function(arr) {
  pw_names <- dimnames(arr)[[3]]
  rows     <- vector("list", length(pw_names))

  for (k in seq_along(pw_names)) {
    slice            <- arr[, , k, drop = TRUE]
    mat_df           <- as.data.frame(as.table(slice), stringsAsFactors = FALSE)
    colnames(mat_df) <- c("source", "target", "prob")
    mat_df$pathway   <- pw_names[k]
    rows[[k]]        <- mat_df
  }

  out      <- do.call(rbind, rows)
  out$prob <- as.numeric(out$prob)
  out
}
