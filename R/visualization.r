#' Apply consistent integer-spaced legend breaks to a FeaturePlot_scCustom plot
#'
#' Modifies the color scale of a `FeaturePlot_scCustom()` output so that legend
#' tick marks are evenly spaced at a fixed step size. The plot's existing color
#' palette, limits, and `na_cutoff` (encoded as the lower limit) are preserved;
#' only the `breaks` argument is changed.
#'
#' @param plot A ggplot or patchwork object returned by `FeaturePlot_scCustom()`.
#' @param break_size Numeric step size between legend ticks. Defaults to `1`.
#' @param colors_use Color palette to use for the gradient. Defaults to
#'   `viridis_plasma_dark_high` to match `FeaturePlot_scCustom`'s default.
#' @param na_color Color for points below the lower limit (i.e. below
#'   `na_cutoff`). Defaults to `"lightgray"` to match `FeaturePlot_scCustom`'s
#'   default.
#'
#' @return A plot of the same type as the input, with updated legend breaks.
#'
#' @importFrom ggplot2 scale_color_gradientn
#' @importFrom scales breaks_width
#'
#' @export
apply_integer_breaks <- function(plot,
                                 break_size = 1,
                                 colors_use = viridis_plasma_dark_high,
                                 na_color = "lightgray") {

    # Helper: pull existing color-scale limits from a single ggplot
    get_color_limits <- function(p) {
        scales <- p$scales$scales
        for (s in scales) {
            if (any(c("colour", "color") %in% s$aesthetics)) {
                if (!is.null(s$limits)) return(s$limits)
            }
        }
        NULL
    }

    # Helper: rebuild scale for one ggplot using its own limits
    rescale_one <- function(p) {
        lims <- get_color_limits(p)
        
        suppressMessages(
            p + ggplot2::scale_color_gradientn(
                colors   = colors_use,
                na.value = na_color,
                limits   = lims,
                breaks   = scales::breaks_width(break_size)
            )
        )
    }

    # List output (combine = FALSE)
    if (is.list(plot) && !inherits(plot, "ggplot") && !inherits(plot, "patchwork")) {
        return(lapply(plot, rescale_one))
    }

    # patchwork: apply per-subplot to preserve each one's own limits
    if (inherits(plot, "patchwork")) {
        plot$patches$plots <- lapply(plot$patches$plots, rescale_one)
        # The "main" plot of a patchwork is also a ggplot — rescale it too
        plot <- rescale_one(plot)
        return(plot)
    }

    # Single ggplot
    rescale_one(plot)
}



#' Create a patchwork of genomic tracks for multiple regions/genes
#'
#' @param seurat A Seurat object with ATAC data
#' @param regions Character vector of regions or gene names to plot.
#'   Regions can be given as coordinate strings ("chr1-1000-2000") or as
#'   gene symbols; gene symbols are resolved to coordinates using the
#'   annotation stored in the object's ATAC assay.
#' @param theme A ggplot2 theme to apply to all plots
#' @param nrow Number of rows in the patchwork grid
#' @param ncol Number of columns in the patchwork grid
#' @param title Optional shared title for the full patchwork
#' @param subtitle Optional shared subtitle
#' @param title_theme Optional theme element for the title
#' @param assay_atac Name of the ATAC/peaks assay (default "ATAC")
#' @param extend_upstream Bases to extend upstream (default 20000)
#' @param extend_downstream Bases to extend downstream (default 20000)
#' @param add_links Logical, whether to add link plots (default FALSE)
#' @param gene_symbol Logical, whether to add expression plots (default FALSE)
#' @param assay Assay to use for expression plot (default "SCT")
#' @param heights Numeric vector of relative heights for tracks
#' @param widths Numeric vector of relative widths (track vs expression)
#'
#' @return A patchwork object
create_tracks_patchwork <- function(seurat,
                                    regions,
                                    theme = NULL,
                                    nrow = NULL,
                                    ncol = NULL,
                                    title = NULL,
                                    subtitle = NULL,
                                    title_theme = NULL,
                                    assay_atac = "ATAC",
                                    extend_upstream = 20000,
                                    extend_downstream = 20000,
                                    add_links = FALSE,
                                    gene_symbol = FALSE,
                                    expression_features = NULL,
                                    assay = "SCT",
                                    heights = c(10, 1, 2, 3),
                                    widths = c(10, 1)) {

  # Determine grid dimensions if not supplied
  n <- length(regions)
  if (is.null(ncol) && is.null(nrow)) {
    ncol <- ceiling(sqrt(n))
    nrow <- ceiling(n / ncol)
  } else if (is.null(ncol)) {
    ncol <- ceiling(n / nrow)
  } else if (is.null(nrow)) {
    nrow <- ceiling(n / ncol)
  }

  # Helper: strip y-axis ticks and labels from a ggplot
  strip_y_ticks_labels <- function(p) {
    if (is.null(p)) return(NULL)
    p + theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank()
    )
  }

  # Helper: fully remove every y-axis element (title, text, ticks, line)
  strip_y_axis_full <- function(p) {
    if (is.null(p)) return(NULL)
    p + theme(
      axis.title.y = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y  = element_blank()
    )
  }

  # Helper: remove legends from a ggplot
  strip_legend <- function(p) {
    if (is.null(p)) return(NULL)
    p + theme(legend.position = "none")
  }

  # Helper: rotate facet strip labels 90° and remove the strip background box.
  # Applied AFTER the user theme so it isn't overridden.
  style_facet_labels <- function(p) {
    if (is.null(p)) return(NULL)
    p + theme(
      strip.text.y        = element_text(angle = 0),
      strip.text.y.left   = element_text(angle = 0),
      strip.text.y.right  = element_text(angle = 0),
      strip.text.x        = element_text(angle = 270, hjust = 0, vjust = 0.5),
      strip.text.x.top    = element_text(angle = 270, hjust = 0, vjust = 0.5),
      strip.text.x.bottom = element_text(angle = 270, hjust = 0, vjust = 0.5),
      strip.background    = element_blank()
    )
  }
  # Helper: remove x and y axis titles, text, and ticks from the expression plot
  strip_expr_axes <- function(p) {
    if (is.null(p)) return(NULL)
    p + theme(
      # axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x.bottom  = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank()
    )
  }
  # Helper: remove facet strip labels entirely (for non-first-column coverage
  # plots and for expression plots). Targets every strip position explicitly
  # so it works regardless of which side the facet strip is on.
  strip_facet_labels <- function(p) {
    if (is.null(p)) return(NULL)
    p + theme(
      strip.text          = element_blank(),
      strip.text.x        = element_blank(),
      strip.text.x.top    = element_blank(),
      strip.text.x.bottom = element_blank(),
      strip.text.y        = element_blank(),
      strip.text.y.left   = element_blank(),
      strip.text.y.right  = element_blank(),
      strip.background    = element_blank()
    )
  }

  # Helper: resolve a gene symbol to a coordinate string.
  resolve_region <- function(region) {
    if (grepl("^[^:\\-]+[:\\-][0-9]+[:\\-][0-9]+$", region)) {
      return(region)
    }
    annotations <- Annotation(object = seurat[[assay_atac]])
    if (is.null(annotations)) {
      stop("No annotation found on assay '", assay_atac,
           "'. Cannot resolve gene symbol '", region, "'.")
    }
    gr <- annotations[annotations$gene_name == region]
    if (length(gr) == 0) {
      stop("Gene '", region, "' not found in annotations of assay '",
           assay_atac, "'.")
    }
    chr   <- as.character(GenomicRanges::seqnames(gr))[1]
    start <- min(GenomicRanges::start(gr))
    end   <- max(GenomicRanges::end(gr))
    paste(chr, start, end, sep = "-")
  }
  

  # Build one combined track plot for a single region
  build_combined_plot <- function(region_input, col_idx, ncol_total, expr_feature = NULL) {

    is_first_col <- col_idx == 1
    coord_region <- resolve_region(region_input)

    cov_plot <- CoveragePlot(
      object = seurat,
      region = coord_region,
      extend.upstream = extend_upstream,
      extend.downstream = extend_downstream,
      annotation = FALSE,
      links = FALSE,
      peaks = FALSE
    ) + theme(axis.title.y = element_text(angle = 90, vjust = 0.5))

    gene_plot <- AnnotationPlot_custom(
      object = seurat,
      region = coord_region,
      extend.upstream = extend_upstream,
      extend.downstream = extend_downstream
    ) + theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

    peak_plot <- PeakPlot(
      object = seurat,
      region = coord_region,
      extend.upstream = extend_upstream,
      extend.downstream = extend_downstream
    ) + theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

    link_plot <- if (add_links) {
    LinkPlot(
      object            = seurat,
      region            = coord_region,
      min.cutoff        = 0,
      sep               = c("-", "-"),
      extend.upstream   = extend_upstream,
      extend.downstream = extend_downstream,
      scale.linewidth   = FALSE,
    ) + scale_color_gradientn(colors = 'orange3') +
       theme(axis.title.y = element_text(angle = 0, vjust = 0.5))
  } else {
    NULL
  }

    expr_plot <- if (gene_symbol) {
      feat <- if (!is.null(expr_feature)) expr_feature else region_input
      avail <- rownames(seurat[[assay]])
      if (feat %in% avail) {
        ExpressionPlot(
          object   = seurat,
          features = feat,
          assay    = assay
        )
      } else {
        warning("Feature '", feat, "' not found in assay '", assay,
                "' — skipping expression plot.")
        NULL
      }
    } else {
      NULL
    }

    # Apply the user-supplied theme FIRST so our later theme() calls
    # (especially strip label rotation) take precedence
    if (!is.null(theme)) {
      cov_plot  <- cov_plot  + theme
      gene_plot <- gene_plot + theme
      peak_plot <- peak_plot + theme
      if (!is.null(link_plot)) link_plot <- link_plot + theme
      if (!is.null(expr_plot)) expr_plot <- expr_plot + theme
    }

    # Handle y-axis handling based on column position
    if (is_first_col) {
      # First column: keep axis, blank ticks/text, override coverage title
      cov_plot  <- strip_y_ticks_labels(cov_plot)
      gene_plot <- strip_y_ticks_labels(gene_plot)
      peak_plot <- strip_y_ticks_labels(peak_plot)
      link_plot <- strip_y_ticks_labels(link_plot)
      expr_plot <- strip_y_ticks_labels(expr_plot)

      cov_plot <- cov_plot + ylab("Normalized signal")
    } else {
      # Non-first column: remove everything y-axis related
      cov_plot  <- strip_y_axis_full(cov_plot)
      gene_plot <- strip_y_axis_full(gene_plot)
      peak_plot <- strip_y_axis_full(peak_plot)
      link_plot <- strip_y_axis_full(link_plot)
      expr_plot <- strip_y_axis_full(expr_plot)
    }

    # Remove legends from every track
    cov_plot  <- strip_legend(cov_plot)
    gene_plot <- strip_legend(gene_plot)
    peak_plot <- strip_legend(peak_plot)
    link_plot <- strip_legend(link_plot)
    expr_plot <- strip_legend(expr_plot)

    # Facet labels: style for first column, remove for others — applied
    # LAST so the user-supplied theme can't reset strip.text back
    if (is_first_col) {
      cov_plot <- style_facet_labels(cov_plot)
    } else {
      cov_plot <- strip_facet_labels(cov_plot)
    }
    expr_plot <- strip_facet_labels(expr_plot)

    # Remove axis titles, text, and ticks from the expression plot
    expr_plot <- strip_expr_axes(expr_plot)

    # Assemble plotlist, dropping NULLs
    plotlist <- list(cov_plot, peak_plot, gene_plot, link_plot)
    keep <- !vapply(plotlist, is.null, logical(1))
    plotlist <- plotlist[keep]
    track_heights <- heights[keep]


    combined <- CombineTracks(
      plotlist        = plotlist,
      expression.plot = expr_plot,
      heights         = track_heights,
      widths          = widths
    )

    # Adjust left/right margins: zero out interior boundaries
    is_leftmost  <- col_idx == 1
    is_rightmost <- col_idx == ncol_total
    l_margin <- if (is_leftmost)  1 else 0
    r_margin <- if (is_rightmost) 1 else 0
    combined <- combined & theme(plot.margin = margin(1, r_margin, 1, l_margin))

    combined
  }

  # Build one combined plot per region, flagging first-column positions
  combined_plots <- lapply(seq_along(regions), function(i) {
    col_idx <- ((i - 1) %% ncol) + 1
    ef <- if (!is.null(expression_features)) expression_features[[i]] else NULL
    build_combined_plot(regions[[i]], col_idx = col_idx, ncol_total = ncol,
                        expr_feature = ef)
  })

  # Gray dashed vertical separator between columns (zero margins)
  sep_line <- wrap_elements(grid::linesGrob(
    x  = grid::unit(c(0.5, 0.5), "npc"),
    y  = grid::unit(c(0, 1), "npc"),
    gp = grid::gpar(col = "gray", lty = "dashed", lwd = 1)
  )) & theme(plot.margin = margin(0, 0, 0, 0))
  sep_width <- 0.005            # very narrow relative width

  # Build each row, interleaving separator lines between columns
  row_patches <- vector("list", nrow)
  for (r in seq_len(nrow)) {
    row_items  <- list()
    row_widths <- c()
    for (c in seq_len(ncol)) {
      idx <- (r - 1) * ncol + c
      if (idx > n) break
      if (c > 1) {
        row_items  <- c(row_items, list(sep_line))
        row_widths <- c(row_widths, sep_width)
      }
      row_items  <- c(row_items, list(combined_plots[[idx]]))
      row_widths <- c(row_widths, 1)
    }
    row_patches[[r]] <- wrap_plots(row_items, nrow = 1) +
      plot_layout(widths = row_widths)
  }

  # Stack rows vertically
  final <- wrap_plots(row_patches, ncol = 1)

  # Add a shared title / subtitle if requested
  if (!is.null(title) || !is.null(subtitle)) {
    default_title_theme <- theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    )
    final <- final + plot_annotation(
      title    = title,
      subtitle = subtitle,
      theme    = if (is.null(title_theme)) default_title_theme else title_theme
    )
  }

  return(final)
}


#' Plot gene annotations
#'
#' Display gene annotations in a given region of the genome.
#'
#' @param object A [SeuratObject::Seurat()] object
#' @param region A genomic region to plot
#' @param assay Name of assay to use. If NULL, use the default assay.
#' @param mode Display mode. Choose either "gene" or "transcript" to determine
#' whether genes or transcripts are plotted.
#' @param sep Separators to use for strings encoding genomic coordinates. First
#' element is used to separate the chromosome from the coordinates, second
#' element is used to separate the start from end coordinate.
#' @param extend.upstream Number of bases to extend the region upstream.
#' @param extend.downstream Number of bases to extend the region downstream.
#' 
#' @return Returns a [ggplot2::ggplot()] object
#' @export
#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomicRanges start end
#' @importFrom GenomeInfoDb seqnames
#' @importFrom ggplot2 theme_classic ylim xlim ylab xlab
#' geom_segment geom_text aes_string scale_color_manual
#' @importFrom grid arrow
#' @importFrom S4Vectors split
#' @importFrom fastmatch fmatch
#' @concept visualization
#' @examples
#' \donttest{
#' AnnotationPlot_custom(object = atac_small, region = c("chr1-29554-39554"))
#' }
AnnotationPlot_custom <- function(
  object,
  region,
  assay = NULL,
  mode = "gene",
  sep = c("-", "-"),
  extend.upstream = 0,
  extend.downstream = 0,
  annotation_font_size = 1.5
) {
  if(mode == "gene") {
    collapse_transcript <- TRUE
    label <- "gene_name"
  } else if (mode == "transcript") {
    collapse_transcript <- FALSE
    label <- "tx_id"
  } else {
    stop("Unknown mode requested, choose either 'gene' or 'transcript'")
  }
  annotation <- Annotation(object = object)
  if (is.null(x = annotation)) {
    return(NULL)
  }
  region <- Signac:::FindRegion(
    object = object,
    region = region,
    sep = sep,
    assay = assay,
    extend.upstream = extend.upstream,
    extend.downstream = extend.downstream
  )
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  chromosome <- seqnames(x = region)

  # get names of genes that overlap region, then subset to include only those
  # genes. This avoids truncating the gene if it runs outside the region
  annotation.subset <- subsetByOverlaps(x = annotation, ranges = region)
  if (mode == "gene") {
    genes.keep <- unique(x = annotation.subset$gene_name)
    annotation.subset <- annotation[
      fastmatch::fmatch(x = annotation$gene_name, table = genes.keep, nomatch = 0L) > 0L
    ]
  } else {
    tx.keep <- unique(x = annotation.subset$tx_id)
    annotation.subset <- annotation[
      fastmatch::fmatch(x = annotation$tx_id, table = tx.keep, nomatch = 0L) > 0L
    ]
  }

  if (length(x = annotation.subset) == 0) {
    # make empty plot
    p <- ggplot(data = data.frame())
    y_limit <- c(0, 1)
  } else {
    annotation_df_list <- Signac:::reformat_annotations(
      annotation = annotation.subset,
      start.pos = start.pos,
      end.pos = end.pos,
      collapse_transcript = collapse_transcript
    )
    p <- ggplot() +
      # exons
      geom_segment(
        data = annotation_df_list$exons,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$exons$dodge,
          xend = "end",
          yend = annotation_df_list$exons$dodge,
          color = "strand"
        ),
        show.legend = FALSE,
        size = 1
      ) +
      # gene body
      geom_segment(
        data = annotation_df_list$labels,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$labels$dodge,
          xend = "end",
          yend = annotation_df_list$labels$dodge,
          color = "strand"
        ),
        show.legend = FALSE,
        size = 1/5
      )
    if (nrow(x = annotation_df_list$plus) > 0) {
      # forward strand arrows
      p <- p + geom_segment(
        data = annotation_df_list$plus,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$plus$dodge,
          xend = "end",
          yend = annotation_df_list$plus$dodge,
          color = "strand"
        ),
        arrow = arrow(
          ends = "last",
          type = "open",
          angle = 45,
          length = unit(x = 0.01, units = "inches")
        ),
        show.legend = FALSE,
        size = 1/6
      )
    }
    if (nrow(x = annotation_df_list$minus) > 0) {
      # reverse strand arrows
      p <- p + geom_segment(
        data = annotation_df_list$minus,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$minus$dodge,
          xend = "end",
          yend = annotation_df_list$minus$dodge,
          color = "strand"
        ),
        arrow = arrow(
          ends = "first",
          type = "open",
          angle = 45,
          length = unit(x = 0.01, units = "inches")
        ),
        show.legend = FALSE,
        size = 1/6
      )
    }
    # label genes
    n_stack <- max(annotation_df_list$labels$dodge)
    annotation_df_list$labels$dodge <- annotation_df_list$labels$dodge + 0.2
    p <- p + geom_text(
      data = annotation_df_list$labels,
      mapping = aes_string(x = "position", y = "dodge", label = label),
      size = annotation_font_size 
    )
    y_limit <- c(0.9, n_stack + 0.4)
  }
  p <- p +
    theme_classic() +
    ylab("Genes") +
    xlab(label = paste0(chromosome, " position (bp)")) +
    xlim(start.pos, end.pos) +
    ylim(y_limit) +
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    ) +
    scale_color_manual(values = c("darkblue", "darkgreen"))
  return(p)
}

#' Custom CoveragePlot with stacked annotation, peak, and link tracks
#'
#' Builds a publication-ready Signac coverage plot for a single genomic region,
#' combining a coverage track with peak, gene annotation, and (optionally)
#' peak-to-gene link tracks, plus an optional expression panel on the right.
#' Y-axis ticks and labels are stripped from all tracks for a clean look,
#' facet strip labels are rotated, and legends are removed. The user can pass
#' an additional `ggplot2` theme that is applied before the function's own
#' styling adjustments so those adjustments always win.
#'
#' @param seurat A `Seurat` object containing a `ChromatinAssay` (typically
#'   named `"ATAC"`) with fragment information and gene annotations attached.
#'   For `gene_symbol = TRUE`, the object must also contain an expression
#'   assay (e.g. `"SCT"` or `"RNA"`) with the feature of interest.
#' @param region Either a coordinate string (e.g. `"chr1-1000-2000"` or
#'   `"chr1:1000-2000"`) or a gene symbol present in the chromatin assay's
#'   annotations. Gene symbols are resolved to their full transcript span via
#'   [resolve_region()].
#' @param theme Optional `ggplot2` theme object applied to every track before
#'   the function's internal styling. Use this to set base font sizes, panel
#'   backgrounds, etc.
#' @param title Optional plot title shown above the assembled tracks.
#' @param subtitle Optional plot subtitle shown beneath the title.
#' @param title_theme Optional `ggplot2` theme controlling the appearance of
#'   the title and subtitle. If `NULL`, a sensible default is used (centered,
#'   bold title at size 16; centered subtitle at size 12).
#' @param assay_atac Name of the chromatin assay in `seurat`. Defaults to
#'   `"ATAC"`. Used both for plotting and for resolving gene symbols against
#'   the assay's annotations.
#' @param extend_upstream Integer. Number of bases to extend the plotted
#'   region upstream of `region`. Defaults to `20000`.
#' @param extend_downstream Integer. Number of bases to extend the plotted
#'   region downstream of `region`. Defaults to `20000`.
#' @param add_links Logical. If `TRUE`, a `LinkPlot` track is added showing
#'   peak-to-gene links stored in `Links(seurat[[assay_atac]])`. Defaults to
#'   `FALSE`.
#' @param gene_symbol Logical. If `TRUE`, an `ExpressionPlot` panel is
#'   appended to the right of the tracks for the feature given by
#'   `expression_feature` (or `region` itself if `expression_feature` is
#'   `NULL`). Defaults to `FALSE`.
#' @param expression_feature Optional character string giving the feature
#'   (typically a gene symbol) to plot in the expression panel. Only used
#'   when `gene_symbol = TRUE`. If `NULL`, the value of `region` is used. A
#'   warning is issued and the panel is dropped if the feature is not found
#'   in `assay`.
#' @param assay Name of the expression assay used for the expression panel.
#'   Defaults to `"SCT"`.
#' @param heights Numeric vector of relative heights for the tracks in this
#'   order: coverage, peak, gene annotation, link. Length 4 by default
#'   (`c(10, 1, 2, 3)`); entries corresponding to omitted tracks (e.g. the
#'   link track when `add_links = FALSE`) are dropped automatically.
#' @param widths Numeric vector of length 2 giving the relative widths of
#'   the track panel and the expression panel, passed to [CombineTracks()].
#'   Defaults to `c(10, 1)`.
#'
#' @return A `patchwork` object combining the coverage and (optionally)
#'   expression panels. Can be printed directly, further composed with other
#'   `patchwork` objects, or saved with [ggplot2::ggsave()].
#'
#' @details
#' Tracks are assembled with [Signac::CombineTracks()] in this order from
#' top to bottom: coverage, peak, gene annotation, link. Any of these can be
#' omitted (links via `add_links = FALSE`; the others are always present).
#'
#' Styling is applied in three passes so user customizations and internal
#' adjustments compose predictably:
#' \enumerate{
#'   \item The user-supplied `theme` (if any) is added to every track.
#'   \item Y-axis ticks/labels are blanked, the coverage y-label is set to
#'         `"Normalized signal"`, and legends are removed.
#'   \item Facet strip labels are rotated on the coverage plot and removed
#'         from the expression panel.
#' }
#'
#' This ordering ensures the function's facet-label rotation and y-axis
#' stripping cannot be overridden by a user theme that resets `strip.text`
#' or `axis.text.y`.
#'
#' @section Required helpers:
#' This function depends on `AnnotationPlot_custom()` (a custom variant of
#' [Signac::AnnotationPlot()]) and a `sequential_palette_dotplot` palette
#' object available in the calling environment.
#'
#' @seealso
#'   [Signac::CoveragePlot()], [Signac::PeakPlot()], [Signac::LinkPlot()],
#'   [Signac::ExpressionPlot()], [Signac::CombineTracks()]
#'
#' @examples
#' \dontrun{
#' # Basic usage with a gene symbol
#' CoveragePlot_custom(seurat_obj, region = "Lypla1")
#'
#' # Coordinate region with link track and expression panel
#' CoveragePlot_custom(
#'   seurat_obj,
#'   region            = "chr1-4768107-4769042",
#'   add_links         = TRUE,
#'   gene_symbol       = TRUE,
#'   expression_feature = "Lypla1",
#'   extend_upstream   = 50000,
#'   extend_downstream = 50000,
#'   title             = "Lypla1 locus",
#'   subtitle          = "scATAC + scRNA"
#' )
#'
#' # Apply a custom base theme
#' CoveragePlot_custom(
#'   seurat_obj,
#'   region = "Tcea1",
#'   theme  = theme_classic(base_size = 11)
#' )
#' }
#'
#' @export
CoveragePlot_custom <- function(seurat,
                                region,
                                theme = NULL,
                                title = NULL,
                                subtitle = NULL,
                                title_theme = NULL,
                                assay_atac = "ATAC",
                                extend_upstream = 20000,
                                extend_downstream = 20000,
                                add_links = FALSE,
                                gene_symbol = FALSE,
                                expression_feature = NULL,
                                assay = "SCT",
                                heights = c(10, 1, 2, 3),
                                widths = c(10, 1)) {

  # Helper: strip y-axis ticks and labels from a ggplot
  strip_y_ticks_labels <- function(p) {
    if (is.null(p)) return(NULL)
    p + theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank()
    )
  }

  # Helper: remove legends from a ggplot
  strip_legend <- function(p) {
    if (is.null(p)) return(NULL)
    p + theme(legend.position = "none")
  }

  # Helper: rotate facet strip labels and remove the strip background box.
  # Applied AFTER the user theme so it isn't overridden.
  style_facet_labels <- function(p) {
    if (is.null(p)) return(NULL)
    p + theme(
      strip.text.y        = element_text(angle = 0),
      strip.text.y.left   = element_text(angle = 0),
      strip.text.y.right  = element_text(angle = 0),
      strip.text.x        = element_text(angle = 270, hjust = 0, vjust = 0.5),
      strip.text.x.top    = element_text(angle = 270, hjust = 0, vjust = 0.5),
      strip.text.x.bottom = element_text(angle = 270, hjust = 0, vjust = 0.5),
      strip.background    = element_blank()
    )
  }

  # Helper: remove x and y axis titles, text, and ticks from the expression plot
  strip_expr_axes <- function(p) {
    if (is.null(p)) return(NULL)
    p + theme(
      axis.title.y       = element_blank(),
      axis.text.x.bottom = element_blank(),
      axis.text.y        = element_blank(),
      axis.ticks.x       = element_blank(),
      axis.ticks.y       = element_blank()
    )
  }

  # Helper: remove facet strip labels entirely (used on expression plot)
  strip_facet_labels <- function(p) {
    if (is.null(p)) return(NULL)
    p + theme(
      strip.text          = element_blank(),
      strip.text.x        = element_blank(),
      strip.text.x.top    = element_blank(),
      strip.text.x.bottom = element_blank(),
      strip.text.y        = element_blank(),
      strip.text.y.left   = element_blank(),
      strip.text.y.right  = element_blank(),
      strip.background    = element_blank()
    )
  }

  # Helper: resolve a gene symbol to a coordinate string
  resolve_region <- function(region) {
    if (grepl("^[^:\\-]+[:\\-][0-9]+[:\\-][0-9]+$", region)) {
      return(region)
    }
    annotations <- Annotation(object = seurat[[assay_atac]])
    if (is.null(annotations)) {
      stop("No annotation found on assay '", assay_atac,
           "'. Cannot resolve gene symbol '", region, "'.")
    }
    gr <- annotations[annotations$gene_name == region]
    if (length(gr) == 0) {
      stop("Gene '", region, "' not found in annotations of assay '",
           assay_atac, "'.")
    }
    chr   <- as.character(GenomicRanges::seqnames(gr))[1]
    start <- min(GenomicRanges::start(gr))
    end   <- max(GenomicRanges::end(gr))
    paste(chr, start, end, sep = "-")
  }

  coord_region <- resolve_region(region)

  # Build the individual track plots
  cov_plot <- CoveragePlot(
    object            = seurat,
    region            = coord_region,
    extend.upstream   = extend_upstream,
    extend.downstream = extend_downstream,
    annotation        = FALSE,
    links = FALSE,
    peaks             = FALSE
  )

  gene_plot <- AnnotationPlot_custom(
    object            = seurat,
    region            = coord_region,
    extend.upstream   = extend_upstream,
    extend.downstream = extend_downstream
  )

  peak_plot <- PeakPlot(
    object            = seurat,
    region            = coord_region,
    extend.upstream   = extend_upstream,
    extend.downstream = extend_downstream
  )

  link_plot <- if (add_links) {
    LinkPlot(
      object            = seurat,
      region            = coord_region,
      min.cutoff        = 0,
      sep               = c("-", "-"),
      extend.upstream   = extend_upstream,
      extend.downstream = extend_downstream,
      scale.linewidth   = FALSE,
    ) + scale_color_gradientn(colors = 'orange3')
  } else {
    NULL
  }

  expr_plot <- if (gene_symbol) {
    feat <- if (!is.null(expression_feature)) expression_feature else region
    avail <- rownames(seurat[[assay]])
    if (feat %in% avail) {
      ExpressionPlot(
        object   = seurat,
        features = feat,
        assay    = assay
      )
    } else {
      warning("Feature '", feat, "' not found in assay '", assay,
              "' — skipping expression plot.")
      NULL
    }
  } else {
    NULL
  }

  # Apply user-supplied theme FIRST so our later theme() calls take precedence
  if (!is.null(theme)) {
    cov_plot  <- cov_plot  + theme
    gene_plot <- gene_plot + theme
    peak_plot <- peak_plot + theme
    if (!is.null(link_plot)) link_plot <- link_plot + theme
    if (!is.null(expr_plot)) expr_plot <- expr_plot + theme
  }

  # Strip y-axis ticks/text and add coverage y-label
  cov_plot  <- strip_y_ticks_labels(cov_plot)
  gene_plot <- strip_y_ticks_labels(gene_plot)
  peak_plot <- strip_y_ticks_labels(peak_plot)
  link_plot <- strip_y_ticks_labels(link_plot)
  expr_plot <- strip_y_ticks_labels(expr_plot)
  cov_plot  <- cov_plot + ylab("Normalized signal")

  # Remove legends
  cov_plot  <- strip_legend(cov_plot)
  gene_plot <- strip_legend(gene_plot)
  peak_plot <- strip_legend(peak_plot)
  link_plot <- strip_legend(link_plot)
  expr_plot <- strip_legend(expr_plot)

  # Facet label styling (rotated for coverage, removed for expression)
  cov_plot  <- style_facet_labels(cov_plot)
  expr_plot <- strip_facet_labels(expr_plot)
  expr_plot <- strip_expr_axes(expr_plot)

  # Assemble plotlist, dropping NULLs
  plotlist <- list(cov_plot, gene_plot, peak_plot, link_plot)
  keep <- !vapply(plotlist, is.null, logical(1))
  plotlist <- plotlist[keep]
  track_heights <- heights[keep]

  combined <- CombineTracks(
    plotlist        = plotlist,
    expression.plot = expr_plot,
    heights         = track_heights,
    widths          = widths
  )

  # Add title / subtitle if requested
  if (!is.null(title) || !is.null(subtitle)) {
    default_title_theme <- theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    )
    combined <- combined + plot_annotation(
      title    = title,
      subtitle = subtitle,
      theme    = if (is.null(title_theme)) default_title_theme else title_theme
    )
  }

  # combined <- combined & theme()
  return(combined)
}

#' Create a patchwork of genomic tracks for multiple regions/genes
#'
#' @param seurat A Seurat object with ATAC data
#' @param regions Character vector of regions or gene names to plot.
#'   Regions can be given as coordinate strings ("chr1-1000-2000") or as
#'   gene symbols; gene symbols are resolved to coordinates using the
#'   annotation stored in the object's ATAC assay.
#' @param theme A ggplot2 theme to apply to all plots
#' @param nrow Number of rows in the patchwork grid
#' @param ncol Number of columns in the patchwork grid
#' @param title Optional shared title for the full patchwork
#' @param subtitle Optional shared subtitle
#' @param title_theme Optional theme element for the title
#' @param assay_atac Name of the ATAC/peaks assay (default "ATAC")
#' @param extend_upstream Bases to extend upstream (default 20000)
#' @param extend_downstream Bases to extend downstream (default 20000)
#' @param add_links Logical, whether to add link plots (default FALSE)
#' @param gene_symbol Logical, whether to add expression plots (default FALSE)
#' @param expression_features Optional vector of features for the expression
#'   plots, parallel to `regions`. If NULL, the region input itself is used.
#' @param assay Assay to use for expression plot (default "SCT")
#' @param heights Numeric vector of relative heights for tracks
#' @param widths Numeric vector of relative widths (track vs expression)
#' @param add_track_titles Logical; if TRUE, add a per-region title above each
#'   track stack. The title is built from the gene symbol (the value passed in
#'   `regions`) combined with `track_title_suffix`. Default FALSE.
#' @param track_title_suffix Optional string or vector to append after each
#'   gene symbol in the per-track title. Accepts:
#'   \itemize{
#'     \item A single string: applied to all regions
#'     \item A named character vector: looked up by region/gene name
#'     \item An unnamed character vector: applied positionally to `regions`
#'   }
#'   Default NULL (gene symbol only).
#' @param track_titles Optional character vector of explicit titles to use
#'   verbatim, length must equal `length(regions)`. Overrides the
#'   `track_title_suffix` logic. Useful when passing coordinate strings as
#'   regions. Default NULL.
#' @param track_title_theme Optional theme for per-track titles. Default is
#'   centered, italic, size 12.
#'
#' @return A patchwork object
create_tracks_patchwork2 <- function(seurat,
                                    regions,
                                    theme = NULL,
                                    nrow = NULL,
                                    ncol = NULL,
                                    title = NULL,
                                    subtitle = NULL,
                                    title_theme = NULL,
                                    assay_atac = "ATAC",
                                    extend_upstream = 20000,
                                    extend_downstream = 20000,
                                    annotation_font_size = 1.5,
                                    add_links = FALSE,
                                    gene_symbol = FALSE,
                                    expression_features = NULL,
                                    assay = "SCT",
                                    heights = c(10, 1, 2, 3),
                                    widths = c(10, 1),
                                    add_track_titles = FALSE,
                                    track_title_suffix = NULL,
                                    track_titles = NULL,
                                    track_title_theme = NULL) {

  # Validate explicit titles length, if provided
  if (!is.null(track_titles) && length(track_titles) != length(regions)) {
    stop("`track_titles` must have the same length as `regions`.")
  }

  # Determine grid dimensions if not supplied
  n <- length(regions)
  if (is.null(ncol) && is.null(nrow)) {
    ncol <- ceiling(sqrt(n))
    nrow <- ceiling(n / ncol)
  } else if (is.null(ncol)) {
    ncol <- ceiling(n / nrow)
  } else if (is.null(nrow)) {
    nrow <- ceiling(n / ncol)
  }

  # Helper: strip y-axis ticks and labels from a ggplot
  strip_y_ticks_labels <- function(p) {
    if (is.null(p)) return(NULL)
    p + theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank()
    )
  }

  # Helper: fully remove every y-axis element (title, text, ticks, line)
  strip_y_axis_full <- function(p) {
    if (is.null(p)) return(NULL)
    p + theme(
      axis.title.y = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y  = element_blank()
    )
  }

  # Helper: remove legends from a ggplot
  strip_legend <- function(p) {
    if (is.null(p)) return(NULL)
    p + theme(legend.position = "none")
  }

  # Helper: rotate facet strip labels 90° and remove the strip background box.
  # Applied AFTER the user theme so it isn't overridden.
  style_facet_labels <- function(p) {
    if (is.null(p)) return(NULL)
    p + theme(
      strip.text.y        = element_text(angle = 0),
      strip.text.y.left   = element_text(angle = 0),
      strip.text.y.right  = element_text(angle = 0),
      strip.text.x        = element_text(angle = 270, hjust = 0, vjust = 0.5),
      strip.text.x.top    = element_text(angle = 270, hjust = 0, vjust = 0.5),
      strip.text.x.bottom = element_text(angle = 270, hjust = 0, vjust = 0.5),
      strip.background    = element_blank()
    )
  }

  # Helper: remove x and y axis titles, text, and ticks from the expression plot
  strip_expr_axes <- function(p) {
    if (is.null(p)) return(NULL)
    p + theme(
      # axis.title.x = element_blank(),
      axis.title.y       = element_blank(),
      axis.text.x.bottom = element_blank(),
      axis.text.y        = element_blank(),
      axis.ticks.x       = element_blank(),
      axis.ticks.y       = element_blank()
    )
  }

  # Helper: remove facet strip labels entirely (for non-first-column coverage
  # plots and for expression plots). Targets every strip position explicitly
  # so it works regardless of which side the facet strip is on.
  strip_facet_labels <- function(p) {
    if (is.null(p)) return(NULL)
    p + theme(
      strip.text          = element_blank(),
      strip.text.x        = element_blank(),
      strip.text.x.top    = element_blank(),
      strip.text.x.bottom = element_blank(),
      strip.text.y        = element_blank(),
      strip.text.y.left   = element_blank(),
      strip.text.y.right  = element_blank(),
      strip.background    = element_blank()
    )
  }

  # Helper: regex test for a coordinate-style region string
  is_coord_string <- function(s) {
    grepl("^[^:\\-]+[:\\-][0-9]+[:\\-][0-9]+$", s)
  }

  # Helper: resolve a gene symbol to a coordinate string.
  resolve_region <- function(region) {
    if (is_coord_string(region)) {
      return(region)
    }
    annotations <- Annotation(object = seurat[[assay_atac]])
    if (is.null(annotations)) {
      stop("No annotation found on assay '", assay_atac,
           "'. Cannot resolve gene symbol '", region, "'.")
    }
    gr <- annotations[annotations$gene_name == region]
    if (length(gr) == 0) {
      stop("Gene '", region, "' not found in annotations of assay '",
           assay_atac, "'.")
    }
    chr   <- as.character(GenomicRanges::seqnames(gr))[1]
    start <- min(GenomicRanges::start(gr))
    end   <- max(GenomicRanges::end(gr))
    paste(chr, start, end, sep = "-")
  }

  # Helper: compute the per-track title for a region
  make_track_title <- function(region_input, region_idx) {
    # Explicit override takes priority
    if (!is.null(track_titles)) {
      return(track_titles[[region_idx]])
    }

    # If region is a coordinate string, just use it verbatim (no suffix)
    if (is_coord_string(region_input)) {
      return(region_input)
    }

    # Resolve suffix for this region
    suffix <- if (is.null(track_title_suffix)) {
      ""
    } else if (length(track_title_suffix) == 1) {
      track_title_suffix
    } else if (!is.null(names(track_title_suffix)) &&
               region_input %in% names(track_title_suffix)) {
      track_title_suffix[[region_input]]
    } else if (is.null(names(track_title_suffix)) &&
               region_idx <= length(track_title_suffix)) {
      track_title_suffix[[region_idx]]
    } else {
      ""
    }

    if (nzchar(suffix)) paste(region_input, suffix) else region_input
  }

  # Build one combined track plot for a single region
  build_combined_plot <- function(region_input, col_idx, ncol_total,
                                  expr_feature = NULL, region_idx = NULL) {

    is_first_col <- col_idx == 1
    coord_region <- resolve_region(region_input)

    cov_plot <- CoveragePlot(
      object            = seurat,
      region            = coord_region,
      extend.upstream   = extend_upstream,
      extend.downstream = extend_downstream,
      annotation        = FALSE,
      links             = FALSE,
      peaks             = FALSE
    ) + theme(axis.title.y = element_text(angle = 90, vjust = 0.5)) 

    gene_plot <- AnnotationPlot_custom(
      object            = seurat,
      region            = coord_region,
      extend.upstream   = extend_upstream,
      extend.downstream = extend_downstream,
      annotation_font_size = annotation_font_size
    ) + theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

    peak_plot <- PeakPlot(
      object            = seurat,
      region            = coord_region,
      extend.upstream   = extend_upstream,
      extend.downstream = extend_downstream
    ) + theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

    link_plot <- if (add_links) {
      LinkPlot(
        object            = seurat,
        region            = coord_region,
        min.cutoff        = 0,
        sep               = c("-", "-"),
        extend.upstream   = extend_upstream,
        extend.downstream = extend_downstream,
        scale.linewidth   = FALSE,
      ) + scale_color_gradientn(colors = 'orange3') +
         theme(axis.title.y = element_text(angle = 0, vjust = 0.5))
    } else {
      NULL
    }

    expr_plot <- if (gene_symbol) {
      feat <- if (!is.null(expr_feature)) expr_feature else region_input
      avail <- rownames(seurat[[assay]])
      if (feat %in% avail) {
        ExpressionPlot(
          object   = seurat,
          features = feat,
          assay    = assay
        )
      } else {
        warning("Feature '", feat, "' not found in assay '", assay,
                "' — skipping expression plot.")
        NULL
      }
    } else {
      NULL
    }

    # Apply the user-supplied theme FIRST so our later theme() calls
    # (especially strip label rotation) take precedence
    if (!is.null(theme)) {
      cov_plot  <- cov_plot  + theme
      gene_plot <- gene_plot + theme
      peak_plot <- peak_plot + theme
      if (!is.null(link_plot)) link_plot <- link_plot + theme
      if (!is.null(expr_plot)) expr_plot <- expr_plot + theme
    }

    # Handle y-axis based on column position
    if (is_first_col) {
      # First column: keep axis, blank ticks/text, override coverage title
      cov_plot  <- strip_y_ticks_labels(cov_plot)
      gene_plot <- strip_y_ticks_labels(gene_plot)
      peak_plot <- strip_y_ticks_labels(peak_plot)
      link_plot <- strip_y_ticks_labels(link_plot)
      expr_plot <- strip_y_ticks_labels(expr_plot)

      cov_plot <- cov_plot + ylab("Normalized signal")
    } else {
      # Non-first column: remove everything y-axis related
      cov_plot  <- strip_y_axis_full(cov_plot)
      gene_plot <- strip_y_axis_full(gene_plot)
      peak_plot <- strip_y_axis_full(peak_plot)
      link_plot <- strip_y_axis_full(link_plot)
      expr_plot <- strip_y_axis_full(expr_plot)
    }

    # Remove legends from every track
    cov_plot  <- strip_legend(cov_plot)
    gene_plot <- strip_legend(gene_plot)
    peak_plot <- strip_legend(peak_plot)
    link_plot <- strip_legend(link_plot)
    expr_plot <- strip_legend(expr_plot)

    # Facet labels: style for first column, remove for others — applied
    # LAST so the user-supplied theme can't reset strip.text back
    if (is_first_col) {
      cov_plot <- style_facet_labels(cov_plot)
    } else {
      cov_plot <- strip_facet_labels(cov_plot)
    }
    expr_plot <- strip_facet_labels(expr_plot)

    # Remove axis titles, text, and ticks from the expression plot
    expr_plot <- strip_expr_axes(expr_plot)

    # Assemble plotlist, dropping NULLs
    plotlist      <- list(cov_plot, peak_plot, gene_plot, link_plot)
    keep          <- !vapply(plotlist, is.null, logical(1))
    plotlist      <- plotlist[keep]
    track_heights <- heights[keep]

    combined <- CombineTracks(
      plotlist        = plotlist,
      expression.plot = expr_plot,
      heights         = track_heights,
      widths          = widths
    )

    # Adjust left/right margins: zero out interior boundaries
    is_leftmost  <- col_idx == 1
    is_rightmost <- col_idx == ncol_total
    l_margin <- if (is_leftmost)  1 else 0
    r_margin <- if (is_rightmost) 1 else 0
    combined <- combined & theme(plot.margin = margin(1, r_margin, 1, l_margin))

    # Add per-track title if requested, then freeze the inner patchwork so
    # the annotation survives nesting inside the outer wrap_plots()
    if (add_track_titles || !is.null(track_titles)) {
      track_title <- make_track_title(region_input, region_idx)

      default_track_title_theme <- fig_theme + theme(
        plot.title = element_text(hjust = 0.5, face = "italic", margin = margin(t = 2, b = 0)))
    
#         default_track_title_theme <- theme(
#   plot.title = element_text(hjust = 0.5, face = "bold", size = 20,
#                             color = "red",
#                             margin = margin(t = 30, b = 30)),
#   plot.background = element_rect(fill = "yellow", color = NA),
#   plot.margin = margin(t = 30, r = 0, b = 0, l = 0)
# )
        
    combined <- combined + plot_annotation(
        title = track_title,
        theme = if (is.null(track_title_theme)) default_track_title_theme else track_title_theme
      )

      # Bake the title into a single grob so the outer patchwork can't strip it
      combined <- patchwork::wrap_elements(full = combined) 
    }

    combined
  }

  # Build one combined plot per region, flagging first-column positions
  combined_plots <- lapply(seq_along(regions), function(i) {
    col_idx <- ((i - 1) %% ncol) + 1
    ef <- if (!is.null(expression_features)) expression_features[[i]] else NULL
    build_combined_plot(regions[[i]], col_idx = col_idx, ncol_total = ncol,
                        expr_feature = ef, region_idx = i)
  })

  # Gray dashed vertical separator between columns (zero margins)
  sep_line <- wrap_elements(grid::linesGrob(
    x  = grid::unit(c(0.5, 0.5), "npc"),
    y  = grid::unit(c(0, 1), "npc"),
    gp = grid::gpar(col = "gray", lty = "dashed", lwd = 1)
  )) & theme(plot.margin = margin(0, 0, 0, 0))
  sep_width <- 0.005            # very narrow relative width

  # Build each row, interleaving separator lines between columns
  row_patches <- vector("list", nrow)
  for (r in seq_len(nrow)) {
    row_items  <- list()
    row_widths <- c()
    for (c in seq_len(ncol)) {
      idx <- (r - 1) * ncol + c
      if (idx > n) break
      if (c > 1) {
        row_items  <- c(row_items, list(sep_line))
        row_widths <- c(row_widths, sep_width)
      }
      row_items  <- c(row_items, list(combined_plots[[idx]]))
      row_widths <- c(row_widths, 1)
    }
    row_patches[[r]] <- wrap_plots(row_items, nrow = 1) +
      plot_layout(widths = row_widths)
  }

  # Stack rows vertically
  final <- wrap_plots(row_patches, ncol = 1)

  # Add a shared title / subtitle if requested
  if (!is.null(title) || !is.null(subtitle)) {
    default_title_theme <- theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    )
    final <- final + plot_annotation(
      title    = title,
      subtitle = subtitle,
      theme    = if (is.null(title_theme)) default_title_theme else title_theme
    )
  }

  return(final)
}

