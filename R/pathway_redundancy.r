## Pathway Redundancy Collapse
## Collapses a redundant over-representation result into non-redundant gene-overlap
## clusters, and visualises them. Similarity is Jaccard on the enriched (leading-edge)
## genes, so it works across the mixed MSigDB collections used by
## pathway_overrepresentation_analysis() (HALLMARK_/REACTOME_/WP_/BIOCARTA_/GOBP_), where
## GO semantic similarity -- and therefore clusterProfiler::simplify() -- is undefined.
## Clustering is Louvain on the thresholded similarity graph.

#' Shorthand Topic Vocabularies for Pathway Redundancy Analysis
#'
#' Returns the built-in keyword vocabularies used by [pathway_topic_terms()] and the
#' \code{topics} argument of [pathway_overrepresentation_analysis()]. A bare keyword like
#' \code{"lipid"} misses \code{Cholesterol efflux} and \code{Acylglycerol homeostasis}, so
#' each entry expands to a synonym set that is OR-ed into a single case-insensitive regex.
#'
#' Anywhere a topic is accepted you may pass one of these names, an arbitrary regex vector,
#' or a named list of regex vectors; the vocabulary is a convenience, not a constraint.
#'
#' @return Named list of character vectors: \code{lipid}, \code{mitochondria},
#'   \code{translation}, \code{cell_cycle}, \code{chromatin}, \code{immune}, \code{stress}.
#'
#' @export
pathway_topics <- function() {
    list(
        lipid        = c('lipid', 'lipoprotein', 'cholesterol', 'sterol', 'fatty acid',
                         'triglyceride', 'acylglycerol', 'acyl', 'lipase', 'adipo', 'eicosanoid'),
        mitochondria = c('mitochondri', 'respiratory', 'oxidative phosphorylation',
                         'electron transport', 'atp synth', 'cellular respiration', 'complex i'),
        translation  = c('translation', 'ribosom', 'trna', 'peptide chain'),
        cell_cycle   = c('cell cycle', 'mitotic', 'mitosis', 'chromosome', 'spindle',
                         'dna replication'),
        chromatin    = c('chromatin', 'histone', 'epigenetic', 'nucleosome', 'telomer',
                         'methylat', 'acetylat', 'dna repair'),
        immune       = c('immune', 'interleukin', 'interferon', 'cytokine', 'antigen',
                         'inflammat', 'antimicrobial', 'defense response'),
        stress       = c('oxidant', 'detoxif', 'glutathione', 'reactive oxygen',
                         'xenobiotic', 'chemical stress', 'unfolded protein')
    )
}

#' Shared Theme for Pathway Redundancy Plots
#'
#' @return A ggplot2 theme
#'
#' @noRd
pathway_redundancy_theme <- function() {
    theme_minimal(base_size = 11) +
        theme(panel.grid.minor = element_blank(), plot.title = element_text(face = 'bold'))
}

#' Is a Pathway Cluster Result Empty?
#'
#' Guard used by every plot function so a thin or absent collapse degrades to a placeholder
#' rather than an error: these plots run inside DE wrappers, after the expensive steps.
#'
#' @param res Object returned by [collapse_pathway_redundancy()], possibly NULL
#'
#' @return Logical
#'
#' @noRd
pathway_clusters_is_empty <- function(res) {
    is.null(res) || !inherits(res, 'pathway_clusters') || nrow(res$clusters) == 0
}

#' Build a ggplot Layer Without the "unknown aesthetics" Warning
#'
#' The \code{text} aesthetic carries plotly hover text. ggplot2 itself ignores it and warns,
#' which would surface in every Quarto render, so only that warning is muffled.
#'
#' @param layer_expression A ggplot2 layer
#'
#' @return The layer
#'
#' @noRd
suppress_unknown_aes_warning <- function(layer_expression) {
    withCallingHandlers(
        layer_expression,
        warning = function(w) {
            if (grepl('unknown aesthetics', conditionMessage(w), fixed = TRUE)) {
                invokeRestart('muffleWarning')
            }
        }
    )
}


#' Sanitize a String for Use in a Filename
#'
#' Same rule as \code{de_contrast_labels()} uses for contrast file tags.
#'
#' @param x Character vector
#'
#' @return Character vector safe to embed in a filename
#'
#' @noRd
sanitize_for_filename <- function(x) {
    gsub('[^A-Za-z0-9._-]+', '.', x)
}

#' Collapse a Redundant Over-Representation Result into Gene-Overlap Clusters
#'
#' Groups the significant terms of an over-representation result into non-redundant clusters
#' by Jaccard similarity of their enriched (leading-edge) genes, followed by Louvain
#' community detection on the thresholded similarity graph. Because similarity is computed
#' on genes rather than ontology structure, it works across mixed MSigDB collections
#' (\code{HALLMARK_}, \code{REACTOME_}, \code{WP_}, \code{BIOCARTA_}, \code{GOBP_}) where
#' \code{clusterProfiler::simplify()} does not apply.
#'
#' @param x A clusterProfiler \code{enrichResult}, a data frame, or a path to an exported
#'   CSV. Needs \code{ID}, \code{p.adjust} and \code{geneID} columns; \code{Description} and
#'   \code{Count} are used when present and derived otherwise. The \code{geneID} column must
#'   be \code{'/'}-separated, as written by \code{clusterProfiler::setReadable()}.
#' @param padj Numeric; adjusted p-value cutoff for including a term. Default 0.05
#' @param min_sim Numeric; Jaccard similarity required for an edge between two terms.
#'   0.2-0.3 is the useful range. Default 0.25
#' @param resolution Numeric; Louvain granularity. Values above 1 split large clusters.
#'   Default 1
#' @param seed Integer; Louvain is stochastic, so the seed is fixed for reproducibility.
#'   Default 1
#' @param name Character or NULL; label used as the default plot title. Default NULL
#'
#' @return An object of class \code{pathway_clusters}: a list with \code{terms} (one row per
#'   term, with its \code{collection}, cleaned \code{label}, \code{genes} list-column,
#'   \code{cluster} and cluster \code{rank}), \code{clusters} (one row per cluster, ranked by
#'   the adjusted p-value of its best member), \code{sim} (the square Jaccard matrix) and
#'   \code{params}. Returns \code{NULL} with a message when fewer than two terms pass
#'   \code{padj}, so the function is safe to call inside a pipeline.
#'
#' @export
collapse_pathway_redundancy <- function(x, padj = 0.05, min_sim = 0.25, resolution = 1,
                                        seed = 1, name = NULL) {

    rlang::check_installed('igraph', reason = 'to cluster pathway terms by gene overlap.')

    enrichment_table <- if (inherits(x, 'enrichResult') || inherits(x, 'gseaResult')) {
        as.data.frame(x)
    } else if (is.character(x) && length(x) == 1) {
        readr::read_csv(x, show_col_types = FALSE, name_repair = 'minimal')
    } else {
        as.data.frame(x)
    }
    enrichment_table <- as_tibble(enrichment_table, .name_repair = 'minimal')

    # The package writes enrichment CSVs with base write.csv and no row.names = FALSE, so
    # they carry a leading index column: unnamed on re-read with utils, '...1' with readr.
    index_columns <- names(enrichment_table) == '' |
        str_detect(names(enrichment_table), '^\\.\\.\\.[0-9]+$')
    enrichment_table <- enrichment_table[, !index_columns, drop = FALSE]

    missing_columns <- setdiff(c('ID', 'p.adjust', 'geneID'), names(enrichment_table))
    if (length(missing_columns)) {
        stop('missing column(s): ', paste(missing_columns, collapse = ', '), call. = FALSE)
    }

    terms <- enrichment_table |>
        filter(.data$p.adjust < padj) |>
        mutate(
            Description = if ('Description' %in% names(enrichment_table)) .data$Description else .data$ID,
            collection  = str_extract(.data$ID, '^[A-Z]+'),
            label       = str_to_sentence(str_replace_all(
                str_remove(.data$ID, '^[A-Z0-9]+_'), '_', ' ')),
            genes       = map(str_split(.data$geneID, fixed('/')),
                              \(genes) unique(genes[!is.na(genes) & nzchar(genes)]))
        ) |>
        mutate(Count = if ('Count' %in% names(enrichment_table)) .data$Count else lengths(.data$genes)) |>
        filter(lengths(.data$genes) > 0) |>
        arrange(.data$p.adjust)

    if (nrow(terms) < 2) {
        message('Pathway redundancy: only ', nrow(terms), ' term passes p.adjust < ', padj,
                ' - nothing to collapse.')
        return(NULL)
    }

    sim_matrix <- pathway_jaccard_matrix(terms)

    ## Louvain on the thresholded similarity graph
    adjacency <- sim_matrix
    adjacency[adjacency < min_sim] <- 0
    diag(adjacency) <- 0
    term_graph <- igraph::graph_from_adjacency_matrix(adjacency, mode = 'undirected',
                                                      weighted = TRUE)
    set.seed(seed)
    communities <- igraph::membership(igraph::cluster_louvain(term_graph, resolution = resolution))
    terms <- terms |> mutate(cluster = as.integer(communities[.data$ID]))

    clusters <- terms |>
        group_by(.data$cluster) |>
        summarise(
            n_terms       = n(),
            best_padj     = min(.data$p.adjust),
            cluster_label = .data$label[which.min(.data$p.adjust)],
            n_genes       = length(unique(unlist(.data$genes))),
            collections   = paste(sort(unique(.data$collection)), collapse = ','),
            .groups       = 'drop'
        ) |>
        arrange(.data$best_padj) |>
        mutate(rank = row_number())

    terms <- terms |>
        left_join(select(clusters, 'cluster', 'cluster_label', 'rank'), by = 'cluster')

    structure(
        list(
            terms    = terms,
            clusters = clusters,
            sim      = sim_matrix,
            params   = list(padj = padj, min_sim = min_sim, resolution = resolution,
                            seed = seed, name = name)
        ),
        class = c('pathway_clusters', 'list')
    )
}

#' Jaccard Similarity Matrix over the Enriched Genes of Each Term
#'
#' Pairwise intersections come from a self-join on the exploded gene-by-term table, which
#' says what Jaccard is more directly than a membership matrix does.
#'
#' @param terms Tibble with an \code{ID} column and a \code{genes} list-column
#'
#' @return Square numeric matrix of Jaccard similarities, rows and columns in the order of
#'   \code{terms$ID}, with 1 on the diagonal
#'
#' @noRd
pathway_jaccard_matrix <- function(terms) {

    genes_long <- terms |>
        select('ID', 'genes') |>
        tidyr::unnest('genes') |>
        distinct(.data$ID, .data$genes)

    set_sizes <- genes_long |> count(.data$ID, name = 'size')

    sim_matrix <- genes_long |>
        inner_join(genes_long, by = 'genes', relationship = 'many-to-many') |>
        count(.data$ID.x, .data$ID.y, name = 'intersection') |>
        left_join(set_sizes, by = c('ID.x' = 'ID')) |>
        left_join(set_sizes, by = c('ID.y' = 'ID'), suffix = c('.x', '.y')) |>
        mutate(jaccard = .data$intersection /
                   (.data$size.x + .data$size.y - .data$intersection)) |>
        select('ID.x', 'ID.y', 'jaccard') |>
        tidyr::pivot_wider(names_from = 'ID.y', values_from = 'jaccard', values_fill = 0) |>
        column_to_rownames('ID.x') |>
        as.matrix()

    sim_matrix[terms$ID, terms$ID, drop = FALSE]
}

#' Print a Pathway Cluster Result
#'
#' One-line summary of the collapse followed by the five most significant clusters.
#'
#' @param x A \code{pathway_clusters} object, as returned by
#'   [collapse_pathway_redundancy()]
#' @param ... Ignored, present for consistency with the \code{print} generic
#'
#' @return \code{x}, invisibly
#'
#' @export
print.pathway_clusters <- function(x, ...) {
    params <- x$params
    cat(sprintf('%s%d terms -> %d clusters (jaccard >= %.2f, louvain, seed %s)\n',
                if (is.null(params$name)) '' else paste0(params$name, ': '),
                nrow(x$terms), nrow(x$clusters), params$min_sim, params$seed))
    top_clusters <- head(x$clusters, 5)
    cat(paste0('  ', top_clusters$rank, '. ', str_trunc(top_clusters$cluster_label, 50),
               ' (', top_clusters$n_terms, ' terms, p.adj ',
               signif(top_clusters$best_padj, 2), ')',
               collapse = '\n'), '\n')
    invisible(x)
}

#' Write Pathway Redundancy Clusters to CSV
#'
#' Writes the per-term assignments to \code{file} and the per-cluster summary alongside it.
#'
#' @param res Object returned by [collapse_pathway_redundancy()]
#' @param file Character; path for the per-term CSV. The summary is written to the same path
#'   with \code{_summary} inserted before the extension.
#'
#' @return \code{res}, invisibly
#'
#' @noRd
write_pathway_clusters <- function(res, file) {
    if (pathway_clusters_is_empty(res)) {
        return(invisible(res))
    }
    res$terms |>
        transmute(
            rank          = .data$rank,
            cluster_label = .data$cluster_label,
            ID            = .data$ID,
            Description   = .data$Description,
            collection    = .data$collection,
            p.adjust      = .data$p.adjust,
            Count         = .data$Count,
            genes         = map_chr(.data$genes, paste, collapse = '/')
        ) |>
        arrange(.data$rank, .data$p.adjust) |>
        readr::write_csv(file)
    readr::write_csv(res$clusters, sub('\\.csv$', '_summary.csv', file))
    invisible(res)
}

#' Bar Plot of Collapsed Pathway Clusters
#'
#' One bar per cluster, labelled by its most significant member, so a block of eight
#' near-duplicate terms reads as one signal instead of eight.
#'
#' @param res Object returned by [collapse_pathway_redundancy()]
#' @param n Integer; number of clusters to show. Default 15
#' @param title Character or NULL; plot title. Defaults to the \code{name} given to
#'   [collapse_pathway_redundancy()]
#'
#' @return A ggplot object, or a placeholder plot when \code{res} is NULL or empty
#'
#' @export
plot_pathway_clusters <- function(res, n = 15, title = res$params$name) {

    if (pathway_clusters_is_empty(res)) {
        return(create_no_data_plot())
    }

    plot_data <- head(res$clusters, n) |>
        mutate(
            log_padj = -log10(.data$best_padj),
            # Rank-prefixed so two clusters whose labels truncate identically stay distinct
            cluster_display = fct_reorder(
                paste0(.data$rank, '. ', str_trunc(.data$cluster_label, 44)),
                .data$log_padj
            )
        )

    ggplot(plot_data, aes(x = .data$log_padj, y = .data$cluster_display, fill = .data$n_genes)) +
        geom_col(width = 0.7) +
        geom_text(aes(label = paste(.data$n_terms, if_else(.data$n_terms == 1, 'term', 'terms'))),
                  hjust = 0, nudge_x = 0.3, size = 2.9, colour = 'grey30') +
        scale_x_continuous(expand = expansion(c(0, 0.16))) +
        scale_fill_viridis_c(option = 'mako', direction = -1, begin = 0.15, end = 0.9,
                             name = 'genes in\ncluster') +
        labs(x = expression(-log[10] ~ italic(p)[adj] ~ 'of best member'), y = '', title = title,
             subtitle = sprintf('%d terms collapsed to %d clusters',
                                nrow(res$terms), nrow(res$clusters))) +
        pathway_redundancy_theme()
}

#' Enrichment Map of Collapsed Pathway Clusters
#'
#' Network of the significant terms: nodes are terms sized by gene count, edges are gene
#' overlap above \code{min_sim}, colour is the cluster. One label per cluster, on its most
#' significant member.
#'
#' @param res Object returned by [collapse_pathway_redundancy()]
#' @param max_clusters Integer; number of clusters to colour. The rest are drawn grey as
#'   \code{'other'}. Default 12
#' @param title Character or NULL; plot title. Defaults to the \code{name} given to
#'   [collapse_pathway_redundancy()]
#' @param seed Integer; seed for the Fruchterman-Reingold layout. Default 1
#' @param label_length Integer; truncation width for node labels. Default 40
#' @param interactive Logical; \code{FALSE} (the default) returns a ggplot, \code{TRUE} a
#'   \code{plotly} htmlwidget whose nodes each carry a hover tooltip with the full term name,
#'   collection, adjusted p-value, gene count and cluster. Worth asking for, because the static
#'   version can only label one term per cluster. Default FALSE, so the return value stays
#'   something \code{ggsave()} accepts
#'
#' @return A ggplot object, a \code{plotly} htmlwidget when \code{interactive = TRUE}, or a
#'   placeholder plot when \code{res} is NULL or empty
#'
#' @details
#' \code{ggplotly()} cannot render repelled labels, so the interactive version places the
#' per-cluster labels without repulsion and leans on the tooltips instead.
#'
#' knitr embeds a widget only when it is the value of a top-level expression, so to see one in a
#' knitted document call this function in its own chunk on the object the pipeline hands back --
#' the pipeline itself only ever builds the static version:
#' \preformatted{
#' de <- pseudobulk(..., run_pathway_enrichment = "ClusterProfiler",
#'                  collapse_redundancy = TRUE)
#' plot_pathway_cluster_map(de$enrichment$redundancy_up, interactive = TRUE)
#' }
#'
#' @export
plot_pathway_cluster_map <- function(res, max_clusters = 12, title = res$params$name,
                                     seed = 1, label_length = 40, interactive = FALSE) {

    rlang::check_installed('igraph', reason = 'to lay out the pathway cluster map.')
    if (interactive) {
        rlang::check_installed('plotly', reason = 'to draw an interactive pathway cluster map.')
    }
    if (pathway_clusters_is_empty(res)) {
        return(create_no_data_plot())
    }

    adjacency <- res$sim
    adjacency[adjacency < res$params$min_sim] <- 0
    diag(adjacency) <- 0
    term_graph <- igraph::graph_from_adjacency_matrix(adjacency, mode = 'undirected',
                                                      weighted = TRUE)

    set.seed(seed)
    layout_coordinates <- igraph::layout_with_fr(term_graph,
                                                 weights = igraph::E(term_graph)$weight)

    top_clusters <- head(res$clusters$cluster, max_clusters)
    # Rank-prefixed so two clusters whose labels truncate identically get distinct colours
    cluster_key <- res$clusters |>
        mutate(cluster_key = paste0(.data$rank, '. ', str_trunc(.data$cluster_label, 30))) |>
        select('cluster', 'cluster_key')

    labelled_terms <- res$terms |>
        filter(.data$cluster %in% top_clusters) |>
        group_by(.data$cluster) |>
        slice_min(.data$p.adjust, n = 1, with_ties = FALSE) |>
        ungroup() |>
        pull('ID')

    nodes <- tibble(ID = igraph::V(term_graph)$name,
                    x = layout_coordinates[, 1],
                    y = layout_coordinates[, 2]) |>
        left_join(res$terms, by = 'ID') |>
        left_join(cluster_key, by = 'cluster') |>
        mutate(
            cluster_group = if_else(.data$cluster %in% top_clusters, .data$cluster_key, 'other'),
            node_label = if_else(.data$ID %in% labelled_terms,
                                 str_trunc(.data$label, label_length), NA_character_),
            hover = paste0(
                '<b>', .data$label, '</b>',
                '<br>', .data$collection,
                '<br>cluster ', .data$rank, ': ', .data$cluster_label,
                '<br>p.adjust = ', signif(.data$p.adjust, 3),
                '<br>', .data$Count, ' genes'
            )
        )

    edge_list <- igraph::as_edgelist(term_graph)
    edges <- tibble(from = edge_list[, 1], to = edge_list[, 2]) |>
        left_join(select(nodes, 'ID', 'x', 'y'), by = c('from' = 'ID')) |>
        left_join(select(nodes, 'ID', 'x', 'y'), by = c('to' = 'ID'),
                  suffix = c('', '_end'))

    cluster_palette <- c(
        setNames(scales::hue_pal()(length(top_clusters)),
                 cluster_key$cluster_key[match(top_clusters, cluster_key$cluster)]),
        other = 'grey85'
    )

    subtitle <- sprintf('%d terms, %d clusters (jaccard >= %.2f)',
                        nrow(res$terms), nrow(res$clusters), res$params$min_sim)

    cluster_map <- ggplot(nodes, aes(x = .data$x, y = .data$y)) +
        geom_segment(data = edges,
                     aes(x = .data$x, y = .data$y, xend = .data$x_end, yend = .data$y_end),
                     colour = 'grey85', alpha = 0.6, inherit.aes = FALSE) +
        # 'text' is inert in ggplot2 and only read by ggplotly, so it warns about an unknown
        # aesthetic; the warning is muffled rather than the tooltip dropped.
        suppress_unknown_aes_warning(
            geom_point(aes(size = .data$Count, fill = .data$cluster_group, text = .data$hover),
                       shape = 21, colour = 'grey30', stroke = 0.3)
        ) +
        scale_size_area(max_size = 8, name = 'genes in term') +
        scale_fill_manual(values = cluster_palette, name = 'cluster') +
        theme_void(base_size = 11) +
        theme(plot.title = element_text(face = 'bold'))

    if (!interactive) {
        return(
            cluster_map +
                geom_text_repel(aes(label = .data$node_label), size = 2.8, fontface = 'bold',
                                box.padding = 0.5, max.overlaps = 30, na.rm = TRUE,
                                min.segment.length = 0, segment.colour = 'grey50') +
                labs(title = title, subtitle = subtitle)
        )
    }

    # ggplotly() drops both the subtitle and repelled labels, so the subtitle moves into the
    # plotly title and the labels are placed unrepelled above each node.
    interactive_map <- cluster_map +
        geom_text(aes(label = .data$node_label), size = 2.8, fontface = 'bold',
                  vjust = -1.1, na.rm = TRUE)

    plotly::ggplotly(interactive_map, tooltip = 'text') |>
        plotly::layout(title = list(
            text = paste0('<b>', title %||% '', '</b><br><sup>', subtitle, '</sup>'),
            x = 0.02, xanchor = 'left'
        ))
}

#' Heatmap of Pathway Term Overlap
#'
#' Jaccard similarity heatmap with the cluster blocks marked out: the evidence for the
#' collapse, and a direct picture of the redundancy itself.
#'
#' @param res Object returned by [collapse_pathway_redundancy()]
#' @param title Character or NULL; heatmap title. Defaults to the \code{name} given to
#'   [collapse_pathway_redundancy()]
#' @param min_label Integer; cluster row titles are dropped for clusters with fewer terms
#'   than this, because they collide. Default 3
#'
#' @return A \code{ComplexHeatmap} Heatmap object, or a placeholder ggplot when \code{res}
#'   is NULL or empty
#'
#' @export
plot_pathway_overlap <- function(res, title = res$params$name, min_label = 3) {

    rlang::check_installed('ComplexHeatmap', reason = 'to draw the pathway overlap heatmap.')
    if (pathway_clusters_is_empty(res)) {
        return(create_no_data_plot())
    }

    terms <- res$terms |> arrange(.data$rank, .data$p.adjust)
    cluster_ranks <- sort(unique(terms$rank))
    cluster_sizes <- as.integer(table(terms$rank)[as.character(cluster_ranks)])
    row_titles <- if_else(cluster_sizes >= min_label, as.character(cluster_ranks), '')

    ComplexHeatmap::Heatmap(
        res$sim[terms$ID, terms$ID],
        name                 = 'jaccard',
        row_split            = factor(terms$rank),
        column_split         = factor(terms$rank),
        show_row_names       = FALSE,
        show_column_names    = FALSE,
        cluster_row_slices   = FALSE,
        cluster_column_slices = FALSE,
        row_title            = row_titles,
        row_title_gp         = grid::gpar(fontsize = 8),
        row_gap              = grid::unit(0.6, 'mm'),
        column_gap           = grid::unit(0.6, 'mm'),
        column_title         = sprintf('%s\n%d terms, %d clusters', title %||% '',
                                       nrow(terms), nrow(res$clusters)),
        col                  = circlize::colorRamp2(c(0, 0.5, 1),
                                                    c('white', '#7fcdbb', '#253494'))
    )
}

#' Resolve a Topic Specification to a Single Regex
#'
#' The pattern is passed through verbatim; case-insensitivity is applied at match time with
#' \code{stringr::regex(ignore_case = TRUE)} rather than by lower-casing the pattern, which
#' would break character classes like \code{[A-Z]} and invert escapes like \code{\\D}.
#'
#' @param pattern Character vector; either a single [pathway_topics()] name or a set of
#'   regexes to OR together
#'
#' @return Character; one regex
#'
#' @noRd
resolve_topic_pattern <- function(pattern) {
    topics <- pathway_topics()
    if (length(pattern) == 1L && pattern %in% names(topics)) {
        pattern <- topics[[pattern]]
    }
    paste(pattern, collapse = '|')
}

#' Pathway Terms Belonging to a Topic
#'
#' Filters a collapsed result to the terms matching a topic, keeping the cluster each term
#' belongs to so a keyword hit can be read against the redundancy structure.
#'
#' @param res Object returned by [collapse_pathway_redundancy()]
#' @param pattern Character vector; either a single [pathway_topics()] name (e.g.
#'   \code{"lipid"}), which expands to that vocabulary, or an arbitrary regex vector, whose
#'   elements are OR-ed into one pattern. Matching is case-insensitive, and the regex is used
#'   verbatim, so metacharacters and character classes work as written
#'   (\code{"^(fatty|lanosterol)"}, \code{"Il[0-9]+"}). A name is only looked up in the
#'   vocabulary when \code{pattern} has length 1, so \code{c("lipid", "immune")} is treated as
#'   the regex \code{lipid|immune}, not as two vocabularies
#' @param where One of \code{"label"} (match the term name, the default), \code{"genes"}
#'   (match the enriched gene symbols, anchored so the topic must equal the whole symbol) or
#'   \code{"both"}
#'
#' @return Tibble of matching rows of \code{res$terms}, ordered by adjusted p-value, with
#'   the resolved regex in an \code{rx} attribute and the total number of terms in an
#'   \code{n_total} attribute. Zero rows (with a warning) when nothing matches.
#'
#' @export
pathway_topic_terms <- function(res, pattern, where = c('label', 'genes', 'both')) {

    where <- match.arg(where)
    if (pathway_clusters_is_empty(res)) {
        return(NULL)
    }

    topic_regex <- resolve_topic_pattern(pattern)
    terms <- res$terms
    label_hits <- str_detect(terms$label, regex(topic_regex, ignore_case = TRUE))
    # Gene matching is anchored: a topic gene must equal the whole symbol, not appear inside it
    gene_regex <- regex(paste0('^(', topic_regex, ')$'), ignore_case = TRUE)
    gene_hits <- map_lgl(terms$genes, ~ any(str_detect(.x, gene_regex)))
    keep <- switch(where, label = label_hits, genes = gene_hits, both = label_hits | gene_hits)

    matching_terms <- terms[keep, ] |> arrange(.data$p.adjust)
    if (!nrow(matching_terms)) {
        warning('no pathway terms matched /', topic_regex, '/', call. = FALSE)
    }
    attr(matching_terms, 'rx') <- topic_regex
    attr(matching_terms, 'n_total') <- nrow(terms)
    matching_terms
}

#' Lollipop Plot of the Pathway Terms in a Topic
#'
#' The topic's terms, faceted by the cluster each belongs to. This is what a flat keyword
#' grep hides: whether the topic is one redundant block or several distinct signals.
#'
#' @param res Object returned by [collapse_pathway_redundancy()]
#' @param pattern Character vector; a [pathway_topics()] name or an arbitrary regex vector
#' @param where One of \code{"label"} (default), \code{"genes"} or \code{"both"}; see
#'   [pathway_topic_terms()]
#' @param n Integer; maximum number of terms to show. Default 30
#' @param title Character or NULL; plot title. Defaults to the resolved regex
#'
#' @return A ggplot object, or a placeholder plot when nothing matches
#'
#' @export
plot_pathway_topic <- function(res, pattern, where = 'label', n = 30, title = NULL) {

    matching_terms <- pathway_topic_terms(res, pattern, where)
    if (is.null(matching_terms) || !nrow(matching_terms)) {
        return(create_no_data_plot())
    }

    plot_data <- head(matching_terms, n) |>
        mutate(
            log_padj = -log10(.data$p.adjust),
            term_display = fct_reorder(make.unique(str_trunc(.data$label, 56), sep = ' '),
                                       .data$log_padj),
            cluster_group = paste0(.data$rank, '. ', str_trunc(.data$cluster_label, 26))
        )

    topic_plot <- ggplot(plot_data, aes(x = .data$log_padj, y = .data$term_display)) +
        geom_segment(aes(x = 0, xend = .data$log_padj, yend = .data$term_display),
                     colour = 'grey80') +
        geom_point(aes(size = .data$Count, colour = .data$collection)) +
        scale_size_area(max_size = 6, name = 'genes in term') +
        scale_colour_brewer(palette = 'Dark2', name = 'collection') +
        scale_x_continuous(expand = expansion(c(0, 0.08))) +
        labs(x = expression(-log[10] ~ italic(p)[adj]), y = '',
             title = title %||% sprintf('Terms matching /%s/',
                                        str_trunc(attr(matching_terms, 'rx'), 50)),
             subtitle = sprintf('%d of %d enriched terms, grouped by cluster',
                                nrow(matching_terms), attr(matching_terms, 'n_total'))) +
        pathway_redundancy_theme() +
        theme(axis.text.y = element_text(size = 7))

    if (n_distinct(plot_data$cluster_group) > 1) {
        topic_plot <- topic_plot +
            facet_grid(cluster_group ~ ., scales = 'free_y', space = 'free_y', switch = 'y',
                       labeller = label_wrap_gen(18)) +
            theme(strip.placement = 'outside',
                  panel.spacing.y = grid::unit(3, 'pt'),
                  strip.background = element_rect(fill = 'grey96', colour = NA),
                  strip.text.y.left = element_text(angle = 0, size = 6.5, hjust = 1))
    }
    topic_plot
}

#' Bar Plot of the Genes Driving a Pathway Topic
#'
#' How many of the topic's terms each gene appears in. A gene present in many terms is
#' precisely why those terms are redundant.
#'
#' @param res Object returned by [collapse_pathway_redundancy()]
#' @param pattern Character vector; a [pathway_topics()] name or an arbitrary regex vector
#' @param where One of \code{"label"} (default), \code{"genes"} or \code{"both"}; see
#'   [pathway_topic_terms()]
#' @param top_genes Integer; number of genes to show. Default 25
#' @param lfc Named numeric or NULL; log2 fold changes named by gene symbol. When supplied,
#'   bars are coloured by direction. Default NULL
#' @param title Character or NULL; plot title
#'
#' @return A ggplot object, or a placeholder plot when nothing matches
#'
#' @export
plot_pathway_topic_genes <- function(res, pattern, where = 'label', top_genes = 25,
                                     lfc = NULL, title = NULL) {

    matching_terms <- pathway_topic_terms(res, pattern, where)
    if (is.null(matching_terms) || !nrow(matching_terms)) {
        return(create_no_data_plot())
    }

    plot_data <- tibble(gene = unlist(matching_terms$genes)) |>
        count(.data$gene, sort = TRUE) |>
        head(top_genes) |>
        mutate(gene = fct_reorder(.data$gene, .data$n))

    gene_plot <- ggplot(plot_data, aes(x = .data$n, y = .data$gene)) +
        scale_x_continuous(expand = expansion(c(0, 0.06))) +
        labs(x = sprintf('appears in N of the %d matching terms', nrow(matching_terms)),
             y = '', title = title %||% 'Genes driving the topic',
             subtitle = 'a gene in many terms is why those terms are redundant') +
        pathway_redundancy_theme()

    if (is.null(lfc)) {
        gene_plot + geom_col(width = 0.7, fill = '#4DBBD5')
    } else {
        gene_plot +
            geom_col(aes(fill = unname(lfc[as.character(.data$gene)])), width = 0.7) +
            scale_fill_gradient2(low = '#3C5488', mid = 'grey92', high = '#E64B35',
                                 midpoint = 0, name = expression(log[2] ~ FC))
    }
}

#' Run and Save the Pathway Redundancy Collapse for One Enrichment Result
#'
#' Collapses an over-representation result and writes the cluster bar plot, the term-overlap
#' network and both cluster tables into \code{local_path}, plus a pair of plots per topic when
#' \code{topics} is supplied. Nothing is printed: the plots are returned so the caller can place
#' them where it wants. Called by [pathway_overrepresentation_analysis()]; the output directory
#' already encodes the group, cluster and contrast, so the filenames need no further tagging.
#'
#' @param enrichment_results clusterProfiler \code{enrichResult}
#' @param local_path Character; directory to write into
#' @param filename Character; prefix for output files
#' @param title Character; plot title, also the \code{name} of the returned object
#' @param topics Character vector or named list of regex vectors; topics to plot. NULL for
#'   none
#' @param padj,min_sim,resolution,seed Passed to [collapse_pathway_redundancy()]
#' @param n_clusters Integer; clusters to show in the bar plot
#' @param max_clusters Integer; clusters to colour in the network
#'
#' @return NULL when the collapse found nothing, otherwise a list of \code{clusters} (the
#'   \code{pathway_clusters} object) and \code{plots}, itself a list of \code{clusters} (the bar
#'   plot), \code{network} and \code{topics} -- the last a list named by topic, each holding
#'   \code{terms} and \code{genes} plots. Topics that matched no term are absent.
#'
#' @noRd
run_pathway_redundancy_analysis <- function(enrichment_results, local_path, filename = '',
                                            title = NULL, topics = NULL, padj = 0.05,
                                            min_sim = 0.25, resolution = 1, seed = 1,
                                            n_clusters = 15, max_clusters = 12) {

    redundancy <- collapse_pathway_redundancy(enrichment_results, padj = padj,
                                              min_sim = min_sim, resolution = resolution,
                                              seed = seed, name = title)
    if (pathway_clusters_is_empty(redundancy)) {
        return(NULL)
    }

    cluster_barplot <- plot_pathway_clusters(redundancy, n = n_clusters)
    # Static: this is what gets saved, and ggsave() cannot take a widget. A caller wanting the
    # interactive network calls plot_pathway_cluster_map(interactive = TRUE) on the result.
    cluster_network <- plot_pathway_cluster_map(redundancy, max_clusters = max_clusters)

    ggsave(plot = cluster_barplot,
           filename = paste0(filename, 'Pathway_redundancy_clusters_barplot.pdf'),
           width = 8, height = 6, path = local_path)
    ggsave(plot = cluster_network,
           filename = paste0(filename, 'Pathway_redundancy_cluster_network.pdf'),
           width = 9, height = 7, path = local_path)
    write_pathway_clusters(redundancy,
                           here::here(local_path,
                                      paste0(filename, 'Pathway_redundancy_terms.csv')))

    # Accept c('lipid', 'immune') as well as list(my_topic = c('Il2', 'Il21'))
    topic_plots <- list()
    if (!is.null(topics) && length(topics)) {
        topic_list <- if (is.list(topics)) topics else setNames(as.list(topics), topics)
        topic_plots <- purrr::imap(topic_list, function(topic_pattern, topic_name) {
            # A topic that matches nothing in this direction writes no file and is left out of
            # the returned list, following the rest of the package: placeholder plots are
            # returned by the plot functions but never saved.
            if (!nrow(pathway_topic_terms(redundancy, topic_pattern))) {
                return(NULL)
            }
            slug <- sanitize_for_filename(topic_name)
            topic_terms_plot <- plot_pathway_topic(redundancy, topic_pattern)
            topic_genes_plot <- plot_pathway_topic_genes(redundancy, topic_pattern)
            ggsave(plot = topic_terms_plot,
                   filename = paste0(filename, 'Pathway_redundancy_topic_', slug, '_terms.pdf'),
                   width = 8, height = 7, path = local_path)
            ggsave(plot = topic_genes_plot,
                   filename = paste0(filename, 'Pathway_redundancy_topic_', slug, '_genes.pdf'),
                   width = 6, height = 6, path = local_path)
            list(terms = topic_terms_plot, genes = topic_genes_plot)
        }) |>
            purrr::compact()
    }

    list(
        clusters = redundancy,
        plots    = list(clusters = cluster_barplot,
                        network  = cluster_network,
                        topics   = topic_plots)
    )
}
