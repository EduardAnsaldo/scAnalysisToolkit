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
set_enrichment_color_scale <- function() {
    color_scale <- viridis::viridis(n = 4, direction = -1)
    options(enrichplot.colours = color_scale)
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
        aes(label = !!sym(label_column), segment.size = 0.3, segment.alpha = 0.8, segment.curvature = 0)
    )
}

#' VDJ (Immune Receptor) Gene Regex Pattern
#'
#' Returns a regular expression matching VDJ genes: immunoglobulin (Igh, Igk, Igl)
#' and T-cell receptor (Trav/Trbv/Trdv/Trgv, Traj/Trbj/Trdj/Trgj, Trac/Trbc/Trdc/Trgc)
#' genes. The pattern is anchored (^) to gene-symbol prefixes to avoid matching
#' unrelated genes (e.g. Trap1), and wrapped in a scoped case-insensitive inline
#' flag group `(?i:...)` so it matches both mouse (Igh, Trav) and human (IGH, TRAV)
#' symbols without affecting any patterns it is combined with.
#'
#' @return Character; a single regex string.
#'
#' @noRd
vdj_gene_pattern <- function() {
    '(?i:^Igh|^Igl|^Igk|^Trav|^Trbv|^Trdv|^Trgv|^Traj|^Trbj|^Trdj|^Trgj|^Trac|^Trbc|^Trdc|^Trgc)'
}

#' Identify Genes to Exclude from Differential Expression
#'
#' Given a vector of gene names, returns a logical vector flagging which genes
#' should be excluded based on user-supplied regular expressions and/or the VDJ
#' gene defaults. User patterns are matched exactly as written (case-sensitive);
#' the VDJ default carries its own case-insensitive flag (see \code{vdj_gene_pattern}).
#'
#' @param genes Character vector of gene names to test.
#' @param genes_to_exclude Character vector of regular expressions; a gene is
#'   excluded if it matches any pattern. Plain gene names work as literal regexes.
#'   Default c()
#' @param exclude_vdj Logical; if TRUE, also exclude VDJ (immunoglobulin and
#'   T-cell receptor) genes. Default FALSE
#'
#' @return Logical vector the same length as \code{genes}; TRUE means exclude.
#'
#' @noRd
match_genes_to_exclude <- function(genes, genes_to_exclude = c(), exclude_vdj = FALSE) {
    patterns <- as.character(genes_to_exclude)
    if (isTRUE(exclude_vdj)) {
        patterns <- c(patterns, vdj_gene_pattern())
    }
    if (length(patterns) == 0) {
        return(rep(FALSE, length(genes)))
    }
    stringr::str_detect(genes, paste(patterns, collapse = '|'))
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

#' Derive Plot Labels from a Fitted DESeq2 Design and Contrast
#'
#' Works out which comparison a DESeq2 \code{contrast} actually encodes, and returns the
#' plot labels that describe it: a subtitle, the volcano x-axis title and direction-arrow
#' labels, and which \code{Avg_} columns the scatterplot should plot.
#'
#' Everything is read off the DESeqDataSet, which is the only authoritative source: the design
#' formula gives the terms, and \code{levels(colData(dds)[[var]])} gives each variable's groups
#' with the DESeq2 reference level first. Both are available before \code{DESeq()} is run.
#'
#' Classification builds a dictionary of the coefficient names DESeq2 will report and looks the
#' contrast up in it. It deliberately does not parse coefficient names in the other direction,
#' because they are sanitized and lossy: levels '12 mo', '12.mo' and '12-mo' all collapse to
#' \code{Age12.mo.conditionCF}, which has no unique split point, so neither the term boundary nor
#' the original level text can be recovered from the name alone.
#'
#' The design matters as much as the contrast. Under an interaction design,
#' \code{c('condition', 'CF', 'CT')} is \emph{not} a main effect: it is the condition effect at the
#' reference level of the interacting variable. So the design's term labels are inspected too.
#'
#' @param dds A DESeqDataSet, before or after \code{DESeq()} has been run.
#' @param contrast Character vector, list, or NULL; the contrast passed to
#'   \code{DESeq2::results}. See \code{?DESeq2::results}.
#' @param condition_var Character; the colData column holding the group assignment. Default
#'   'condition', as built by \code{pseudobulk_de}.
#'
#' @return A list with elements \code{type}, \code{numerator}, \code{denominator},
#'   \code{stratum}, \code{covariate}, \code{subtitle}, \code{up_label}, \code{down_label},
#'   \code{volcano_x}, \code{x_column}, \code{y_column}, \code{x_group_label},
#'   \code{y_group_label}, \code{scatter_applicable} and \code{file_tag}.
#'
#' @noRd
de_contrast_labels <- function(dds, contrast = NULL, condition_var = 'condition') {

    col_data <- as.data.frame(SummarizedExperiment::colData(dds))
    design <- DESeq2::design(dds)
    term_labels <- attr(stats::terms(design), 'term.labels')

    # Only factors get treatment-contrast coefficients worth labelling
    factor_vars <- Filter(function(v) is.factor(col_data[[v]]),
                          intersect(all.vars(design), names(col_data)))

    if (!condition_var %in% factor_vars) {
        stop("'", condition_var, "' is not a factor in the design of this DESeqDataSet",
             call. = FALSE)
    }
    condition_levels <- levels(col_data[[condition_var]])
    g1 <- condition_levels[1]                      # DESeq2 reference level
    g2 <- if (length(condition_levels) > 1) condition_levels[2] else NA_character_

    # Which variables interact with the condition in this design? Those are the ones that make a
    # condition contrast stratum-specific rather than a main effect.
    interaction_terms <- lapply(grep(':', term_labels, fixed = TRUE, value = TRUE),
                                function(x) strsplit(x, ':', fixed = TRUE)[[1]])
    condition_partners <- unique(unlist(lapply(interaction_terms, function(vs) {
        if (condition_var %in% vs && length(vs) == 2) setdiff(vs, condition_var) else NULL
    })))
    # Everything else in the design is being controlled for
    covariates <- setdiff(factor_vars, c(condition_var, condition_partners))

    # Build the dictionary of coefficient names DESeq2 will report. make.names() is what turns the
    # model matrix spelling into the reported one -- notably the interaction ':' into '.' -- and
    # mirrors DESeq2:::renameModelMatrixColumns, which applies the identical call.
    coef_map <- list()
    for (v in factor_vars) {
        lv <- levels(col_data[[v]])
        for (l in lv[-1]) {
            coef_map[[make.names(paste0(v, '_', l, '_vs_', lv[1]))]] <-
                list(kind = 'main', var = v, num = l, den = lv[1])
        }
    }
    for (vs in interaction_terms) {
        if (length(vs) != 2 || !all(vs %in% factor_vars)) next
        for (l1 in levels(col_data[[vs[1]]])[-1]) {
            for (l2 in levels(col_data[[vs[2]]])[-1]) {
                # the model matrix column order follows the formula, so accept either spelling
                for (nm in make.names(c(paste0(vs[1], l1, ':', vs[2], l2),
                                        paste0(vs[2], l2, ':', vs[1], l1)))) {
                    coef_map[[nm]] <- list(kind = 'interaction',
                                           var1 = vs[1], lev1 = l1, var2 = vs[2], lev2 = l2)
                }
            }
        }
    }

    type <- 'unknown'
    numerator <- NA_character_
    denominator <- NA_character_
    stratum <- NA_character_
    stratum_var <- NA_character_
    covariate <- NA_character_
    contrast_echo <- NULL

    covariate_text <- function(vars) if (length(vars)) paste(vars, collapse = ' and ') else NA_character_

    # A contrast on the condition means different things depending on the design: the effect within
    # the reference level of the interacting variable, or a main effect across it.
    classify_main <- function(v, num, den) {
        partners <- unique(unlist(lapply(interaction_terms, function(vs) {
            if (v %in% vs && length(vs) == 2) setdiff(vs, v) else NULL
        })))
        if (length(partners)) {
            # stratified at the reference level of the (first) interacting variable
            sv <- partners[1]
            list(if (v == condition_var) 'condition_in_stratum' else 'other_in_stratum',
                 num, den, levels(col_data[[sv]])[1], sv, NA_character_)
        } else {
            list(if (v == condition_var) 'main_condition' else 'main_other',
                 num, den, NA_character_, NA_character_,
                 covariate_text(setdiff(factor_vars, v)))
        }
    }
    assign_class <- function(x) {
        type <<- x[[1]]; numerator <<- x[[2]]; denominator <<- x[[3]]
        stratum <<- x[[4]]; stratum_var <<- x[[5]]; covariate <<- x[[6]]
    }

    if (is.null(contrast)) {
        type <- 'default'
        numerator <- g2
        denominator <- g1

    } else if (is.character(contrast) && length(contrast) == 3L) {
        v <- contrast[[1]]; num <- contrast[[2]]; den <- contrast[[3]]
        contrast_echo <- paste0(v, ': ', num, ' vs ', den)
        if (v %in% factor_vars && all(c(num, den) %in% levels(col_data[[v]]))) {
            # honour a user-flipped contrast rather than assuming the level order
            assign_class(classify_main(v, num, den))
            if (type == 'main_condition' && is.na(covariate)) type <- 'default'
        }

    } else if (is.list(contrast) && length(contrast) == 1L) {
        nm <- as.character(contrast[[1]])
        contrast_echo <- paste(nm, collapse = ' + ')
        entries <- coef_map[nm[nm %in% names(coef_map)]]
        if (length(entries) == length(nm) && length(nm) %in% c(1L, 2L)) {
            kinds <- vapply(entries, function(e) e$kind, character(1))
            if (length(nm) == 1L && kinds == 'main') {
                e <- entries[[1]]
                assign_class(classify_main(e$var, e$num, e$den))
                if (type == 'main_condition' && is.na(covariate)) type <- 'default'
            } else if (length(nm) == 1L && kinds == 'interaction') {
                e <- entries[[1]]
                # the interaction is the difference in the condition effect between the two levels
                # of the other variable; orient it so 'condition' is the effect being compared
                if (e$var1 == condition_var) {
                    assign_class(list('interaction', e$lev2, levels(col_data[[e$var2]])[1],
                                      NA_character_, e$var2, NA_character_))
                } else {
                    assign_class(list('interaction', e$lev1, levels(col_data[[e$var1]])[1],
                                      NA_character_, e$var1, NA_character_))
                }
            } else if (length(nm) == 2L && setequal(kinds, c('main', 'interaction'))) {
                main <- entries[[which(kinds == 'main')]]
                inter <- entries[[which(kinds == 'interaction')]]
                # main effect + interaction = the main variable's effect in the non-reference
                # level of the interacting variable
                other <- if (inter$var1 == main$var) inter$var2 else inter$var1
                other_lev <- if (inter$var1 == main$var) inter$lev2 else inter$lev1
                if (main$var %in% c(inter$var1, inter$var2)) {
                    assign_class(list(
                        if (main$var == condition_var) 'condition_in_stratum' else 'other_in_stratum',
                        main$num, main$den, other_lev, other, NA_character_))
                }
            }
        }

    } else if (is.numeric(contrast)) {
        contrast_echo <- paste0('numeric contrast of length ', length(contrast))
    } else if (is.list(contrast) && length(contrast) >= 2L) {
        # DESeq2's listValues form: sum of the first set minus sum of the second
        contrast_echo <- paste0(paste(as.character(contrast[[1]]), collapse = ' + '),
                                ' minus ', paste(as.character(contrast[[2]]), collapse = ' + '))
    } else {
        contrast_echo <- paste0(class(contrast)[1], ' of length ', length(contrast))
    }

    if (type == 'unknown') {
        warning('Could not classify the DESeq2 contrast, so plot labels fall back to neutral, ',
                'direction-agnostic wording rather than guessing which group a positive log2 ',
                'fold change favours. Contrast: ', contrast_echo, call. = FALSE)
    }

    in_stratum <- if (is.na(stratum)) '' else paste0(' in ', stratum)
    paren_stratum <- if (is.na(stratum)) '' else paste0(' (', stratum, ')')
    controlling <- if (is.na(covariate)) '' else paste0(' (controlling for ', covariate, ')')

    two_group <- function(tag) list(
        subtitle = if (identical(tag, 'default')) NULL else
            paste0(numerator, ' vs ', denominator, controlling, in_stratum),
        volcano_x = paste0('Average log2 FC (', numerator, '/', denominator, ')', in_stratum),
        up_label = paste0('UP in ', numerator, paren_stratum),
        down_label = paste0('UP in ', denominator, paren_stratum),
        x_group_label = paste0(denominator, paren_stratum),
        y_group_label = paste0(numerator, paren_stratum),
        file_tag = if (identical(tag, 'default')) '' else
            paste0('_', numerator, '_vs_', denominator,
                   if (is.na(stratum)) '' else paste0('_in_', stratum))
    )

    labels <- switch(type,
        default = two_group('default'),
        main_condition = two_group('tagged'),
        main_other = two_group('tagged'),
        condition_in_stratum = two_group('tagged'),
        other_in_stratum = two_group('tagged'),
        interaction = list(
            subtitle = paste0(stratum_var, '-dependent ', g2, ' response: (', g2, '/', g1,
                              ') in ', numerator, ' vs in ', denominator),
            volcano_x = paste0('Difference in average log2 FC (', g2, '/', g1, '), ',
                               numerator, ' - ', denominator),
            # Deliberately about the ratio, not "more upregulated": a positive interaction
            # coefficient can equally mean a smaller decrease in the numerator level.
            up_label = paste0(g2, '/', g1, ' higher in ', numerator),
            down_label = paste0(g2, '/', g1, ' higher in ', denominator),
            x_group_label = g1, y_group_label = g2,
            file_tag = paste0('_', g2, 'vs', g1, '_', numerator, '_vs_', denominator,
                              '_interaction')),
        # unknown: never guess a direction
        list(
            subtitle = paste0('Custom contrast: ', contrast_echo),
            volcano_x = 'log2 Fold Change',
            up_label = 'Positive log2 FC',
            down_label = 'Negative log2 FC',
            x_group_label = g1, y_group_label = g2,
            file_tag = '_custom_contrast')
    )

    # Which Avg_ columns should the scatterplot draw? These are created by pseudobulk_de.
    avg <- function(...) paste0('Avg_', paste(c(...), collapse = '_'))
    cols <- if (type %in% c('default', 'main_condition', 'main_other')) {
        # keyed off numerator/denominator, not the level order, so a flipped contrast keeps the
        # axes and their labels in agreement
        c(avg(denominator), avg(numerator))
    } else if (type == 'condition_in_stratum' && !is.na(stratum)) {
        c(avg(denominator, stratum), avg(numerator, stratum))
    } else if (type == 'other_in_stratum' && !is.na(stratum)) {
        c(avg(stratum, denominator), avg(stratum, numerator))
    } else {
        # interaction / unknown: no meaningful pair, so fall back to the pooled columns
        c(avg(g1), avg(g2))
    }

    list(
        type = type,
        numerator = numerator,
        denominator = denominator,
        stratum = stratum,
        stratum_variable = stratum_var,
        covariate = covariate,
        subtitle = labels$subtitle,
        up_label = labels$up_label,
        down_label = labels$down_label,
        volcano_x = labels$volcano_x,
        x_column = cols[1],
        y_column = cols[2],
        x_group_label = labels$x_group_label,
        y_group_label = labels$y_group_label,
        scatter_applicable = type != 'interaction',
        file_tag = gsub('[^A-Za-z0-9._-]+', '.', labels$file_tag)
    )
}

#' Check a DESeq2 Contrast Names Real Coefficients
#'
#' Validates a \code{contrast} against a fitted DESeqDataSet before it reaches
#' \code{DESeq2::results}, so a mistyped coefficient name produces a message naming the available
#' coefficients instead of a bare DESeq2 error.
#'
#' @param contrast Character vector, list, numeric vector, or NULL; the contrast to check.
#' @param dds A DESeqDataSet that \code{DESeq()} has already been run on.
#'
#' @return NULL, invisibly. Called for its side effect of erroring on an invalid contrast.
#'
#' @noRd
validate_contrast_coefs <- function(contrast, dds) {
    if (is.null(contrast)) return(invisible(NULL))
    available <- DESeq2::resultsNames(dds)

    if (is.list(contrast)) {
        named <- unlist(lapply(contrast, as.character), use.names = FALSE)
        missing_coefs <- setdiff(named, available)
        if (length(missing_coefs)) {
            stop('contrast names ', paste0("'", missing_coefs, "'", collapse = ', '),
                 ', which ', if (length(missing_coefs) == 1) 'is not a coefficient' else
                     'are not coefficients',
                 ' of this fit. Available coefficients: ',
                 paste(available, collapse = ', '), call. = FALSE)
        }
    } else if (is.character(contrast) && length(contrast) == 3L) {
        col_data <- as.data.frame(SummarizedExperiment::colData(dds))
        v <- contrast[[1]]
        if (!v %in% names(col_data)) {
            stop("contrast refers to variable '", v, "', which is not a column of colData(dds). ",
                 'Available columns: ', paste(names(col_data), collapse = ', '), call. = FALSE)
        }
        lv <- levels(as.factor(col_data[[v]]))
        missing_levels <- setdiff(contrast[-1], lv)
        if (length(missing_levels)) {
            stop("contrast refers to level(s) ", paste0("'", missing_levels, "'", collapse = ', '),
                 " of '", v, "', which ", if (length(missing_levels) == 1) 'is' else 'are',
                 ' not present. Available levels: ', paste(lv, collapse = ', '), call. = FALSE)
        }
    }
    invisible(NULL)
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
plot_gene_lists <- function(results, gene_lists_to_plot, group1, group2,
                           cluster, local_figures_path,
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
                   cluster = cluster, 
                   local_figures_path = local_figures_path,
                   FC_threshold = FC_threshold, p_value_threshold = p_value_threshold,
                   max_overlaps = max_overlaps, label_size = label_size,
                   label_threshold = label_threshold,
                   distance_from_diagonal_threshold = distance_from_diagonal_threshold,
                   test_type = test_type, ...)

        volcano_plot(genes_to_plot = genes_to_plot, gene_list_name = gene_list,
                    results = results, group1 = group1, group2 = group2,
                    cluster = cluster, 
                    local_figures_path = local_figures_path,
                    FC_threshold = FC_threshold, p_value_threshold = p_value_threshold,
                    max_overlaps = max_overlaps, label_size = label_size,
                    label_threshold = label_threshold,
                    test_type = test_type, ...)
    }
}
