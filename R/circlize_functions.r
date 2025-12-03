## TCR Circos Visualization Functions
## Functions for creating circular chord diagrams of TCR repertoire data

#' Create Color Palette for Grouping Variables
#'
#' Generates a named color vector for samples based on their group membership.
#' Automatically creates color shades for multiple samples within the same group.
#'
#' @param samples Character vector; sample names
#' @param groups Character vector; group names (same length as samples or unique groups)
#' @param base_colors Character vector; base colors to use for groups. Default colorblind-friendly palette
#'
#' @return Named character vector with colors for each sample
create_group_palette <- function(samples, groups, base_colors = c(
                "skyblue", "darkorange", "wheat", "darkolivegreen", "orchid",
                "seagreen", "gold", "lavenderblush", "tomato")) {
                # Helper to get shades for a base color
                get_shades <- function(color_name, n) {
                        shades <- c(                            
                                color_name,
                                paste0(color_name, "1"),
                                paste0(color_name, "2"),
                                paste0(color_name, "3"),
                                paste0(color_name, "4")
                        )
                        # If more than 5, use grDevices::colors() to find similar colors
                        if (n > 5) {
                                all_colors <- grDevices::colors()
                                similar <- grep(color_name, all_colors, value = TRUE)
                                # Remove already used shades
                                similar <- setdiff(similar, shades)
                                shades <- c(shades, head(similar, n - 5))
                        }
                        head(shades, n)
                }
                palette <- setNames(rep(NA, length(samples)), samples)
                for (grp in groups) {
                        grp_samples <- samples[grepl(grp, samples)]
                        n <- length(grp_samples)
                        # Pick base color for group (from base_colors, in order of groups)
                        idx <- which(groups == grp)
                        base_color <- base_colors[idx]
                        shades <- get_shades(base_color, n)
                        palette[grp_samples] <- shades[seq_along(grp_samples)]
                }
                return(palette)
        } 


plot_circos_clonotypes <- function(
    clonotype_data_plot,
    clonotype_data_plot_distinct,
    figures_path,
    grouping_variable,
    variables_to_color_by = NULL,
    file_type = "pdf",
    circle_margin = 0.2,
    major_ticks = 500,
    cex = 0.6,
    samples,
    groups,
    color_palette = NULL,
    alpha_col = 0.8,
    color_across_groups = NULL,
    base_colors = c(
        "skyblue", "darkorange", "wheat", "darkolivegreen", "orchid", 
        "seagreen", "gold", "lavenderblush", "tomato")
) {
    cluster <- 'all'
    if (is.null(color_palette)) {
        grid_cols <- create_group_palette(samples, groups, base_colors)
    } else {
        grid_cols <- color_palette
    }

    file_name <- paste0(figures_path, 'Circos_clonotypes_per_', grouping_variable, '.', file_type)
    if (file_type == "pdf") {
        grDevices::pdf(file = file_name)
    } else if (file_type == "png") {
        grDevices::png(file = file_name, width = 15, height = 15, units = 'cm', res = 300)
    } else {
        stop("file_type must be either 'pdf' or 'png'")
    }
    circlize::circos.clear()
    circlize::circos.par(gap.degree = 2, track.height = 0.1, cell.padding = c(0, 0, 0, 0), circle.margin = circle_margin)
    circlize::circos.initialize(xlim = clonotype_data_plot_distinct)

    circlize::circos.track(ylim = c(0,1),
        panel.fun = function(x, y) {
            print(CELL_META$xrange[[1]])
            if (CELL_META$cell.width < 45) {
                circlize::circos.text(
                    CELL_META$xcenter, 
                    CELL_META$cell.ylim[2] + circlize::mm_y(9), 
                    adj = c(0 , 1),
                    CELL_META$sector.index,
                    facing = 'clockwise', 
                    niceFacing = TRUE, 
                    cex  = cex
                )
            } else {
                circlize::circos.text(
                    CELL_META$xcenter, 
                    CELL_META$cell.ylim[2] + circlize::mm_y(11),
                    CELL_META$sector.index,
                    facing = 'bending.inside', 
                    niceFacing = TRUE,
                    cex = cex + 0.2
                )
            }
            if (CELL_META$xrange[[1]] > major_ticks * 2) {
                circlize::circos.axis(
                    labels.cex = 0.5, 
                    minor.ticks =  0,
                    major.tick = 1,
                    labels.facing = 'clockwise',
                    major.at = seq(major_ticks, CELL_META$xrange[[1]], by = major_ticks)
                )
            }
            circlize::highlight.sector(CELL_META$sector.index, col = grid_cols[CELL_META$sector.index])
        }
    )

    done <- c()
    for (origin in levels(clonotype_data_plot$group)) {
        for (target in levels(clonotype_data_plot$group)) {
            if (origin == target | target %in% done) {
                next
            } else {
                table_one <- clonotype_data_plot |>
                    filter(group == origin & group_counts != 0)  |>
                    dplyr::select(all_of(c('clonotype', 'group', 'coordinates')))
                table_two <- clonotype_data_plot |>
                    filter(group == target & group_counts != 0) |>
                    dplyr::select(all_of(c('clonotype', 'group', 'coordinates')))
                link_table  <-  inner_join(table_one, table_two, by = 'clonotype') |> column_to_rownames(var = 'clonotype')
                for (clonotype1 in rownames(link_table)) {
                    group1 <- link_table[[clonotype1, 'group.x']]
                    group2 <- link_table[[clonotype1, 'group.y']]
                    coordinates1 <- as.vector(link_table[[clonotype1, 'coordinates.x']])
                    coordinates2 <- as.vector(link_table[[clonotype1, 'coordinates.y']])
                    if (!is.null(variables_to_color_by)) {
                        for (variable in variables_to_color_by) {
                            if (str_detect(origin, variable) | str_detect(target, variable)) {
                                color = alpha(grid_cols[variable], alpha_col)
                                break
                            } else {
                                color = alpha("gray", 0.25)                                                               
                            }                        
                        }
                    } else if (!is.null(color_across_groups)) {
                        across_groups = F
                        for (combination in color_across_groups) {
                            if ((str_detect(origin, combination[1]) & str_detect(target, combination[2])) | (str_detect(origin, combination[2]) & str_detect(target, combination[1])) ) {
                                across_groups = T
                            }
                        }
                        if (across_groups) {
                            color = alpha("purple", alpha_col)                                                               
                        } else {
                            color = alpha(grid_cols[origin], alpha_col)                                                               
                        }
                        
                    }else {
                        color = alpha(grid_cols[origin], alpha_col)                                                               
                    }
                    circos.link(
                        origin, 
                        coordinates1,
                        target,
                        coordinates2,
                        col = color
                    )
                }
            }
        }
        done <- c(done, origin)
    }
    title(paste0('Clonotype Overlap ', cluster, ' per ', grouping_variable))
    grDevices::dev.off()
    circlize::circos.clear()
    return()
}

#' Create TCR Clonotype Overlap Circos Diagram and Summary Tables
#'
#' Generates circular chord diagrams showing TCR clonotype overlap across samples/groups
#' and creates summary tables of clonotype distribution. Uses circlize package to
#' visualize shared clonotypes between conditions with colored links.
#'
#' @param seurat Seurat object with scRepertoire TCR data
#' @param grouping_variable Character; metadata column for grouping samples. Default 'Samples'
#' @param results_path Character; base results directory. Default './'
#' @param variables_to_color_by Character vector; specific variables to highlight with color. Default NULL
#' @param cell_types_column Character; metadata column with cell type labels. Default 'cell_types'
#' @param write_table Logical; whether to write summary table. Default TRUE
#' @param figures_path Character; path to save circos plots
#' @param tables_path Character; path to save summary tables
#' @param circle_margin Numeric; margin around circle. Default 0.2
#' @param major_ticks Integer; spacing for axis ticks. Default 500
#' @param cex Numeric; text size for labels. Default 0.6
#' @param samples Character vector; sample names
#' @param groups Character vector; group assignments for samples
#' @param color_palette Named character vector; custom color palette. If NULL, creates automatically. Default NULL
#' @param alpha_col Numeric; transparency for links (0-1). Default 0.8
#' @param color_across_groups List; group pairs to highlight with special color. Default NULL
#' @param base_colors Character vector; base colors for automatic palette. Default colorblind-friendly palette
#'
#' @return NULL (invisibly); creates PDF and PNG circos plots and CSV table
#'
#' @examples
#' overlap_circos_and_tables(seurat_tcr, grouping_variable = "Samples",
#'                          figures_path = "./figures", tables_path = "./tables",
#'                          samples = c("WT1", "WT2", "KO1", "KO2"),
#'                          groups = c("WT", "WT", "KO", "KO"))
#' @export
overlap_circos_and_tables <- function(
    seurat,
    grouping_variable = 'Samples',
    results_path = './',
    variables_to_color_by = NULL,
    cell_types_column = 'cell_types',
    write_table = TRUE,
    figures_path,
    tables_path,
    circle_margin = 0.2,
    major_ticks = 500,
    cex = 0.6,
    samples,
    groups,
    color_palette = NULL,
    alpha_col = 0.8,
    color_across_groups = NULL,
    base_colors = c(
        "skyblue", "darkorange", "wheat", "darkolivegreen", "orchid",
        "seagreen", "gold", "lavenderblush", "tomato")
) {
    Mode <- function(x) {
        ux <- unique(x)
        ux[which.max(tabulate(match(x, ux)))]
    }
    
    combined2 <- scRepertoire:::.expression2List(seurat, split.by = 'orig.ident')
    
    clonotype_data <- combined2[[1]] |>
        as_tibble() |>
        dplyr::select(all_of(c('CTaa', 'CTgene', grouping_variable, cell_types_column))) |>
        add_count(CTaa, !!as.name(grouping_variable), sort = TRUE, name = 'counts_per_condition') |>
        group_by(CTaa, !!as.name(grouping_variable)) |>
        summarize(across(everything(), Mode), .groups = 'drop') |>
        arrange(desc(counts_per_condition)) |>
        pivot_wider(
            id_cols = c('CTaa'),
            names_from = !!as.name(grouping_variable),
            values_from = counts_per_condition,
            unused_fn = Mode,
            values_fill = 0
        ) |>
        ungroup() |>
        mutate(clonotype = paste0('clonotype ', as.character(row_number())))
    
    if (write_table) {
        readr::write_csv(clonotype_data, file = here::here(tables_path, 'TCR_data_per_group.csv'))
    }
    
    group_levels <- levels(pull(seurat@meta.data, grouping_variable))
    
    clonotype_data_plot <- clonotype_data |>
        pivot_longer(
            c(group_levels),
            names_to = 'group',
            values_to = 'group_counts'
        ) |>
        group_by(group) |>
        mutate(sequence_count_by_grouping_variable = sum(group_counts)) |>
        arrange(desc(group_counts)) |>
        mutate(clonotype_position_on_circos_by_grouping_variable = cumsum(group_counts)) |>
        mutate(clonotype_start_position_by_grouping_variable = c(
            0,
            clonotype_position_on_circos_by_grouping_variable[-length(clonotype_position_on_circos_by_grouping_variable)]
        )) |>
        ungroup() |>
        rowwise() |>
        mutate(coordinates = list(c(
            clonotype_start_position_by_grouping_variable,
            clonotype_position_on_circos_by_grouping_variable
        ))) |>
        ungroup() |>
        mutate(group = fct(group, levels = group_levels))
    
    clonotype_data_plot_distinct <- clonotype_data_plot |>
        dplyr::select(all_of(c('group', 'sequence_count_by_grouping_variable'))) |>
        distinct() |>
        arrange(group) |>
        mutate(origin = 0) |>
        column_to_rownames('group') |>
        relocate(origin) |>
        as.matrix()
    
    plot_circos_clonotypes(
        clonotype_data_plot = clonotype_data_plot,
        clonotype_data_plot_distinct = clonotype_data_plot_distinct,
        figures_path = figures_path,
        grouping_variable = grouping_variable,
        variables_to_color_by = variables_to_color_by,
        file_type = 'pdf',
        circle_margin = circle_margin,
        major_ticks = major_ticks,
        cex = cex,
        samples = samples,
        groups = groups,
        alpha_col = alpha_col,
        base_colors = base_colors,
        color_palette = color_palette,
        color_across_groups = color_across_groups
    )

    plot_circos_clonotypes(
        clonotype_data_plot = clonotype_data_plot,
        clonotype_data_plot_distinct = clonotype_data_plot_distinct,
        figures_path = figures_path,
        grouping_variable = grouping_variable,
        variables_to_color_by = variables_to_color_by,
        file_type = 'png',
        circle_margin = circle_margin,
        major_ticks = major_ticks,
        cex = cex,
        samples = samples,
        groups = groups,
        base_colors = base_colors,
        color_palette = color_palette,
        alpha_col = alpha_col,
        color_across_groups = color_across_groups
    )

    return()
}

# Helper function to render chord diagram
render_chord_diagram <- function(trbv_j_processed, grid_colors, colored_function,
                                  trbv_j_transparency, chain_a, chain_b,
                                  group_name, cex) {
    # Create chord diagram
    circlize::chordDiagram(trbv_j_processed,
        grid.col = grid_colors,
        col = colored_function,
        link.zindex = rank(trbv_j_processed$interaction),
        transparency = trbv_j_transparency$transparency,
        annotationTrack = c("grid", 'axis'),
        preAllocateTracks = 1)

    # Add labels to sectors
    circlize::circos.track(track.index = 1,
        bg.border = NA,
        panel.fun = function(x, y) {
            if (CELL_META$cell.width < 365) {
                circlize::circos.text(
                    CELL_META$xcenter,
                    CELL_META$ylim[1] + circlize::mm_y(3),
                    adj = c(0, 1),
                    CELL_META$sector.index,
                    facing = 'clockwise',
                    niceFacing = T,
                    cex = cex
                )
            } else {
                circlize::circos.text(
                    CELL_META$xcenter,
                    CELL_META$ylim[1] + circlize::mm_y(6),
                    CELL_META$sector.index,
                    facing = 'bending.inside',
                    niceFacing = T,
                    cex = cex + 0.2
                )
            }
        }
    )

    title(paste0(chain_a, ' - ', chain_b, ' gene pairings in ', group_name))
}

#' Create TCR V-J Gene Pairing Chord Diagram
#'
#' Generates a circular chord diagram showing the pairing frequency between TCR
#' V-region and J-region genes. Link color and transparency represent interaction
#' strength. Commonly used to visualize TRBV-TRBJ or TRAV-TRAJ pairings.
#'
#' @param trbv_j_table Data frame or matrix; V-J gene pairing frequency table with
#'   rownames as "V_J" pairs and column with interaction counts
#' @param figures_path Character; path to save chord diagram PDF
#' @param chain_a Character; name of first chain (V-region). Default 'TRBV'
#' @param chain_b Character; name of second chain (J-region). Default 'TRBJ'
#' @param group_name Character; sample/group identifier for plot title and filename
#' @param trbv_palette Character; HCL color palette name for V genes. Default "Teal"
#' @param trbj_palette Character; HCL color palette name for J genes. Default "Emrl"
#' @param link_palette Character; HCL color palette name for links. Default "YlOrRd"
#' @param cex Numeric; text size for gene labels. Default 0.7
#' @param major_ticks Integer; spacing for axis ticks. Default 10
#'
#' @return Recorded plot object of the chord diagram
#'
#' @examples
#' create_V_J_chord_diagram(vj_table, figures_path = "./figures",
#'                         chain_a = "TRBV", chain_b = "TRBJ",
#'                         group_name = "WT_sample")
#' @export
create_V_J_chord_diagram <- function(trbv_j_table,
                                        figures_path,
                                        chain_a = 'TRBV',
                                        chain_b = 'TRBJ',
                                        group_name,
                                        trbv_palette = "Teal",
                                        trbj_palette = "Emrl",
                                        link_palette = "YlOrRd",
                                        cex = 0.7,
                                        major_ticks = 10) {


    chain_pair <- paste0(chain_a, '_', chain_b)

    # Prepare the data
    trbv_j_processed <- trbv_j_table |>
        as.data.frame() |>
        rownames_to_column(var = chain_pair) |>
        select(c(all_of(chain_pair), all_of(group_name))) |>
        rename(interaction = all_of(group_name)) |>
        separate_wider_delim(all_of(chain_pair), delim = '_', names = c(chain_a, chain_b))

    # Calculate transparency
    trbv_j_transparency <- trbv_j_processed |>
        mutate(transparency = (interaction - min(interaction)) / (max(interaction) - min(interaction))) |>
        mutate(transparency = 1 - transparency * 1.5)

    # Generate colors
    trbv_colors <- grDevices::hcl.colors(length(unique(trbv_j_processed[[chain_a]])), trbv_palette, rev = TRUE)
    trbj_colors <- grDevices::hcl.colors(length(unique(trbv_j_processed[[chain_b]])), trbj_palette, rev = TRUE)
    grid_colors <- c(trbv_colors, trbj_colors)

    # Create color function for links
    colored_function <- colorRamp2::colorRamp2(range(trbv_j_processed$interaction), hcl_palette = link_palette, reverse = TRUE)
    
    # Save to PDF
    grDevices::pdf(file = here::here(figures_path, paste0(chain_pair, '_chord_diagram_', group_name, '.pdf')))
    render_chord_diagram(trbv_j_processed, grid_colors, colored_function,
                        trbv_j_transparency, chain_a, chain_b, group_name, cex)
    grDevices::dev.off()
    circlize::circos.clear()

    # Create plot for return
    render_chord_diagram(trbv_j_processed, grid_colors, colored_function,
                        trbv_j_transparency, chain_a, chain_b, group_name, cex)
    p1 <- recordPlot()
    circlize::circos.clear()
    return(p1)
}
