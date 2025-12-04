## TCR Repertoire Analysis Functions
## Functions for analyzing T cell receptor repertoire diversity and clonality

#' Visualize TCR Clonal Repertoire Using Circle Packing
#'
#' Creates a circle packing visualization where each circle represents a unique
#' T cell clonotype, sized by clonal frequency. Colors indicate clone size category
#' (small, medium, large, hyperexpanded). Requires scRepertoire-processed TCR data
#' with CTaa and clonalFrequency metadata.
#'
#' @param seurat Seurat object with TCR repertoire metadata (CTaa, clonalFrequency,
#'   cloneSize, Groups)
#' @param group_name Character; value in 'Groups' metadata column to visualize
#' @param figures_path Character; path to save the circle packing plot
#' @param viridis_option Character; viridis color palette option ('A', 'B', 'C', 'D',
#'   'E'). Default 'D'
#' @param label_size_range Numeric vector of length 2; range for scaling clonotype
#'   frequency labels. Default c(3, 10)
#'
#' @return ggplot object with circle packing visualization
#'
#' @export
plot_clonal_repertoire <- function(seurat, group_name, figures_path, viridis_option = 'D', label_size_range = c(3,10)) {
  # Prepare repertoire data
  repertoire <- seurat@meta.data |>
    filter(Groups == group_name & !is.na(CTaa)) |>
    arrange(desc(clonalFrequency))

  # Determine unique clone sizes and generate appropriate number of colors
  repertoire2 <- repertoire |> arrange(cloneSize)
  clone_sizes <- unique(repertoire2$cloneSize[!is.na(repertoire$cloneSize)]) |> as.character()
  n_colors <- length(clone_sizes)
  viridis_colors <- viridis::viridis(n = n_colors, option = viridis_option, direction = 1)
  print(clone_sizes)

  # Create color mapping
  repertoire <- repertoire |>
    mutate(color = case_when(
      str_equal(cloneSize, clone_sizes[1]) ~ viridis_colors[1],
      str_equal(cloneSize, clone_sizes[min(2, n_colors)]) ~ viridis_colors[min(2, n_colors)],
      str_equal(cloneSize, clone_sizes[min(3, n_colors)]) ~ viridis_colors[min(3, n_colors)],
      str_equal(cloneSize, clone_sizes[min(4, n_colors)]) ~ viridis_colors[min(4, n_colors)],
      str_equal(cloneSize, clone_sizes[min(5, n_colors)]) ~ viridis_colors[n_colors],
      is.na(cloneSize) ~ 'grey90'
    )) |>
    select(clonalFrequency, CTaa, color) |>
    distinct()

  # Calculate circle packing layout
  packing <- packcircles::circleProgressiveLayout(repertoire$clonalFrequency, sizetype = 'area')
  repertoire <- repertoire |>
    mutate(x = packing$x,
           y = packing$y)
  ggplot_data_circles <- packcircles::circleLayoutVertices(packing, npoints = 50)

  # Create plot
  plot <- ggplot(data = ggplot_data_circles) +
    geom_polygon(aes(x, y, group = id, fill = factor(id)),
                 colour = "black", linewidth = 0.1, alpha = 0.7, show.legend = FALSE) +
    scale_fill_manual(values = repertoire$color, labels = repertoire$cloneSize) +
    geom_text(data = repertoire,
              aes(x, y, size = clonalFrequency, label = clonalFrequency)) +
    guides(fill = guide_legend(title = 'Clonotype size'), size = "none") +
    scale_size_continuous(range = label_size_range) +
    labs(title = paste0("Clonal repertoire for group: ", group_name)) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 20)) +
    coord_equal()
  ggsave(plot = plot, filename = paste0('clonal_repertoire_', group_name, '.pdf'), width = 8, height = 8, path = figures_path)

#print(plot)
return(plot)
}


#' Calculate D50 TCR Diversity Metric
#'
#' Computes the D50 diversity index for T cell receptor repertoires across cell types
#' and replicates. D50 measures the proportion of unique clonotypes needed to account
#' for the most abundant 50% of cells. Higher D50 indicates greater clonal diversity.
#' Requires at least 20 cells per group for reliable calculation.
#'
#' @param seurat Seurat object with TCR repertoire metadata from scRepertoire
#' @param cell_grouping_var Unquoted column name; cell type or cluster grouping variable
#' @param replicate_var Unquoted column name; biological replicate identifier (e.g., sample ID)
#' @param replicate_group_var Unquoted column name; grouping of replicates for plotting
#'   (e.g., experimental condition). If NULL, uses replicate_var. Default NULL
#' @param results_path Character; path for saving DESeq2::results (currently unused)
#' @param figures_path Character; path to save D50 bar plot
#' @param tables_path Character; path to save D50 table CSV
#' @param colors Named character vector; colors for replicate_group_var levels. If NULL,
#'   uses ggplot2 defaults. Default NULL
#'
#' @return ggplot object showing D50 values per cell type and replicate group with
#'   mean bars and error bars
#'
#' @export
calculate_D50 <- function (seurat, cell_grouping_var, replicate_var, replicate_group_var = NULL, results_path, figures_path, tables_path, colors  = NULL) {

    #Extracting TCR data for clusters of interest
    Idents(seurat) <- englue("{{cell_grouping_var}}")
    combined2 <- scRepertoire:::.expression2List(seurat, split.by ='ident')
    combined3 <- scRepertoire:::.expression2List(seurat, split.by ='orig.ident')
    grouping_var_levels <- levels(seurat@meta.data |> pull({{cell_grouping_var}}))
    replicate_var_levels <- levels(seurat@meta.data |> pull({{replicate_var}}))

    #Initiating results data frame
    results <- as.data.frame(matrix(nrow = 0,ncol = length(grouping_var_levels)))
    # colnames(results)
    rnames <- c()

    #Calculate D50
    for (HTO in replicate_var_levels) {

        result <- c()

        cell_type_HTO_data <- combined3[[1]]|>
                                        filter({{replicate_var}} == HTO) |>
                                        add_count(CTaa, sort=TRUE)
            #Calculating D50
            if (nrow(cell_type_HTO_data) < 20) {
                D50 <- NA
            } else {
                L50 <- floor(nrow(cell_type_HTO_data)/2)
                number_unique_50 <- cell_type_HTO_data[1:L50,] |> summarise(n_distinct(CTaa)) |> as.numeric()
                number_unique_total <- cell_type_HTO_data[] |> summarise(n_distinct(CTaa)) |> as.numeric()
                D50 <- number_unique_50/number_unique_total
            }
            result <- c(result, D50)


        for (cell_type in grouping_var_levels) {
            #Extracting data for cell type and HTO
            cell_type_HTO_data <- combined2[[cell_type]]|>
                                        filter({{replicate_var}} == HTO) |>
                                        add_count(CTaa, sort=TRUE)

            #Calculating D50
            if (nrow(cell_type_HTO_data) < 20) {
                D50 <- NA
            } else {
                L50 <- floor(nrow(cell_type_HTO_data)/2)
                number_unique_50 <- cell_type_HTO_data[1:L50,] |> summarise(n_distinct(CTaa)) |> as.numeric()
                number_unique_total <- cell_type_HTO_data[] |> summarise(n_distinct(CTaa)) |> as.numeric()
                D50 <- number_unique_50/number_unique_total
            }
            result <- c(result, D50)
        }
        results <- rbind(results, result)
        rnames <- c(rnames, HTO)
        # print(HTO)

    }
    colnames(results) <- paste0(c('All', grouping_var_levels), '_D50')
    results <- results |> mutate({{replicate_var}} := rnames) |> arrange({{replicate_var}}) |> relocate({{replicate_var}})

    write.csv(results, file = englue("{tables_path}/D50_per_{{cell_grouping_var}}.csv"), row.names=FALSE)
    print(head(results))

    # Convert results to long format for plotting
    results_long <- results |>
        pivot_longer(-{{replicate_var}}, names_to = englue("{{cell_grouping_var}}"), values_to = "D50") |>
        mutate({{cell_grouping_var}} := fct_inorder({{cell_grouping_var}}))



# Add replicate grouping variable information
    if (!(englue('{{replicate_var}}') == englue('{{replicate_group_var}}'))) {
        replicate_grouping_var_info <- seurat@meta.data |>
            dplyr::select({{replicate_var}}, {{replicate_group_var}})
        replicate_grouping_var_info <- replicate_grouping_var_info[!duplicated(replicate_grouping_var_info), ]
        results_long <- results_long |>
            left_join(replicate_grouping_var_info, by = join_by({{replicate_var}} == {{replicate_var}}))
    }


    # Remove NA values for plotting
    results_long <- results_long |> filter(!is.na(D50))

    # Plot: x axis is {{replicate_var}}, show only dots (no columns)
    plot1 <- ggplot(results_long, aes(x = {{replicate_group_var}}, y = D50, color = {{replicate_group_var}})) +
        geom_quasirandom(size = 1.5) +
        stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.7), width = 0.5, aes(fill = {{replicate_group_var}}, alpha = 0.5)) +
        stat_summary(fun.data = mean_se, fun.args = list(mult = 1),
            geom = "errorbar", width = 0.2, position = position_dodge(width = 0.7)) +
        scale_y_continuous(limits = c(0, 0.5), expand = expansion(mult = c(0, 0.05))) +
        theme_classic() +
        labs(x = englue("{{replicate_group_var}}"), y = "D50 Diversity", title = englue("D50 per {{replicate_group_var}} by {{cell_grouping_var}}")) +
        guides(alpha = 'none')+
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
            text = element_text(size = 14),
            axis.line = element_line(colour = "black"))
    if (!is.null(colors)) {
        plot1 <- plot1 + scale_color_manual(values = colors) + scale_fill_manual(values = colors)
    }
    ggsave(plot = plot1, filename = paste0(englue('D50_per_{{cell_grouping_var}}.pdf')), width = length(levels(results_long |> pull({{cell_grouping_var}})))+2, height = 6, path = figures_path)
    print(plot1)
    return(plot1)
}
