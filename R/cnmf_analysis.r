## cNMF Analysis Functions
## Functions for consensus Non-negative Matrix Factorization analysis

#' Run cNMF Gene Expression Program Analysis
#'
#' Performs comprehensive consensus Non-negative Matrix Factorization (cNMF) analysis
#' to identify gene expression programs (GEPs) in single-cell data. Optionally runs
#' cNMF consensus command, loads results, normalizes program usage scores, creates
#' visualizations (UMAP feature plots, violin plots, heatmaps), and performs pathway
#' enrichment analysis on top genes per program.
#'
#' @param seurat Seurat object with gene expression data
#' @param data_dir Character; path to directory containing cNMF run folder
#' @param runname Character; cNMF run name identifier
#' @param k_used Integer; number of gene expression programs (components) to use. Default 6
#' @param local_density_threshold Numeric; threshold for consensus clustering. Default 0.1
#' @param do_cnmf Logical; whether to run system cNMF consensus command. Default TRUE
#' @param sequential_palette Character vector; color palette for FeaturePlot. Default NULL
#' @param sequential_palette_dotplot Character vector; color palette for DotPlot. Default NULL
#' @param color_palette Character vector; color palette for gene loading heatmap. Default NULL
#' @param output_path Character; path to save output files. Default here::here('results', 'cNMF', runname, 'outputs')
#' @param object_annotations Character; string to append to output file names. Default ''
#' @param top_n Integer; number of top genes per program to extract. Default 100
#' @param topn_plot Integer; number of top genes per program for heatmap visualization. Default 6
#' @param run_pathway_enrichment Character vector; enrichment methods ('Metascape',
#'   'ClusterProfiler', or NULL). Default NULL
#' @param ... Additional arguments passed to pathway enrichment functions
#'
#' @return List with elements:
#'   \item{seurat}{Seurat object with cNMF program scores added to metadata}
#'   \item{top_colnames_long}{Long-format data frame of top genes with annotations}
#'   \item{top_colnames}{Wide-format data frame of top genes per program}
#'   \item{ORA_results}{Pathway overrepresentation analysis results or NULL}
#'   \item{genes_score_table}{Data frame of gene loadings for visualization}
#'   \item{spectra_score}{Matrix of gene spectra scores}
#'
#' @export
run_cnmf_results <- function (
    seurat,                    # Seurat object
    data_dir,                  # path to directory containing cNMF run folder (string)
    runname,                   # cNMF run name (string)
    k_used = 6,
    local_density_threshold = 0.1,
    do_cnmf = TRUE,            # whether to run system cnmf command
    sequential_palette = NULL, # palette for FeaturePlot_scCustom
    sequential_palette_dotplot = NULL, # palette for DotPlot_scCustom
    color_palette = NULL, # heatmap color palette
    output_path = here::here('results', 'cNMF', runname, 'outputs'),
    object_annotations = '',
    top_n = 100,                # top genes per program to extract
    topn_plot = 6,             # top genes used in dotplot
    run_pathway_enrichment = NULL,       # run Metascape_overrepresentation_analysis
    ...
) {

    dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

    # 1) optionally run cNMF consensus
    if (isTRUE(do_cnmf)) {
        cmd <- paste("cnmf consensus --output-dir", data_dir,
                                 "--name", runname,
                                 "--components", k_used,
                                 "--local-density-threshold", local_density_threshold,
                                 "--show-clustering", sep = " ")
        system(cmd)
    }

    dt_str <- gsub("\\.", "_", as.character(local_density_threshold))
    # 2) build file paths (match your original filenames)
    usage_file <- here::here(data_dir, runname, paste0(runname, ".usages", ".k_", k_used, ".dt_",  dt_str , ".consensus", ".txt"))
    spectra_score_file <- here::here(data_dir, runname, paste0(runname, ".gene_spectra_score", ".k_", k_used, ".dt_",  dt_str , ".txt"))
    spectra_tpm_file <- here::here(data_dir, runname, paste0(runname, ".gene_spectra_tpm", ".k_", k_used, ".dt_",  dt_str , ".txt"))

    # 3) read tables
    usage <- utils::read.table(usage_file, sep = "\t", row.names = 1, header = TRUE, check.names = FALSE)
    spectra_score <- utils::read.table(spectra_score_file, sep = "\t", row.names = 1, header = TRUE, check.names = FALSE)
    spectra_tpm <- utils::read.table(spectra_tpm_file, sep = "\t", row.names = 1, header = TRUE, check.names = FALSE)

    # 4) normalize usages (per-cell sums to 1)
    usage_norm <- as.data.frame(t(apply(usage, 1, function(x) x / sum(x))))

    # 5) attach normalized usages to seurat meta.data (preserve existing meta columns except X*)
    barcodes <- Seurat::Cells(seurat)
    seurat@meta.data <- seurat@meta.data %>%
        dplyr::select(-dplyr::starts_with('X')) %>%
        dplyr::mutate(barcode = barcodes) %>%
        dplyr::left_join(tibble::rownames_to_column(usage_norm, 'barcode'), by = 'barcode') %>%
        tibble::column_to_rownames('barcode')

    # 4) normalize usages (per-cell sums to 1)
    usage_norm <- as.data.frame(t(apply(usage, 1, function(x) x / sum(x))))
    colnames(usage_norm) <- paste0('cNMF_program_', seq_len(ncol(usage_norm)))


    barcodes <- Cells(seurat)

    seurat@meta.data <- seurat@meta.data |>
        select(-starts_with('cNMF_program_')) |>
        # select(-barcode) |>
        mutate(barcode = barcodes) |>
        left_join(usage_norm |> rownames_to_column('barcode'), by  = 'barcode') |>
        column_to_rownames('barcode')


    # 5) Feature plot for programs (relies on FeaturePlot_scCustom existing in environment)
    if (!is.null(sequential_palette)) {
        p1 <- FeaturePlot_scCustom(seurat, features = colnames(usage_norm), colors_use = sequential_palette, num_columns = 3)
    } else {
        p1 <- FeaturePlot_scCustom(seurat, features = colnames(usage_norm), num_columns = 3)
    }
    print(p1)
    ggsave(plot = p1, filename = paste0('FeaturePlot_cNMF_k_', k_used, '_', object_annotations, '.pdf'), path = output_path, width = 15, height = 10)


    signature_violin_plot <- function (signature) {
        plot1 <- VlnPlot(seurat, features = paste0(signature, ''), group.by = 'Groups', pt.size = 0) + labs(title = signature) + NoLegend()
        # print(plot1)
        ggsave(plot = plot1, filename = paste0('VlnPlot_', k_used, signature, '_by_group_', object_annotations, '.pdf'), path = output_path, width = 5, height = 4)
        return(plot1)
    }
    plot_list <- map(colnames(usage_norm), signature_violin_plot)
    plots <- wrap_plots(plot_list)
    print(plots)

    #### Interpreting Top genes per gene expression program

    # 8) Extracting top_n genes per program
    top_colnames <- spectra_score |>
        tibble::rownames_to_column("program") |>
        tidyr::pivot_longer(-program, names_to = "gene", values_to = "score") |>
        mutate(program = fct_inseq(program)) |>
        arrange(program) |>
        dplyr::group_by(program) |>
        filter(score > 0) |>
        dplyr::slice_max(order_by = score, n = top_n, with_ties = FALSE) |>
        dplyr::arrange(desc(score), .by_group = TRUE) |>
        dplyr::select(program, gene) |>
        dplyr::mutate(rank = dplyr::row_number()) |>
        ungroup() |>
        tidyr::pivot_wider(names_from = program, values_from = gene) |>
        dplyr::select(-rank)

    # spectra_score |>
    #     tibble::rownames_to_column("program") |>
    #     tidyr::pivot_longer(-program, names_to = "gene", values_to = "score") |>
    #     mutate(program = fct_inseq(program)) |>
    #     arrange(program) |>
    #     dplyr::group_by(program) |>
    #     filter(score > 0) |>
    #     select(program, gene) |>
    #     summarize(n = n()) |>
    #     print()


    #Add gene annotations:
    top_colnames_long <- top_colnames |>
        pivot_longer(everything(), names_to = 'gene_expression_program', values_to = 'gene')  |>
        mutate(gene_expression_program = fct_inseq(gene_expression_program)) |>
        arrange(gene_expression_program) |>
        left_join(y= unique(annotations[,c('gene_name', 'description')]), by = c('gene' = 'gene_name'))
    write.table(top_colnames_long,file=here::here(output_path, paste0('top50_per_gene_program_k_', k_used, ".tsv")), sep="\t",row.names = FALSE)

    top_colnames <- top_colnames |>
        select(levels(top_colnames_long$gene_expression_program))

        # 8.2) Extracting the top genes per program for heatmap plotting
    genes_score <- spectra_score |>
        as.matrix() |>
        t() |>
        as.data.frame() |>
        rownames_to_column('gene')  |>
        pivot_longer(-gene, names_to = 'GEP', values_to = 'loading')  |>
        mutate(GEP = fct_inseq(GEP)) |>
        group_by(GEP) |>
        # Z-score normalize the loadings within each GEP
        mutate(loading_zscore = loading) |>
        arrange(desc(loading_zscore), .by_group = TRUE) #|>
        # slice_max(order_by = loading_zscore, n = topn_plot, with_ties = FALSE)

    genes_to_plot <- genes_score |>
        group_by(GEP) |>
        slice_max(order_by = loading_zscore, n = topn_plot, with_ties = FALSE) |>
        ungroup() |>
        distinct(gene) |>
        pull(gene)

    genes_score <- genes_score  |>
        filter(gene %in% genes_to_plot)  |>
        mutate(gene = factor(gene, levels = genes_to_plot))

    # Calculate symmetric limits for the color scale based on z-scores
    max_abs <- max(abs(genes_score$loading_zscore), na.rm = TRUE)
    heatmap_plot <- ggplot(genes_score, aes(x = gene, y = GEP, fill = loading_zscore)) +
        geom_tile(linewidth = 0) +
        scale_fill_gradientn(
            colors = color_palette,
            name = "z-score<br>loading",
            limits = c(-max_abs, max_abs)
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
            axis.text.y = element_text(),
            plot.title = element_text(hjust = 0.5, size = 14)
        ) +
        labs(x = NULL, y = NULL, title = 'Top genes per\nGene Expression Program')
    print(heatmap_plot)

    # 10)  overrepresentation analysis per program
        local_path_pathway_enrichment <- here::here(output_path, paste0('pathway_enrichment_k_', k_used))
        dir.create(local_path_pathway_enrichment, recursive = TRUE, showWarnings = FALSE)
        over_representation_results <- NULL
        # Run pathway enrichment analysis
    if (!is.null(run_pathway_enrichment)) {
        if ('Metascape' %in% run_pathway_enrichment) {
            plot_list <- list()
            for (gene_expression_program in colnames(top_colnames)) {
                local_path_2 <- here::here(local_path_pathway_enrichment, gene_expression_program)
                dir.create(local_path_2, recursive = TRUE, showWarnings = FALSE)
                plot_list[[gene_expression_program]] <- Metascape_overrepresentation_analysis(top_colnames |> dplyr::pull(gene_expression_program),
                            local_path = local_path_2,
                            group = gene_expression_program,
                            filename = paste0(gene_expression_program, '_'),
                            grouping_var = 'GEP',
                            nterms_to_plot_metascape = 20, ...)
                }
            }
        plots <- wrap_plots(plot_list) &
            theme(plot.title = element_text(size = 9, hjust = 0.5),
                axis.text.y = element_text(size = 6),
                axis.text.x = element_text(size = 6),
                axis.title.x  = element_text(size = 7),
            axis.title.y  = element_text(size = 7))
        print(plots)

        if ('ClusterProfiler' %in% run_pathway_enrichment) {
            all_genes <- Features(seurat[['RNA']]) |>unique()
            gene_list <- top_colnames |> as.list() |>  map(~ .x[!is.na(.x)])
            over_representation_results <- GO_overrepresentation_analysis_multiple_lists(gene_list, all_genes, local_path_pathway_enrichment, ontology = 'ALL', minGSSize = 5, maxGSSize = 400, filename = '',  drop_levels = F, levels_to_drop = c(), simplify_function = min, simplify_by = 'p.adjust', simplify_terms = F, run_network = F, network_n_terms = 100, nterms_to_plot = 5, font_size = 8, ...)
        }
    }

    return(list(seurat = seurat, top_colnames_long = top_colnames_long, top_colnames=top_colnames,ORA_results = over_representation_results, genes_score_table = genes_score, spectra_score = spectra_score) )
    }
