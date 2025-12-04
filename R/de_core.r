## Core Differential Expression Analysis Functions
## Functions that perform the statistical analysis without visualization

#' Pseudobulk Differential Expression Analysis
#'
#' Performs differential expression analysis using DESeq2 on aggregated pseudobulk
#' counts from single-cell data. Aggregates counts by biological replicate and
#' uses DESeq2's negative binomial model for statistical testing. Includes PCA
#' visualization for quality control.
#'
#' @param scRNAseq Seurat object containing single-cell RNA-seq data
#' @param comparison Character; metadata column name for comparison groups
#' @param group1 Character; pattern to match for first group (control)
#' @param group2 Character; pattern to match for second group (treatment)
#' @param cluster Character; cluster or cell type identifier. Default 'all_clusters'
#' @param path Character; output directory path. Default './'
#' @param FC_threshold Numeric; log2 fold change threshold for filtering. Default 0.3
#' @param p_value_threshold Numeric; adjusted p-value threshold. Default 0.05
#' @param expression_threshold_for_gene_list Numeric; minimum average normalized counts
#'   for a gene to be included in filtered gene lists. Default 20
#' @param minimum_cell_number Integer; minimum cells required per group. Default 10
#' @param genes_to_exclude Character vector; genes to exclude from analysis. Default c()
#' @param ... Additional arguments (unused)
#'
#' @return List with elements:
#'   \item{all_count}{Integer; number of significant DEGs}
#'   \item{UP_count}{Integer; number of upregulated genes}
#'   \item{DOWN_count}{Integer; number of downregulated genes}
#'   \item{results}{Data frame; complete DE results with gene annotations and counts}
#'
#' @export
pseudobulk_de <- function(scRNAseq, comparison, group1, group2, cluster = 'all_clusters',
                          path = './', FC_threshold = 0.3, p_value_threshold = 0.05,
                          expression_threshold_for_gene_list = 20, minimum_cell_number = 10,
                          genes_to_exclude = c(), ...) {

    # Subset seurat object
    scRNAseq <- subset(scRNAseq, subset = (str_detect(!!as.name(comparison), group1) | str_detect(!!as.name(comparison), group2)))

    # Set Paths
    gene_lists_path <- here::here(path, 'gene_lists')
    dir.create(gene_lists_path, showWarnings = FALSE, recursive = TRUE)

    print(paste('Cluster', cluster))

    group1 <- stringr::fixed(group1)
    group2 <- stringr::fixed(group2)

    Idents(scRNAseq) <- comparison

    # Check there are enough cells
    cell_validation <- validate_cell_counts(scRNAseq, comparison, group1, group2, minimum_cell_number)

    print('number of cells in group 1')
    print(cell_validation$n_group1)
    print('number of cells in group 2')
    print(cell_validation$n_group2)

    if (!cell_validation$valid) {
        return(list(
            all_count = cell_validation$message,
            UP_count = cell_validation$message,
            DOWN_count = cell_validation$message,
            results = NULL
        ))
    }

    # Aggregate counts
    counts <- AggregateExpression(scRNAseq, group.by = c(comparison),
                                 assays = 'RNA',
                                 slot = 'counts',
                                 return.seurat = FALSE)

    counts <- counts$RNA |>
        as.data.frame() |>
        rownames_to_column('genes') |>
        as_tibble() |>
        dplyr::select(-any_of(genes_to_exclude)) |>
        column_to_rownames('genes')

    # Generate sample level metadata
    colData <- data.frame(samples = colnames(counts)) |>
        mutate(condition = ifelse(grepl(group1, samples), group1, group2))

    # Filter
    counts <- counts |>
        mutate(row_sums = rowSums(counts)) |>
        filter(row_sums >= 10) |>
        dplyr::select(-row_sums)

    print('Group 1 Length')
    print(nrow(colData |> filter(condition == group1)))
    print('Group 2 Length')
    print(nrow(colData |> filter(condition == group2)))

    # Check for sufficient replicates
    if ((length(unique(colData$condition)) != 2) |
        (nrow(colData |> filter(condition == group1)) < 2) |
        (nrow(colData |> filter(condition == group2)) < 2)) {
        return(list(
            all_count = 'Not enough biological replicates per group',
            UP_count = 'Not enough biological replicates per group',
            DOWN_count = 'Not enough biological replicates per group',
            results = NULL
        ))
    }

    # Create DESeq2 object
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                  colData = colData,
                                  design = ~condition)
    dds$condition <- factor(dds$condition, levels = c(group1, group2))

    # DESeq2 QC
    rld <- DESeq2::rlog(dds, blind = TRUE)
    local_figures_path <- here::here(path, 'figures')
    dir.create(local_figures_path, showWarnings = FALSE, recursive = TRUE)

    DESeq2::plotPCA(rld, ntop = 500, intgroup = 'condition')
    ggsave(filename = paste0('Pseudobulk_PCA_', cluster, '.pdf'), path = local_figures_path)

    PCA_table <- DESeq2::plotPCA(rld, ntop = 500, intgroup = 'condition', returnData = TRUE)
    write.csv(PCA_table, file = here::here(path, paste('PCA_pseudobulk', cluster, group2, 'vs', group1, '.csv', sep = '_')))

    # Run DESeq2
    dds <- DESeq2::DESeq(dds)
    DESeq2::resultsNames(dds)

    # Generate results object
    results <- DESeq2::results(dds) |> as.data.frame()

    # Get Normalized Counts
    normalized_counts <- DESeq2::counts(dds, normalized = TRUE)
    normalized_counts <- normalized_counts |>
        as.data.frame() |>
        rownames_to_column('genes') |>
        as_tibble() |>
        rowwise() |>
        mutate(
            !!paste0('Avg_', group2) := mean(c_across(contains(group2))),
            !!paste0('Avg_', group1) := mean(c_across(contains(group1)))
        ) |>
        ungroup()

    # Add gene annotations
    results <- results |>
        rownames_to_column('genes') |>
        left_join(y = unique(annotations[, c('gene_name', 'description')]),
                 by = c('genes' = 'gene_name')) |>
        left_join(y = normalized_counts, by = c('genes' = 'genes'))

    # Filter results
    results_filtered <- filter(
        results,
        padj < p_value_threshold &
            ((!!sym(paste0('Avg_', group2)) > expression_threshold_for_gene_list) |
             (!!sym(paste0('Avg_', group1)) > expression_threshold_for_gene_list)) &
            (log2FoldChange >= FC_threshold | log2FoldChange <= -1 * FC_threshold)
    ) |>
        arrange(padj)

    results_filtered_UP <- filter(results_filtered, log2FoldChange >= FC_threshold)
    results_filtered_DOWN <- filter(results_filtered, log2FoldChange <= -1 * FC_threshold)

    # Write results to CSV files
    write.csv(results |> arrange(padj),
             file = here::here(gene_lists_path, paste('ALL_GENES_DEG_Analysis', cluster, 'pseudobulk', group2, 'vs', group1, '.csv', sep = '_')))
    write.csv(results_filtered_UP |> arrange(desc(log2FoldChange)),
             file = here::here(gene_lists_path, paste('DEG_UP_in', group2, cluster, 'pseudobulk', group2, 'vs', group1, '.csv', sep = '_')))
    write.csv(results_filtered_DOWN |> arrange(log2FoldChange),
             file = here::here(gene_lists_path, paste('DEG_DOWN_in', group2, cluster, 'pseudobulk', group2, 'vs', group1, '.csv', sep = '_')))

    # Return counts and results
    list(
        all_count = nrow(results_filtered),
        UP_count = nrow(results_filtered_UP),
        DOWN_count = nrow(results_filtered_DOWN),
        results = results
    )
}


#' Wilcoxon Differential Expression Analysis for RNA Assay
#'
#' Performs differential expression analysis using Seurat's FindMarkers function
#' with Wilcoxon rank-sum test on the RNA assay. Calculates counts per million (CPM)
#' for visualization and filtering. Suitable for non-normalized or log-normalized data.
#'
#' @param scRNAseq Seurat object containing single-cell RNA-seq data
#' @param comparison Character; metadata column name for comparison groups
#' @param group1 Character; pattern to match for first group (control)
#' @param group2 Character; pattern to match for second group (treatment)
#' @param cluster Character; cluster or cell type identifier. Default 'all_clusters'
#' @param path Character; output directory path. Default './'
#' @param FC_threshold Numeric; log2 fold change threshold for filtering. Default 0.3
#' @param p_value_threshold Numeric; adjusted p-value threshold. Default 0.05
#' @param minimum_cell_number Integer; minimum cells required per group. Default 30
#' @param ... Additional arguments (unused)
#'
#' @return List with elements:
#'   \item{all_count}{Integer; number of significant DEGs}
#'   \item{UP_count}{Integer; number of upregulated genes}
#'   \item{DOWN_count}{Integer; number of downregulated genes}
#'   \item{results}{Data frame; complete DE results with gene annotations and CPMs}
#'
#' @export
DEG_FindMarkers_RNA_assay_de <- function(scRNAseq, comparison, group1, group2,
                                         cluster = 'all_clusters', path = './',
                                         FC_threshold = 0.3, p_value_threshold = 0.05,
                                         minimum_cell_number = 30, ...) {

    # Set Paths
    gene_lists_path <- here::here(path, 'gene_lists')
    dir.create(gene_lists_path, showWarnings = FALSE, recursive = TRUE)

    print(paste('Cluster', cluster))

    group1 <- stringr::fixed(group1)
    group2 <- stringr::fixed(group2)

    Idents(scRNAseq) <- comparison

    # Check there are enough cells
    cell_validation <- validate_cell_counts(scRNAseq, comparison, group1, group2, minimum_cell_number)

    print('number of cells in group 1')
    print(cell_validation$n_group1)
    print('number of cells in group 2')
    print(cell_validation$n_group2)

    if (!cell_validation$valid) {
        return(list(
            all_count = cell_validation$message,
            UP_count = cell_validation$message,
            DOWN_count = cell_validation$message,
            results = NULL
        ))
    }

    # Run FindMarkers
    results <- FindMarkers(object = scRNAseq, ident.1 = group1, ident.2 = group2,
                          assay = 'RNA', slot = 'data', test.use = 'wilcox')

    # Calculate CPMs
    scRNAseq_CPM <- scRNAseq |> AggregateExpression(group.by = c(comparison),
                                        assays = 'RNA',
                                        return.seurat = TRUE,
                                        normalization.method = 'RC',
                                        scale.factor = 1e6)

    counts_CPM <- scRNAseq_CPM |> GetAssayData(assay = 'RNA', layer = 'data') |>
                            as.data.frame() |>
                            rownames_to_column(var = 'gene') |>
                            mutate(
                                !!paste0('Avg_', group2) := !!as.name(group2),
                                !!paste0('Avg_', group1) := !!as.name(group1)
                            )

    # Add gene annotations
    results <- results |>
                    rownames_to_column('genes') |>
                    rename(
                        log2FoldChange = avg_log2FC,
                        padj = p_val_adj,
                        pvalue = p_val
                    ) |>
                    mutate(log2FoldChange = log2FoldChange * -1) |>
                    left_join(y = unique(annotations[, c('gene_name', 'description')]),
                             by = c('genes' = 'gene_name')) |>
                    left_join(y = counts_CPM, by = c('genes' = 'gene'))

    # Filter results
    results_filtered <- filter(results,
                               padj < p_value_threshold &
                               (log2FoldChange >= FC_threshold | log2FoldChange <= -1 * FC_threshold)) |>
                        arrange(padj)

    results_filtered_UP <- filter(results_filtered, log2FoldChange >= FC_threshold)
    results_filtered_DOWN <- filter(results_filtered, log2FoldChange <= -1 * FC_threshold)

    # Write results to CSV files
    write.csv(results_filtered |> arrange(padj),
             file = here::here(gene_lists_path, paste('ALL_GENES_DEG_Analysis', cluster, 'Wilcox', group2, 'vs', group1, '.csv', sep = '_')))
    write.csv(results_filtered_UP |> arrange(desc(log2FoldChange)),
             file = here::here(gene_lists_path, paste('DEG_UP_in', group2, cluster, 'Wilcox', group2, 'vs', group1, '.csv', sep = '_')))
    write.csv(results_filtered_DOWN |> arrange(log2FoldChange),
             file = here::here(gene_lists_path, paste('DEG_DOWN_in', group2, cluster, 'Wilcox', group2, 'vs', group1, '.csv', sep = '_')))

    # Return counts and results
    list(
        all_count = nrow(results_filtered),
        UP_count = nrow(results_filtered_UP),
        DOWN_count = nrow(results_filtered_DOWN),
        results = results
    )
}


#' Wilcoxon Differential Expression Analysis for SCT Assay
#'
#' Performs differential expression analysis using Seurat's FindMarkers function
#' with Wilcoxon rank-sum test on the SCTransform-normalized assay. Calculates CPMs
#' from SCT counts for visualization. Handles integrated datasets with optional
#' UMI recorrection.
#'
#' @param scRNAseq Seurat object containing SCT-normalized single-cell RNA-seq data
#' @param comparison Character; metadata column name for comparison groups
#' @param group1 Character; pattern to match for first group (control)
#' @param group2 Character; pattern to match for second group (treatment)
#' @param is_integrated_subset Logical; if TRUE, skips UMI recorrection for integrated
#'   subset objects. Default FALSE
#' @param cluster Character; cluster or cell type identifier. Default 'all_clusters'
#' @param path Character; output directory path. Default './'
#' @param FC_threshold Numeric; log2 fold change threshold for filtering. Default 0.3
#' @param p_value_threshold Numeric; adjusted p-value threshold. Default 0.05
#' @param min_fraction Numeric; minimum fraction of cells expressing a gene. Default 0.01
#' @param expression_threshold_for_gene_list Numeric; minimum average CPM for filtering.
#'   Default 20
#' @param minimum_cell_number Integer; minimum cells required per group. Default 30
#' @param ... Additional arguments (unused)
#'
#' @return List with elements:
#'   \item{all_count}{Integer; number of significant DEGs}
#'   \item{UP_count}{Integer; number of upregulated genes}
#'   \item{DOWN_count}{Integer; number of downregulated genes}
#'   \item{results}{Data frame; complete DE results with gene annotations and CPMs}
#'
#' @export
DEG_FindMarkers_SCT_assay_de <- function(scRNAseq, comparison, group1, group2, is_integrated_subset = FALSE,
                                         cluster = 'all_clusters', path = './', FC_threshold = 0.3,
                                         p_value_threshold = 0.05, min_fraction = 0.01, expression_threshold_for_gene_list = 20,
                                         minimum_cell_number = 30, ...) {

    # Set Paths
    gene_lists_path <- here::here(path, 'gene_lists')
    dir.create(gene_lists_path, showWarnings = FALSE, recursive = TRUE)

    print(paste('Cluster', cluster))

    group1 <- stringr::fixed(group1)
    group2 <- stringr::fixed(group2)

    Idents(scRNAseq) <- comparison

    print('number of cells in group 1')
    print(scRNAseq@meta.data |> filter(str_detect(!!as.name(comparison), group1)) |> nrow())
    print('number of cells in group 2')
    print(scRNAseq@meta.data |> filter(str_detect(!!as.name(comparison), group2)) |> nrow())

    # Check there are enough cells
    cell_validation <- validate_cell_counts(scRNAseq, comparison, group1, group2, minimum_cell_number)

    if (!cell_validation$valid) {
        return(list(
            all_count = cell_validation$message,
            UP_count = cell_validation$message,
            DOWN_count = cell_validation$message,
            results = NULL
        ))
    }

    # Run Wilcox test
    results <- FindMarkers(object = scRNAseq, ident.1 = group1, ident.2 = group2,
                          assay = 'SCT', slot = 'data', test.use = 'wilcox',
                          recorrect_umi = !is_integrated_subset, min.pct = min_fraction)

    # Aggregate expression and calculate CPMs
    scRNAseq_CPM <- scRNAseq |> AggregateExpression(group.by = c(comparison),
                                                     assays = 'SCT',
                                                     slot = 'counts')
    # Calculating CPMs
    counts_CPM <- scRNAseq_CPM$SCT |>
        as.data.frame() |>
        rename_with( ~ gsub("-", "_", .x, fixed = TRUE)) |>
        mutate(across(where(is.numeric), ~ .x / sum(.x) * 1e6)) |>
        rownames_to_column(var = 'gene') |>
        mutate(
            !!paste0('Avg_', group2) := !!as.name(group2),
            !!paste0('Avg_', group1) := !!as.name(group1)
        )
    counts_CPM |> head() |> print()
    # Add gene annotations
    results <- results |>
        rownames_to_column('genes') |>
        rename(
            log2FoldChange = avg_log2FC,
            padj = p_val_adj,
            pvalue = p_val
        ) |>
        mutate(log2FoldChange = log2FoldChange * -1) |>
        left_join(y = unique(annotations[, c('gene_name', 'description')]),
                 by = c('genes' = 'gene_name')) |>
        left_join(y = counts_CPM, by = c('genes' = 'gene'))

    # Filter results
    results_filtered <- filter(
        results,
        padj < p_value_threshold &
            (log2FoldChange >= FC_threshold | log2FoldChange <= -1 * FC_threshold)
    ) |>
        arrange(padj)

    results_filtered_UP <- filter(results_filtered, log2FoldChange >= FC_threshold)
    results_filtered_DOWN <- filter(results_filtered, log2FoldChange <= -1 * FC_threshold)

    # Write results to CSV files
    write.csv(results |> arrange(padj),
             file = here::here(gene_lists_path, paste('ALL_GENES_DEG_Analysis', cluster, 'Wilcox', group2, 'vs', group1, '.csv', sep = '_')))
    write.csv(results_filtered_UP |> arrange(desc(log2FoldChange)),
             file = here::here(gene_lists_path, paste('DEG_UP_in', group2, cluster, 'Wilcox', group2, 'vs', group1, '.csv', sep = '_')))
    write.csv(results_filtered_DOWN |> arrange(log2FoldChange),
             file = here::here(gene_lists_path, paste('DEG_DOWN_in', group2, cluster, 'Wilcox', group2, 'vs', group1, '.csv', sep = '_')))

    # Return counts and results
    list(
        all_count = nrow(results_filtered),
        UP_count = nrow(results_filtered_UP),
        DOWN_count = nrow(results_filtered_DOWN),
        results = results
    )
}


#' Bulk RNA-seq Differential Expression Analysis
#'
#' Performs differential expression analysis using DESeq2 on bulk RNA-seq count data.
#' Uses DESeq2's negative binomial model for statistical testing and includes PCA
#' visualization for quality control. Requires at least 2 biological replicates per group.
#'
#' @param counts_table Data frame; count matrix with genes as rows and samples as columns,
#'   must include a 'genes' column
#' @param comparison Character; column name in metadata for comparison groups. Default 'Groups'
#' @param group1 Character; pattern to match for first group (control)
#' @param group2 Character; pattern to match for second group (treatment)
#' @param cluster Character; identifier for the dataset. Default ''
#' @param path Character; output directory path. Default './'
#' @param FC_threshold Numeric; log2 fold change threshold for filtering. Default 0.3
#' @param p_value_threshold Numeric; adjusted p-value threshold. Default 0.05
#' @param expression_threshold_for_gene_list Numeric; minimum average normalized counts
#'   for a gene to be included in filtered gene lists. Default 20
#' @param ... Additional arguments (unused)
#'
#' @return List with elements:
#'   \item{all_count}{Integer; number of significant DEGs}
#'   \item{UP_count}{Integer; number of upregulated genes}
#'   \item{DOWN_count}{Integer; number of downregulated genes}
#'   \item{results}{Data frame; complete DE results with gene annotations and counts}
#'
#' @export
bulk_analysis_de <- function(counts_table, comparison = 'Groups', group1, group2,
                             cluster = '', path = './', FC_threshold = 0.3,
                             p_value_threshold = 0.05,
                             expression_threshold_for_gene_list = 20, ...) {

    # Set Paths
    gene_lists_path <- here::here(path, 'gene_lists')
    dir.create(gene_lists_path, showWarnings = FALSE, recursive = TRUE)

    print(paste('Cluster', cluster))

    group1 <- stringr::fixed(group1)
    group2 <- stringr::fixed(group2)

    counts <- tibble(counts_table) |> column_to_rownames('genes')

    # Generate sample level metadata
    colData <- data.frame(samples = colnames(counts)) |>
                mutate(condition = ifelse(grepl(group1, samples), group1, group2))

    # Filter
    counts <- counts |> mutate(row_sums = rowSums(counts)) |>
              filter(row_sums >= 10) |> dplyr::select(-row_sums)

    print('Group 1 Length')
    print(nrow(colData |> filter(condition == group1)))
    print('Group 2 Length')
    print(nrow(colData |> filter(condition == group2)))

    # Check for sufficient replicates
    if ((length(unique(colData$condition)) != 2) |
        (nrow(colData |> filter(condition == group1)) < 2) |
        (nrow(colData |> filter(condition == group2)) < 2)) {
        return(list(
            all_count = 'Not enough biological replicates per group',
            UP_count = 'Not enough biological replicates per group',
            DOWN_count = 'Not enough biological replicates per group',
            results = NULL
        ))
    }

    # Create DESeq2 object
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
        colData = colData,
        design = ~condition)
    dds$condition <- factor(dds$condition, levels = c(group1, group2))

    # DESeq2 QC
    rld <- DESeq2::rlog(dds, blind = TRUE)
    local_figures_path <- here::here(path, 'figures')
    dir.create(local_figures_path, showWarnings = FALSE, recursive = TRUE)

    DESeq2::plotPCA(rld, ntop = 500, intgroup = 'condition')
    ggsave(filename = paste0('Bulk_PCA_', cluster, '.pdf'), path = local_figures_path)

    PCA_table <- DESeq2::plotPCA(rld, ntop = 500, intgroup = 'condition', returnData = TRUE)
    write.csv(PCA_table, file = here::here(path, paste('PCA_bulk', cluster, group2, 'vs', group1, '.csv', sep = '_')))

    # Run DESeq2
    dds <- DESeq2::DESeq(dds)
    DESeq2::resultsNames(dds)

    # Generate results object
    results <- DESeq2::results(dds) |> as.data.frame()

    # Get Normalized Counts
    normalized_counts <- DESeq2::counts(dds, normalized = TRUE)
    normalized_counts <- normalized_counts |>
        as.data.frame() |>
        rownames_to_column('genes') |>
        as_tibble() |>
        rowwise() |>
        mutate(
            !!paste0('Avg_', group2) := mean(c_across(contains(group2))),
            !!paste0('Avg_', group1) := mean(c_across(contains(group1)))
        ) |>
        ungroup()

    # Add gene annotations
    results <- results |>
                    rownames_to_column('genes') |>
                    left_join(y = unique(annotations[, c('gene_name', 'description')]),
                             by = c('genes' = 'gene_name')) |>
                    left_join(y = normalized_counts, by = c('genes' = 'genes'))

    # Filter results
    results_filtered <- filter(
        results,
        padj < p_value_threshold &
            ((!!sym(paste0('Avg_', group2)) > expression_threshold_for_gene_list) |
             (!!sym(paste0('Avg_', group1)) > expression_threshold_for_gene_list)) &
            (log2FoldChange >= FC_threshold | log2FoldChange <= -1 * FC_threshold)
    ) |>
        arrange(padj)

    results_filtered_UP <- filter(results_filtered, log2FoldChange >= FC_threshold)
    results_filtered_DOWN <- filter(results_filtered, log2FoldChange <= -1 * FC_threshold)

    # Write results to CSV files
    write.csv(results |> arrange(padj),
             file = here::here(gene_lists_path, paste('ALL_GENES_DEG_Analysis', cluster, 'bulk', group2, 'vs', group1, '.csv', sep = '_')))
    write.csv(results_filtered_UP |> arrange(desc(log2FoldChange)),
             file = here::here(gene_lists_path, paste('DEG_UP_in', group2, cluster, 'bulk', group2, 'vs', group1, '.csv', sep = '_')))
    write.csv(results_filtered_DOWN |> arrange(log2FoldChange),
             file = here::here(gene_lists_path, paste('DEG_DOWN_in', group2, cluster, 'bulk', group2, 'vs', group1, '.csv', sep = '_')))

    # Return counts and results
    list(
        all_count = nrow(results_filtered),
        UP_count = nrow(results_filtered_UP),
        DOWN_count = nrow(results_filtered_DOWN),
        results = results
    )
}
