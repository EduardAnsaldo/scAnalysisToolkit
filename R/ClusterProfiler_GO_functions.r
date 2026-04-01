## ClusterProfiler Gene Ontology (GO) Enrichment Functions
## Functions for GO enrichment analysis using clusterProfiler

#' Gene Ontology Overrepresentation Analysis
#'
#' Performs Gene Ontology (GO) overrepresentation analysis using clusterProfiler's
#' enrichGO function. Tests which GO terms are statistically enriched in a gene list
#' compared to a background. Optionally simplifies redundant terms and creates network
#' visualizations.
#'
#' @param significant_genes Character vector; gene symbols to test for enrichment
#' @param all_genes Character vector; background/universe gene symbols
#' @param local_path Character; path to save results and plots
#' @param ontology Character; GO ontology ('BP', 'MF', 'CC', or 'ALL'). Default 'ALL'
#' @param minGSSize Integer; minimum gene set size. Default 5
#' @param maxGSSize Integer; maximum gene set size. Default 400
#' @param filename Character; prefix for output files. Default ''
#' @param group Character; group identifier for plot titles. Default ''
#' @param drop_levels Logical; whether to drop specific GO levels. Default FALSE
#' @param levels_to_drop Integer vector; GO levels to exclude if drop_levels=TRUE. Default c()
#' @param simplify_function Function; function for simplifying redundant terms (min, max, mean). Default min
#' @param simplify_by Character; metric for simplification ('p.adjust', 'pvalue'). Default 'p.adjust'
#' @param simplify_terms Logical; whether to remove redundant GO terms. Default TRUE
#' @param run_network Logical; whether to create enrichment network plots. Default FALSE
#' @param network_n_terms Integer; number of terms to include in network. Default 100
#' @param nterms_to_plot Integer; number of terms to show in dot plot. Default 50
#' @param font_size Numeric; font size for plots. Default 8
#' @param ... Additional arguments (unused)
#'
#' @return Invisible NULL; creates CSV files and PDF plots
#'
#' @export
GO_overrepresentation_analysis <- function (significant_genes, all_genes, local_path, ontology = 'ALL', minGSSize = 5, maxGSSize = 400, filename = '', group = '', drop_levels = F, levels_to_drop = c(), simplify_function = min, simplify_by = 'p.adjust', simplify_terms = T, run_network = F, network_n_terms = 100, nterms_to_plot = 50, font_size = 8,  ...)  {

     set_enrichment_color_scale()

     if (length(significant_genes) <= 4) {
          p1 <- create_no_data_plot()
          print(p1)
          return(invisible(NULL))
     }

     enrichment_results <- clusterProfiler::enrichGO(gene = significant_genes,
                    universe = all_genes,
                    keyType = "SYMBOL",
                    OrgDb = org.Mm.eg.db,
                    ont = ontology,
                    pAdjustMethod = "BH",
                    minGSSize    = minGSSize,
                    maxGSSize    = maxGSSize,
                    qvalueCutoff = 0.25)

     if (drop_levels  == T ) {
            enrichment_results <- clusterProfiler::dropGO(enrichment_results, level = levels_to_drop)

     }

     enrichment_results_table <- as_tibble(enrichment_results)
     write.csv(enrichment_results_table, here::here(local_path, paste0(filename,'GO_OverRepresentation_analysis_results_', ontology, '.csv')))

     if (nrow(enrichment_results_table) > 1) {
       if (simplify_terms == T) {
            ## Add similarity matrix to the termsim slot of enrichment result
          enrichment_results <- enrichplot::pairwise_termsim(enrichment_results, showCategory = dim(enrichment_results)[1])
          enrichment_results_unfiltered <-  enrichment_results
          enrichment_results <- clusterProfiler::simplify(enrichment_results, cutoff=0.7, by=simplify_by, select_fun=simplify_function)
       }

          write.csv(as_tibble(enrichment_results), here::here(local_path, paste0(filename,'GO_OverRepresentation_analysis_results_filtered_', ontology, '.csv')))
          p1 <- enrichplot::dotplot(enrichment_results,
               showCategory=nterms_to_plot,
               title = paste0(filename,'GO ORA UP in ', group, ' - ', ontology),
            #    x = 'p.adjust',
               label_format = 60,
               font.size = font_size)
          print(p1)
          ggsave(plot = p1, filename = paste0(filename, 'GO overrepresentation_analysis_dotplot_', ontology,'.pdf'), width = 10, height = 18, path = local_path)

          ## =  Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
          if (run_network) {
               p2 <- try(enrichplot::emapplot(enrichment_results_unfiltered, showCategory = network_n_terms) + ggtitle(paste0(filename, 'GO Overrepresentation analysis ', ontology)))
               print(p2)
               p2 <- try(enrichplot::emapplot(enrichment_results, showCategory = network_n_terms) + ggtitle(paste0(filename, 'GO Overrepresentation analysis filtered terms', ontology)))
               print(p2)
            #    ggsave(plot = p2, filename = paste0(filename, 'GO_overrepresentation_analysis_network_', ontology,'.pdf'), width = 14, height = 18, path = local_path)

          }

          # ggsave(plot = p2, filename = paste0(filename, 'GO_overrepresentation_analysis_network_', ontology,'.pdf'), width = 14, height = 18, path = local_path)
     } else {
          p1 <- create_no_data_plot()
          print(p1)
     }
}

#' Gene Ontology Gene Set Enrichment Analysis
#'
#' Performs Gene Ontology (GO) Gene Set Enrichment Analysis (GSEA) using clusterProfiler's
#' gseGO function. Tests whether genes ranked by fold change show enrichment for specific
#' GO terms. Creates dot plots and saves results to CSV files.
#'
#' @param results Data frame; differential expression results with columns: genes,
#'   log2FoldChange, padj
#' @param local_path Character; path to save results and plots
#' @param ontology Character; GO ontology ('BP', 'MF', 'CC', or 'ALL'). Default 'ALL'
#' @param group Character; group identifier for plot titles
#' @param minGSSize Integer; minimum gene set size. Default 5
#' @param maxGSSize Integer; maximum gene set size. Default 500
#' @param pvalueCutoff Numeric; p-value cutoff for significance. Default 0.05
#' @param ... Additional arguments (unused)
#'
#' @return Invisible NULL; creates CSV files and PDF plots
#'
#' @details
#' Genes are ranked by log2 fold change. The analysis uses gseGO with mouse
#' annotation database (org.Mm.eg.db) and gene symbols as identifiers.
#'
#' @export
GO_GSEA_analysis <- function (results, local_path, ontology = 'ALL', group, minGSSize = 5, maxGSSize = 500, pvalueCutoff = 0.05, ...) {

     set_enrichment_color_scale()

          #### GSEA ####

     fold_changes <- results |> arrange(padj, desc(log2FoldChange)) |> pull(log2FoldChange)
     names(fold_changes) <- results |> arrange(padj, desc(log2FoldChange)) |> pull(genes)

     gsea_results <- clusterProfiler::gseGO(geneList     = fold_changes,
               OrgDb        = org.Mm.eg.db,
               ont          = ontology,
               keyType = "SYMBOL",
               minGSSize    = minGSSize,
               maxGSSize    = maxGSSize,
               pvalueCutoff = pvalueCutoff,
               verbose      = FALSE,
            eps = 0)

     gsea_results_table <- as_tibble(gsea_results)
     write.csv(gsea_results_table, here::here(local_path, paste0('Gene_Set_Enrichment_Analysis_results_', ontology,'.csv')))

     if (nrow(gsea_results_table) > 1) {
          p3 <- enrichplot::dotplot(gsea_results,
               showCategory=50,
               title = paste0('GSEA analysis ', group, ' - ', ontology),
               label_format = 60)
          print(p3)
          ggsave(plot = p3, filename = paste0('GSEA_dotplot_', ontology,'.pdf'), width = 10, height = 18, path = local_path)

          ## Add similarity matrix to thenes,  termsim slot of enrichment result
          gsea_results <- enrichplot::pairwise_termsim(gsea_results, showCategory = dim(gsea_results)[1])

        #   ## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
        # try(enrichplot::emapplot(gsea_results) + ggtitle(paste0('GSEA analysis ', ontology)))
        #   ggsave(paste0('GSEA_network_', ontology,'.pdf'), width = 14, height = 18, path = local_path)
           }
}

#' GO Functional Analysis on Differential Expression Results
#'
#' Wrapper function that performs comprehensive GO enrichment analysis on differential
#' expression results. Runs overrepresentation analysis (ORA) separately for upregulated
#' and downregulated genes, with optional Gene Set Enrichment Analysis (GSEA). Creates
#' organized output directories with plots and tables.
#'
#' @param results Data frame; differential expression results with columns: genes,
#'   log2FoldChange, padj
#' @param grouping_var Character; identifier for cluster or group being analyzed
#' @param path Character; base output directory path. Default './'
#' @param FC_threshold Numeric; log2 fold change threshold for filtering genes. Default 0.3
#' @param p_value_threshold Numeric; adjusted p-value threshold. Default 0.05
#' @param group1 Character; name of control group. Default ''
#' @param group2 Character; name of treatment group. Default ''
#' @param run_GSEA Logical; whether to run GSEA in addition to ORA. Default FALSE
#' @param top_n_genes Integer or NULL; number of top genes to use for ORA. If NULL, uses all significant genes. Default NULL
#' @param ... Additional arguments passed to GO_overrepresentation_analysis
#'
#' @return Invisible NULL; creates output directories with enrichment results
#'
#' @export
GO_functional_analysis <- function (results,  grouping_var, path='./', FC_threshold = 0.3, p_value_threshold = 0.05, group1 = '', group2 = '', run_GSEA = FALSE, top_n_genes = NULL, ...) {

cluster <- grouping_var

####################################### UP ########################################

      local_path <- create_analysis_directory(here::here(path, paste0('GO_functional_analysis_UP_in_', group2, '_', cluster)))

      significant_genes <- results |>
            filter((padj < p_value_threshold) & (log2FoldChange > FC_threshold)) |>
            arrange(padj)

      if (!is.null(top_n_genes)) {
            significant_genes <- significant_genes |> slice_head(n = top_n_genes)
      }

      significant_genes <- significant_genes |> pull(genes)
      all_genes <- results |> arrange(padj) |> pull(genes)

      ######################## ORA ########################

      GO_overrepresentation_analysis(significant_genes, all_genes, local_path, group = group2, ...)

      #################### GSEA ####################
     if (run_GSEA) {
      GO_GSEA_analysis(results, local_path, group = group2, ...)
     }


######################################## DOWN ########################################

      local_path <- create_analysis_directory(here::here(path, paste0('GO_functional_analysis_UP_in_', group1, '_', cluster)))

      significant_genes <- results |>
            filter((padj < p_value_threshold) & (log2FoldChange < -1*FC_threshold)) |>
            arrange(padj)

      if (!is.null(top_n_genes)) {
            significant_genes <- significant_genes |> slice_head(n = top_n_genes)
      }

      significant_genes <- significant_genes |> pull(genes)
      all_genes <- results |> arrange(padj) |> pull(genes)

      ######################## ORA ########################

      GO_overrepresentation_analysis(significant_genes, all_genes, local_path, group = group1, ...)


      #################### GSEA ####################
      if (run_GSEA) {
      GO_GSEA_analysis(results, local_path, group1, ...)
      }

      return()
}

#' GO Functional Analysis for Cluster Identification
#'
#' Performs Gene Ontology (GO) enrichment analysis on marker genes from multiple
#' clusters to characterize cluster identities. Extracts top marker genes per cluster,
#' runs GO overrepresentation analysis using clusterProfiler's compareCluster, and
#' generates a dot plot visualization comparing enriched terms across clusters.
#'
#' @param scRNAseq Seurat object with RNA assay containing gene expression data
#' @param results Data frame; marker gene results from FindAllMarkers with columns:
#'   cluster, gene, avg_log2FC, p_val_adj
#' @param path Character; base output directory path. Default './'
#' @param object_annotations Character; string to append to output directory name.
#'   Default ''
#' @param top_gene_number Integer; number of top marker genes per cluster to use
#'   for enrichment analysis. Default 50
#' @param ... Additional arguments passed to GO_overrepresentation_analysis_multiple_lists
#'   (e.g., ontology, minGSSize, maxGSSize, simplify_terms)
#'
#' @return ggplot object; dot plot showing enriched GO terms across clusters
#'
#' @details
#' This function:
#' 1. Extracts the top N marker genes per cluster ranked by log2 fold change
#' 2. Creates gene lists for each cluster
#' 3. Runs GO overrepresentation analysis comparing all clusters simultaneously
#' 4. Creates an output directory: 'Cluster_identification_functional_analysis_GO_annotations'
#' 5. Saves enrichment results to CSV files
#' 6. Returns a dot plot for visualization
#'
#' The analysis uses all genes in the RNA assay as the background/universe for
#' statistical testing.
#'
#' @export
GO_functional_analysis_cluster_identification <- function (scRNAseq, results, path='./', object_annotations = '', top_gene_number = 50, ...) {

    set_enrichment_color_scale()
     local_path <- create_analysis_directory(here::here(path, paste0('Cluster_identification_functional_analysis_GO_', object_annotations)))

     all_genes <- Features(scRNAseq[['RNA']]) |>unique()

    gene_lists <- results |>
        group_by(cluster) |>
        arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) |>
        slice_head(n = top_gene_number)  |>
        ungroup() |>
        select(cluster, gene)  |>
        group_by(cluster) |>
        mutate(id = row_number()) |>
        ungroup() |>
        pivot_wider(names_from = cluster, values_from = gene)  |>
        select(-id) |>
        as.list() |>
        map(~ .x[!is.na(.x)])

    overrepresentation_results <- GO_overrepresentation_analysis_multiple_lists(gene_lists, all_genes = all_genes, local_path = local_path, grouping_var = 'Clusters', ...)

    return(overrepresentation_results$plot)
}

#' GO Overrepresentation Analysis for Multiple Gene Lists
#'
#' Performs Gene Ontology (GO) overrepresentation analysis on multiple gene lists using
#' clusterProfiler's compareCluster with enrichGO. Compares enriched GO terms across
#' multiple groups/clusters simultaneously and generates a comparative dot plot visualization.
#'
#' @param gene_list Named list; each element contains gene symbols for a group/cluster
#' @param all_genes Character vector; background/universe gene symbols
#' @param local_path Character; path to save results and plots
#' @param ontology Character; GO ontology ('BP', 'MF', 'CC', or 'ALL'). Default 'ALL'
#' @param minGSSize Integer; minimum gene set size. Default 5
#' @param maxGSSize Integer; maximum gene set size. Default 400
#' @param filename Character; prefix for output files. Default ''
#' @param drop_levels Logical; whether to drop specific GO levels. Default FALSE
#' @param levels_to_drop Integer vector; GO levels to exclude if drop_levels=TRUE. Default c()
#' @param simplify_function Function; function for simplifying redundant terms (min, max, mean). Default min
#' @param simplify_by Character; metric for simplification ('p.adjust', 'pvalue'). Default 'p.adjust'
#' @param simplify_terms Logical; whether to remove redundant GO terms. Default FALSE
#' @param run_network Logical; whether to create enrichment network plots. Default FALSE
#' @param network_n_terms Integer; number of terms to include in network. Default 100
#' @param nterms_to_plot Integer; number of terms to show per cluster in dot plot. Default 5
#' @param font_size Numeric; font size for plots. Default 8
#' @param grouping_var Character; label for x-axis grouping variable. Default ''
#' @param ... Additional arguments (unused)
#'
#' @return List with two elements: plot (ggplot object) and results (enrichment results)
#'
#' @details
#' Uses mouse annotation database (org.Mm.eg.db) with gene symbols. Creates a dot plot
#' comparing enriched GO terms across all provided gene lists. Automatically orders
#' clusters numerically if possible.
#'
#' @export
GO_overrepresentation_analysis_multiple_lists <- function (gene_list, all_genes, local_path, ontology = 'ALL', minGSSize = 5, maxGSSize = 400, filename = '',  drop_levels = F, levels_to_drop = c(), simplify_function = min, simplify_by = 'p.adjust', simplify_terms = F, run_network = F, network_n_terms = 100, nterms_to_plot = 5, font_size = 8, grouping_var = '', ...)  {

     set_enrichment_color_scale()

     enrichment_results <- clusterProfiler::compareCluster(
          geneCluster = gene_list,
          fun = "enrichGO",
          universe = all_genes,
          keyType = "SYMBOL",
          OrgDb = org.Mm.eg.db,
          ont = ontology,
          pAdjustMethod = "BH",
          minGSSize    = minGSSize,
          maxGSSize    = maxGSSize,
          qvalueCutoff = 0.2,
          pvalueCutoff = 0.05
     )

     if (drop_levels  == T ) {
            enrichment_results <- clusterProfiler::dropGO(enrichment_results, level = levels_to_drop)

     }

     enrichment_results <- enrichment_results |>
          mutate(Cluster = if (all(!is.na(suppressWarnings(as.numeric(levels(factor(Cluster))))))) {
               fct_inseq(Cluster)
          } else {
               factor(Cluster)
          })

     enrichment_results_table <- as_tibble(enrichment_results)
     write.csv(enrichment_results_table, here::here(local_path, paste0(filename,'GO_OverRepresentation_analysis_gene_lists_results_', ontology, '.csv')))

     if (nrow(enrichment_results_table) > 1) {
       if (simplify_terms == T) {
            ## Add similarity matrix to the termsim slot of enrichment result
          enrichment_results <- enrichplot::pairwise_termsim(enrichment_results, showCategory = dim(enrichment_results)[1])
          enrichment_results_unfiltered <-  enrichment_results
          enrichment_results <- clusterProfiler::simplify(enrichment_results, cutoff=0.7, by=simplify_by, select_fun=simplify_function)
       }

          write.csv(as_tibble(enrichment_results), here::here(local_path, paste0(filename,'GO_OverRepresentation_analysis_gene_lists_results__filtered_', ontology, '.csv')))
          p1 <- enrichplot::dotplot(enrichment_results,
               showCategory=nterms_to_plot,
               title = paste0(filename,'GO Overrepresentation analysis ', ' - ', ontology),
            #    x = 'p.adjust',
               label_format = 60,
               font.size = font_size) +
               labs(x = grouping_var) +
               theme(axis.text.x = element_text(angle = 45, hjust = 1))
#           print(p1)
          ggsave(plot = p1, filename = paste0(filename, 'GO overrepresentation_analysis_dotplot_', ontology,'.pdf'), width = 10, height = 18, path = local_path)
     }
     return(list(plot = p1, results = enrichment_results))
}
