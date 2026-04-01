## ClusterProfiler MSigDB Pathway Enrichment Functions
## Functions for pathway enrichment analysis using MSigDB collections with clusterProfiler

#' Pathway Overrepresentation Analysis using MSigDB
#'
#' Performs pathway overrepresentation analysis using clusterProfiler's enricher
#' function with MSigDB pathway collections. Tests which pathways are statistically
#' enriched in a gene list compared to a background.
#'
#' @param significant_genes Character vector; gene symbols to test for enrichment
#' @param all_genes Character vector; background/universe gene symbols
#' @param local_path Character; path to save results and plots
#' @param minGSSize Integer; minimum gene set size. Default 5
#' @param maxGSSize Integer; maximum gene set size. Default 400
#' @param filename Character; prefix for output files. Default ''
#' @param group Character; group identifier for plot titles. Default ''
#' @param nterms_to_plot Integer; number of terms to show in dot plot. Default 50
#' @param font_size Numeric; font size for plots. Default 8
#' @param ... Additional arguments (unused)
#'
#' @return Invisible NULL; creates CSV files and PDF plots
#'
#' @export
pathway_overrepresentation_analysis <- function (significant_genes, all_genes, local_path, minGSSize = 5, maxGSSize = 400, filename = '', group = '', nterms_to_plot = 30, font_size = 10,  ...)  {

     set_enrichment_color_scale()

     if (length(significant_genes) <= 4) {
          p1 <- create_no_data_plot()
          print(p1)
          return(invisible(NULL))
     }

     # Convert gene symbols to Entrez IDs
     significant_genes_ids <- clusterProfiler::bitr(significant_genes, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
     all_genes_ids <- clusterProfiler::bitr(all_genes, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")

     significant_genes_entrez <- significant_genes_ids$ENTREZID
     all_genes_entrez <- all_genes_ids$ENTREZID

     if (length(significant_genes_entrez) <= 4) {
          p1 <- create_no_data_plot()
          print(p1)
          return(invisible(NULL))
     }

     # Create TERM2GENE database from MSigDB pathway collections
     pathway_collections <- list(c('MH'), c('M2', 'CP:BIOCARTA'), c('M2', 'CP:REACTOME'), c('M2', 'CP:WIKIPATHWAYS'), c('M5', 'GO:BP'))
     msigdbr_pathways <- purrr::map_dfr(
          pathway_collections,
          function(pc) {
               if (length(pc) >= 2 && !is.na(pc[2]) && nzchar(pc[2])) {
                    msigdbr::msigdbr(db_species = "MM", species = "mouse", collection = pc[1], subcategory = pc[2])
               } else {
                    msigdbr::msigdbr(db_species = "MM", species = "mouse", collection = pc[1])
               }
          }
     ) |> select(gs_name, ncbi_gene)

     enrichment_results <- clusterProfiler::enricher(gene = significant_genes_entrez,
                    universe = all_genes_entrez,
                    pAdjustMethod = "BH",
                    minGSSize    = minGSSize,
                    maxGSSize    = maxGSSize,
                    qvalueCutoff = 0.25,
                    TERM2GENE = msigdbr_pathways)

     enrichment_results_table <- as_tibble(enrichment_results)
     write.csv(enrichment_results_table, here::here(local_path, paste0(filename,'Pathway_OverRepresentation_analysis_results.csv')))

     colnames(enrichment_results_table) |> print()

     # if (nrow(enrichment_results_table) > 1) {
     #      p1 <- enrichplot::dotplot(enrichment_results,
     #           showCategory=nterms_to_plot,
     #           title = paste0(filename,'Pathways UP in ', group),
     #           label_format = 60,
     #           font.size = font_size)
     #      print(p1)
     #      ggsave(plot = p1, filename = paste0(filename, 'Pathway_overrepresentation_analysis_dotplot.pdf'), width = 10, height = 18, path = local_path)
     # } else {
     #      p1 <- create_no_data_plot()
     #      print(p1)
     # }

     if (nrow(enrichment_results_table) > 1) {
          p1 <- enrichment_results_table |>       
               mutate(log_q_value = -log10(qvalue)) |>
               slice_max(order_by = log_q_value, n = nterms_to_plot) |>
          ggplot(aes(x = log_q_value, y = fct_reorder(Description, log_q_value , .desc = F), fill = log_q_value)) +
               geom_col(width = 0.7) +
               scale_fill_binned(palette = grDevices::hcl.colors(n = 6, 'YlOrRd', rev = T), breaks = c(2, 4, 6, 10, 20)) +
               labs(x = '-log10(p-value)', y = '', title = paste0(filename,'Pathways UP in ', group)) +
               theme_minimal() +
               scale_y_discrete(labels = function(x) stringr::str_trunc(x, 60)) + 
               theme(axis.text.y = element_text(size = 9), title = element_text(size = 16), plot.title.position = 'plot', legend.position = 'none', axis.text.x = element_text(size = 10))
     #    print(plot2)
     ggsave(plot = p1, filename = paste0(filename, '_Pathway_enrichment_analysis_metascape', '.pdf'), width = 6, height = 6, path = local_path)
     print(p1)
     } else {
          p1 <- create_no_data_plot()
          print(p1)
     }
}

#' Pathway Gene Set Enrichment Analysis using MSigDB
#'
#' Performs Gene Set Enrichment Analysis (GSEA) using clusterProfiler's GSEA function
#' with MSigDB pathway collections. Tests whether genes ranked by fold change show
#' enrichment for specific pathways. Creates dot plots and saves results to CSV files.
#'
#' @param results Data frame; differential expression results with columns: genes,
#'   log2FoldChange, padj
#' @param local_path Character; path to save results and plots
#' @param group Character; group identifier for plot titles
#' @param minGSSize Integer; minimum gene set size. Default 5
#' @param maxGSSize Integer; maximum gene set size. Default 500
#' @param pvalueCutoff Numeric; p-value cutoff for significance. Default 0.05
#' @param ... Additional arguments (unused)
#'
#' @return Invisible NULL; creates CSV files and PDF plots
#'
#' @details
#' Uses MSigDB pathway collections: Hallmark (MH), BIOCARTA, REACTOME, WIKIPATHWAYS,
#' GO:BP, and Immunologic signatures (M7). Genes are ranked by log2 fold change.
#'
#' @export
pathway_GSEA_analysis <- function (results, local_path, group, minGSSize = 5, maxGSSize = 500, pvalueCutoff = 0.05, ...) {

     set_enrichment_color_scale()

     #### GSEA ####

     # Convert gene symbols to Entrez IDs
     gene_ids <- clusterProfiler::bitr(results$genes, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")

     # Merge with results to keep fold changes
     results_with_ids <- results |>
          inner_join(gene_ids, by = c("genes" = "SYMBOL")) |>
          arrange(padj, desc(log2FoldChange))

     # Create TERM2GENE database from MSigDB pathway collections
     pathway_collections <- list(c('MH'), c('M2', 'CP:BIOCARTA'), c('M2', 'CP:REACTOME'), c('M2', 'CP:WIKIPATHWAYS'), c('M5', 'GO:BP'))
     msigdbr_pathways <- purrr::map_dfr(
          pathway_collections,
          function(pc) {
               if (length(pc) >= 2 && !is.na(pc[2]) && nzchar(pc[2])) {
                    msigdbr::msigdbr(db_species = "MM", species = "mouse", collection = pc[1], subcategory = pc[2])
               } else {
                    msigdbr::msigdbr(db_species = "MM", species = "mouse", collection = pc[1])
               }
          }
     ) |> select(gs_name, ncbi_gene)

     fold_changes <- results_with_ids |> pull(log2FoldChange)
     names(fold_changes) <- results_with_ids |> pull(ENTREZID)

     gsea_results <- clusterProfiler::GSEA(geneList     = fold_changes,
               minGSSize    = minGSSize,
               maxGSSize    = maxGSSize,
               pvalueCutoff = pvalueCutoff,
               verbose      = FALSE,
               eps = 0,
               TERM2GENE = msigdbr_pathways)

     gsea_results_table <- as_tibble(gsea_results)
     write.csv(gsea_results_table, here::here(local_path, paste0('Pathway_GSEA_analysis_results.csv')))

     if (nrow(gsea_results_table) > 1) {
          p3 <- enrichplot::dotplot(gsea_results,
               showCategory=50,
               title = paste0('Pathway GSEA analysis ', group),
               label_format = 60)
          print(p3)
          ggsave(plot = p3, filename = paste0('Pathway_GSEA_dotplot.pdf'), width = 10, height = 18, path = local_path)
           }
}

#' Pathway Functional Analysis on Differential Expression Results using MSigDB
#'
#' Wrapper function that performs comprehensive pathway enrichment analysis on differential
#' expression results using MSigDB pathway collections. Runs overrepresentation analysis (ORA)
#' separately for upregulated and downregulated genes, with optional Gene Set Enrichment
#' Analysis (GSEA). Creates organized output directories with plots and tables.
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
#' @param ... Additional arguments passed to pathway_overrepresentation_analysis
#'
#' @return Invisible NULL; creates output directories with enrichment results
#'
#' @details
#' Uses MSigDB pathway collections: Hallmark (MH), BIOCARTA, REACTOME, WIKIPATHWAYS,
#' GO:BP, and Immunologic signatures (M7). Creates separate analyses for genes
#' upregulated in group2 vs group1 and downregulated genes (upregulated in group1).
#'
#' @export
pathway_functional_analysis <- function (results,  grouping_var, path='./', FC_threshold = 0.3, p_value_threshold = 0.05, group1 = '', group2 = '', run_GSEA = FALSE, top_n_genes = NULL, ...) {

cluster <- grouping_var

####################################### UP ########################################

      local_path <- create_analysis_directory(here::here(path, paste0('Pathway_functional_analysis_UP_in_', group2, '_', cluster)))

      significant_genes <- results |>
            filter((padj < p_value_threshold) & (log2FoldChange > FC_threshold)) |>
            arrange(padj)

      if (!is.null(top_n_genes)) {
            significant_genes <- significant_genes |> slice_head(n = top_n_genes)
      }

      significant_genes <- significant_genes |> pull(genes)
      all_genes <- results |> arrange(padj) |> pull(genes)

      ######################## ORA ########################

      pathway_overrepresentation_analysis(significant_genes, all_genes, local_path, group = group2, ...)

      #################### GSEA ####################
     if (run_GSEA) {
      pathway_GSEA_analysis(results, local_path, group = group2, ...)
     }


######################################## DOWN ########################################

      local_path <- create_analysis_directory(here::here(path, paste0('Pathway_functional_analysis_UP_in_', group1, '_', cluster)))

      significant_genes <- results |>
            filter((padj < p_value_threshold) & (log2FoldChange < -1*FC_threshold)) |>
            arrange(padj)

      if (!is.null(top_n_genes)) {
            significant_genes <- significant_genes |> slice_head(n = top_n_genes)
      }

      significant_genes <- significant_genes |> pull(genes)
      all_genes <- results |> arrange(padj) |> pull(genes)

      ######################## ORA ########################

      pathway_overrepresentation_analysis(significant_genes, all_genes, local_path, group = group1, ...)


      #################### GSEA ####################
      if (run_GSEA) {
      pathway_GSEA_analysis(results, local_path, group1, ...)
      }

      return()
}

#' Pathway Functional Analysis for Cluster Identification using MSigDB
#'
#' Performs pathway enrichment analysis using MSigDB collections on marker genes from
#' multiple clusters to characterize cluster identities. Extracts top marker genes per
#' cluster, runs pathway overrepresentation analysis using clusterProfiler's compareCluster
#' with enricher, and generates a dot plot visualization comparing enriched pathways
#' across clusters.
#'
#' @param scRNAseq Seurat object with RNA assay containing gene expression data
#' @param results Data frame; marker gene results from FindAllMarkers with columns:
#'   cluster, gene, avg_log2FC, p_val_adj
#' @param path Character; base output directory path. Default './'
#' @param object_annotations Character; string to append to output directory name.
#'   Default ''
#' @param top_gene_number Integer; number of top marker genes per cluster to use
#'   for enrichment analysis. Default 50
#' @param ... Additional arguments passed to pathway_overrepresentation_analysis_multiple_lists
#'   (e.g., minGSSize, maxGSSize)
#'
#' @return ggplot object; dot plot showing enriched pathways across clusters
#'
#' @details
#' This function:
#' 1. Extracts the top N marker genes per cluster ranked by log2 fold change
#' 2. Creates gene lists for each cluster
#' 3. Runs pathway overrepresentation analysis comparing all clusters simultaneously
#' 4. Creates an output directory: 'Cluster_identification_functional_analysis_Pathway_annotations'
#' 5. Saves enrichment results to CSV files
#' 6. Returns a dot plot for visualization
#'
#' Uses MSigDB pathway collections: Hallmark (MH), BIOCARTA, REACTOME, WIKIPATHWAYS,
#' GO:BP, and Immunologic signatures (M7). The analysis uses all genes in the RNA
#' assay as the background/universe for statistical testing.
#'
#' @export
pathway_functional_analysis_cluster_identification <- function (scRNAseq, results, path='./', object_annotations = '', top_gene_number = 50, ...) {

    set_enrichment_color_scale()
     local_path <- create_analysis_directory(here::here(path, paste0('Cluster_identification_functional_analysis_Pathway_', object_annotations)))

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

    overrepresentation_results <- pathway_overrepresentation_analysis_multiple_lists(gene_lists, all_genes = all_genes, local_path = local_path, grouping_var = 'Clusters', ...)

    return(overrepresentation_results$plot)
}

#' Pathway Overrepresentation Analysis for Multiple Gene Lists using MSigDB
#'
#' Performs pathway overrepresentation analysis on multiple gene lists using
#' clusterProfiler's compareCluster with enricher and MSigDB pathway collections.
#' Compares enriched pathways across multiple groups/clusters simultaneously and
#' generates a comparative dot plot visualization.
#'
#' @param gene_list Named list; each element contains gene symbols for a group/cluster
#' @param all_genes Character vector; background/universe gene symbols
#' @param local_path Character; path to save results and plots
#' @param minGSSize Integer; minimum gene set size. Default 5
#' @param maxGSSize Integer; maximum gene set size. Default 400
#' @param filename Character; prefix for output files. Default ''
#' @param nterms_to_plot Integer; number of terms to show per cluster in dot plot. Default 5
#' @param font_size Numeric; font size for plots. Default 8
#' @param grouping_var Character; label for x-axis grouping variable. Default ''
#' @param ... Additional arguments (unused)
#'
#' @return List with two elements: plot (ggplot object) and results (enrichment results)
#'
#' @details
#' Uses MSigDB pathway collections: Hallmark (MH), BIOCARTA, REACTOME, WIKIPATHWAYS,
#' GO:BP, and Immunologic signatures (M7). Creates a dot plot comparing enriched
#' pathways across all provided gene lists.
#'
#' @export
pathway_overrepresentation_analysis_multiple_lists <- function (gene_list, all_genes, local_path, minGSSize = 5, maxGSSize = 400, filename = '', nterms_to_plot = 5, font_size = 8, grouping_var = '', ...)  {

     set_enrichment_color_scale()

     # Convert gene symbols to Entrez IDs for each gene list
     gene_list_entrez <- lapply(gene_list, function(genes) {
          ids <- clusterProfiler::bitr(genes, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
          return(ids$ENTREZID)
     })

     # Convert background genes to Entrez IDs
     all_genes_ids <- clusterProfiler::bitr(all_genes, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
     all_genes_entrez <- all_genes_ids$ENTREZID

     # Create TERM2GENE database from MSigDB pathway collections
     pathway_collections <- list(c('MH'), c('M2', 'CP:BIOCARTA'), c('M2', 'CP:REACTOME'), c('M2', 'CP:WIKIPATHWAYS'), c('M5', 'GO:BP'))
     msigdbr_pathways <- purrr::map_dfr(
          pathway_collections,
          function(pc) {
               if (length(pc) >= 2 && !is.na(pc[2]) && nzchar(pc[2])) {
                    msigdbr::msigdbr(db_species = "MM", species = "mouse", collection = pc[1], subcategory = pc[2])
               } else {
                    msigdbr::msigdbr(db_species = "MM", species = "mouse", collection = pc[1])
               }
          }
     ) |> select(gs_name, ncbi_gene)

     enrichment_results <- clusterProfiler::compareCluster(
          geneCluster = gene_list_entrez,
          fun = "enricher",
          universe = all_genes_entrez,
          pAdjustMethod = "BH",
          minGSSize    = minGSSize,
          maxGSSize    = maxGSSize,
          qvalueCutoff = 0.2,
          pvalueCutoff = 0.05,
          TERM2GENE = msigdbr_pathways
     )

     enrichment_results <- enrichment_results |>
          mutate(Cluster = if (all(!is.na(suppressWarnings(as.numeric(levels(factor(Cluster))))))) {
               fct_inseq(Cluster)
          } else {
               factor(Cluster)
          })

     enrichment_results_table <- as_tibble(enrichment_results)
     write.csv(enrichment_results_table, here::here(local_path, paste0(filename,'Pathway_OverRepresentation_analysis_gene_lists_results.csv')))

     if (nrow(enrichment_results_table) > 1) {
          p1 <- enrichplot::dotplot(enrichment_results,
               showCategory=nterms_to_plot,
               title = paste0(filename,'Pathway Overrepresentation analysis'),
               label_format = 60,
               font.size = font_size) +
               labs(x = grouping_var) +
               theme(axis.text.x = element_text(angle = 45, hjust = 1))
          ggsave(plot = p1, filename = paste0(filename, 'Pathway_overrepresentation_analysis_dotplot.pdf'), width = 10, height = 18, path = local_path)
     }
     return(list(plot = p1, results = enrichment_results))
}
