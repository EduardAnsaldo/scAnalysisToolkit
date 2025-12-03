## Metascape Pathway Enrichment Functions
## Functions for web-based pathway enrichment using Metascape API

#' Metascape Overrepresentation Analysis
#'
#' Performs pathway overrepresentation analysis using the Metascape web-based tool
#' (msbio command-line interface). Submits gene list to Metascape, retrieves results,
#' and creates a bar plot of top enriched pathways.
#'
#' @param significant_genes_FC_ordered Character vector; gene symbols ordered by fold change
#' @param local_path Character; path to save Metascape results and plots
#' @param group Character; group identifier for plot title
#' @param grouping_var Character; grouping variable name for plot title and file naming
#' @param nterms_to_plot_metascape Integer; number of top pathways to plot. Default 12
#' @param filename Character; prefix for output files. Default ''
#' @param ... Additional arguments (unused)
#'
#' @return ggplot object with pathway enrichment bar plot, or "N/A" plot if no results
#'
#' @examples
#' Metascape_overrepresentation_analysis(deg_genes, local_path = "./enrichment",
#'                                      group = "treatment", grouping_var = "cluster_1",
#'                                      nterms_to_plot_metascape = 10)
#' @export
Metascape_overrepresentation_analysis <- function (significant_genes_FC_ordered, local_path, group, grouping_var, nterms_to_plot_metascape = 12 , filename = '', ...) {

#     # Create local folder
#     local_path <- here::here(local_path, 'Metascape_overrepresentation_analysis')
#     unlink(local_path, recursive = T)
#     dir.create(local_path, recursive = T, showWarnings = F)
    
    group_kebab_case <- str_replace_all(group, ' ', '_')
    grouping_var_kebab_case <- str_replace_all(grouping_var, ' ', '_')
    filename <- str_c(filename, '_')
    
    # Create folder in msbio_v3-2
    project_name <- here::here() |> basename()
    msbio_path <- paste0('~/msbio_v3-2/', 'data/', project_name, '/', basename(local_path))
    msbio_path_short <- paste0('data/', project_name, '/', basename(local_path))
    dir.create(msbio_path, recursive = T)
    msbio_file <- paste0(msbio_path, '/', filename, grouping_var_kebab_case, '_UP_in_', group_kebab_case, '_gene_list.txt')
    msbio_file_path_short <- paste0(msbio_path_short, '/', filename, grouping_var_kebab_case, '_UP_in_', group_kebab_case, '_gene_list.txt')

    # Write gene list to file in msbio folder
    writeLines(significant_genes_FC_ordered, msbio_file)

    # Run bash command
    bash_command <- paste0("cd ~/msbio_v3-2/; ls; bin/ms.sh -up -t Symbol --source_tax_id 9606 --target_tax_id 9606 -o ", msbio_path_short, " ", msbio_file_path_short)
#     print(bash_command)
    system(bash_command)
    
    # Copy output back to current working directory
    # Note: You may need to adjust this based on what output files are created
     output_files <- list.files(path = msbio_path, full.names = TRUE)
     # Copy output back to current working directory
     copy_success <- file.copy(output_files, local_path, recursive = TRUE)
     if (all(copy_success)) {
          message("Successfully copied ", length(output_files), " file(s) to ", local_path)
     } else {
          warning("Failed to copy some files. Check permissions and paths.")
     }
     unlink(msbio_path, recursive = TRUE) 

   if (!file.exists(here::here(local_path, 'metascape_result.xlsx'))) {
     p1 <- create_no_data_plot()
     return(p1)
   }       
     
  # Read results
  enrichment_results <- readxl::read_excel(here::here(local_path, 'metascape_result.xlsx'), sheet = 'Enrichment')

# Check if there are enrichment results to plot
if (nrow(enrichment_results) > 1) {
     # Filter and prepare for plotting
     enrichment_results <- enrichment_results |>
          filter(str_detect(GroupID, 'Summary') ) |>
          #select(Description, `Log(q-value)`) |>
          select(Description, LogP) |>
          rename(term = Description, log_q_value = LogP)  |>
          mutate(log_q_value = -1*log_q_value) |> 
          # rename(term = Description, log_q_value = LogP)  |>
          head(n = nterms_to_plot_metascape)

     if (nrow(enrichment_results) > 1) {
          
     plot2 <- enrichment_results |>       
          ggplot(aes(x = log_q_value, y = fct_reorder(term, log_q_value , .desc = F), fill = log_q_value)) +
               geom_col(width = 0.7) +
               scale_fill_binned(palette = grDevices::hcl.colors(n = 6, 'YlOrRd', rev = T), breaks = c(2, 4, 6, 10, 20)) +
               labs(x = '-log10(p-value)', y = '', title = paste0('UP in ', group, ' - ', grouping_var)) +
               theme_minimal() +
               scale_y_discrete(labels = scales::label_wrap(22))+
               theme(axis.text.y = element_text(size = 9), title = element_text(size = 16), plot.title.position = 'plot', legend.position = 'none', axis.text.x = element_text(size = 10))
     #    print(plot2)
     ggsave(plot = plot2, filename = paste0(filename, '_Pathway_enrichment_analysis_metascape', '.pdf'), width = 12, height = 6, path = local_path)
     return(plot2)
     } else {
          return(create_no_data_plot())
     }
     } else {
     return(create_no_data_plot())
}}


#' Metascape Functional Analysis on Differential Expression Results
#'
#' Wrapper function that performs Metascape pathway enrichment analysis on differential
#' expression results. Runs enrichment separately for upregulated and downregulated genes,
#' creating organized output directories for each direction.
#'
#' @param results Data frame; differential expression results with columns: genes,
#'   log2FoldChange, padj
#' @param grouping_var Character; identifier for cluster or group being analyzed
#' @param group2 Character; name of treatment group (upregulated genes)
#' @param group1 Character; name of control group (downregulated genes)
#' @param path Character; base output directory path. Default './'
#' @param FC_threshold Numeric; log2 fold change threshold for filtering genes
#' @param p_value_threshold Numeric; adjusted p-value threshold. Default 0.05
#' @param ... Additional arguments passed to Metascape_overrepresentation_analysis
#'
#' @return Invisible NULL; creates output directories with Metascape enrichment results
#'
#' @examples
#' Metascape_functional_analysis(de_results, grouping_var = "cluster_1",
#'                              group2 = "KO", group1 = "WT",
#'                              path = "./results", FC_threshold = 0.5)
#' @export
Metascape_functional_analysis <- function (results, grouping_var, group2, group1, path='./', FC_threshold, p_value_threshold = 0.05, ...) {

    grouping_var_kebab_case <- str_replace_all(grouping_var, ' ', '_')
    group2_kebab_case <- str_replace_all(group2, ' ', '_')
    group1_kebab_case <- str_replace_all(group1, ' ', '_')

     results <- results[which(duplicated(results$genes) == F),]

####################################### UP ########################################

     local_path <- create_analysis_directory(here::here(path, paste0('Functional_analysis_metascape_UP_in_', group2_kebab_case, '_', grouping_var_kebab_case)))
     print(local_path)

     significant_genes <- results |> filter((padj < p_value_threshold) & (log2FoldChange > FC_threshold)) |> arrange(padj) |> arrange(desc(log2FoldChange)) |> pull(genes)
     
     if (length(significant_genes) > 2) {
          p1 <- Metascape_overrepresentation_analysis(significant_genes, local_path =  local_path , group = group2, grouping_var = grouping_var, filename = '', ...)
          print(p1)

     }else {
          p1 <- create_no_data_plot()
          print(p1)
     }

######################################## DOWN ########################################

     local_path <- create_analysis_directory(here::here(path, paste0('Functional_analysis_metascape_DOWN_in', group1_kebab_case, '_', grouping_var_kebab_case)))
     print(local_path)

     significant_genes <- results |> filter((padj < p_value_threshold) & (log2FoldChange < -1*FC_threshold)) |> arrange(padj) |> arrange(log2FoldChange) |> pull(genes)

     if (length(significant_genes) > 2) {
          p1 <- Metascape_overrepresentation_analysis(significant_genes, local_path = local_path , group = group1, grouping_var = grouping_var, filename = '', ...)
          print(p1)
     }
     else {
          p1 <- create_no_data_plot()
          print(p1)
     }
} 



#' Metascape Functional Analysis for Cluster Marker Genes
#'
#' Performs Metascape pathway enrichment analysis on marker genes for each cluster
#' from Seurat FindAllMarkers results. Creates individual enrichment plots for each
#' cluster and combines them into a summary figure.
#'
#' @param seurat Seurat object with cluster assignments
#' @param results Data frame; marker gene results from FindAllMarkers with columns:
#'   gene, cluster
#' @param identities Character; metadata column name for cluster identity.
#'   Default 'seurat_clusters'
#' @param path Character; base output directory path. Default './'
#' @param object_annotations Character; string to append to output directory name. Default ''
#' @param ... Additional arguments passed to Metascape_overrepresentation_analysis
#'
#' @return Patchwork object with combined enrichment plots for all clusters
#'
#' @examples
#' Metascape_functional_analysis_cluster_identification(seurat, marker_results,
#'                                                     identities = "seurat_clusters",
#'                                                     path = "./results")
#' @export
Metascape_functional_analysis_cluster_identification <- function (seurat, results, identities = 'seurat_clusters', path='./', object_annotations = '', ...) {

     local_path <- create_analysis_directory(here::here(path, paste0(object_annotations, '_Cluster_identification_functional_analysis')))
     results <- results |> 
          mutate(cluster = fct_inseq(cluster))

     plot_list <- list()
     for (cluster in levels(results$cluster)) {

          name <- paste0('Cluster_', cluster)
          local_path_cluster <- create_analysis_directory(here::here(local_path, name))

          significant_results_cluster <- results |> 
               filter(cluster == {{cluster}}) |>
               pull(gene) |>
               unique()
     

plot_list[[cluster]] <- Metascape_overrepresentation_analysis(significant_genes_FC_ordered = significant_results_cluster, local_path = local_path_cluster, nterms_to_plot_metascape = 6, group = '', filename =  name,  grouping_var = cluster) 

  }
  plots <- wrap_plots(plot_list) + plot_annotation(title = 'Cluster identification - Metascape functional analysis') & 
    theme(plot.title = element_text(size = 9, hjust = 0.5),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          axis.title.x  = element_text(size = 8),
     axis.title.y  = element_text(size = 8))
  print(plots)
  return(plots)
}
