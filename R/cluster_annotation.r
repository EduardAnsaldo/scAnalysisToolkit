## Cluster Annotation Functions
## Functions for annotating clusters and identifying marker genes

#' Annotate Seurat Object with SingleR
#'
#' Performs automated cell type annotation using SingleR reference databases.
#' Supports both cluster-level and cell-level annotation with fine or coarse
#' granularity. Generates UMAP visualizations of annotations and optionally
#' splits by experimental groups.
#'
#' @param seurat Seurat object to annotate
#' @param local_path Character; path to save annotation plots
#' @param database Character; reference database to use. Currently only 'ImmGen'
#'   is supported. Default c("ImmGen")
#' @param annotation_basis Character; annotation strategy. Options:
#'   \itemize{
#'     \item "cluster_fine": Fine-grained annotation at cluster level
#'     \item "cluster_coarse": Coarse annotation at cluster level
#'     \item "cell_coarse": Coarse annotation at single-cell level
#'     \item "cell_fine": Fine-grained annotation at single-cell level
#'   }
#'   Default c("cluster_fine", "cell_coarse", "cell_fine")
#' @param split_by_groups Logical; if TRUE, creates additional plots split by
#'   'Groups' metadata column. Default TRUE
#'
#' @return Seurat object with annotations added to metadata. Annotation columns:
#'   \itemize{
#'     \item labels_per_cluster_fine: Fine cluster annotations
#'     \item labels_per_cluster_coarse: Coarse cluster annotations
#'     \item labels_per_cell_coarse: Coarse cell-level annotations
#'     \item labels_per_cell_fine: Fine cell-level annotations
#'   }
#'
#' @export
annotate_seurat_with_SingleR_Eduard <- function(
    seurat,
    local_path,
    database = c("ImmGen"),
    annotation_basis = c("cluster_fine", "cell_coarse", "cell_fine"),
    split_by_groups = TRUE
) {
    # Select database
    if (database == "ImmGen") {
        ref <- celldex::ImmGenData(ensembl = FALSE)
    } else {
        stop("Currently only 'ImmGen' database is supported.")
    }

    DefaultAssay(seurat) <- 'RNA'

    # Annotation logic
    if (annotation_basis == "cluster_fine") {
        predictions <- SingleR::SingleR(
            test  =  as.SingleCellExperiment(seurat),
            assay.type.test = 1,
            ref = ref,
            labels = ref$label.fine,
            cluster = seurat$seurat_clusters
        )
        row.names <- rownames(predictions)
        predictions_tbl <- predictions |>
            as_tibble() |>
            dplyr::select(labels)
        predictions_tbl$cluster <- row.names
        annotations <- seurat@meta.data |>
            left_join(predictions_tbl, by = join_by('seurat_clusters' == 'cluster')) |>
            pull(labels)
        seurat$labels_per_cluster_fine <- annotations
        Idents(seurat) <- 'labels_per_cluster_fine'
        p <- DimPlot(seurat, label = FALSE, label.size = 2.5)
        print(p)
        ggsave(plot = p, filename = paste0('UMAP_cluster_SingleR_annotations_fine','.pdf'), path = local_path, width = 8, height = 5)
        if (split_by_groups) {
            p1 <- DimPlot(seurat, label = TRUE, label.size = 2.5, split.by = 'Groups')
            ggsave(plot = p1, filename = paste0('UMAP_cluster_SingleR_annotations_fine_by_group','.pdf'), path = local_path, width = 10, height = 5)
        }
    } else if (annotation_basis == "cluster_coarse") {
        predictions <- SingleR::SingleR(
            test = as.SingleCellExperiment(seurat),
            assay.type.test = 1,
            ref = ref,
            labels = ref$label.main,
            cluster = seurat$seurat_clusters
        )
        row.names <- rownames(predictions)
        predictions_tbl <- predictions |>
            as_tibble() |>
            dplyr::select(labels)
        predictions_tbl$cluster <- row.names
        annotations <- seurat@meta.data |>
            left_join(predictions_tbl, by = join_by('seurat_clusters' == 'cluster')) |>
            pull(labels)
        seurat$labels_per_cluster_coarse <- annotations
        Idents(seurat) <- 'labels_per_cluster_coarse'
        p <- DimPlot(seurat, label = TRUE, label.size = 2.5)
        print(p)
        ggsave(plot = p, filename = paste0('UMAP_cluster_SingleR_annotations_coarse','.pdf'), path = local_path, width = 8, height = 5)
        if (split_by_groups) {
            p1 <- DimPlot(seurat, label = TRUE, label.size = 2.5, split.by = 'Groups')
            ggsave(plot = p1, filename = paste0('UMAP_cluster_SingleR_annotations_coarse_by_group','.pdf'), path = local_path, width = 10, height = 5)
        }
    } else if (annotation_basis == "cell_coarse") {
        predictions <- SingleR::SingleR(
            test = as.SingleCellExperiment(seurat),
            assay.type.test = 1,
            ref = ref,
            labels = ref$label.main
        )
        predictions_tbl <- predictions |>
            as_tibble() |>
            dplyr::select(labels) |>
            rename(labels_per_cell_coarse = labels)
        seurat$labels_per_cell_coarse <- predictions_tbl |> pull(labels_per_cell_coarse)
        Idents(seurat) <- 'labels_per_cell_coarse'
        p <- DimPlot_scCustom(seurat, label = FALSE)
        print(p)
        ggsave(plot = p, filename = paste0('UMAP_cell_SingleR_annotations_coarse','.pdf'), path = local_path, width = 5, height = 5)
        if (split_by_groups) {
            p1 <- DimPlot_scCustom(seurat, label = FALSE, split.by = 'Groups')
            ggsave(plot = p1, filename = paste0('UMAP_cell_SingleR_annotations_coarse_by_group','.pdf'), path = local_path, width = 6, height = 5)
        }
    } else if (annotation_basis == "cell_fine") {
        predictions <- SingleR::SingleR(
            test = as.SingleCellExperiment(seurat),
            assay.type.test = 1,
            ref = ref,
            labels = ref$label.fine
        )
        predictions_tbl <- predictions |>
            as_tibble() |>
            dplyr::select(labels) |>
            rename(labels_per_cell_fine = labels)
        seurat$labels_per_cell_fine <- predictions_tbl |> pull(labels_per_cell_fine)
        Idents(seurat) <- 'labels_per_cell_fine'
        p <- DimPlot_scCustom(seurat, label = FALSE)
        print(p)
        ggsave(plot = p, filename = paste0('UMAP_cell_SingleR_annotations_fine','.pdf'), path = local_path, width = 26, height = 5)
        if (split_by_groups) {
            p1 <- DimPlot_scCustom(seurat, label = FALSE, split.by = 'Groups')
            ggsave(plot = p1, filename = paste0('UMAP_cell_SingleR_annotations_fine_by_group','.pdf'), path = local_path, width = 30, height = 5)
        }
    } else {
        stop("annotation_basis must be one of 'cluster_fine', 'cluster_coarse', 'cell_coarse', or 'cell_fine'.")
    }

    DefaultAssay(seurat) <- 'SCT'
    return(seurat)
}


#' Find and Visualize Top Marker Genes Per Cluster
#'
#' Identifies marker genes for each cluster using Seurat's FindAllMarkers and
#' generates a dot plot visualization. Saves ranked gene lists (top 10, 50, 100)
#' to TSV files and optionally performs pathway enrichment analysis using
#' Metascape and/or ClusterProfiler.
#'
#' @param seurat Seurat object with cluster assignments
#' @param n_genes_to_plot Integer; number of top genes per cluster to include in
#'   dot plot. Default 3
#' @param grouping_var Character; metadata column name for cluster identity.
#'   Default 'seurat_clusters'
#' @param object_annotations Character; string to append to output file names.
#'   Default ''
#' @param tables_path Character; path to save output tables. Default 'results/tables/'
#' @param figures_path Character; path to save figures. Default 'results/figures/'
#' @param results_path Character; path for pathway enrichment results. Default 'results/'
#' @param run_pathway_enrichment Character vector; enrichment methods to run
#'   ('Metascape', 'ClusterProfiler', or NULL). Default NULL
#' @param filter_ig Logical; whether to filter out immunoglobulin genes (Igh, Igl, Igk)
#'   from FindAllMarkers results. Default FALSE
#' @param filter_tcr Logical; whether to filter out T cell receptor genes (Trav, Trbv,
#'   Trdv, Trgv, Traj, Trbj, Trdj, Trgj, Trac, Trbc, Trdc, Trgc) from FindAllMarkers
#'   results. Default FALSE
#' @param ... Additional arguments passed to ClusterProfiler analysis
#'
#' @return List with elements:
#'   \item{plot}{ggplot; dot plot of top marker genes}
#'   \item{ClusterProfiler_results}{ClusterProfiler enrichment results or NULL}
#'   \item{metascape_results}{Metascape enrichment results or NULL}
#'   \item{topn}{Data frame; top n markers per cluster}
#'   \item{top100}{Data frame; top 100 markers per cluster}
#'
#' @export
top_genes_per_cluster <- function (seurat, n_genes_to_plot = 3, grouping_var = 'seurat_clusters', object_annotations = '', tables_path = 'results/tables/', figures_path = 'results/figures/', results_path = 'results/', run_pathway_enrichment = NULL, filter_ig = FALSE, filter_tcr = FALSE, ...) {

    sequential_palette_dotplot <- grDevices::hcl.colors(n = 20,'YlGn',rev = T)

    # Set the identity class for clustering
    Idents(seurat) <- grouping_var

    seurat.markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

    # Filter out Ig genes if requested
    if (filter_ig) {
        message("Filtering out Ig genes from FindAllMarkers results...")
        n_before <- nrow(seurat.markers)
        seurat.markers <- seurat.markers |>
            filter(!str_detect(gene, "^Igh|^Igl|^Igk"))
        n_after <- nrow(seurat.markers)
        message(sprintf("Removed %d Ig genes (%d -> %d markers)", n_before - n_after, n_before, n_after))
    }

    # Filter out TCR genes if requested
    if (filter_tcr) {
        message("Filtering out TCR genes from FindAllMarkers results...")
        n_before <- nrow(seurat.markers)
        seurat.markers <- seurat.markers |>
            filter(!str_detect(gene, "^Trav|^Trbv|^Trdv|^Trgv|^Traj|^Trbj|^Trdj|^Trgj|^Trac|^Trbc|^Trdc|^Trgc"))
        n_after <- nrow(seurat.markers)
        message(sprintf("Removed %d TCR genes (%d -> %d markers)", n_before - n_after, n_before, n_after))
    }

    # saveRDS(seurat.markers, file = 'seurat.markers.rds')



    #Add gene annotations:
    seurat.markers <- seurat.markers |>
                    left_join(y= unique(annotations[,c('gene_name', 'description')]),
                        by = c('gene' = 'gene_name')) 

    
    # Only convert to ordered factor if all cluster levels can be coerced to numeric
    if (all(!is.na(suppressWarnings(as.numeric(as.character(unique(seurat.markers$cluster))))))) {
        seurat.markers <- seurat.markers |>
            mutate(cluster = fct_inseq(cluster))
    }

    #Top10 markers
    seurat.markers |>
        group_by(cluster) |>
        arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) |>
        slice_head(n = 10) -> top10

    #Top25 markers
    seurat.markers |>
        group_by(cluster) |>
        arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) |>
        slice_head(n = 50) -> top50

    #Top100 markers
    seurat.markers |>
        group_by(cluster) |>
        arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) |>
        slice_head(n = 100) -> top100

    #Topn markers
    seurat.markers |>
        group_by(cluster) |>
        arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) |>
        slice_head(n = n_genes_to_plot) -> topn

    # Save the top markers to files
    write.table(top100,file=here::here(tables_path, paste0('top100', '_',object_annotations, ".tsv")), sep="\t",row.names = FALSE)
    # write.table(top25,file=here::here(path,'top25',object_annotations, ".tsv"), sep="\t",row.names = FALSE)
    write.table(top10,file=here::here(tables_path, paste0('top10', '_',object_annotations, ".tsv")), sep="\t",row.names = FALSE)

    top100_genes_per_cluster <- top100 |>
        group_by(cluster) |>
        summarise(genes = str_flatten_comma(gene))

    write.table(top100_genes_per_cluster,
                file = here::here(tables_path, paste0('top100_gene_names_per_cluster_', object_annotations, ".tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)

    gene_list_plot <- topn |> pull(gene)

    gene_list_plot <- gene_list_plot |> unique() |> rev()
    plot1 <- DotPlot_scCustom(seurat,
                    features = gene_list_plot,
                    colors_use = hcl.colors(12, palette = "RdBu", rev = TRUE),
                    flip_axes = T,
                    dot.scale = 8,
                    dot.min = 0,
                    scale.min = 0,
                    scale.max = 80,
                    x_lab_rotate = T,
                    y_lab_rotate = F) +
        theme(axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            legend.title = element_text(size = 18))

    metascape_results <- NULL
    ClusterProfiler_results <- NULL
        # Run pathway enrichment analysis
    if (!is.null(run_pathway_enrichment)) {
        if ('Metascape' %in% run_pathway_enrichment) {
            metascape_results <- Metascape_functional_analysis_cluster_identification(seurat, top100, identities = 'seurat_clusters', path=results_path, object_annotations = object_annotations)
        }
        if ('ClusterProfiler' %in% run_pathway_enrichment) {
            ClusterProfiler_results <- GO_functional_analysis_cluster_identification(seurat, seurat.markers, path=results_path, object_annotations = object_annotations, top_gene_number = 100, ...)
        }

    }

    return(list(plot = plot1, ClusterProfiler_results = ClusterProfiler_results, metascape_results = metascape_results, topn = topn, top100 = top100, seurat.markers = seurat.markers) )

}
