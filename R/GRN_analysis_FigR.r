#' Convert JASPAR Motif ID to Gene Symbol
#'
#' Retrieves the gene symbol name from a JASPAR motif ID by querying the JASPAR2024 database.
#' The returned name is converted to title case for consistency.
#'
#' @param motif_id Character string containing a JASPAR motif ID (e.g., "MA0003.1")
#'
#' @return Character string with the gene symbol in title case, or NA if motif_id is NA
#'
#' @details
#' This function:
#' \itemize{
#'   \item Connects to the JASPAR2024 SQLite database
#'   \item Retrieves the motif matrix associated with the given ID
#'   \item Extracts and formats the gene name
#'   \item Automatically disconnects from the database
#' }
#'
#' @examples
#' \dontrun{
  #' ID_to_symbol("MA0003.1")  # Returns gene symbol for this motif ID
  #' ID_to_symbol(NA)           # Returns NA
#' }
#'
#' @seealso \code{\link{gene_to_jaspar}} for the reverse operation
ID_to_symbol <- function(motif_id){
  
  # Handle missing input
  if (is.na(motif_id)) {
    return(NA_character_)
  }
  
  # Initialize JASPAR database connection
  JASPAR2024 <- JASPAR2024::JASPAR2024()
  JASPARConnect <- RSQLite::dbConnect(RSQLite::SQLite(), db(JASPAR2024))
  
  # Retrieve motif matrix by ID
  motif_matrix <- TFBSTools::getMatrixByID(
    x = JASPARConnect, 
    ID = motif_id   
  )
  
  # Close database connection
  RSQLite::dbDisconnect(JASPARConnect)
  
  # Extract and format gene name
  return(motif_matrix@name |> stringr::str_to_title())
}


#' Convert Gene Symbol to JASPAR2024 Motif ID(s)
#'
#' Queries the JASPAR2024 database to find transcription factor binding site motifs
#' associated with a given gene symbol. Returns the matrix IDs of matching entries.
#'
#' @param gene_symbol Character string with the gene symbol (e.g., "STAT3", "TP53")
#' @param species Character string specifying the species (default: "Homo sapiens")
#' @param collection Character string specifying JASPAR collection (default: "CORE")
#'
#' @return Character vector of JASPAR matrix IDs (e.g., "MA0003.1"), or NA if no matches found
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Connects to JASPAR2024 SQLite database
#'   \item Searches for motifs matching the gene symbol in vertebrates
#'   \item If no results, retries with uppercase gene symbol
#'   \item Extracts detailed information including species, class, and alternative symbols
#'   \item Returns only the matrix IDs from matching entries
#' }
#'
#' The search is restricted to:
#' \itemize{
#'   \item Taxonomic group: vertebrates
#'   \item Latest versions only (all_versions = FALSE)
#'   \item Specified collection (default: CORE)
#' }
#'
#' @examples
#' \dontrun{
#' # Find JASPAR IDs for STAT3
#' gene_to_jaspar("STAT3")
#' 
#' # Search with specific parameters
#' gene_to_jaspar("TP53", species = "Homo sapiens", collection = "CORE")
#' }
#'
#' @seealso \code{\link{ID_to_symbol}} for converting JASPAR IDs back to gene symbols
gene_to_jaspar <- function(gene_symbol, species = "Homo sapiens", collection = "CORE") {
  
  # Initialize database connection
  JASPAR2024 <- JASPAR2024::JASPAR2024()
  JASPARConnect <- RSQLite::dbConnect(RSQLite::SQLite(), db(JASPAR2024))
  
  # Query JASPAR database with original gene symbol
  pfm_list <- TFBSTools::getMatrixSet(
  JASPARConnect,
  opts = list(
    name = gene_symbol,
    tax_group = 'vertebrates',
    collection = collection,
    all_versions = FALSE
  )
  )
  
  # Retry with uppercase if initial search yields no results
  if (length(pfm_list) == 0) {
  pfm_list <- TFBSTools::getMatrixSet(
    JASPARConnect,
    opts = list(
    name = toupper(gene_symbol),
    tax_group = 'vertebrates',
    collection = collection,
    all_versions = FALSE
    )
  )
  }
  
  # Close database connection
  RSQLite::dbDisconnect(JASPARConnect)
  
  # Return NA if no matches found
  if (length(pfm_list) == 0) {
  message(paste("No JASPAR entry found for gene symbol:", gene_symbol))
  return(NA_character_)
  }
  
  # Extract motif information with safe handling of potentially missing fields
  results <- data.frame(
  jaspar_id = sapply(pfm_list, function(x) TFBSTools::name(x)),
  matrix_id = names(pfm_list),
  gene_symbol = sapply(pfm_list, function(x) {
    sym <- TFBSTools::tags(x)$symbol
    if (is.null(sym) || length(sym) == 0) NA_character_ else as.character(sym)
  }),
  species = sapply(pfm_list, function(x) {
    sp <- TFBSTools::tags(x)$species
    if (is.null(sp) || length(sp) == 0) NA_character_ else as.character(sp)
  }),
  class = sapply(pfm_list, function(x) {
    cls <- TFBSTools::tags(x)$class
    if (is.null(cls) || length(cls) == 0) NA_character_ else as.character(cls)
  }),
  stringsAsFactors = FALSE
  )
  
  # Return only the matrix IDs
  return(results$matrix_id)
}

#' Run FigR DORC Analysis on Seurat Object. Function written by hand, refactored using Claude 
#'
#' @param seurat Seurat object with ATAC assay
#' @param rna_assay Name of RNA assay (default: "RNA")
#' @param atac_assay Name of ATAC assay (default: "ATAC")
#' @param reduction Reduction to use for neighbors (default: "lsi")
#' @param dims Dimensions to use for neighbor finding (default: 2:12)
#' @param n_neighbors Number of neighbors for smoothing (default: 20)
#' @param genome Genome assembly (default: "mm10")
#' @param n_cores Number of cores for parallel processing (default: 8)
#' @param p_cutoff P-value cutoff for filtering correlations (default: 0.05)
#' @param dorc_cutoff Cutoff for DORC selection in dorcJPlot (default: 7)
#' @param dorc_k_pct Percentage of cells for dorcK parameter (default: 0.03)
#' @param n_bg Number of background peaks for correlation (default: 100)
#' @return List containing filtered correlations, DORC matrix, smoothed matrices, and FigR results
run_figr_analysis <- function(seurat,
                 rna_assay = "RNA",
                 atac_assay = "ATAC",
                 reduction = "lsi",
                 dims = 2:12,
                 n_neighbors = 20,
                 genome = "mm10",
                 n_cores = 8,
                 p_cutoff = 0.05,
                 dorc_cutoff = 7,
                 dorc_k_pct = 0.03,
                 n_bg = 100) {
  
  # Step 1: Filter to standard chromosomes
  message("Filtering to standard chromosomes...")
  DefaultAssay(seurat) <- atac_assay
  features.keep <- as.character(GenomeInfoDb::seqnames(granges(seurat))) %in% GenomeInfoDb::standardChromosomes(granges(seurat))
  
  seurat@assays[[atac_assay]]@counts <- seurat@assays[[atac_assay]]@counts[features.keep, ]
  seurat@assays[[atac_assay]]@data <- seurat@assays[[atac_assay]]@data[features.keep, ]
  seurat@assays[[atac_assay]]@meta.features <- seurat@assays[[atac_assay]]@meta.features[features.keep, ]
  
  # Get standard chromosomes as a vector
  standard_chroms <- GenomeInfoDb::standardChromosomes(granges(seurat))

  # Create a logical vector for subsetting
  gr_seqnames <- as.character(GenomeInfoDb::seqnames(gr))
  keep_idx <- gr_seqnames %in% standard_chroms

  # Subset using the logical vector
  gr <- gr[keep_idx]

  seurat@assays[[atac_assay]]@annotation <- gr

  # # Update annotation
  # gr <- seurat@assays[[atac_assay]]@annotation
  # seurat@assays[[atac_assay]]@annotation <- gr[GenomeInfoDb::seqnames(gr) %in% 
  #                         GenomeInfoDb::standardChromosomes(granges(seurat))]
  
  # Fix seqinfo
  gr <- seurat@assays[[atac_assay]]@seqinfo
  standard_chroms <- GenomeInfoDb::standardChromosomes(granges(seurat)) |> as.vector()
  standard_chroms_short <- stringr::str_remove(standard_chroms, 'chr')
  GenomeInfoDb::seqlevels(gr) <- standard_chroms_short
  seurat@assays[[atac_assay]]@seqinfo <- gr[standard_chroms_short]
  
  # Update ranges
  gr <- seurat@assays[[atac_assay]]@ranges
  # gr <- gr[GenomeInfoDb::seqnames(gr) %in% GenomeInfoDb::standardChromosomes(granges(seurat))]
  
  # Create a logical vector for subsetting
  gr_seqnames <- as.character(GenomeInfoDb::seqnames(gr))
  keep_idx <- gr_seqnames %in% standard_chroms

  # Subset using the logical vector
  gr <- gr[keep_idx]
  
  sequence_levels <- GenomeInfoDb::seqlevels(gr)
  newStyle <- GenomeInfoDb::mapSeqlevels(sequence_levels, style = "UCSC")
  newStyle <- newStyle[complete.cases(newStyle)]
  gr <- GenomeInfoDb::renameSeqlevels(gr, newStyle)
  gr <- GenomeInfoDb::keepStandardChromosomes(gr)
  seurat@assays[[atac_assay]]@ranges <- gr
  
  # Update variable features
  VariableFeatures(seurat) <- VariableFeatures(seurat)[features.keep]
  
  # Step 2: Convert to SingleCellExperiment
  message("Converting to SingleCellExperiment...")
  chromatin_assay <- seurat[[atac_assay]]
  peak_ranges <- granges(chromatin_assay)
  ATAC_se <- Seurat::as.SingleCellExperiment(seurat, assay = atac_assay)
  SummarizedExperiment::rowRanges(ATAC_se) <- peak_ranges
  
  # Preserve genomic information
  sequence_levels <- GenomeInfoDb::seqlevels(GenomeInfoDb::seqinfo(chromatin_assay)) |> 
  GenomeInfoDb::mapSeqlevels(style = "UCSC")
  sequence_information <- GenomeInfoDb::seqinfo(chromatin_assay) |> 
  GenomeInfoDb::renameSeqlevels(sequence_levels)
  GenomeInfoDb::genome(sequence_information) <- genome
  GenomeInfoDb::seqinfo(SummarizedExperiment::rowRanges(ATAC_se)) <- sequence_information
  S4Vectors::metadata(ATAC_se)$annotation <- Signac::Annotation(chromatin_assay)
  S4Vectors::metadata(ATAC_se)$fragments <- Signac::Fragments(chromatin_assay)
  
  # Step 3: Prepare RNA data
  message("Preparing RNA data...")
  seurat <- JoinLayers(seurat, assay = rna_assay)
  seurat <- NormalizeData(seurat, assay = rna_assay, 
              normalization.method = 'LogNormalize', 
              scale.factor = 10000)
  RNA_matrix <- LayerData(seurat, assay = rna_assay, layer = 'data')
  
  # Step 4: Find neighbors
  message("Finding neighbors...")
  seurat_neighbors <- FindNeighbors(seurat, reduction = reduction, dims = dims,
                  assay = atac_assay, graph.name = 'ATAC_figr',
                  return.neighbor = TRUE)
  cellkNN <- seurat_neighbors@neighbors$ATAC_figr@nn.idx
  rownames(cellkNN) <- seurat_neighbors@neighbors$ATAC_figr@cell.names
  
  # Step 5: Run peak-gene correlation
  message("Running peak-gene correlation...")
  cisCorrelation <- FigR::runGenePeakcorr(
    ATAC.se = ATAC_se,
    RNAmat = RNA_matrix,
    genome = genome,
    nCores = n_cores,
    p.cut = NULL,
    n_bg = n_bg
  )
  
  # Step 6: Filter correlations and identify DORCs
  message("Identifying DORCs...")
  cisCorrelation_filtered <- cisCorrelation |> dplyr::filter(pvalZ < p_cutoff)
  dorcGenes <- FigR::dorcJPlot(cisCorrelation_filtered,
                 cutoff = dorc_cutoff,
                 labelTop = 25,
                 returnGeneList = TRUE)
  
  # Step 7: Calculate DORC scores
  message("Calculating DORC scores...")
  dorcMat <- FigR::getDORCScores(ATAC.se = ATAC_se,
                 dorcTab = cisCorrelation_filtered,
                 geneList = dorcGenes,
                 nCores = n_cores)
  
  # Step 8: Smooth scores
  message("Smoothing scores...")
  dorcMat.s <- FigR::smoothScoresNN(NNmat = cellkNN[, 1:n_neighbors], 
                  mat = dorcMat, 
                  nCores = n_cores)
  RNAmat.s <- FigR::smoothScoresNN(NNmat = cellkNN[, 1:n_neighbors], 
                   mat = RNA_matrix, 
                   nCores = n_cores)
  
  # Step 9: Run FigR GRN
  message("Running FigR GRN...")
  kNN_K <- ceiling(nrow(dorcMat) * dorc_k_pct)
  figR.d <- FigR::runFigRGRN(ATAC.se = ATAC_se,
               dorcK = kNN_K,
               dorcTab = cisCorrelation_filtered,
               genome = genome,
               dorcMat = dorcMat.s,
               rnaMat = RNAmat.s,
               nCores = n_cores)
  
  message("Analysis complete!")
  
  # Return results
  list(
    seurat = seurat,
    ATAC_se = ATAC_se,
    cisCorrelation = cisCorrelation_filtered,
    dorcGenes = dorcGenes,
    dorcMat = dorcMat,
    dorcMat_smoothed = dorcMat.s,
    RNAmat_smoothed = RNAmat.s,
    figR_results = figR.d
  )
}

#' Plot DORC and RNA expression for genes
#'
#' @param results Output from run_figr_analysis()
#' @param genes Character vector of genes to plot (default: top 15 valid DORCs)
#' @param umap_reduction Name of UMAP reduction in ATAC_se (default: "UMAP.WNN")
#' @return List of combined ggplot objects
plot_dorc_genes <- function(results, 
              genes = NULL, 
              umap_reduction = "UMAP.WNN") {
  
  # Get UMAP coordinates
  umap.d <- SingleCellExperiment::reducedDims(results$ATAC_se)[[umap_reduction]] |> as.data.frame()
  
  # Select genes if not provided
  if (is.null(genes)) {
    valid_genes <- results$dorcGenes[!grepl("^[0-9]", results$dorcGenes)]
    genes <- head(valid_genes, 15)
  }
  
  # Create plots for each gene
  plots <- purrr::map(genes, function(gene) {
  dorc_plot <- FigR::plotMarker2D(
    umap.d,
    results$dorcMat_smoothed,
    markers = gene,
    maxCutoff = "q0.99",
    colorPalette = "brewer_heat"
  ) + ggplot2::ggtitle(paste(gene, "DORC"))
  
  rna_plot <- FigR::plotMarker2D(
    umap.d,
    results$RNAmat_smoothed,
    markers = gene,
    maxCutoff = "q0.99",
    colorPalette = "brewer_purple"
  ) + ggplot2::ggtitle(paste(gene, "RNA"))
  
  dorc_plot + rna_plot
  })
  
  plots
}
#' Analyze and visualize FigR results
#'
#' @param output_dir Directory to save plots (default: "figr_plots")
#'
#' @param figR.d FigR results data frame
#' @param seurat Seurat object with ATAC, RNA/SCT, and chromvar assays
#' @param dorcMat.s DORC matrix (sparse matrix)
#' @param RNAmat.s RNA matrix (sparse matrix)
#' @param umap.d UMAP coordinates data frame
#' @param score_cuts Named list of score cutoffs for different plots (default: list(scatter = 1, heatmap1 = 2, heatmap2 = 1.5, network = 1.6))
#' @param drivers Character vector of driver TF gene names to analyze
#' @param genes Character vector of genes to create DORC plots for
#' @param top_dorcs_genes Character vector of top DORC
#' @param umap_reduction Name of UMAP reduction to use (default: "umap.atac")
#' @param species Species name for JASPAR lookup (default: "Mus musculus")
#'
#' @return List containing plots and results
#' @export
analyze_figr_results <- function(figR.d,
                  seurat,
                  dorcMat.s,
                  RNAmat.s,
                  umap.d,
                  top_dorcs_genes,
                  output_dir = here("results/figr_plots"),
                  score_cuts = list(scatter = 1, heatmap1 = 2, heatmap2 = 1.5, network = 1.6),
                  drivers = NULL,
                  genes = NULL,
                  cluster_ident = "seurat_clusters_atac",
                  umap_reduction = "umap.atac",
                  species = "Mus musculus") {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  }

  # Create a text-only plot showing the cluster identifier
  text_plot <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = paste('Results for', cluster_ident), 
             size = 8, fontface = "bold") +
    theme_void() +
    theme(plot.margin = margin(20, 20, 20, 20))
  
  print(text_plot)
  ggsave(here(output_dir, "cluster_info.png"), text_plot, width = 8, height = 2, dpi = 300)

  
  results <- list()
  
  # 1. Scatter plot of FigR results
  tryCatch({
  results$scatter_plot <- figR.d |>
    ggplot(aes(Corr.log10P, Enrichment.log10P, color = Score)) +
    ggrastr::geom_point_rast(size = 0.01, shape = 16) +
    theme_classic() +
    scale_color_gradientn(
    colours = jdb_palette("solar_extra"),
    limits = c(-3, 3),
    oob = scales::squish,
    breaks = scales::breaks_pretty(n = 3)
    )
  print(results$scatter_plot)
  ggsave(here(output_dir, "scatter_plot.png"), results$scatter_plot, width = 8, height = 6, dpi = 300)
  message("✓ Scatter plot completed")
  }, error = function(e) {
  warning("Failed to create scatter plot: ", e$message)
  results$scatter_plot <<- NULL
  })
  
  # 2. Rank drivers
  tryCatch({
  results$rank_by_mean <- rankDrivers(figR.d, rankBy = "meanScore", interactive = FALSE)
  print(results$rank_by_mean)
  ggsave(here(output_dir, "rank_by_mean.png"), results$rank_by_mean, width = 10, height = 8, dpi = 300)
  message("✓ Rank by mean completed")
  }, error = function(e) {
  warning("Failed to rank drivers by mean: ", e$message)
  results$rank_by_mean <<- NULL
  })
  
  tryCatch({
  results$rank_by_targets <- rankDrivers(figR.d, score.cut = score_cuts$scatter, rankBy = "nTargets", interactive = TRUE)
  print(results$rank_by_targets)
  message("✓ Rank by targets completed")
  }, error = function(e) {
  warning("Failed to rank drivers by targets: ", e$message)
  results$rank_by_targets <<- NULL
  })
  
  # 3. Heatmaps
  tryCatch({
  results$heatmap1 <- plotfigRHeatmap_EA(
    figR.d = figR.d,
    score.cut = score_cuts$heatmap1,
    column_names_gp = gpar(fontsize = 12),
    row_names_gp = gpar(fontsize = 12),
    show_row_dend = FALSE
  )
  png(here(output_dir, "heatmap1.png"), width = 10, height = 8, units = "in", res = 300)
  print(results$heatmap1)
  dev.off()
  message("✓ Heatmap 1 completed")
  }, error = function(e) {
  warning("Failed to create heatmap1: ", e$message)
  results$heatmap1 <<- NULL
  if (dev.cur() > 1) dev.off()
  })
  
  tryCatch({
  results$heatmap2 <- plotfigRHeatmap_EA(
    figR.d = figR.d,
    score.cut = score_cuts$heatmap2,
    column_names_gp = gpar(fontsize = 10),
    row_names_gp = gpar(fontsize = 10),
    show_row_dend = FALSE
  )
  png(here(output_dir, "heatmap2.png"), width = 10, height = 8, units = "in", res = 300)
  print(results$heatmap2)
  dev.off()
  message("✓ Heatmap 2 completed")
  }, error = function(e) {
  warning("Failed to create heatmap2: ", e$message)
  results$heatmap2 <<- NULL
  if (dev.cur() > 1) dev.off()
  })
  
  tryCatch({
  results$heatmap_clustered <- plotfigRHeatmap_EA(
    figR.d = figR.d,
    score.cut = score_cuts$scatter,
    TFs = NULL,
    DORCs = NULL,
    column_names_gp = gpar(fontsize = 8.5),
    row_names_gp = gpar(fontsize = 8.5),
    show_row_dend = FALSE,
    cluster_columns = TRUE,
    cluster_rows = TRUE
  )
  png(here(output_dir, "heatmap_clustered.png"), width = 10, height = 8, units = "in", res = 300)
  print(results$heatmap_clustered)
  dev.off()
  message("✓ Clustered heatmap completed")
  }, error = function(e) {
  warning("Failed to create clustered heatmap: ", e$message)
  results$heatmap_clustered <<- NULL
  if (dev.cur() > 1) dev.off()
  })
  
  # 4. Network view
  tryCatch({
  results$network <- plotfigRNetwork(
    figR.d,
    score.cut = score_cuts$network,
    TFs = NULL,
    weight.edges = TRUE
  )
  png(here(output_dir, "network.png"), width = 10, height = 8, units = "in", res = 300)
  print(results$network)
  dev.off()
  message("✓ Network plot completed")
  }, error = function(e) {
  warning("Failed to create network plot: ", e$message)
  results$network <<- NULL
  if (dev.cur() > 1) dev.off()
  })
  
  # 5. DORC plots for specified genes
  if (!is.null(genes)) {
  tryCatch({
    results$dorc_plots <- map(genes, function(gene) {
    tryCatch({
      dorc_plot <- plotMarker2D(
      umap.d,
      dorcMat.s,
      markers = gene,
      maxCutoff = "q0.99",
      colorPalette = "brewer_heat"
      ) + ggtitle(paste(gene, "DORC"))
      
      rna_plot <- plotMarker2D(
      umap.d,
      RNAmat.s,
      markers = gene,
      maxCutoff = "q0.99",
      colorPalette = "brewer_purple"
      ) + ggtitle(paste(gene, "RNA"))
      
      combined_plot <- dorc_plot + rna_plot
      print(combined_plot)
      ggsave(here(output_dir, paste0("dorc_", gene, ".png")), combined_plot, width = 12, height = 5, dpi = 300)
      return(combined_plot)
    }, error = function(e) {
      warning("Failed to create DORC plot for gene ", gene, ": ", e$message)
      return(NULL)
    })
    })
    message("✓ DORC plots completed")
  }, error = function(e) {
    warning("Failed to create DORC plots: ", e$message)
    results$dorc_plots <<- NULL
  })
  }
  
  # 6. TF Chromvar visualization
  tryCatch({
  figR.summ <- figR.d %>% 
    group_by(Motif) %>%
    dplyr::summarise(Score = mean(Score)) %>%
    arrange(desc(Score)) %>%
    mutate(Motif = factor(Motif, levels = as.character(Motif)))
  
  drivers_UP <- figR.summ |> pull(Motif) |> head(8)
  drivers_DOWN <- figR.summ |> pull(Motif) |> tail(4)
  
  drivers_names <- c(drivers_UP, drivers_DOWN)
  drivers_ID <- map(drivers_names, ~ gene_to_jaspar(.x, species = species, collection = "CORE"))
  drivers_ID <- list_c(drivers_ID)
  drivers_ID <- drivers_ID[!is.na(drivers_ID)]
  drivers_names <- map(drivers_ID, ID_to_symbol)
  
  results$driver_plots <- map(drivers_ID, function(.x) {
    tryCatch({
    DefaultAssay(seurat) <- "chromvar"
    p1 <- FeaturePlot_scCustom(seurat, features = .x, reduction = 'umap.wnn', colors_use = sequential_palette) +
      labs(title = paste0(drivers_names[[which(drivers_ID == .x)]], ' Chromvar score')) + NoAxes()
    
    DefaultAssay(seurat) <- "SCT"
    gene_name <- drivers_names[[which(drivers_ID == .x)]]
    if (gene_name %in% rownames(seurat[["SCT"]])) {
      p2 <- FeaturePlot_scCustom(seurat, features = gene_name, reduction = 'umap.wnn', colors_use = sequential_palette) +
      labs(title = paste0(gene_name, ' expression')) + NoAxes()
    } else {
      p2 <- ggplot() + theme_void()
    }
    
    combined_plot <- p1 + p2
    print(combined_plot)
    ggsave(here(output_dir, paste0("driver_", gene_name, ".png")), combined_plot, width = 12, height = 5, dpi = 300)
    return(combined_plot)
    }, error = function(e) {
    warning("Failed to create driver plot for ", .x, ": ", e$message)
    return(NULL)
    })
  })
  message("✓ Driver plots completed")
  }, error = function(e) {
  warning("Failed to create driver plots: ", e$message)
  results$driver_plots <<- NULL
  })
  
  # 7. Add DORCs to seurat and find markers
  tryCatch({
  seurat[['DORCs']] <- CreateAssayObject(counts = dorcMat.s)
  Idents(seurat) <- cluster_ident
  message("✓ DORCs added to Seurat object")
  }, error = function(e) {
  warning("Failed to add DORCs to Seurat object: ", e$message)
  })
  
  # 8. Create detailed DORC plots with coverage
  tryCatch({
  available_genes <- Annotation(seurat[['ATAC']])$gene_name
  top_dorcs_genes_filtered <- top_dorcs_genes[top_dorcs_genes %in% available_genes]
  
  missing_genes <- setdiff(top_dorcs_genes, top_dorcs_genes_filtered)
  if (length(missing_genes) > 0) {
    message("Genes not found in ATAC annotation: ", paste(missing_genes, collapse = ", "))
  }
  
  top_dorcs_genes_filtered <- head(top_dorcs_genes_filtered, 6)
  
  results$detailed_dorc_plots <- map(top_dorcs_genes_filtered, function(gene) {
    tryCatch({
    DefaultAssay(seurat) <- 'DORCs'
    p1 <- FeaturePlot_scCustom(seurat, feature = gene, reduction = umap_reduction, colors_use = sequential_palette, pt.size = 1.3) +
      labs(title = paste0('DORC: ', gene)) + NoAxes() + NoLegend()
    
    DefaultAssay(seurat) <- 'SCT'
    p2 <- FeaturePlot_scCustom(seurat, feature = gene, reduction = umap_reduction, colors_use = sequential_palette, pt.size = 1.3) +
      labs(title = paste0('RNA: ', gene)) + NoAxes() + NoLegend()
    
    DefaultAssay(seurat) <- 'ATAC'
    if (gene %in% top_dorcs_genes_filtered) {
      p3 <- CoveragePlot(
      seurat,
      region = gene,
      group.by = 'seurat_clusters_atac',
      extend.upstream = 50000,
      extend.downstream = 50000
      )
    } else {
      p3 <- ggplot() + theme_void()
    }
    
    combined_plot <- (p1 | p2) / p3 + plot_layout(heights = c(1, 1.3)) + plot_annotation(title = paste0('DORC: ', gene))
    print(combined_plot)
    ggsave(here(output_dir, paste0("detailed_dorc_", gene, ".png")), combined_plot, width = 12, height = 10, dpi = 300)
    return(combined_plot)
    }, error = function(e) {
    warning("Failed to create detailed DORC plot for gene ", gene, ": ", e$message)
    return(NULL)
    })
  })
  message("✓ Detailed DORC plots completed")
  }, error = function(e) {
  warning("Failed to create detailed DORC plots: ", e$message)
  results$detailed_dorc_plots <<- NULL
  })
  
  message("\n=== Analysis complete ===")
  return(results)
}


#' Plot FigR heatmap with modifications
#'
#' Heatmap visualization of TF-DORC associations based on the regulation scores inferred by FigR
#'@param figR.d data.frame of results returned by \code{\link[FigR]{runFigRGRN}}).
#'@param score.cut numeric specifying the absolute regulation score to threshold TF-DORC connections on. Default is 1
#'@param DORCs character specifying valid DORC gene symbols to subset heatmap to. Default is NULL (no subsetting)
#'@param TFs character specifying valid TF gene symbols to subset heatmap to. Default is NULL (no subsetting)
#'@param ... additional parameters passed to the \code{\link[ComplexHeatmap]{Heatmap}})
#'@return a TF-DORC filtered Heatmap generatd using \code{\link[ComplexHeatmap]{Heatmap}})
#'@export
plotfigRHeatmap_EA <- function(figR.d,
                 score.cut=1,
                 DORCs=NULL,
                 TFs=NULL,
                 ... # Additional params passed to ComplexHeatmap
){
  
  
  message("Using absolute score cut-off of: ",score.cut," ..\n")
  
  DORCsToKeep <- figR.d %>% dplyr::filter(abs(Score) >= score.cut) %>% pull(DORC) %>% unique()
  TFsToKeep <- figR.d %>% dplyr::filter(abs(Score) >= score.cut) %>% pull(Motif) %>% unique()
  
  
  if(!is.null(DORCs)){
  if(!all(DORCs %in% figR.d$DORC))
    stop("One or more DORCs specified is not a valid DORC symbol found in the data.frame")
  DORCsToKeep <- intersect(DORCsToKeep,DORCs)
  TFsToKeep <- figR.d %>% dplyr::filter(abs(Score) >= score.cut & DORC %in% DORCsToKeep) %>% pull(Motif) %>% unique()
  }
  
  
  if(!is.null(TFs)){
  if(!all(TFs %in% figR.d$Motif))
    stop("One or more TFs specified is not a valid TF symbol found in the data.frame")
  TFsToKeep <- intersect(TFsToKeep,TFs)
  DORCsToKeep <- figR.d %>% dplyr::filter(abs(Score) >= score.cut & Motif %in% TFsToKeep) %>% pull(DORC) %>% unique()
  }
  
  
  net.d <- figR.d %>% dplyr::filter(DORC %in% DORCsToKeep & Motif %in% TFsToKeep) %>%
  reshape2::dcast(DORC ~ Motif) %>%
  tibble::column_to_rownames("DORC") %>% as.matrix()
  
  message("Plotting ",nrow(net.d)," DORCs x ",ncol(net.d), "TFs\n")
  
  # Heatmap view
  
  myCols <- circlize::colorRamp2(seq(-2,2,length.out = 9),colors = BuenColors::jdb_palette("solar_flare"))
  myHeat <- ComplexHeatmap::Heatmap(net.d,
                  col=myCols,
                  clustering_distance_rows = "pearson",
                  clustering_distance_columns = "pearson",
                  name="Score",
                  border = TRUE,
                  ...)
  
  myHeat
  return(myHeat)
  
}


#' Compare DORC Scores Between Two Clusters in a Pairwise Manner
#'
#' This function performs differential activity analysis between two specified clusters
#' in a Seurat object and generates visualization plots for top differentially expressed
#' genes (DORCs). It identifies genes with significantly higher or lower activity in
#' one cluster compared to another and creates combined DORC and RNA expression plots.
#'
#' @param seurat A Seurat object containing DORCs assay with chromatin accessibility data
#' @param dorcMat.s Matrix containing DORC (Domain of Regulatory Chromatin) activity scores
#' @param umap.d Data frame or matrix containing UMAP coordinates for dimensionality reduction
#' @param cluster_id Character or numeric. The primary cluster identifier (ident.1) for comparison
#' @param other_cluster_id Character or numeric. The comparison cluster identifier (ident.2)
#' @param figures_path Character. Directory path where output plots will be saved. 
#'   Should exist before function call
#' @param identities Character. Column name in metadata to use for cell identities. 
#'   Default is 'seurat_clusters_atac'
#' @param topn Numeric. Number of top differentially active genes to visualize for both
#'   up and downregulated genes. Default is 9
#' @param ... Additional arguments passed to volcano_plot function
#'
#' @return A data frame containing differential activity results with the following columns:
#'   \itemize{
#'     \item Gene names as row names
#'     \item p_val: Raw p-values from differential testing
#'     \item avg_diff: Average difference in activity between clusters
#'     \item pct.1: Percentage of cells in cluster_id with detected activity
#'     \item pct.2: Percentage of cells in other_cluster_id with detected activity
#'     \item p_val_adj: Bonferroni-corrected p-values
#'   }
#'
#' 
#' The function performs the following workflow:
#' \enumerate{
#'   \item Sets the default assay to 'DORCs' and identities to specified column
#'   \item Runs FindMarkers using Wilcoxon rank sum test on the DORCs assay
#'   \item Identifies top upregulated genes (avg_diff > 0, padj < 0.05)
#'   \item Identifies top downregulated genes (avg_diff < 0, padj < 0.05)
#'   \item Generates a volcano plot of differential activity (saved to figures_path)
#'   \item Creates side-by-side DORC and RNA expression plots for top upregulated genes
#'   \item Creates side-by-side DORC and RNA expression plots for top downregulated genes
#'   \item Saves all individual gene plots as PNG files to figures_path
#' }
#'
#' # View top differentially active genes
#' head(results)
#' }
#'
#' @seealso 
#' \code{\link[Seurat]{FindMarkers}} for differential expression testing
#' 
#' @export
compare_DORC_scores_pairwise <- function(seurat, dorcMat.s, umap.d, cluster_id, other_cluster_id, figures_path, identities='seurat_clusters_atac', topn=9,  ...) {
  
  # Set the default assay to  DORCs  for motif activity analysis
  # and update cell identities to the specified metadata column
  DefaultAssay(seurat) <- 'DORCs'
  Idents(seurat) <- identities  
  
  # Perform differential activity analysis between the two specified clusters
  # using the DORCs assay to identify genes with significantly different activity
  differential_activity <- FindMarkers(
    object = seurat,
    ident.1 = cluster_id,
    ident.2 = other_cluster_id,
    assay = 'DORCs',
    slot = 'data',
    only.pos = FALSE,
    fc.name = 'avg_diff'
  )

  # Extract top upregulated genes (positive avg_diff) with significant adjusted p-values
  # Filter for positive fold change, significant p-value, sort by significance,
  # and select the top n genes
  top_motif_activities_up <- differential_activity |> 
    filter(avg_diff > 0) |> 
    filter(p_val_adj < 0.05) |>
    arrange(p_val_adj) |>
    slice_head(n = topn) |>
    rownames_to_column(var = 'gene') |>
    pull(gene) |> 
    unique()

  # Extract top downregulated genes (negative avg_diff) with significant adjusted p-values
  # Filter for negative fold change, significant p-value, sort by significance,
  # and select the top n genes
  top_motif_activities_down <- differential_activity |> 
    filter(avg_diff < 0) |> 
    filter(p_val_adj < 0.05) |>
    arrange(p_val_adj) |>
    slice_head(n = topn) |>
    rownames_to_column(var = 'gene') |>
    pull(gene) |> 
    unique()
  
  # Reformat differential activity results for volcano plot
  # Rename columns to match expected input format
  differential_activity_cluster <- differential_activity |> 
    rownames_to_column(var = 'genes') |>
    rename(log2FoldChange = avg_diff, pvalue = p_val, padj = p_val_adj)

  # Generate volcano plot showing differential activity between clusters
  plot <- volcano_plot(differential_activity_cluster, group1 = other_cluster_id, group2 = cluster_id, cluster = cluster_id, local_figures_path = figures_path, FC_threshold = 0.25, p_value_threshold = 0.05, test_type = 'Wilcox_ATAC', ...)

  #################### UP ####################

  # Store upregulated genes for plotting
  dorcGenes_up <- top_motif_activities_up

  # Filter out genes that start with numbers as they create invalid R variable names
  valid_genes <- dorcGenes_up[!grepl("^[0-9]", dorcGenes_up)]

  # Create side-by-side DORC and RNA expression plots for each upregulated gene
  # Each plot shows both DORC activity and RNA expression on UMAP embeddings
  dorc_plots <- map(valid_genes, function(gene) {
  # Create DORC activity plot with heat color palette
  dorc_plot <- plotMarker2D(
    umap.d, 
    dorcMat.s, 
    markers = gene,
    maxCutoff = "q0.99",
    colorPalette = "brewer_heat"
  ) + ggtitle(paste(gene, "DORC"))
  
  # Create RNA expression plot with purple color palette
  rna_plot <- plotMarker2D(
    umap.d,
    RNAmat.s,
    markers = gene,
    maxCutoff = "q0.99",
    colorPalette = "brewer_purple"
  ) + ggtitle(paste(gene, "RNA"))
  
  # Display combined plot and save to file
  print(dorc_plot + rna_plot)
  ggsave(filename = paste0(figures_path, "/", gene, "_DORC_RNA_plot_up.png"), plot = dorc_plot + rna_plot, width = 8, height = 5)
  })

  #################### DOWN ####################

  # Store downregulated genes for plotting
  dorcGenes_down <- top_motif_activities_down

  # BUG: This line should use dorcGenes_down instead of dorcGenes
  # Filter out genes that start with numbers as they create invalid R variable names
  valid_genes <- dorcGenes_down[!grepl("^[0-9]", dorcGenes_down)]

  # Create side-by-side DORC and RNA expression plots for each downregulated gene
  # Each plot shows both DORC activity and RNA expression on UMAP embeddings
  dorc_plots <- map(valid_genes, function(gene) {
  # Create DORC activity plot with heat color palette
  dorc_plot <- plotMarker2D(
    umap.d, 
    dorcMat.s, 
    markers = gene,
    maxCutoff = "q0.99",
    colorPalette = "brewer_heat"
  ) + ggtitle(paste(gene, "DORC"))
  
  # Create RNA expression plot with purple color palette
  rna_plot <- plotMarker2D(
    umap.d,
    RNAmat.s,
    markers = gene,
    maxCutoff = "q0.99",
    colorPalette = "brewer_purple"
  ) + ggtitle(paste(gene, "RNA"))
  
  # Display combined plot and save to file
  print(dorc_plot + rna_plot)
  ggsave(filename = paste0(figures_path, "/", gene, "_DORC_RNA_plot_down.png"), plot = dorc_plot + rna_plot, width = 8, height = 5)
  })

  # Return the complete differential activity results table
  return(differential_activity)
}