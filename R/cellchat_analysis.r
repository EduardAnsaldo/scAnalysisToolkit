## CellChat Analysis Functions
## Functions for cell-cell communication analysis

#' Run CellChat Analysis on Seurat Object
#'
#' @param seurat A Seurat object with cell type annotations
#' @param group_by Character string specifying the metadata column for cell grouping
#' @param species Character string, either "mouse" or "human"
#' @param smooth_data Logical, whether to smooth data using PPI network
#' @param subset_cell_groups Character vector of cell groups to subset, or NULL for all
#' @param min_cells Integer, minimum number of cells required for communication
#' @param use_parallel Logical, whether to use parallel processing
#' @param n_workers Integer, number of parallel workers (only used if use_parallel = TRUE)
#' @param max_memory_mb Numeric, maximum memory for future in MB (only used if use_parallel = TRUE)
#' @param population_size Logical, whether to consider cell population size effect
#' @param nboot Integer, number of bootstrap samples for computing p-values
#' @param output_file Character string, path to save RDS file, or NULL to skip saving
#'
#' @return A CellChat object
#' @export
run_cellchat_analysis <- function(seurat,
           group_by = "major_cell_types",
           species = "mouse",
           smooth_data = FALSE,
           subset_cell_groups = NULL,
           min_cells = 10,
           use_parallel = TRUE,
           n_workers = 10,
           max_memory_mb = 4000,
           population_size = FALSE,
           nboot = 100,
           output_file = NULL) {

  # Check for required packages
  if (!requireNamespace("CellChat", quietly = TRUE)) {
    stop("Package 'CellChat' is required but not installed. Please install it with:\n",
         "  devtools::install_github('jinworks/CellChat')",
         call. = FALSE)
  }

  # Part I: Setup and preprocessing
  message("Starting CellChat analysis...")

  ## LogNormalize data
  message("Normalizing data...")
  seurat <- NormalizeData(seurat, assay = "RNA", normalization.method = "LogNormalize")

  ## Extract expression data and metadata from the Seurat object
  if ("Samples" %in% colnames(seurat@meta.data)) {
  seurat@meta.data <- seurat@meta.data |> mutate(samples = Samples)
  }

  message("Creating CellChat object...")
  cellchat_object <- CellChat::createCellChat(object = seurat, group.by = group_by, assay = "RNA")

  ## Set the CellChat database
  cellchat_object@DB <- if (species == "mouse") {
  CellChatDB.mouse
  } else if (species == "human") {
  CellChatDB.human
  } else {
  stop("Species must be either 'mouse' or 'human'")
  }

  ## Preprocessing the data for CellChat analysis
  message("Preprocessing data...")
  cellchat_object <- subsetData(cellchat_object)

  # Set up parallel processing if requested
  if (use_parallel) {
  message(sprintf("Setting up parallel processing with %d workers...", n_workers))
  future::plan("multisession", workers = n_workers)
  options(future.globals.maxSize = max_memory_mb * 1024^2)
  } else {
  message("Running in sequential mode...")
  future::plan("sequential")
  }

  cellchat_object <- CellChat::identifyOverExpressedGenes(cellchat_object)

  message("Identifying over-expressed interactions...")
  ptm <- Sys.time()
  cellchat_object <- CellChat::identifyOverExpressedInteractions(cellchat_object)
  execution_time <- Sys.time() - ptm
  message(sprintf("Execution time: %.2f mins", as.numeric(execution_time, units = "mins")))

  ## Optional: smooth data using PPI network
  if (smooth_data) {
  message("Smoothing data using PPI network...")
  ppi_network <- if (species == "mouse") PPI.mouse else PPI.human
  cellchat_object <- smoothData(cellchat_object, adj = ppi_network)
  }

  # Part II: Compute communication probabilities
  message(sprintf("Computing communication probabilities (nboot = %d)...", nboot))
  ptm <- Sys.time()
  cellchat_object <- CellChat::computeCommunProb(
  cellchat_object,
  type = "triMean",
  raw.use = !smooth_data,
  population.size = population_size,
  seed.use = 3514L,
  nboot = nboot
  )
  execution_time <- Sys.time() - ptm
  message(sprintf("Execution time: %.2f mins", as.numeric(execution_time, units = "mins")))

  ## Filter out cell-cell communication
  message(sprintf("Filtering communication (min.cells = %d)...", min_cells))
  cellchat_object <- CellChat::filterCommunication(cellchat_object, min.cells = min_cells)

  ## Infer signaling pathway level communication

  message("Computing communication probabilities at pathway level...")
  ptm <- Sys.time()
  cellchat_object <- CellChat::computeCommunProbPathway(cellchat_object)
  execution_time <- Sys.time() - ptm
  message(sprintf("Execution time: %.2f mins", as.numeric(execution_time, units = "mins")))

  ## Calculate aggregated network
  message("Aggregating cell-cell communication network...")
  if (!is.null(subset_cell_groups)) {
  cellchat_object <- CellChat::aggregateNet(
  cellchat_object,
  sources.use = subset_cell_groups,
  targets.use = subset_cell_groups
  )
  } else {
  cellchat_object <- CellChat::aggregateNet(cellchat_object)
  }

  # Save if output file specified
  if (!is.null(output_file)) {
  message(sprintf("Saving CellChat object to: %s", output_file))
  saveRDS(cellchat_object, file = output_file)
  }

  future::plan("sequential")  # Reset to sequential plan
  message("CellChat analysis complete!")
  return(invisible(NULL))
}

# Example usage:
# seurat <- readRDS(here::here('data/SPF_GF_annotated.rds'))
# cellchat_object <- run_cellchat_analysis(
#   seurat,
#   group_by = "major_cell_types",
#   species = "mouse",
#   use_parallel = TRUE,
#   output_file = here::here('data/CellChat_SPF_GF_major_cell_types.rds')
# )
