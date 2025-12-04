## Pathway Analysis Dispatcher Functions
## Dispatcher wrappers for pathway analyses that support multiple methods

#' Run Differential Expression Gene Functional Analysis
#'
#' Dispatcher function that performs functional pathway analysis on differential
#' expression results using one or more methods. Supports ClusterProfiler (GO/KEGG),
#' Metascape web-based analysis, and gProfiler2. Automatically calls the appropriate
#' analysis function based on method selection.
#'
#' @param results Data frame; differential expression results containing gene names
#'   and statistics (typically from DE analysis functions)
#' @param method Character vector; analysis method(s) to run. Options: 'ClusterProfiler',
#'   'Metascape', 'gProfiler2'. Can specify up to 3 methods. If NULL or empty,
#'   returns invisible(NULL). Default NULL
#' @param ... Additional arguments passed to the specific analysis functions
#'   (e.g., grouping_var, path, FC_threshold, pathways_of_interest)
#'
#' @return If single method: results from that method. If multiple methods: named list
#'   with results from each method. If method is NULL: invisible(NULL)
#'
#' @export
run_DEG_functional_analysis <- function(results,
                                        method = NULL,
                                        ...) {
  if (is.null(method) || length(method) == 0) return(invisible(NULL))

  method <- match.arg(method, choices = c("ClusterProfiler", "Metascape", "gProfiler2"), several.ok = TRUE)
  if (length(method) > 3) stop("up to 3 methods allowed")

  res <- lapply(method, function(m) {
    switch(m,
           ClusterProfiler = {
             if (exists("GO_functional_analysis", mode = "function")) {
               GO_functional_analysis(results = results, ...)
             } else stop("GO_functional_analysis not found")
           },
           Metascape = {
             if (exists("Metascape_functional_analysis", mode = "function")) {
               Metascape_functional_analysis(results = results, ...)
             } else stop("Metascape_functional_analysis not found")
           },
           gProfiler2 = {
             if (exists("gProfiler2_functional_analysis", mode = "function")) {
               gProfiler2_functional_analysis(results = results, ...)
             } else stop("gProfiler2_functional_analysis not found")
           },
           stop("Unknown DEG functional method: ", m))
  })
  names(res) <- method
  if (length(res) == 1) res[[1]] else res
}


#' Run Gene List Overrepresentation Analysis
#'
#' Dispatcher function that performs overrepresentation analysis (enrichment) on a
#' gene list using one or more methods. Supports ClusterProfiler (GO/KEGG), gProfiler2,
#' and Metascape web-based analysis. Tests which biological pathways/processes are
#' statistically enriched in the gene list compared to a background.
#'
#' @param genes Character vector; gene symbols to test for enrichment
#' @param method Character vector; analysis method(s) to run. Options: 'ClusterProfiler',
#'   'gProfiler2', 'Metascape'. Can specify up to 3 methods. If NULL or empty,
#'   returns invisible(NULL). Default NULL
#' @param ... Additional arguments passed to the specific analysis functions
#'   (e.g., all_genes, local_path, ontology, organism)
#'
#' @return If single method: results from that method. If multiple methods: named list
#'   with results from each method. If method is NULL: invisible(NULL)
#'
#' @export
run_overrepresentation_analysis <- function(genes,
                                           method = NULL,
                                           ...) {
  if (is.null(method) || length(method) == 0) return(invisible(NULL))

  method <- match.arg(method, choices = c("ClusterProfiler", "gProfiler2", "Metascape"), several.ok = TRUE)
  if (length(method) > 3) stop("up to 3 methods allowed")

  res <- lapply(method, function(m) {
    switch(m,
           ClusterProfiler = {
             if (exists("GO_overrepresentation_analysis", mode = "function")) {
               GO_overrepresentation_analysis(genes, ...)
             } else stop("GO_overrepresentation_analysis_overrepresentation_analysis not found")
           },
           gProfiler2 = {
             if (exists("gProfiler2_overrepresentation_analysis", mode = "function")) {
               gProfiler2_overrepresentation_analysis(genes, ...)
             } else stop("gProfiler2_overrepresentation_analysis not found")
           },
           Metascape = {
             if (exists("Metascape_overrepresentation_analysis", mode = "function")) {
               Metascape_overrepresentation_analysis(genes, ...)
             } else stop("Metascape_overrepresentation_analysis not found")
           },
           stop("Unknown overrepresentation method: ", m))
  })
  names(res) <- method
  if (length(res) == 1) res[[1]] else res
}
