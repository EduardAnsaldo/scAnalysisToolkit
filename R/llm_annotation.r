## LLM Cluster Annotation Functions
## Context-aware cell type annotation of Seurat clusters using the Anthropic API
## (Claude), grounded with each cluster's top marker genes, the tissue /
## experimental-design context, and optionally retrieved PubMed abstracts.

# ----------------------------------------------------------------------------
# Internal helpers
# ----------------------------------------------------------------------------

# System prompt + JSON spec (kept verbatim from the benchmarked context-aware
# approach so behaviour matches).
.llm_system_prompt <- function() {
    paste(
        "You are an expert single-cell transcriptomics annotator. Given marker genes",
        "for one cell cluster, identify the most likely cell type. Respond with ONLY a",
        "JSON object, no prose, no markdown fences."
    )
}

.llm_json_spec <- function() {
    paste0(
        'Return exactly: {"cell_type": "<specific cell type>", ',
        '"confidence": <0.0-1.0>, "rationale": "<=25 words>"}'
    )
}

# Build the context-aware prompt for one cluster.
.llm_build_context_prompt <- function(markers, tissue, design, abstract_block) {
    paste0(
        "Tissue: ", tissue, "\n",
        "Experimental design: ", design, "\n\n",
        "Top marker genes (high to low): ", paste(markers, collapse = ", "), "\n\n",
        "Relevant published literature:\n", abstract_block, "\n\n",
        "Using the tissue context and the literature above, name the most specific ",
        "defensible cell type for this cluster.\n", .llm_json_spec()
    )
}

# Best-effort extraction of the first JSON object from a model response.
.llm_parse_json <- function(text) {
    empty <- list(cell_type = "unknown", confidence = NA_real_, rationale = "")
    if (is.null(text) || !nzchar(trimws(text))) return(empty)

    cleaned <- text |>
        trimws() |>
        stringr::str_remove_all("(?m)^```(?:json)?") |>
        stringr::str_remove_all("(?m)```$") |>
        trimws()

    parsed <- tryCatch(jsonlite::fromJSON(cleaned), error = function(e) NULL)
    if (is.null(parsed)) {
        m <- stringr::str_extract(cleaned, "(?s)\\{.*\\}")
        if (!is.na(m)) parsed <- tryCatch(jsonlite::fromJSON(m), error = function(e) NULL)
    }
    if (is.null(parsed) || is.null(parsed$cell_type)) return(empty)

    list(
        cell_type  = as.character(parsed$cell_type %||% "unknown"),
        confidence = suppressWarnings(as.numeric(parsed$confidence %||% NA)),
        rationale  = as.character(parsed$rationale %||% "")
    )
}

# Retrieve up to `n` PubMed abstracts relevant to a tissue + marker set, with an
# on-disk cache. Best-effort: network/parse failures yield no abstracts.
.llm_pubmed_abstracts <- function(tissue, markers, cache_dir, n = 3,
                                  max_chars = 1200, email = NULL) {
    esearch <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    efetch  <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    top <- utils::head(markers, 6)
    query <- paste0(
        "(", tissue, ") AND (", paste(top, collapse = " OR "), ") ",
        "AND (marker OR cell type OR single cell)"
    )

    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    key <- substr(rlang::hash(paste0(query, "|", n)), 1, 16)
    cache_file <- file.path(cache_dir, paste0(key, ".json"))
    if (file.exists(cache_file)) {
        return(tibble::as_tibble(jsonlite::fromJSON(cache_file)))
    }

    common <- list(db = "pubmed")
    if (!is.null(email)) common$email <- email

    out <- tryCatch(
        {
            Sys.sleep(0.4)  # respect NCBI rate limit (~3 req/s without a key)
            ids_json <- httr2::request(esearch) |>
                httr2::req_url_query(
                    !!!common, term = query, retmax = n,
                    sort = "relevance", retmode = "json"
                ) |>
                httr2::req_timeout(30) |>
                httr2::req_perform() |>
                httr2::resp_body_json()
            ids <- ids_json$esearchresult$idlist %||% list()
            if (length(ids) == 0) {
                tibble::tibble(pmid = character(), text = character())
            } else {
                Sys.sleep(0.4)
                raw <- httr2::request(efetch) |>
                    httr2::req_url_query(
                        !!!common, id = paste(unlist(ids), collapse = ","),
                        rettype = "abstract", retmode = "text"
                    ) |>
                    httr2::req_timeout(30) |>
                    httr2::req_perform() |>
                    httr2::resp_body_string()
                .llm_parse_efetch_text(raw, unlist(ids), max_chars)
            }
        },
        error = function(e) tibble::tibble(pmid = character(), text = character())
    )

    jsonlite::write_json(out, cache_file, auto_unbox = TRUE)
    out
}

# efetch text mode separates records with blank lines; pair each with its PMID.
.llm_parse_efetch_text <- function(raw, ids, max_chars) {
    chunks <- raw |> stringr::str_split("\\n\\n\\n+") |> unlist() |> trimws()
    chunks <- chunks[nzchar(chunks)]
    if (length(chunks) == 0) {
        return(tibble::tibble(pmid = character(), text = character()))
    }
    chunks <- utils::head(chunks, length(ids))
    tibble::tibble(
        pmid = as.character(ids[seq_along(chunks)]),
        text = substr(stringr::str_squish(chunks), 1, max_chars)
    )
}

.llm_format_abstracts <- function(abstracts) {
    if (is.null(abstracts) || nrow(abstracts) == 0) {
        return("(no relevant abstracts retrieved)")
    }
    paste(paste0("[PMID ", abstracts$pmid, "] ", abstracts$text), collapse = "\n\n")
}

# Single-turn Anthropic Messages API call. The shared system prompt is sent with
# cache_control (prompt caching) to amortise cost across the per-cluster calls.
.llm_anthropic_message <- function(prompt, system, model, api_key,
                                   max_tokens = 1024, temperature = 0,
                                   anthropic_version = "2023-06-01") {
    body <- list(
        model = model,
        max_tokens = max_tokens,
        temperature = temperature,
        system = list(list(
            type = "text", text = system,
            cache_control = list(type = "ephemeral")
        )),
        messages = list(list(role = "user", content = prompt))
    )

    t0 <- Sys.time()
    resp <- tryCatch(
        httr2::request("https://api.anthropic.com/v1/messages") |>
            httr2::req_headers(
                "x-api-key" = api_key,
                "anthropic-version" = anthropic_version,
                "content-type" = "application/json"
            ) |>
            httr2::req_body_json(body) |>
            httr2::req_retry(
                max_tries = 5,
                is_transient = function(r) {
                    httr2::resp_status(r) %in% c(408, 429, 500, 502, 503, 529)
                }
            ) |>
            httr2::req_timeout(300) |>
            httr2::req_perform(),
        error = function(e) e
    )
    latency <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    if (inherits(resp, "error")) {
        return(list(text = "", input_tokens = 0L, output_tokens = 0L,
                    latency_s = latency, error = conditionMessage(resp)))
    }
    payload <- httr2::resp_body_json(resp)
    text <- vapply(
        payload$content,
        function(b) if (!is.null(b$text)) b$text else "",
        character(1)
    ) |> paste(collapse = "") |> trimws()
    usage <- payload$usage %||% list()
    list(
        text = text,
        input_tokens  = as.integer(usage$input_tokens  %||% 0),
        output_tokens = as.integer(usage$output_tokens %||% 0),
        latency_s = latency,
        error = NULL
    )
}

# Resolve the Anthropic API key from (in priority order) an explicit argument,
# the system keyring (macOS Keychain etc.), or the ANTHROPIC_API_KEY env var.
# Returns "" when none is available.
.llm_resolve_api_key <- function(api_key = NULL, keyring_service = "anthropic") {
    if (!is.null(api_key) && nzchar(api_key)) return(api_key)
    if (rlang::is_installed("keyring")) {
        key <- tryCatch(keyring::key_get(keyring_service), error = function(e) NULL)
        if (!is.null(key) && nzchar(key)) return(key)
    }
    env <- Sys.getenv("ANTHROPIC_API_KEY")
    if (nzchar(env)) return(env)
    ""
}

# Top-N markers per cluster as a tibble with a list-column of genes.
.llm_cluster_markers <- function(seurat, cluster_col, n, assay, markers) {
    meta <- seurat@meta.data
    cluster_ids <- as.character(meta[[cluster_col]])
    n_cells <- tibble::tibble(cluster = cluster_ids) |>
        count(cluster, name = "n_cells")

    if (is.null(markers)) {
        Idents(seurat) <- factor(cluster_ids)
        if (is.null(assay)) assay <- DefaultAssay(seurat)
        markers <- FindAllMarkers(
            seurat, only.pos = TRUE, min.pct = 0.1,
            logfc.threshold = 0.25, assay = assay, verbose = FALSE
        )
    }
    if (!all(c("cluster", "gene", "avg_log2FC") %in% colnames(markers))) {
        stop("`markers` must contain columns 'cluster', 'gene', and 'avg_log2FC'.")
    }

    arrange_expr <- if ("p_val_adj" %in% colnames(markers)) {
        function(df) arrange(df, p_val_adj, desc(avg_log2FC), .by_group = TRUE)
    } else {
        function(df) arrange(df, desc(avg_log2FC), .by_group = TRUE)
    }

    top <- markers |>
        tibble::as_tibble() |>
        mutate(cluster = as.character(cluster)) |>
        group_by(cluster) |>
        arrange_expr() |>
        slice_head(n = n) |>
        summarise(markers = list(gene), .groups = "drop")

    n_cells |>
        left_join(top, by = "cluster") |>
        mutate(markers = lapply(markers, function(m) if (is.null(m)) character() else m)) |>
        arrange(cluster)
}

# ----------------------------------------------------------------------------
# Exported function
# ----------------------------------------------------------------------------

#' Annotate Seurat Clusters with an LLM (Claude)
#'
#' Context-aware, marker-based cell type annotation of a clustered Seurat object
#' using the Anthropic API. For each cluster the top marker genes (computed with
#' [Seurat::FindAllMarkers()] or supplied) are sent to Claude together with the
#' tissue and experimental-design context and, optionally, abstracts retrieved
#' from PubMed, and the model returns the most specific defensible cell type. The
#' per-cluster label is written back into the object metadata (broadcast to every
#' cell in the cluster).
#'
#' Requires an Anthropic API key and the `httr2` and `jsonlite` packages. The key
#' is resolved (in priority order) from the `api_key` argument, the system keyring
#' (recommended — store it once with `keyring::key_set("anthropic")`), then the
#' `ANTHROPIC_API_KEY` environment variable. Use `dry_run = TRUE` to assemble and
#' inspect the prompts without calling the API (no key needed).
#'
#' @param seurat Seurat object with cluster assignments.
#' @param tissue Character; tissue / sample description used as context, e.g.
#'   "human peripheral blood mononuclear cells (PBMC)".
#' @param design Character; experimental-design context, e.g.
#'   "healthy donor, 10x Genomics 3' scRNA-seq". Default is a generic placeholder.
#' @param cluster_col Character; metadata column holding cluster ids.
#'   Default 'seurat_clusters'.
#' @param label_col Character; metadata column to write predicted labels into.
#'   A companion `<label_col>_confidence` column is also added. Default 'llm_celltype'.
#' @param model Character; Anthropic model id. Default 'claude-opus-4-8'.
#' @param markers Optional data frame of precomputed markers (a
#'   [Seurat::FindAllMarkers()] result, or any data frame with columns 'cluster',
#'   'gene', and 'avg_log2FC'); skips marker computation. Default NULL.
#' @param n_markers Integer; number of top markers per cluster sent to the model.
#'   Default 25.
#' @param assay Character; assay used for marker computation. Default NULL
#'   (uses [Seurat::DefaultAssay()]).
#' @param use_pubmed Logical; retrieve PubMed abstracts for grounding. Default TRUE.
#' @param pubmed_n Integer; abstracts to retrieve per cluster. Default 10.
#' @param pubmed_email Character; contact email for NCBI E-utilities etiquette.
#'   Default NULL.
#' @param max_tokens Integer; max output tokens per call. Default 1024.
#' @param temperature Numeric; sampling temperature. Default 0.
#' @param api_key Character; Anthropic API key. Default NULL, in which case the
#'   key is taken from the system keyring (`keyring_service`) or the
#'   `ANTHROPIC_API_KEY` environment variable.
#' @param keyring_service Character; keyring service name under which the API key
#'   is stored (see `keyring::key_set()`). Default 'anthropic'.
#' @param cache_dir Character; directory for the PubMed abstract cache.
#'   Default a session temp directory.
#' @param dry_run Logical; if TRUE, assemble prompts only (no API call) and store
#'   them in `seurat@misc$llm_annotation_prompts`. Default FALSE.
#' @param verbose Logical; print progress messages. Default TRUE.
#'
#' @return The Seurat object with added metadata columns `label_col` and
#'   `<label_col>_confidence`, and the per-cluster results tibble stored in
#'   `seurat@misc$llm_annotation` (columns: cluster, n_cells, top_markers,
#'   cell_type, confidence, rationale, n_refs, input_tokens, output_tokens, raw,
#'   error). With `dry_run = TRUE`, metadata is unchanged and the assembled
#'   prompts are stored in `seurat@misc$llm_annotation_prompts` instead.
#'
#' @examples
#' \dontrun{
#' # One-time: store the key securely in the OS keyring (prompts for the value).
#' keyring::key_set("anthropic")
#'
#' seurat <- annotate_seurat_with_LLM(
#'     seurat,
#'     tissue = "human peripheral blood mononuclear cells (PBMC)",
#'     design = "healthy donor, 10x Genomics 3' scRNA-seq"
#' )
#' seurat@misc$llm_annotation        # per-cluster table
#'
#' # Preview prompts without calling the API:
#' seurat <- annotate_seurat_with_LLM(seurat, tissue = "human PBMC", dry_run = TRUE)
#' seurat@misc$llm_annotation_prompts
#' }
#'
#' @importFrom utils head
#' @export
annotate_seurat_with_LLM <- function(
    seurat,
    tissue,
    design = "Single-cell RNA-seq; experimental design unspecified.",
    cluster_col = "seurat_clusters",
    label_col = "llm_celltype",
    model = "claude-opus-4-8",
    markers = NULL,
    n_markers = 25,
    assay = NULL,
    use_pubmed = TRUE,
    pubmed_n = 10,
    pubmed_email = NULL,
    max_tokens = 1024,
    temperature = 0,
    api_key = NULL,
    keyring_service = "anthropic",
    cache_dir = file.path(tempdir(), "llm_pubmed_cache"),
    dry_run = FALSE,
    verbose = TRUE
) {
    if (missing(tissue) || !nzchar(tissue)) {
        stop("`tissue` is required (it provides the annotation context).")
    }
    if (!cluster_col %in% colnames(seurat@meta.data)) {
        stop(sprintf(
            "Cluster column '%s' not found in metadata. Available columns: %s",
            cluster_col, paste(colnames(seurat@meta.data), collapse = ", ")
        ))
    }
    # Dependencies needed for any network / parsing work.
    if (use_pubmed || !dry_run) {
        rlang::check_installed(c("httr2", "jsonlite"),
                               reason = "to call the Anthropic API and query PubMed.")
    }
    if (!dry_run) {
        api_key <- .llm_resolve_api_key(api_key, keyring_service)
        if (!nzchar(api_key)) {
            stop("No Anthropic API key found. Store it once with ",
                 "`keyring::key_set(\"", keyring_service, "\")`, or set the ",
                 "ANTHROPIC_API_KEY environment variable, or pass `api_key=`. ",
                 "Use `dry_run = TRUE` to assemble prompts without a key.")
        }
    }

    mk <- .llm_cluster_markers(seurat, cluster_col, n_markers, assay, markers)
    if (verbose) message(sprintf("Prepared markers for %d cluster(s).", nrow(mk)))

    system_prompt <- .llm_system_prompt()
    rows <- vector("list", nrow(mk))

    for (i in seq_len(nrow(mk))) {
        cl <- mk$cluster[i]
        genes <- mk$markers[[i]]
        abstracts <- NULL
        if (use_pubmed && length(genes) > 0) {
            abstracts <- .llm_pubmed_abstracts(
                tissue, genes, cache_dir, n = pubmed_n, email = pubmed_email
            )
        }
        n_refs <- if (is.null(abstracts)) 0L else nrow(abstracts)
        prompt <- .llm_build_context_prompt(
            genes, tissue, design, .llm_format_abstracts(abstracts)
        )

        if (dry_run) {
            rows[[i]] <- tibble::tibble(
                cluster = cl, n_cells = mk$n_cells[i],
                top_markers = paste(genes, collapse = ", "),
                n_refs = n_refs, prompt = prompt
            )
            next
        }

        resp <- .llm_anthropic_message(
            prompt, system_prompt, model = model, api_key = api_key,
            max_tokens = max_tokens, temperature = temperature
        )
        parsed <- .llm_parse_json(resp$text)
        if (verbose) message(sprintf("  cluster %s -> %s", cl, parsed$cell_type))
        rows[[i]] <- tibble::tibble(
            cluster = cl, n_cells = mk$n_cells[i],
            top_markers = paste(genes, collapse = ", "),
            cell_type = parsed$cell_type, confidence = parsed$confidence,
            rationale = parsed$rationale, n_refs = n_refs,
            input_tokens = resp$input_tokens, output_tokens = resp$output_tokens,
            raw = resp$text, error = resp$error %||% NA_character_
        )
    }

    results <- bind_rows(rows)

    if (dry_run) {
        seurat@misc$llm_annotation_prompts <- results
        if (verbose) {
            message(sprintf(
                "DRY RUN: assembled %d prompt(s) in seurat@misc$llm_annotation_prompts (no API call).",
                nrow(results)
            ))
        }
        return(seurat)
    }

    # Broadcast each cluster's label to its cells.
    idx <- match(as.character(seurat@meta.data[[cluster_col]]), results$cluster)
    seurat[[label_col]] <- results$cell_type[idx]
    seurat[[paste0(label_col, "_confidence")]] <- results$confidence[idx]
    seurat@misc$llm_annotation <- results

    if (verbose) {
        message(sprintf(
            "Done. Wrote '%s' to metadata. Tokens: %d in / %d out.",
            label_col, sum(results$input_tokens, na.rm = TRUE),
            sum(results$output_tokens, na.rm = TRUE)
        ))
    }
    return(seurat)
}
