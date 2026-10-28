# Clustering and annotation

**Shared contract for this file.** `perform_seurat_clustering()` and `process_multimodal_seurat()`
return **lists, not Seurat objects** — take `$seurat`, and check `$stage` first. `extract_cell_counts()`
is one of only two functions in the package taking a **bare** column name. None of these functions
create their output directory.

---

## `perform_seurat_clustering()`

```r
perform_seurat_clustering(seurat, npcs = 100, dimensions = 30, k_param = 20,
    n_neighbors = 30, resolutions = c(0.05, 0.1, 0.15, 0.25, 0.375, 0.5),
    filter_variable = NULL, keep_values = NULL, plot_elbow = TRUE, seed = 3514L,
    verbose = FALSE, selected_resolution = NULL, object_annotations = 'full_object',
    figures_path = NULL, quiet_vdj = FALSE, integrate_rna = FALSE,
    integration_column = 'Origin', vars_to_regress = NULL, regress_cell_cycle = FALSE,
    cell_cycle_species = 'mouse', plot_cell_cycle_umap = TRUE,
    algorithm = 4, cluster_random_seed = seed)
```

SCTransform → PCA → UMAP → `FindNeighbors` → `FindClusters` at every resolution. Use this rather than
hand-rolling normalisation.

**Returns a 7-element list:** `seurat`, `resolution_plot` (a `gtable`), `detailed_plots` (`NULL`, or
a list of `cluster_plot`, `group_split_plot`, `groups_plot`), `dimensions`, `resolutions`,
`selected_resolution`, `cell_cycle_plot`.

⚠️ **`seurat <- perform_seurat_clustering(seurat)` silently replaces your object with a list.**

```r
res <- perform_seurat_clustering(seurat, dimensions = 30, n_neighbors = 30,
                                 selected_resolution = 0.25, figures_path = figures_path)
seurat <- res$seurat
```

**Staged early return.** With `regress_cell_cycle = TRUE` **and** `dimensions = NULL` it returns only
`list(seurat, stage = "select_dimensions")` after printing the regressed elbow, so you can pick
dimensions and call again. Note `dimensions = NULL` only short-circuits in that combination — with
`regress_cell_cycle = FALSE` it fails at `dims = 1:NULL`.

**Intended two-pass use.** Call once to see the elbow and the per-resolution UMAP grid, then set
`dimensions`, `n_neighbors` and `selected_resolution` and call again.

**Writes** — only when `figures_path` is non-`NULL` (the default `NULL` skips all writes), via
`ggsave(path=)` so relative to `getwd()`, into a directory that must already exist:

- `UMAP_cell_cycle_before_after_<object_annotations>.pdf` (12×10) — only with `regress_cell_cycle`
- `UMAP_clusters_R_<selected_resolution>_<object_annotations>.pdf`
- `UMAP_clusters_by_group_R_<selected_resolution>_<object_annotations>.pdf` (5×5)

**Preconditions.** `RNA` assay with counts. `Groups` is required **only when `selected_resolution` is
non-`NULL`** (used as `split.by` / `group.by`). `integration_column` (default `'Origin'`) only when
`integrate_rna = TRUE`. `algorithm = 4` is Leiden and needs `leidenbase` — it is in `Imports` but never
referenced by name, so a broken install shows up only here.

**Metadata written.** `SCT_snn_res.<resolution>` for every resolution. `seurat_clusters` is set **only**
when `selected_resolution` is supplied.

**Gotchas.**
- The roxygen says `filter_variable` defaults to `FALSE`; the real default is `NULL` and the guard is
  `!is.null()`. Passing `FALSE` enters the branch and stops with *"keep_values must be provided"*.
- `quiet_vdj = TRUE` drops variable features matching the **unanchored** pattern
  `"Igh|Igl|Igk|Tra|Trb|Trd|Trg"`, so it also removes genes like `Trap1` and `Trak2`.
- `gridExtra::grid.arrange()` draws immediately as a side effect as well as returning the gtable, so
  the resolution grid appears even if you discard it.

## `process_multimodal_seurat()`

```r
process_multimodal_seurat(seurat, figures_path = NULL, object_annotations = NULL,
    dimensions_RNA = NULL, dimensions_ATAC = NULL,
    resolution_RNA = NULL, resolution_ATAC = NULL, resolution_WNN = NULL,
    resolutions_RNA  = c(0.375, 0.4, 0.475, 0.5, 0.55, 0.6, 0.625, 0.75, 1),
    resolutions_ATAC = c(0.1, 0.25, 0.3, 0.325, 0.375, 0.5, 0.625, 0.75, 1),
    resolutions_WNN  = c(0.1, 0.25, 0.375, 0.5, 0.625, 0.75, 1, 1.25, 1.5),
    integrate_rna = FALSE, integration_column = 'Origin', npcs = 50, ndims_lsi = 40,
    k_param = 20, n_neighbors = 30, random_seed = 3514, remove_ig = FALSE,
    genome = 'GRCm39', generate_plots = TRUE, skip_preprocessing = FALSE,
    skip_umap_and_neighbors = FALSE)
```

RNA + ATAC + WNN, as a **three-stage interactive pipeline** with two deliberate early returns:

| returned | when |
|---|---|
| `list(seurat, stage = "dimensions")` — 2 elements, no `missing` | `dimensions_RNA` or `dimensions_ATAC` is `NULL`; prints both elbows and `DepthCor` |
| `list(seurat, stage = "resolutions", missing = <chr>)` | any of the three `resolution_*` is `NULL`; prints the three resolution grids |
| `list(seurat, stage = "complete")` | everything supplied |

Intended sequence: call bare → inspect elbows → recall with dimensions and `skip_preprocessing = TRUE`
→ inspect grids → recall with resolutions, `skip_preprocessing = TRUE` and
`skip_umap_and_neighbors = TRUE`.

**Hard-coded names.** Assays `RNA` and `ATAC` (SCTransform creates `SCT`); graphs `SCT_snn`,
`ATAC_snn`, `wsnn`; reductions `pca` (or `harmony`), `lsi`, `umap.rna`, `umap.atac`, `umap.wnn`.
Writes `seurat_clusters_rna` / `_atac` / `_wnn` and deletes the per-resolution columns once a final
resolution is chosen. ATAC always uses dims `2:dimensions_ATAC` (dropping LSI 1).

**Requires** `Groups` (set as `Idents` before RNA clustering) and, with `generate_plots = TRUE`, also
`Samples`.

**Mouse only.** `genome` must be `'GRCm38'` or `'GRCm39'`, used only when `remove_ig = TRUE`, with
hard-coded Ig coordinates and a mouse-symbol `"Igh|Igl|Igk"` regex.

⚠️ **`figures_path` must end in a slash** — filenames are built with plain
`paste0(figures_path, object_annotations, ...)` with no separator. Nine PDFs are written when
`generate_plots = TRUE` and both `figures_path` and `object_annotations` are non-`NULL`.

⚠️ Harmony integration calls `IntegrateLayers(method = HarmonyIntegration)`, which needs `harmony` —
**not declared in `Suggests` and not guarded**.

⚠️ `skip_umap_and_neighbors = TRUE` still re-runs `FindClusters` over the whole resolution vectors and
re-prints the grids, and requires the three graphs to already exist.

---

## `extract_cell_counts()`

```r
extract_cell_counts(seurat, grouping_var, figures_path, tables_path,
                    object_annotations = '', plot_grouping = 'Samples', title = NULL)
```

⚠️ **`grouping_var` is embraced — pass a bare, unquoted column name.** `plot_grouping` is a plain
string. This is one of only two functions in the package with that convention.

```r
extract_cell_counts(seurat, grouping_var = seurat_clusters,
                    figures_path = figures_path, tables_path = tables_path,
                    plot_grouping = 'Groups')
```

**Requires `Samples` *and* `Groups`** in all three modes — the first line is an unconditional
`FetchData(vars = c(<grouping_var>, "Samples", "Groups"))`.

`plot_grouping` must be `'Samples'`, `'Groups'` or `'grouping_var'`; anything else stops. Call it once
per grouping you want. Note the modes differ in more than layout: `'Samples'` plots the aggregated
per-sample tibble (mean ± SE), while the other two plot per-cell stacked proportions from
`meta.data` directly.

The per-sample table is zero-filled (since 0.12.1): a sample with no cells of a level gets a 0 row, so
the `'Samples'` mean ± SE averages over every sample in the group and the CSVs show 0, not `NA`. Only
levels present in the object are filled — empty factor levels are dropped first.

`title` replaces the auto-generated plot title (`"Frequency per <grouping_var>"` and its equivalents)
in all three modes; `NULL` keeps the default. Set it on subsets, where the column name does not say
what the denominator is.

Always prints the plot and invisibly returns it.

**Writes** — tables via `here::here()` (**project-root** relative), figures via `ggsave(path=)`
(**`getwd()`** relative). If those differ, your tables and figures land in different trees. Filenames
contain **literal spaces**:

- `"frequency per <grouping_var> per condition <object_annotations>.csv"`
- `"counts per <grouping_var> per condition <object_annotations>.csv"`
- `frequency_per_<grouping_var>_per_sample <object_annotations>.pdf` (12×6) for `'Samples'`
- `composition_of_groups_per_<grouping_var> <object_annotations>.pdf` (10×6) for `'Groups'`
- `composition_of_<grouping_var>_per_group <object_annotations>.pdf` (10×6) for `'grouping_var'`

---

## `top_genes_per_cluster()`

```r
top_genes_per_cluster(seurat, n_genes_to_plot = 3, grouping_var = 'seurat_clusters',
    object_annotations = '', tables_path = 'results/tables/',
    figures_path = 'results/figures/', results_path = 'results/',
    run_pathway_enrichment = NULL, n_genes_for_enrichment = 100, filter_ig = FALSE,
    filter_tcr = FALSE, flip_axes = TRUE, dot_scale = 8, ...)
```

⚠️ **`grouping_var` here is a quoted string** — the opposite convention from `extract_cell_counts()`.

**Returns a 6-element list:** `plot` (the dotplot), `ClusterProfiler_results`, `metascape_results`,
`topn`, `top100`, `seurat.markers`.

**Writes** to `tables_path` only, via `here::here()`: `top100_<object_annotations>.tsv`,
`top10_<object_annotations>.tsv`, `top100_gene_names_per_cluster_<object_annotations>.tsv`.

⚠️ **`figures_path` is accepted and documented but never used** — the dotplot is returned, not saved.
`ggsave()` it yourself.

**Preconditions.** `FindAllMarkers` runs on the **current `DefaultAssay` as-is** — no assay is set and
**`PrepSCTFindMarkers()` is not called**, so run it first when grouping by a custom or subclustered
variable on an SCT object with several models. Joins the internal **mouse** `annotations` table for
gene descriptions.

**Gotchas.**
- `run_pathway_enrichment` here is a character vector matched with `%in%` (e.g.
  `c('Metascape','ClusterProfiler')`), not the `method`-style dispatcher argument used elsewhere.
- Its Metascape branch **hard-codes `identities = 'seurat_clusters'`**, ignoring your `grouping_var`.
- `flip_axes = FALSE` mutates the metadata column with `fct_rev()`, so that column must be a factor.
- `results_path` is what reaches the enrichment functions, and therefore the destructive
  `create_analysis_directory()`.

---

## `annotate_seurat_with_SingleR_Eduard()`

```r
annotate_seurat_with_SingleR_Eduard(seurat, local_path, database = c('ImmGen'),
    annotation_basis = c('cluster_fine', 'cell_coarse', 'cell_fine'),
    split_by_groups = TRUE)
```

⚠️ **The defaults are unusable.** `annotation_basis` defaults to a length-3 vector and the dispatch is
a bare `if (annotation_basis == "cluster_fine")`, with no `match.arg()`. On R ≥ 4.2 that errors with
*"the condition has length > 1"*. **Always pass exactly one** of `'cluster_coarse'`, `'cluster_fine'`,
`'cell_coarse'`, `'cell_fine'`.

**One call adds one column.** To get all four, call it four times:

```r
seurat <- annotate_seurat_with_SingleR_Eduard(seurat, local_path = figures_path,
                                              annotation_basis = 'cluster_coarse')
seurat <- annotate_seurat_with_SingleR_Eduard(seurat, local_path = figures_path,
                                              annotation_basis = 'cell_coarse')
```

Adds `labels_per_cluster_{fine,coarse}` / `labels_per_cell_{coarse,fine}` and **sets `Idents()` to
it**. Returns the Seurat object.

**Preconditions.** `seurat_clusters` (hard-coded, both `cluster_*` branches); `Groups` when
`split_by_groups = TRUE`; an `RNA` assay **and an `SCT` assay** — it ends with an unconditional
`DefaultAssay(seurat) <- 'SCT'`, which errors without one. A `umap` reduction. Needs `SingleR` and
`celldex`.

**Mouse immune only** — `celldex::ImmGenData()`. For human, call `SingleR`/`celldex` directly with
`HumanPrimaryCellAtlasData` or `MonacoImmuneData`.

**Writes** one or two PDFs into `local_path` per call, e.g.
`UMAP_cluster_SingleR_annotations_coarse.pdf` (8×5) and `..._by_group.pdf` (10×5). The `cell_fine`
variants are 26×5 and 30×5. ⚠️ **Filenames carry no object tag**, so repeated calls on different
subsets overwrite each other — move or rename between calls.

---

## `annotate_seurat_with_LLM()`

```r
annotate_seurat_with_LLM(seurat, tissue,
    design = 'Single-cell RNA-seq; experimental design unspecified.',
    cluster_col = 'seurat_clusters', label_col = 'llm_celltype',
    model = 'claude-opus-4-8', markers = NULL, n_markers = 25, assay = NULL,
    use_pubmed = TRUE, pubmed_n = 10, pubmed_email = NULL, max_tokens = 1024,
    temperature = 0, api_key = NULL, keyring_service = 'anthropic',
    cache_dir = file.path(tempdir(), 'llm_pubmed_cache'), dry_run = FALSE, verbose = TRUE)
```

⚠️ **The default `model` is not a valid model id.** Every call fails at the API unless you override it
with a current one. Worse, API errors are swallowed into a per-row `error` column rather than
stopping — you get `cell_type = "unknown"` for every cluster and no exception. Check
`seurat@misc$llm_annotation$error`.

Returns the Seurat object with `<label_col>` and `<label_col>_confidence` added, and the full results
tibble at `seurat@misc$llm_annotation` (`cluster`, `n_cells`, `top_markers`, `cell_type`,
`confidence`, `rationale`, `n_refs`, `input_tokens`, `output_tokens`, `raw`, `error`).

`dry_run = TRUE` leaves metadata untouched and stores the prompts at
`seurat@misc$llm_annotation_prompts` — but **still hits PubMed** when `use_pubmed = TRUE`.

`cluster_col`, `label_col` and `assay` are quoted strings. If `markers` is supplied it needs `cluster`,
`gene` and `avg_log2FC`; otherwise `FindAllMarkers` runs. Needs `httr2 (>= 1.0.0)` and `jsonlite`.

Cost note: one PubMed esearch + efetch **and** one API call per cluster, with two `Sys.sleep(0.4)`
pauses each. Only the PubMed responses are cached (under `cache_dir`).
