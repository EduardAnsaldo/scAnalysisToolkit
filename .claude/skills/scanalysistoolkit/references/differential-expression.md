# Differential expression

`pseudobulk()` is the default. Use per-cell Wilcoxon only when the design has no biological
replicates. All arguments here are **quoted strings** — nothing in this domain uses `rlang` embrace.

**Shared contract for this file.** `path` is a *root*: every wrapper creates `<path>/figures/` and
`<path>/gene_lists/` under it. Enrichment output also lands under `path`, via
`create_analysis_directory()`, which **calls `unlink()` first** — so two contrasts written into one
`path` with the same cluster and group names delete each other. Vary `cluster`, or pass `file_suffix`.
`log2FoldChange` is positive when a gene is **up in `group2`** throughout, including the Wilcoxon
paths, where Seurat's `avg_log2FC` is sign-flipped to match.

---

## `pseudobulk()` — the default path

```r
pseudobulk(scRNAseq, comparison, group1, group2, cluster = 'all_clusters', path = './',
    shrink_LFC = FALSE, design_formula = NULL, parameter_2 = NULL,
    group1_parameter_2 = NULL, group2_parameter_2 = NULL, contrast = NULL,
    shrink_type = NULL, FC_threshold = 0.3, p_value_threshold = 0.05,
    max_overlaps = 15, expression_threshold_for_gene_list = 20,
    minimum_cell_number = 10, run_pathway_enrichment = NULL,
    collapse_redundancy = FALSE, topics = NULL,
    genes_to_exclude = c(), exclude_vdj = FALSE, ...)
```

Aggregates RNA counts per level of `comparison` — pass `'Samples'`, one column per replicate — then
runs DESeq2, so the **sample** is the unit of replication rather than the cell.

**Returns a 9-element list:** `all_count`, `UP_count`, `DOWN_count`, `results`, `dds`,
`contrast_labels`, `scatterplot`, `volcanoplot`, `enrichment`.

**On failure it returns only 3 elements** (`all_count`/`UP_count`/`DOWN_count` holding an error
string, no `results`, no plots). Check the length before indexing.

**Preconditions.** An `RNA` assay with a counts layer. ≥ `minimum_cell_number` cells and **≥ 2
pseudobulk samples per group**; below that you get `'Not enough biological replicates per group'`.
`condition` is derived by substring-matching `group1` then `group2` against the aggregated sample
names, so each sample name must contain its group label and **neither label may be a substring of the
other** (`CT` vs `CTKO` labels everything `CT`). Gene pre-filter: `rowSums(counts >= 10) >= min(n1, n2)`.

**Prints** its scatter and volcano — *unless* `collapse_redundancy = TRUE`, which suppresses all
enrichment printing (see below).

**Writes** under `path`: `gene_lists/ALL_GENES_DEG_Analysis_<cluster>_pseudobulk_<g2>_vs_<g1><tag>_.csv`
plus `DEG_UP_in_…` and `DEG_DOWN_in_…`; `figures/Pseudobulk_PCA_<cluster><tag>.pdf`,
`figures/Pseudobulk_scatter_DEG_in_<cluster><tag>.pdf`,
`figures/Pseudobulk_volcano_in_<cluster><tag>.pdf`; and `PCA_pseudobulk_<cluster>_<g2>_vs_<g1><tag>_.csv`.

### Designs and contrasts

`design_formula = NULL` and `contrast = NULL` fit `~condition` and compare `group2` to `group1`.
For a second factor:

- `parameter_2` names it; `group1_parameter_2` / `group2_parameter_2` give its levels.
  **Reference levels come from these arguments, not from the metadata factor order.**
- A mistyped `contrast` errors after `DESeq()` with the available `resultsNames()` listed.
- `shrink_LFC = TRUE` needs `shrink_type = 'apeglm'` (single coefficient, no `contrast`) or `'ashr'`.
  `shrink_type = 'apeglm'` with a non-`NULL` `contrast` is an explicit `stop()`. Other values fail —
  the guard tries to `check_installed()` a package named after the string itself.

`contrast_labels` reports how the contrast was classified in `$type`: `default`, `main_condition`,
`main_other`, `condition_in_stratum`, `other_in_stratum`, `interaction`, `unknown`. It also carries
`file_tag` (appended to every filename so several contrasts coexist in one `path`) and
`scatter_applicable`, which is `FALSE` for an interaction — the "scatterplot" is then a labelled
placeholder, because a difference of differences has no meaningful pair of group means to plot.

**Do not pass** `plot_subtitle`, `x_axis_label`, `up_label`, `down_label`, `x_column`, `y_column`,
`x_group_label`, `y_group_label`, `scatter_applicable` or `file_suffix` through `...` — `pseudobulk()`
sets all of them from `contrast_labels` and you will get *"formal argument matched by multiple actual
arguments"*.

### Enrichment

`run_pathway_enrichment` is a **string**: `NULL`, `'ClusterProfiler'`, `'Metascape'`, or both.
`'ClusterProfiler'` routes to the **MSigDB** family, not the `GO_*` family.

`collapse_redundancy = TRUE` (≥ 0.11.0) additionally collapses enriched terms into gene-overlap
clusters — and **suppresses all printing**, including the plain bar plots. Everything comes back on
`$enrichment`:

```r
de <- pseudobulk(seurat, comparison = 'Samples', group1 = 'WT', group2 = 'KO',
                 cluster = 'Tcells', path = out,
                 run_pathway_enrichment = 'ClusterProfiler',
                 collapse_redundancy = TRUE, topics = c('lipid', 'immune'))

de$enrichment$plot_up                          # over-representation bar plot
de$enrichment$redundancy_up                    # a pathway_clusters object
de$enrichment$redundancy_plots_up$clusters     # bar plot of clusters
de$enrichment$redundancy_plots_up$network      # term-overlap network
de$enrichment$redundancy_plots_up$topics$lipid$terms
de$enrichment$redundancy_plots_up$topics$lipid$genes
```

`$enrichment` is `NULL` when enrichment did not run; the four `redundancy*` elements are `NULL`
unless `collapse_redundancy = TRUE`. `collapse_redundancy` and `topics` do nothing under `'Metascape'`.
`minGSSize` passes through `...` and reaches `pathway_overrepresentation_analysis()`, whose default is
**4** (the GSEA and multiple-lists entry points use 5).

---

## `DEG_FindMarkers_SCT_assay()` — Wilcoxon fallback

```r
DEG_FindMarkers_SCT_assay(scRNAseq, comparison, group1, group2, is_integrated_subset = FALSE,
    cluster = 'all_clusters', path = './', FC_threshold = 0.3, p_value_threshold = 0.05,
    min_fraction = 0.01, max_overlaps = 15, label_size = 5, pathways_of_interest = NULL,
    label_threshold = 100000, distance_from_diagonal_threshold = 0.7,
    gene_lists_to_plot = NULL, expression_threshold_for_gene_list = 20,
    minimum_cell_number = 30, run_pathway_enrichment = NULL, collapse_redundancy = FALSE,
    topics = NULL, exclude_sex_chr_genes = FALSE, genes_to_exclude = c(),
    exclude_vdj = FALSE, ...)
```

Tests **cells**, not replicates, so prefer `pseudobulk()` whenever the design allows. Set `Idents()`
to `comparison` first, subset the object to one cell type, and pass `cluster = <celltype>`.

**Returns** `all_count`, `UP_count`, `DOWN_count`, `results`, `scatterplot`, `volcanoplot`,
`enrichment`. Both plot calls are wrapped in `tryCatch` — on failure the element is `NULL` with a
warning, but the counts still come back. On DE failure: a 3-element list.

**Preconditions.** An `SCT` assay; ≥ `minimum_cell_number` (default **30**) cells per group. Set
`is_integrated_subset = TRUE` for a subset of an SCT-integrated object — it flips `recorrect_umi`.

`group1`/`group2` must be **exact level values**, not patterns: they go straight to
`FindMarkers(ident.1=, ident.2=)`. Only the cell-count pre-check treats them as patterns.

**Writes** `gene_lists/*_Wilcox_<g2>_vs_<g1>_.csv` and `figures/Wilcox_{scatter_DEG,volcano}_in_<cluster>.pdf`.
`exclude_sex_chr_genes = TRUE` drops a fixed 12-gene mouse list (`Xist`, `Ddx3y`, `Uty`, …).

⚠️ **No `file_suffix` reaches enrichment here.** Two comparisons sharing `path`, `cluster` and one
group name will delete each other's enrichment directory.

⚠️ `gene_lists_to_plot` reuses the same output filenames for every list, so each overwrites the last.

`pathways_of_interest` and `expression_threshold_for_gene_list` are accepted but unused.

## `DEG_FindMarkers_RNA_assay()`

Same shape against the RNA `data` layer, with two differences worth knowing: it returns **no**
`scatterplot`/`volcanoplot` elements, and on failure returns a named **character vector** rather than
a list. Its `ALL_GENES_…` CSV is written from the *filtered* table — the full table only comes back
in `$results`.

## `bulk_analysis()` — true bulk

```r
bulk_analysis(counts_table, comparison = 'Groups', group1, group2, cluster = '', path = './', ...)
```

Takes a **data frame with a literal `genes` column**, not a Seurat object. Sample column names must
literally contain `group1` / `group2` (no `_`→`-` normalisation here). Needs ≥ 2 samples per group;
row filter `rowSums >= 10`. Returns `all_count`, `UP_count`, `DOWN_count`, `results`, `enrichment` —
the roxygen "named vector with 3 elements" describes only the failure path. `comparison` and
`minimum_cell_number` are declared but never used.

---

## The statistics cores

`pseudobulk_de()`, `DEG_FindMarkers_SCT_assay_de()`, `DEG_FindMarkers_RNA_assay_de()` and
`bulk_analysis_de()` are what the wrappers call. They return the counts and `results`, write the
`gene_lists/` CSVs and (for pseudobulk/bulk) the PCA, but produce **no scatter, no volcano and no
enrichment**. Call one only when you want the numbers alone.

`pseudobulk_de()` additionally returns `dds` and `contrast_labels`.

---

## Plotting

### `volcano_plot()`

```r
volcano_plot(results, group1, group2, cluster, local_figures_path, FC_threshold,
    p_value_threshold, max_overlaps = 15, label_size = 5,
    my_colors = c('#023FA5', '#8E063B', 'gray'),
    test_type = c('Wilcox','Pseudobulk','Bulk','Wilcox_ATAC','Wilcox_DORC',
                  'Wilcox_ATAC_closest_genes','Wilcox_ChromVar_motif'),
    genes_to_plot = NULL, pt_size = 1.5, pval_cap = 300, y_axis_max = NULL,
    nudge_x = NULL, nudge_y = NULL, label_line_width = 0.5, dashed_line_width = 0.5,
    plot_subtitle = NULL, x_axis_label = NULL, up_label = NULL, down_label = NULL,
    file_suffix = NULL, ...)
```

First seven arguments required. `test_type` uses `match.arg()`, so it defaults to `'Wilcox'` safely.
**Always prints and saves**, and returns the ggplot. Writes
`<local_figures_path>/<test_type>_volcano_in_<cluster><file_suffix>.pdf`; the directory must exist.

Silently drops non-significant genes outside the significant range (±3, or ±2 for
`Wilcox_ChromVar_motif`) and messages how many. Keep `pval_cap` below ~307 or `10^-pval_cap`
underflows.

### `scatterplot()`

```r
scatterplot(results, group1, group2, local_figures_path, FC_threshold, p_value_threshold,
    cluster = 'all_clusters', my_colors = c('#023FA5','#8E063B','gray'), max_overlaps = 15,
    label_size = 5, label_threshold = 10000, distance_from_diagonal_threshold = 0.5,
    test_type = c('Wilcox','Pseudobulk','Bulk'), genes_to_plot = NULL, pt_size = 1.3,
    x_column = NULL, y_column = NULL, x_group_label = NULL, y_group_label = NULL,
    plot_subtitle = NULL, file_suffix = NULL, scatter_applicable = TRUE, ...)
```

⚠️ **`test_type` has no `match.arg()`** — leaving it at the default length-3 vector errors with
*"the condition has length > 1"*. Always pass one.

Always prints and saves. Axes default to `Avg_<group1>` / `Avg_<group2>`, log10 with
`limits = c(0.5, max*6)`. `scatter_applicable = FALSE` short-circuits to a saved placeholder.

**Two `Avg_` conventions:** pseudobulk and bulk give DESeq2 **normalised counts**; the Wilcoxon paths
give **CPM**. The axis label follows `test_type`.

### `plot_deg_counts()`

```r
plot_deg_counts(deg_counts_df, figures_path, group2 = NULL, group1 = NULL, title = NULL,
    data_type = c('genes','peaks','motifs','dorcs'), width = 10, height = 6,
    up_color = '#8E063B', down_color = '#023FA5', file_suffix = NULL)
```

`deg_counts_df` must have exactly `cell_type`, `UP_count`, `DOWN_count` or it stops. Always prints
and saves, via raw `paste0()` — so `figures_path` is relative to `getwd()`, not the project root, and
must exist. Legend labels become "Up in `group2`" / "Up in `group1`" only when **both** are supplied.
(The roxygen text quotes different default colours than the code; the values above are the real ones.)

### `plot_gene_lists()`

```r
plot_gene_lists(results, gene_lists_to_plot, group1, group2, cluster, local_figures_path,
    FC_threshold, p_value_threshold, test_type, max_overlaps = 15, label_size = 5,
    label_threshold = 1e+05, distance_from_diagonal_threshold = 0.7, ...)
```

Scatter + volcano restricted to named gene lists. `gene_lists_to_plot` is a **named list** of symbol
vectors. Returns `NULL` invisibly; the output is the files.

⚠️ Every list writes to the **same filenames** — the list name never reaches the filename — so each
call overwrites the previous. Move or rename between calls.

### `add_gene_labels_layer()`

```r
add_gene_labels_layer(label_column, diffexpressed, label_size = 5, max_overlaps = 15,
    nudge_x = NULL, nudge_y = NULL)
```

Returns the `geom_text_repel` layer the DE plots use, for adding consistent gene labels to a plot you
build yourself.

---

## `validate_cell_counts()`

```r
validate_cell_counts(seurat, comparison, group1, group2, minimum_cell_number)
```

Returns `list(valid, n_group1, n_group2, message)`. Use it to decide whether a comparison is worth
running before paying for the DE, especially inside a loop over cell types where some subsets will be
too small.
