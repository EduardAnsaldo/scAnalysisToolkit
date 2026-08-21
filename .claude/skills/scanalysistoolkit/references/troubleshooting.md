# Troubleshooting

Symptom → cause → fix, for the failures this package actually produces. Several entries are
pre-existing defects rather than usage errors; those are marked **[bug]** and the workaround is what
you want, not a fix to the call.

---

## Version and export problems

**`could not find function "plot_pathway_topic_genes"`** (or any `plot_pathway_*`,
`collapse_pathway_redundancy`, `pathway_topics`, `pathway_topic_terms`)
→ The installed build predates 0.11.0. Those eight exports do not exist.
→ `packageVersion('scAnalysisToolkit')`, then reinstall from GitHub.

**Redundancy figures never appear, and there is no error.**
→ `collapse_redundancy = TRUE` and `topics =` were absorbed by `...` on a pre-0.11.0 build. Silent by
construction — the DE wrapper's `...` accepts anything.
→ `'collapse_redundancy' %in% names(formals(pseudobulk))`. If `FALSE`, upgrade.

**`'gene_to_jaspar' is not an exported object from 'namespace:scAnalysisToolkit'`**
→ Pre-0.12.0. Twelve helpers were internal before that release.
→ Upgrade, or call with `:::` as a stopgap.

**`'run_figr_analysis' is not an exported object`**
→ Correct — it has a `man/` page but no `@export`, and it is **[bug]** broken anyway (`gr` read before
assignment; Seurat-v3 slot writes that fail on a v5 object).
→ Do not call it. Run the FigR steps directly and hand the results to `analyze_figr_results()`; see
`atac-regulatory.md`.

---

## Argument-passing problems

**A column "not found" whose name appears doubly quoted in the error, e.g. `"\"seurat_clusters\""`.**
→ A quoted string was passed to an embraced argument. `rlang::englue()` wrapped it in literal quotes.
→ Pass a **bare** name. Only `extract_cell_counts(grouping_var=)` and
`calculate_D50(cell_grouping_var=, replicate_var=, replicate_group_var=)` are embraced; everything
else wants a string.

**`the condition has length > 1`, from `annotate_seurat_with_SingleR_Eduard()`.** **[bug]**
→ `annotation_basis` defaults to a length-3 vector and the dispatch is a bare `if (x == "cluster_fine")`.
There is no `match.arg()`.
→ Always pass exactly one of `'cluster_coarse'`, `'cluster_fine'`, `'cell_coarse'`, `'cell_fine'`.

**`the condition has length > 1`, from `scatterplot()`.** **[bug]**
→ Same shape: `test_type` defaults to `c('Wilcox','Pseudobulk','Bulk')` with no `match.arg()`.
→ Always pass one explicitly. (`volcano_plot()` *does* use `match.arg`, so it is fine.)

**`formal argument "plot_subtitle" matched by multiple actual arguments`, from `pseudobulk()`.**
→ `pseudobulk()` derives `plot_subtitle`, `x_axis_label`, `up_label`, `down_label`, `x_column`,
`y_column`, `x_group_label`, `y_group_label`, `scatter_applicable` and `file_suffix` from
`contrast_labels` and passes them to the plots itself.
→ Do not pass any of those through `...`.

**`GO_GSEA_analysis()` behaves as if `group` were an ontology.**
→ `group` is required but sits *after* `ontology`, so a positional third argument binds to `ontology`.
→ Name it: `GO_GSEA_analysis(results, local_path, group = 'KO')`.

**`run_pathway_enrichment = TRUE` errors in `match.arg`.**
→ It is a string, not a logical.
→ Use `NULL`, `'ClusterProfiler'`, `'Metascape'`, or `c('ClusterProfiler','Metascape')`.

---

## Empty or wrong results

**`$UP_count` holds a message string instead of a number, and `$results` is missing.**
→ An early return. `'Not enough cells'` (fewer than `minimum_cell_number` in a group) or
`'Not enough biological replicates per group'` (fewer than 2 pseudobulk samples per group).
→ Check the returned list's length before indexing — on failure the DE wrappers return a **3-element**
list, or in `DEG_FindMarkers_RNA_assay`'s case a named **character vector**. For the replicate case,
switch to `DEG_FindMarkers_SCT_assay()`.

**Every sample is assigned to the same group in `pseudobulk()`.**
→ `condition` is derived by substring-matching `group1` then `group2` against the aggregated sample
names. `group1` is tested first, so with labels `CT` and `CTKO` every sample matches `CT`.
→ Rename to disjoint labels. Also ensure each sample name *contains* its group label.

**`pseudobulk()` groups do not match, and the names contain `-` where you wrote `_`.**
→ Seurat's `AggregateExpression()` rewrites `_` to `-` in group names. The code compensates, but with
`parameter_2` it also splits the aggregated name on `_` and assumes exactly one.
→ Avoid `_` inside group and sample *labels* when using `parameter_2`.

**`collapse_pathway_redundancy()` returns `NULL` with a message.**
→ Fewer than two terms passed `p.adjust < padj`. Expected, not an error.
→ Raise `padj`, or handle `NULL` before plotting — the `plot_pathway_*` functions return a
placeholder rather than failing.

**`plot_pathway_topic_genes()` bars are all grey when you passed `padj =`.**
→ `padj` is a **named numeric vector keyed by gene symbol**, and a grey bar means that symbol was not
found in it. Enrichment round-trips symbols through Entrez, so the casing may differ.
→ Check `names(padj)` against `res$terms$genes`.

**`object 'p1' not found` from `GO_overrepresentation_analysis_multiple_lists()`.** **[bug]**
→ When the comparison yields ≤1 term the plotting block is skipped but `p1` is still returned. There
is no no-data fallback (the MSigDB twin has one).
→ Use `pathway_overrepresentation_analysis_multiple_lists()` instead, or pre-check that the gene lists
enrich.

**`calculate_D50()` produces an empty or unlabelled plot.** **[bug]**
→ Two causes. `replicate_group_var = NULL` does **not** fall back to `replicate_var` as documented —
`englue()` on the default yields the literal `"NULL"`, so the plot is built with `aes(x = NULL)`. And
`cell_grouping_var` / `replicate_var` must be **factors**; on character columns `levels()` is `NULL`
and the loops silently produce nothing.
→ Always pass `replicate_group_var` explicitly, and `mutate(across(..., as.factor))` first.

**`plot_dorcs_heatmap()` errors in `fct_inseq()`, or group labels lose a letter.** **[bug]**
→ It does `str_replace(Sample, 'g', '')` then `fct_inseq()`, assuming `AggregateExpression`'s
`g`-prefixed **numeric** group names. With names like `"SPF_1"` it strips the first literal `g` and
then fails because the levels are not numeric.
→ Only usable when the `grouping_var` values are numbers. Otherwise build the heatmap directly.

**Metascape returns nothing useful for mouse data.**
→ The msbio command is hard-coded to `--source_tax_id 9606 --target_tax_id 9606` (human), while the
rest of the package is mouse.
→ Use the `pathway_*` (MSigDB) or `GO_*` families for mouse.

---

## Return-shape surprises

**Downstream code sees a list where a Seurat object should be.**
→ `perform_seurat_clustering()` and `process_multimodal_seurat()` return **lists**.
`seurat <- perform_seurat_clustering(seurat)` replaces your object with a list.
→ `res <- perform_seurat_clustering(seurat); seurat <- res$seurat`.

**Nothing happens after the elbow plot; the returned list has only two elements.**
→ A staged early return. `perform_seurat_clustering(regress_cell_cycle = TRUE, dimensions = NULL)`
returns `list(seurat, stage = "select_dimensions")`. `process_multimodal_seurat()` has three stages
(`"dimensions"`, `"resolutions"`, `"complete"`) keyed on which of five arguments are `NULL`.
→ Branch on `$stage`, pick your dimensions/resolutions from the printed plots, and call again.

**`cc <- run_cellchat_analysis(seurat)` gives `NULL` after a long run.** **[bug]**
→ The function ends in `return(invisible(NULL))` despite documenting a CellChat object.
→ Pass `output_file = 'path/to/cellchat.rds'` and `readRDS()` it.

**Attributes vanish from a `pathway_overrepresentation_analysis()` result.**
→ It returns a bare ggplot carrying `enrichment_results`, `pathway_clusters` and
`pathway_redundancy_plots` as **attributes**. Adding any layer (`p + theme(...)`) drops them.
→ Read the attributes first, then modify the plot.

**`ggsave()` fails on a plot from `plot_pathway_overlap()` or `plot_pathway_cluster_map(interactive = TRUE)`.**
→ The first returns a `ComplexHeatmap::Heatmap`, the second a plotly htmlwidget. Neither is a ggplot.
→ Use `ComplexHeatmap::draw()` into a device for the heatmap. In Quarto, a widget must be the
top-level value of its own chunk.

---

## Nothing printed, nothing saved

**`pseudobulk(collapse_redundancy = TRUE)` prints no figures at all.**
→ Intended. In that mode the ClusterProfiler path prints nothing — not even the plain
over-representation bar plots.
→ Place them yourself from `$enrichment$plot_up`, `$plot_down`,
`$redundancy_plots_up$clusters` / `$network` / `$topics`.

**A figure argument seems to be ignored.**
→ `top_genes_per_cluster()` accepts `figures_path` and never uses it — its dot plot is returned, not
saved. `run_cnmf_results()` prints its gene-loading heatmap but never saves it.
→ `ggsave()` the returned plot yourself.

**`cannot open file` / `No such file or directory` when writing a figure or table.**
→ Most functions do not create their output directory. See the "not created at all" list in
`contracts.md`.
→ `dir.create(figures_path, recursive = TRUE)` first.

**Figures and tables from one call land in different directories.**
→ `extract_cell_counts()` writes tables via `here::here()` (project-root relative) and figures via
`ggsave(path=)` (`getwd()` relative).
→ Pass absolute paths, or run with `getwd()` equal to the project root.

**Output lands beside the intended directory with a mangled name.**
→ Raw `paste0()` joining. `overlap_circos_and_tables`, `compare_DORC_scores_pairwise` and
`process_multimodal_seurat` need a **trailing `/`** on their path argument.

**A previous contrast's enrichment directory has disappeared.**
→ `create_analysis_directory()` calls `unlink()` first, and enrichment directory names are built from
one group name plus the cluster.
→ Vary `cluster`, or pass `file_suffix` (only `pathway_functional_analysis()` and
`Metascape_functional_analysis()` accept it; `pseudobulk()` supplies one automatically).

---

## Dependency and environment problems

**An interactive `check_installed()` prompt, or `The package "X" is required`.**
→ A `Suggests` dependency is missing. See the matrix in `contracts.md`.
→ Note `create_V_J_chord_diagram()` demands `scRepertoire` but never uses it — install it anyway.

**`FindClusters()` fails with `algorithm = 4`.**
→ `leidenbase` is the runtime backend. It is in `Imports` but never referenced by name, so a broken
install is not obvious.
→ Reinstall `leidenbase`, or pass `algorithm = 1` to `perform_seurat_clustering()`.

**`annotate_seurat_with_LLM()` returns `cell_type = "unknown"` for every cluster.** **[bug]**
→ API errors are swallowed into a per-row `error` column rather than stopping. The default
`model = "claude-opus-4-8"` is **not a valid model id**, so every call fails.
→ Inspect `seurat@misc$llm_annotation$error`, and pass a current `model` explicitly.

**`AnnotationPlot_custom()` or `calculate_D50()` breaks after upgrading Signac or scRepertoire.**
→ They call `Signac:::FindRegion`, `Signac:::reformat_annotations` and
`scRepertoire:::.expression2List` — unexported internals.
→ Pin the dependency, or avoid those two functions.
