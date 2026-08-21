# Cross-cutting contracts

Per-function lookup tables. Open this when you have a specific function in hand and need to know what
its object must contain, where its output lands, or what must be installed.

---

## 1. Bare vs quoted column arguments

**Exactly two exported functions use `rlang` embrace (`{{ }}`), covering four arguments. Everything
else takes quoted strings.** Verified by mapping every `{{` and `englue()` occurrence in `R/` to its
enclosing function.

```r
extract_cell_counts(seurat, grouping_var = seurat_clusters,          # BARE
                    figures_path = figures_path, tables_path = tables_path,
                    plot_grouping = 'Groups')                        # ... but this one is a string

calculate_D50(seurat, cell_grouping_var = major_cell_types,          # BARE
              replicate_var = Samples, replicate_group_var = Groups, # BARE
              results_path, figures_path, tables_path)
```

A quoted string passed to an embraced argument **does not error**. `englue()` turns `'seurat_clusters'`
into the literal `"\"seurat_clusters\""`, which matches no column; you get an empty result or a
confusing downstream failure.

### Every other column-name argument (all quoted)

| function | argument(s) |
|---|---|
| `top_genes_per_cluster` | `grouping_var` (default `'seurat_clusters'`) |
| `top_peaks_per_cluster` | `grouping_var` (default `'seurat_clusters_atac'`) |
| `perform_seurat_clustering` | `filter_variable`, `integration_column`, `vars_to_regress` |
| `process_multimodal_seurat` | `integration_column` |
| `annotate_seurat_with_LLM` | `cluster_col`, `label_col`, `assay` |
| `overlap_circos_and_tables` | `grouping_variable`, `cell_types_column` |
| `plot_pathways_heatmap`, `plot_pathways_heatmap2`, `plot_chromvar_heatmap`, `plot_dorcs_heatmap` | `grouping_var` |
| `pseudobulk(_de)`, `DEG_FindMarkers_*(_de)`, `bulk_analysis(_de)` | `comparison`, `parameter_2` |
| `run_differential_accessibility(_FindMarkers)` | `comparison` |
| `compare_motive_activity_pairwise`, `compare_DORC_scores_pairwise` | `identities` |
| `analyze_figr_results` | `cluster_ident` |
| `Metascape_functional_analysis_cluster_identification` | `identities` |
| `run_cellchat_analysis` | `group_by` |
| `plot_pseudotime_heatmap` | `pseudotime_col` |
| `scatterplot` | `x_column`, `y_column` |
| `annotate_seurat_with_SingleR_Eduard` | *none* — `seurat_clusters` and `Groups` are hard-coded |

### ⚠️ The same argument name means three different things

`grouping_var` is:

- an **embraced bare column** in `extract_cell_counts`
- a **quoted column name** in `top_genes_per_cluster`, `top_peaks_per_cluster` and the four heatmaps
- a **free-text label**, never looked up as a column, in `pathway_functional_analysis`,
  `GO_functional_analysis`, `Metascape_functional_analysis` and
  `pathway_overrepresentation_analysis_multiple_lists` — it is only pasted into output directory and
  file names

Likewise `group_name` is a **value** of the `Groups` column in `plot_clonal_repertoire`, but a
**column name of the input table** in `create_V_J_chord_diagram`.

---

## 2. Path arguments and the disk-write contract

### `path` (a root) vs everything else (a leaf)

`path` is a **root** the function creates a subtree beneath — `<path>/figures/`, `<path>/gene_lists/`,
plus enrichment directories. `figures_path`, `tables_path`, `results_path`, `local_path`,
`local_figures_path`, `output_path`, `output_dir` are **leaves**, written into directly.

| argument | default | functions |
|---|---|---|
| `path` | `'./'` | `pseudobulk(_de)`, `DEG_FindMarkers_*(_de)`, `bulk_analysis(_de)`, `run_differential_accessibility(_FindMarkers)`, `GO_functional_analysis(_cluster_identification)`, `pathway_functional_analysis(_cluster_identification)`, `Metascape_functional_analysis(_cluster_identification)` |
| `figures_path` | **required** | `extract_cell_counts`, `plot_clonal_repertoire`, `calculate_D50`, `create_V_J_chord_diagram`, `overlap_circos_and_tables`, `plot_deg_counts`, `compare_motive_activity_pairwise`, `compare_DORC_scores_pairwise` |
| | `NULL` → all writes skipped | `perform_seurat_clustering`, `process_multimodal_seurat` |
| | `'results/figures/'` | `top_peaks_per_cluster`; `top_genes_per_cluster` (**accepted but never used**) |
| `tables_path` | **required** | `extract_cell_counts`, `calculate_D50`, `overlap_circos_and_tables` |
| | `'results/tables/'` | `top_genes_per_cluster`, `top_peaks_per_cluster` |
| `results_path` | `'results/'` | `top_genes_per_cluster`, `top_peaks_per_cluster` (passed on as enrichment `path`) |
| | **required but unused** | `calculate_D50` |
| | `'./'`, unused | `overlap_circos_and_tables` |
| `local_path` | **required** | `annotate_seurat_with_SingleR_Eduard`, all six `*_overrepresentation_analysis*` |
| `local_figures_path` | **required** | `scatterplot`, `volcano_plot` |
| `output_path` | `here::here('results','cNMF',runname,'outputs')` | `run_cnmf_results` |
| `output_dir` | `here::here('results/figr_plots')` | `analyze_figr_results` |
| `output_file` | `NULL` | `run_cellchat_analysis` — a **full file path**, not a directory |
| `cache_dir` | `file.path(tempdir(), 'llm_pubmed_cache')` | `annotate_seurat_with_LLM` |

### ⚠️ Three incompatible joining conventions

1. **`here::here(path, file)`** — resolves against the **project root** (the directory holding
   `.here` / the `.Rproj`), *not* `getwd()`. Pass relative paths. Used for every `.tsv`/`.csv` in
   `extract_cell_counts` and `top_genes_per_cluster`, the `overlap_circos_and_tables` CSV, the
   `create_V_J_chord_diagram` PDF, and every `create_analysis_directory()` target.
2. **`ggsave(filename = ..., path = ...)`** — resolves against **`getwd()`**. Used for every figure in
   `perform_seurat_clustering`, `extract_cell_counts`, `annotate_seurat_with_SingleR_Eduard`,
   `plot_clonal_repertoire`, `calculate_D50`.
3. **Raw `paste0()`** — `plot_deg_counts`, the circos writer reached from `overlap_circos_and_tables`,
   `calculate_D50`'s CSV, `compare_DORC_scores_pairwise`, and `process_multimodal_seurat`. The last
   three **require a trailing `/`** on the path argument (`process_multimodal_seurat` concatenates
   with no separator at all).

`extract_cell_counts()` uses conventions 1 and 2 in a single call. If your project root and `getwd()`
differ, its tables and figures land in different trees.

### Directory creation

**Created for you** (`dir.create(recursive = TRUE)`, additive and safe): `<path>/gene_lists` and
`<path>/figures` in every DE wrapper and core; `run_differential_accessibility*`;
`run_cnmf_results` (`output_path`); `analyze_figr_results` (`output_dir`); `annotate_seurat_with_LLM`
(`cache_dir`).

**Created destructively** — see §3.

**Not created at all** — the directory must already exist or the call fails:
`perform_seurat_clustering`, `extract_cell_counts`, `top_genes_per_cluster`, `top_peaks_per_cluster`,
`annotate_seurat_with_SingleR_Eduard`, `plot_clonal_repertoire`, `calculate_D50`,
`create_V_J_chord_diagram`, `overlap_circos_and_tables`, `plot_deg_counts`, `scatterplot`,
`volcano_plot`, `compare_motive_activity_pairwise`, `compare_DORC_scores_pairwise`, and all six
`*_overrepresentation_analysis*` (`local_path`).

---

## 3. ⚠️ `create_analysis_directory()` deletes

```r
create_analysis_directory <- function(path, recursive = TRUE, showWarnings = FALSE) {
    unlink(path, recursive = TRUE)      # <- the previous run's output is gone
    dir.create(path, recursive = recursive, showWarnings = showWarnings)
    return(path)
}
```

**Consequence:** run a DE wrapper twice into the same `path` and the second run deletes the first's
enrichment output instead of merging. Enrichment directory names are built from one group name plus
the cluster, so two contrasts sharing those collide.

**Reached directly by** the six `*_functional_analysis*` functions:

| function | directories it recreates |
|---|---|
| `pathway_functional_analysis` | `<path>/Pathway_functional_analysis_UP_in_<group2>_<grouping_var><file_suffix>` and the `group1` twin |
| `GO_functional_analysis` | `<path>/GO_functional_analysis_UP_in_<group2>_<grouping_var>` and the `group1` twin — **no `file_suffix` argument at all** |
| `Metascape_functional_analysis` | `<path>/Functional_analysis_metascape_UP_in_<g2>_<gv><file_suffix>` and `..._DOWN_in<g1>_<gv><file_suffix>` (note the missing `_` after `DOWN_in`) |
| `pathway_functional_analysis_cluster_identification` | `<path>/Cluster_identification_functional_analysis_Pathway_<object_annotations>` + `Cluster_<level>/` |
| `GO_functional_analysis_cluster_identification` | `<path>/Cluster_identification_functional_analysis_GO_<object_annotations>` |
| `Metascape_functional_analysis_cluster_identification` | `<path>/<object_annotations>_Cluster_identification_functional_analysis` + `Cluster_<level>/` |

**Reached indirectly** whenever `run_pathway_enrichment` is non-`NULL`: `pseudobulk`,
`DEG_FindMarkers_SCT_assay`, `DEG_FindMarkers_RNA_assay`, `bulk_analysis`,
`run_differential_accessibility`, `top_genes_per_cluster`, `top_peaks_per_cluster`.

**Mitigation.** Only `pathway_functional_analysis()` and `Metascape_functional_analysis()` accept
`file_suffix`, and only `pseudobulk()` supplies one automatically (from `contrast_labels$file_tag`).
For every other looped comparison, vary `cluster` or pass `file_suffix` yourself.

---

## 4. Required Seurat metadata columns

Hard-coded by name — these are not configurable unless a column argument is listed.

| column | required by | notes |
|---|---|---|
| `Samples` | `extract_cell_counts` | **mandatory in all three modes**, alongside `Groups` |
| | `plot_pathways_heatmap(2)`, `plot_chromvar_heatmap`, `plot_dorcs_heatmap` | default `grouping_var` — overridable |
| | `overlap_circos_and_tables` | default `grouping_variable` — overridable |
| | `process_multimodal_seurat` | when `generate_plots = TRUE` |
| | `run_cellchat_analysis` | *optional*; copied to a lowercase `samples` column if present |
| `Groups` | `extract_cell_counts` | **mandatory** |
| | `plot_clonal_repertoire` | **mandatory** — `group_name` is a value of it |
| | `perform_seurat_clustering` | only when `selected_resolution` is non-`NULL` |
| | `annotate_seurat_with_SingleR_Eduard` | only when `split_by_groups = TRUE` |
| | `process_multimodal_seurat`, `run_cnmf_results` | hard-coded `group.by` / `Idents` |
| `orig.ident` | `calculate_D50`, `overlap_circos_and_tables` | the `scRepertoire` split variable |
| `seurat_clusters` | `annotate_seurat_with_SingleR_Eduard` | **hard-coded**, both `cluster_*` branches |
| | `top_genes_per_cluster` | default `grouping_var`; **also hard-coded** in its Metascape branch |
| | `annotate_seurat_with_LLM` | default `cluster_col`, validated with a clear error |
| `seurat_clusters_atac` | `top_peaks_per_cluster`, `compare_motive_activity_pairwise`, `compare_DORC_scores_pairwise`, `analyze_figr_results` | defaults |
| `CTaa` | `plot_clonal_repertoire`, `calculate_D50`, `overlap_circos_and_tables` | scRepertoire output |
| `CTgene` | `overlap_circos_and_tables` | |
| `clonalFrequency`, `cloneSize` | `plot_clonal_repertoire` | from `combineExpression()` |
| `nCount_ATAC` | `top_peaks_per_cluster`, `run_differential_accessibility_FindMarkers` | hard-coded `latent.vars` |
| `Origin` | `perform_seurat_clustering`, `process_multimodal_seurat` | default `integration_column`, only when `integrate_rna = TRUE` |
| `pseudotime` | `plot_pseudotime_heatmap` | default `pseudotime_col`, validated |
| `major_cell_types` | `run_cellchat_analysis` | default `group_by` — overridable |

**Columns that must be factors, not character** (the code reads `levels()` and silently produces
nothing otherwise): `calculate_D50`'s `cell_grouping_var` and `replicate_var`;
`overlap_circos_and_tables`'s `grouping_variable`; `plot_chromvar_heatmap`'s `grouping_var`;
`top_genes_per_cluster`'s `grouping_var` when `flip_axes = FALSE`.

### Required assays and object state

| requirement | functions |
|---|---|
| `RNA` assay with counts | `perform_seurat_clustering`, `process_multimodal_seurat`, `annotate_seurat_with_SingleR_Eduard`, `plot_pathways_heatmap(2)`, `run_cellchat_analysis`, `pathway_functional_analysis_cluster_identification`, all four `*_de` cores |
| `SCT` assay | `DEG_FindMarkers_SCT_assay(_de)`; and `annotate_seurat_with_SingleR_Eduard`, which ends with an unconditional `DefaultAssay(seurat) <- 'SCT'` |
| `ATAC` ChromatinAssay | `top_peaks_per_cluster`, `run_differential_accessibility*`, `process_multimodal_seurat`, `analyze_figr_results`, `CoveragePlot_custom` (`assay_atac` default) |
| `peaks` ChromatinAssay | `SE_from_Signac` — **its `assayName` default is `"peaks"`, unlike everything else** |
| `chromvar` assay | `plot_chromvar_heatmap`, `compare_motive_activity_pairwise`, `analyze_figr_results` |
| `DORCs` assay | `plot_dorcs_heatmap`, `compare_DORC_scores_pairwise` |
| motifs attached (`Signac::AddMotifs`) | `find_enriched_motifs`, `compare_motive_activity_pairwise` |
| fragment files readable | `CoveragePlot_custom`, `create_tracks_patchwork(2)`, `analyze_figr_results` |
| `umap` reduction | `annotate_seurat_with_SingleR_Eduard` |
| `umap.wnn` reduction | `compare_motive_activity_pairwise` (default), `analyze_figr_results` (hard-coded in its driver block) |
| mouse gene symbols | anything joining the internal `annotations` table: `top_genes_per_cluster`, `top_peaks_per_cluster`, `run_cnmf_results`, the four `*_de` cores, `run_differential_accessibility_FindMarkers` |

---

## 5. Suggests matrix

The package loads without these; each is guarded by `rlang::check_installed()` and prompts at the
point of use.

| function | needs installed |
|---|---|
| `annotate_seurat_with_SingleR_Eduard` | `SingleR`, `celldex` |
| `annotate_seurat_with_LLM` | `httr2 (>= 1.0.0)`, `jsonlite`; `keyring` soft-checked |
| `run_cellchat_analysis` | `CellChat` |
| `plot_clonal_repertoire` | `packcircles` |
| `calculate_D50`, `overlap_circos_and_tables`, `create_V_J_chord_diagram` | `scRepertoire` |
| `pseudobulk(_de)` with `shrink_LFC = TRUE` | `apeglm` (no contrast) or `ashr` (with a contrast) |
| `collapse_pathway_redundancy` | `igraph` |
| `plot_pathway_cluster_map` | `igraph`; **plus `plotly`** when `interactive = TRUE` |
| `plot_pathway_overlap` | `ComplexHeatmap` |
| `top_peaks_per_cluster`, `run_differential_accessibility_FindMarkers`, `find_enriched_motifs`, `process_multimodal_seurat` | `Signac` |
| `CoveragePlot_custom`, `create_tracks_patchwork(2)` | `Signac`, `GenomicRanges` |
| `AnnotationPlot_custom` | `Signac`, `GenomicRanges`, `GenomeInfoDb`, `IRanges`, `fastmatch` |
| `SE_from_Signac` | `Signac`, `S4Vectors`, `SingleCellExperiment` |
| `gene_to_jaspar`, `ID_to_symbol` | `JASPAR2024`, `RSQLite`, `TFBSTools` |
| `compare_motive_activity_pairwise` | `JASPAR2024`, `RSQLite`, `TFBSTools`, `Signac` |
| `analyze_figr_results` | `ggrastr`, `Signac`, `FigR`, `BuenColors` (+ `ComplexHeatmap` via `plotfigRHeatmap_EA`) |
| `plotfigRHeatmap_EA` | `ComplexHeatmap`, `BuenColors` |
| `compare_DORC_scores_pairwise` | `FigR` |

**Needing nothing extra:** `perform_seurat_clustering`, `extract_cell_counts`,
`top_genes_per_cluster`, `plot_pathways_heatmap(2)`, `plotfigRHeatmap_EA_gg`, `plot_chromvar_heatmap`,
`plot_dorcs_heatmap`, `plot_pseudotime_heatmap`, `apply_integer_breaks`, `netVisual_bubble_pathways`,
`run_cnmf_results`, and the whole DE and enrichment path apart from the entries above.

### Anomalies worth knowing

- **`leidenbase`** sits in `Imports` but is never called by name. It is the runtime backend for
  `FindClusters(algorithm = 4)` — the default in `perform_seurat_clustering` and hard-coded in
  `process_multimodal_seurat` — so it must be installed even though nothing references it.
- **`create_V_J_chord_diagram` requires `scRepertoire`** via `check_installed()` but never uses it.
  A false dependency; install it or the call refuses to start.
- **`run_cnmf_results` shells out** to `conda run -n cNMF cnmf consensus …` when `do_cnmf = TRUE`,
  requiring a conda environment literally named `cNMF`. The exit status is never checked.
- **The Metascape family shells out** to `~/msbio_v3-2/bin/ms.sh`, writes into
  `~/msbio_v3-2/data/…`, copies results back, then deletes that directory. It needs a local msbio
  v3-2 install and is hard-coded to **human tax id 9606** in both `--source_tax_id` and
  `--target_tax_id`, while the rest of the package is mouse.
- **`:::` into other packages** — `Signac:::FindRegion`, `Signac:::reformat_annotations`
  (`AnnotationPlot_custom`) and `scRepertoire:::.expression2List` (`calculate_D50`,
  `overlap_circos_and_tables`). The most version-fragile lines in the package; an upgrade of either
  dependency can break them.

---

## 6. Version gates

Source and installed versions drift. Check `formals()` before relying on anything below.

| feature | minimum version | how it fails on an older build |
|---|---|---|
| `collapse_pathway_redundancy`, `plot_pathway_clusters`, `plot_pathway_cluster_map`, `plot_pathway_overlap`, `plot_pathway_topic`, `plot_pathway_topic_genes`, `pathway_topics`, `pathway_topic_terms` | 0.11.0 | `could not find function` |
| `collapse_redundancy =`, `topics =` on the DE wrappers | 0.11.0 | **silently absorbed by `...`** — no error, no redundancy output |
| `label_by =` on `collapse_pathway_redundancy` | 0.11.0 | as above |
| `interactive =` on `plot_pathway_cluster_map` | 0.11.0 | as above |
| `gene_to_jaspar`, `ID_to_symbol`, `plot_dorc_genes`, `create_tracks_patchwork(2)`, `create_group_palette`, `create_analysis_directory`, `validate_cell_counts`, `plot_gene_lists`, `add_gene_labels_layer`, `create_no_data_plot`, `set_enrichment_color_scale` | 0.12.0 | `object is not exported by 'namespace:scAnalysisToolkit'` |

```r
'collapse_redundancy' %in% names(formals(pseudobulk))              # >= 0.11.0
'gene_to_jaspar' %in% getNamespaceExports('scAnalysisToolkit')     # >= 0.12.0
packageVersion('scAnalysisToolkit')
```
