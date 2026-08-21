# Cell–cell communication, cNMF, pseudotime and utilities

**Shared contract for this file.** All arguments are **quoted strings**. Two of these functions have
surprising return values — `run_cellchat_analysis()` returns `NULL`, and `create_analysis_directory()`
deletes before it creates. Read those entries before use.

---

## CellChat

### `run_cellchat_analysis()`

```r
run_cellchat_analysis(seurat, group_by = 'major_cell_types', species = 'mouse',
    smooth_data = FALSE, subset_cell_groups = NULL, subset_db = NULL, min_cells = 10,
    method = 'triMean', trim = NULL, use_parallel = TRUE, n_workers = 10,
    max_memory_mb = 4000, population_size = FALSE, nboot = 100, output_file = NULL)
```

⚠️ **It returns `invisible(NULL)`**, despite documenting a CellChat object.
`cc <- run_cellchat_analysis(seurat)` gives you `NULL` after a multi-minute run. **The only way to get
the object is `output_file` + `readRDS()`:**

```r
run_cellchat_analysis(seurat, group_by = 'major_cell_types', species = 'mouse',
                      output_file = here('data/cellchat_KO.rds'))
cc <- readRDS(here('data/cellchat_KO.rds'))
```

`output_file` is a **full file path**, not a directory, and its parent must already exist.

**Setting up a call.**

- The object needs an assay literally named **`RNA`**. The function runs `NormalizeData()` itself, so
  pass raw counts. It does not join layers, so a v5 split-layer RNA assay may misbehave in
  `createCellChat` — `JoinLayers()` first.
- `group_by` is the quoted metadata column holding cell labels.
- If a `Samples` column exists it is copied to a lowercase `samples` column, which is what CellChat
  expects for multi-sample workflows. If your column is named otherwise, add `samples` yourself.
- `species` must be `'mouse'` or `'human'` — anything else stops. It selects `CellChatDB.mouse` /
  `CellChatDB.human` and, with `smooth_data = TRUE`, `PPI.mouse` / `PPI.human`. There is no way to
  pass a custom database object.
- `subset_db` goes to `subsetDB(search = subset_db, key = 'annotation')`. Valid values are the
  database's annotation levels: `"Secreted Signaling"`, `"ECM-Receptor"`, `"Cell-Cell Contact"`, and
  in newer builds `"Non-protein Signaling"`.
- `method` is `'triMean'` (default) or `'truncatedMean'`; `trim` (0–0.25) only matters for the latter.
- `subset_cell_groups` is passed as **both** `sources.use` and `targets.use` to `aggregateNet()` only
  — it does not restrict the earlier probability computation.

Pipeline: `subsetData` → `identifyOverExpressedGenes` → `identifyOverExpressedInteractions` →
(`smoothData`) → `computeCommunProb` → `filterCommunication(min.cells)` → `computeCommunProbPathway`
→ `aggregateNet`. `seed.use = 3514L`.

**Side effects beyond the file.** With `use_parallel = TRUE` it sets `future::plan("multisession")`
and raises `future.globals.maxSize`, then resets the plan to sequential at the end — but **does not
restore the option**. Needs `CellChat`.

### `netVisual_bubble_pathways()`

```r
netVisual_bubble_pathways(object, sources.use = NULL, targets.use = NULL, signaling = NULL,
    top.n = NULL, top.n.by = c('mean','max','sum'), remove.isolate = FALSE, min.prob = 0,
    sort.by.target = FALSE, do.flip = FALSE, show.lr.pairs = FALSE, color.use = NULL,
    color.heatmap = 'Spectral', min.quantile = 0, max.quantile = 1, show.legend = TRUE,
    scale.size = TRUE, size.range = c(1, 10), dot.size = 4, font.size = 10,
    font.size.title = 10, angle.x = 90, vjust.x = NULL, hjust.x = NULL,
    title.name = NULL, return.data = FALSE)
```

Pathway-level communication bubble plot. Returns a ggplot, or with `return.data = TRUE` a list of
`communication` (a tidy data frame with `source`, `target`, `pathway`, `prob`, `pair`) and `gg.obj`.
Never prints, never writes.

**Requires `object@netP$prob`** — i.e. `computeCommunProbPathway()` must have run — otherwise it stops
with an explicit message. `show.lr.pairs = TRUE` additionally needs `object@LR$LRsig`.

`sources.use` / `targets.use` accept **either** integer indices **or** character group names;
unknown names and out-of-range indices warn and are dropped. Unknown `signaling` names warn unless
*all* are unknown, which stops.

`color.use` must be a **named** vector (pathway → colour) to get categorical colouring; otherwise the
continuous `color.heatmap` ramp is used.

Convenient property: it only reads S4 slots and never calls CellChat functions, so it works on a saved
CellChat object without CellChat attached — provided the class definition can be loaded.

---

## `run_cnmf_results()`

```r
run_cnmf_results(seurat, data_dir, runname, k_used = 6, local_density_threshold = 0.1,
    do_cnmf = TRUE, sequential_palette = NULL, sequential_palette_dotplot = NULL,
    color_palette = NULL,
    output_path = here::here('results', 'cNMF', runname, 'outputs'),
    object_annotations = '', top_n = 100, topn_plot = 6,
    run_pathway_enrichment = NULL, ...)
```

Loads a finished cNMF run onto a Seurat object and produces program plots.

⚠️ **`do_cnmf = TRUE` shells out**: `conda run -n cNMF cnmf consensus --output-dir <data_dir> …`. That
requires a conda environment literally named **`cNMF`** with the `cnmf` CLI on `PATH`, and the exit
status is never checked. Set `do_cnmf = FALSE` when consensus has already been run.

**Files it reads**, under `here::here(data_dir, runname, …)`, with the threshold rendered as
`0.1` → `"0_1"`:

- `<runname>.usages.k_<k_used>.dt_<dt>.consensus.txt` — cells × K
- `<runname>.gene_spectra_score.k_<k_used>.dt_<dt>.txt` — K × genes
- `<runname>.gene_spectra_tpm.k_<k_used>.dt_<dt>.txt` — read and then **never used**, but must exist

⚠️ **Barcode alignment is unvalidated.** The usage file's row names must match `Cells(seurat)`
exactly; a mismatch silently produces `NA` program scores.

⚠️ **Metadata side effects.** The join happens twice, so `meta.data` ends up with **both**
numerically-named columns (`"1"`, `"2"`, …) and `cNMF_program_1..K`. Any pre-existing column starting
with a capital `X` is dropped.

**Returns** `seurat`, `top_colnames_long`, `top_colnames`, `ORA_results`, `genes_score_table`,
`spectra_score`.

**Requires** an `RNA` assay and a metadata column literally named **`Groups`** (hard-coded violin
`group.by`). Joins the internal mouse `annotations` table.

**Writes** into `output_path` (created): `FeaturePlot_cNMF_k_<k>_<object_annotations>.pdf`,
one `VlnPlot_<k><program>_by_group_<object_annotations>.pdf` per program, and
`top50_per_gene_program_k_<k>.tsv` — **misnamed**, it holds `top_n` genes (default 100). The gene-loading
heatmap is printed but **never saved**.

**Gotchas.** The `'Metascape'` branch is entirely commented out, so that value runs nothing and
returns `ORA_results = NULL`. `sequential_palette_dotplot` is in the signature but never used.

---

## `plot_pseudotime_heatmap()`

```r
plot_pseudotime_heatmap(seurat_obj, genes, pseudotime_col = 'pseudotime', assay = NULL,
    slot = 'data', normalize = c('minmax', 'minmax_pct', 'zscore', 'none'),
    pct_low = 0.01, pct_high = 0.99, smooth_window = NULL, loess_smooth = FALSE,
    loess_span = 0.1, gene_order = c('input', 'peak', 'hclust'), cell_subset = NULL,
    show_pseudotime_bar = TRUE, palette = NULL, palette_direction = 1,
    pseudotime_palette = 'magma', legend_title = NULL)
```

Cells × genes heatmap ordered by pseudotime. `normalize` and `gene_order` use `match.arg()`, so the
vector defaults are safe — they resolve to `'minmax'` and `'input'`.

Returns a ggplot, or by default a `patchwork` of the pseudotime bar above the heatmap. **No printing,
no disk output.** Species-agnostic, no extra dependencies.

`pseudotime_col` must exist or it stops. Genes missing from the assay warn and are dropped; all
missing stops. Cells with `NA` pseudotime are dropped.

`palette = NULL` resolves to a blue-white-red diverging ramp for `'zscore'` and viridis otherwise. A
palette may be a viridis option name, a Brewer name, a vector of ≥2 colours, or a ggplot `Scale`;
anything else stops.

---

## Utilities

### `apply_integer_breaks()`

```r
apply_integer_breaks(plot, break_size = 1,
    colors_use = scCustomize::viridis_plasma_dark_high, na_color = 'lightgray')
```

Re-applies a `FeaturePlot_scCustom()` colour scale with integer legend breaks. Accepts a ggplot, a
patchwork, or a plain list (from `combine = FALSE`) and returns the same type. No disk output.

It reads the existing scale's `limits`; if the input plot has no colour scale with explicit limits,
the limits reset to the data range. Only the **colour** aesthetic is touched, not `fill`.

⚠️ The patchwork branch **mutates in place** before returning.

### `create_analysis_directory()`

```r
create_analysis_directory(path, recursive = TRUE, showWarnings = FALSE)
```

⚠️ **Destructive.** It calls `unlink(path, recursive = TRUE)` *before* `dir.create()`, so an existing
directory of that name is deleted rather than merged into. Returns the path.

This is what the six `*_functional_analysis*` functions use for their output directories, and it is
why running two contrasts into the same `path` with the same cluster and group names destroys the
first one's enrichment output. Use it deliberately when you want a guaranteed-clean directory; use
`dir.create(recursive = TRUE)` when you do not.
