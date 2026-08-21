# ATAC, motifs and regulatory networks

**Shared contract for this file.** Every argument here is a **quoted string** — nothing in this domain
uses `rlang` embrace. Assay names are mostly **hard-coded and not configurable**: `'ATAC'` for
differential accessibility, `'chromvar'` for motif activity, `'DORCs'` for DORC scores.
`SE_from_Signac()` is the odd one out, defaulting to `'peaks'`. Most functions here do **not** create
their output directory. All of them need `Signac`, and several need `FigR`, `JASPAR2024` +
`TFBSTools`, `ComplexHeatmap` or `BuenColors`.

---

## Differential accessibility

### `run_differential_accessibility()`

```r
run_differential_accessibility(seurat, comparison, group1, group2,
    is_integrated_subset = FALSE, cluster = 'all_clusters', path = './',
    FC_threshold = 0.3, p_value_threshold = 0.05, min_fraction = 0.05,
    max_overlaps = 15, label_size = 5, pathways_of_interest = NULL,
    label_threshold = 100000, gene_lists_to_plot = NULL,
    colors = c('green4', 'darkorchid4'), minimum_cell_number = 30,
    run_pathway_enrichment = NULL, peak_to_gene = NULL,
    closest_feature_fallback = FALSE, ...)
```

Returns `all_count`, `UP_count`, `DOWN_count`, `results`.

⚠️ **On too few cells it returns a 3-element list with no `results` element at all** — the core
function's early return propagates. Handle both shapes.

**Preconditions.** An assay literally named `'ATAC'` (`DefaultAssay` is set to it) and a metadata
column `nCount_ATAC`, used as the hard-coded `latent.vars`.

**Writes.** Creates `<path>/figures/` and `<path>/gene_lists/`. Two volcano PDFs —
`Wilcox_ATAC_volcano_in_<cluster>.pdf` (peak-coordinate labels, `max_overlaps` forced to 1) and
`Wilcox_ATAC_closest_genes_volcano_in_<cluster>.pdf` (gene-name labels). Three CSVs from the core.
`run_pathway_enrichment` output lands under `path` via the destructive
`create_analysis_directory()`.

**Prints** both volcanoes; the plot objects are discarded rather than returned.

`is_integrated_subset` is documented but genuinely unused. `colors` is length 2 (DOWN, UP);
`"gray"` is appended for non-significant.

### `run_differential_accessibility_FindMarkers()`

The statistics core — same preconditions, no plots, no enrichment. Returns `all_count`, `UP_count`,
`DOWN_count`, `results`.

⚠️ On too few cells the three counts hold the **character string** `'Not enough cells'`, not numbers,
and `results` is `NULL`.

⚠️ `group1`/`group2` are used two ways: `str_detect()` for the cell-count check but
`FindMarkers(ident.1=, ident.2=)` for the test. Despite the docs calling them patterns, they must be
**exact level values**.

`log2FoldChange` is **sign-flipped** relative to Seurat's `avg_log2FC`, so positive means up in
`group2`, matching the RNA path. Results carry an internal `has_correlated_gene` column that this
function does not strip (the wrapper does).

Writes into `<path>/gene_lists/`: `ALL_PEAKS_DA_Analysis_…`, `PEAKS_UP_in_<group2>_…`, and
`PEAKS_UP_in_DOWN_in_<group2>_…` — the last name is as odd as it looks.

`peak_to_gene` needs columns `PeakRanges` (colon form `chr1:1-2`), `Gene` and `rObs`.

### `top_peaks_per_cluster()`

```r
top_peaks_per_cluster(seurat, n_genes_to_plot = 3, grouping_var = 'seurat_clusters_atac',
    object_annotations = '', tables_path = 'results/tables/',
    figures_path = 'results/figures/', results_path = 'results/',
    run_pathway_enrichment = NULL, peak_to_gene = NULL,
    closest_feature_fallback = FALSE, ...)
```

Returns `ClusterProfiler_results`, `metascape_results`, `topn`, `da_results`. The `da_results` element
is exactly what `find_enriched_motifs()` wants.

Test is hard-coded `FindAllMarkers(test.use = 'LR', only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)`.

⚠️ `figures_path` is accepted but **never used** — no figures are written. `tables_path` is **not
created**. With the default `object_annotations = ''` filenames end `_..tsv` (e.g. `top100_peaks_.tsv`).

⚠️ Its ClusterProfiler branch passes `identities = grouping_var` to a function that has no such
parameter, so it lands silently in `...` and your `grouping_var` is ignored there.

---

## Motifs and chromVAR

### `find_enriched_motifs()`

```r
find_enriched_motifs(seurat_object, da_results, cluster_id,
    pct_threshold = 0.05, log2fc_threshold = 1, n_motifs_to_plot = 10)
```

⚠️ The roxygen prose says `pct_threshold` defaults to 0.1; the real default is **0.05**.

Returns `list(enriched_motifs, plot, n_peaks)`, or **`NULL`** if no peaks survive filtering.
**Prints** the motif-logo plot as well as returning it. No disk output.

`da_results` must have `peak`, `cluster`, `pct.1`, `avg_log2FC`, `p_val_adj` — i.e.
`top_peaks_per_cluster()$da_results`. The object's **`DefaultAssay` must already be the
motif-annotated ChromatinAssay** (the function does not set it), with `Signac::AddMotifs()` run and a
GC background available.

### `plot_chromvar_heatmap()`

```r
plot_chromvar_heatmap(motifs_to_plot, seurat, title,
    color_palette = grDevices::hcl.colors(n = 20, 'RdBu', rev = TRUE),
    grouping_var = 'Samples', motif_labels = FALSE)
```

Returns a ggplot. No printing, no disk output, no extra dependencies. Needs an assay named `chromvar`.

⚠️ **`grouping_var` must be a factor.** The aggregation is built from `levels(group_vec)`; on a
character column `levels()` is `NULL` and the function collapses.

Stops on duplicated `motifs_to_plot`, or if `motif_labels = TRUE` and the vector is unnamed.

### `compare_motive_activity_pairwise()`

```r
compare_motive_activity_pairwise(seurat, cluster_id, other_cluster_id,
    cluster_title = NULL, figures_path, identities = 'seurat_clusters_atac',
    reduction = 'umap.wnn', topn = 9, sequential_palette, ...)
```

⚠️ `figures_path` and `sequential_palette` have **no defaults** and are required, and `figures_path`
sits after an argument that has one — pass by name.

Returns the `FindMarkers` data frame (rownames = motif IDs) plus a `motif_name` column, or **`NULL`
with a warning** if either group has fewer than 3 cells. **Prints** heavily and also saves.

**Preconditions.** A `chromvar` assay **and** an `ATAC` assay with motifs attached — it reads
`seurat@assays$ATAC@motifs@data` directly. A reduction named by `reduction`. A live JASPAR2024
connection. `figures_path` must exist.

Writes six PDFs per call (`motif_sequences_{UP,DOWN}_…`, `motif_activities_{UP,DOWN}_…`,
`motif_activities_violin_{UP,DOWN}_…`) plus the volcano.

Note the FeaturePlots are drawn on the **full** object while the violins use the subset.

### `gene_to_jaspar()` / `ID_to_symbol()`

```r
gene_to_jaspar(gene_symbol, species = 'Homo sapiens', collection = 'CORE')
ID_to_symbol(motif_id)
```

Gene symbol ↔ JASPAR2024 motif ID. Need `JASPAR2024`, `RSQLite`, `TFBSTools`.

⚠️ `gene_to_jaspar()`'s `species` argument is **inert** — it is never placed into the `getMatrixSet`
`opts`, so it does not filter. The default also disagrees with the rest of the package, which is mouse.

---

## FigR / DORCs

### The FigR gap

`run_figr_analysis()` has a `man/` page but **no `@export`**, and it is **broken**:
`R/GRN_analysis_FigR.r:300` reads `gr` before it is assigned, and it uses Seurat-v3 slot writes
(`@counts <-`, `@data <-`) that fail on a v5 `Assay5` object. Do not call it, even via `:::`.

Drive the upstream steps yourself and hand the outputs to the exported `analyze_figr_results()`:

```r
atac_se <- SE_from_Signac(seurat, assayName = 'ATAC')     # note: default is 'peaks'
cis     <- FigR::runGenePeakcorr(atac_se, RNAmat, genome = 'mm10')
dorcs   <- FigR::dorcJPlot(cis, cutoff = 7)
dorcMat <- FigR::getDORCScores(atac_se, dorcTab = cis, geneList = dorcs)
dorcMat.s <- FigR::smoothScoresNN(NNmat, dorcMat)
RNAmat.s  <- FigR::smoothScoresNN(NNmat, RNAmat)
figR.d  <- FigR::runFigRGRN(ATAC.se = atac_se, dorcTab = cis, dorcMat = dorcMat.s,
                            rnaMat = RNAmat.s, genome = 'mm10')

analyze_figr_results(figR.d, seurat, dorcMat.s, RNAmat.s, umap.d, top_dorcs_genes)
```

### `SE_from_Signac()`

```r
SE_from_Signac(obj, assayName = 'peaks', fetchRawCounts = TRUE, reductions = NULL)
```

Returns a `SingleCellExperiment` — `assay(se, "counts")` is peaks × cells, `rowRanges()` the peak
ranges, `colData()` the metadata, `reducedDims()` the requested embeddings (`NULL` takes all).

⚠️ The default `assayName` is **`"peaks"`**, inconsistent with the `'ATAC'` hard-coded everywhere else.
Stops if the named assay is absent or is not a `ChromatinAssay`.

With `fetchRawCounts = FALSE` you must also set `normalizeATACmat = FALSE` in the FigR calls — the
function only warns.

### `analyze_figr_results()`

```r
analyze_figr_results(figR.d, seurat, dorcMat.s, RNAmat.s, umap.d, top_dorcs_genes,
    output_dir = here::here('results/figr_plots'),
    score_cuts = list(scatter = 1, heatmap1 = 2, heatmap2 = 1.5, network = 1.6),
    drivers = NULL, genes = NULL, cluster_ident = 'seurat_clusters_atac',
    umap_reduction = 'umap.atac', species = 'Mus musculus')
```

The first **six** arguments are required.

⚠️ **Every step is wrapped in `tryCatch` and only warns**, so the function never stops and the returned
list can be missing keys or hold `NULL`s. Possible elements: `scatter_plot`, `rank_by_mean`,
`rank_by_targets`, `heatmap1`, `heatmap2`, `heatmap_clustered`, `network`, `dorc_plots`,
`driver_plots`, `detailed_dorc_plots`. Check with `names()` before indexing.

**Preconditions.** `figR.d` needs `Corr.log10P`, `Enrichment.log10P`, `Score`, `Motif`, `DORC`. The
object needs assays `chromvar`, `SCT` and `ATAC`, a reduction named **`umap.wnn`** (hard-coded in the
driver block — `umap_reduction` only affects the detailed-DORC block), and readable fragment files.

⚠️ It builds a `DORCs` assay internally but **on a local copy that is never returned**, so you do not
get it back. Create it yourself if you want `plot_dorcs_heatmap()` or
`compare_DORC_scores_pairwise()` afterwards.

⚠️ `drivers` is documented but **never used** — drivers are always recomputed as the top-8/bottom-4
motifs by mean Score. `species` is effectively inert (see `gene_to_jaspar()` above).

Creates `output_dir` and writes `cluster_info.png`, `scatter_plot.png`, `rank_by_mean.png`,
`heatmap1.png`, `heatmap2.png`, `heatmap_clustered.png`, `network.png`, plus `dorc_<gene>.png`,
`driver_<gene>.png` and `detailed_dorc_<gene>.png`. Prints nearly everything.

Needs `ggrastr`, `Signac`, `FigR`, `BuenColors` (and `ComplexHeatmap` via the heatmap helper).

### `plot_dorc_genes()`

```r
plot_dorc_genes(results, genes = NULL, umap_reduction = 'UMAP.WNN')
```

`results` is the list shape `run_figr_analysis()` would have returned — it reads `$ATAC_se`,
`$dorcGenes`, `$dorcMat_smoothed`, `$RNAmat_smoothed`. Returns a list of patchworks (DORC | RNA); does
not print, does not save. `genes = NULL` takes the first 15 DORC genes not starting with a digit.

⚠️ The default reduction name is **`"UMAP.WNN"`** (uppercase), unlike `analyze_figr_results()`'s
`"umap.atac"`. It has no `check_installed()` guard despite calling `FigR::`.

### `plot_dorcs_heatmap()`

```r
plot_dorcs_heatmap(dorcs_to_plot, seurat, title,
    color_palette = grDevices::hcl.colors(n = 20, 'RdBu', rev = TRUE),
    grouping_var = 'Samples')
```

Returns a ggplot; no printing, no disk output. Needs an assay named `DORCs`.

⚠️ **Only works when the `grouping_var` values are numbers.** It does
`str_replace(Sample, 'g', '')` then `fct_inseq()`, assuming `AggregateExpression`'s `g`-prefixed
numeric group names. With `"SPF_1"` it strips the first literal `g` and then errors because the levels
are not numeric.

### `compare_DORC_scores_pairwise()`

```r
compare_DORC_scores_pairwise(seurat, dorcMat.s, RNAmat.s, umap.d, cluster_id,
    other_cluster_id, figures_path, identities = 'seurat_clusters_atac',
    topn = 9, cluster = NULL, ...)
```

First seven required. Returns the raw `FindMarkers` data frame. ⚠️ **No cell-count guard** — with too
few cells `FindMarkers` itself errors.

Needs an assay literally named `'DORCs'`. `figures_path` must exist and is joined with raw `paste0()`.
Prints one combined DORC/RNA plot per gene and writes `<gene>_DORC_RNA_plot_{up,down}.png` plus the
volcano. Genes whose names start with a digit are silently dropped from the plots but stay in the
returned table. Needs `FigR`.

### `plotfigRHeatmap_EA()` / `plotfigRHeatmap_EA_gg()`

```r
plotfigRHeatmap_EA(figR.d, score.cut = 1, DORCs = NULL, TFs = NULL, ...)

plotfigRHeatmap_EA_gg(figR.d, score.cut = 1, DORCs = NULL, TFs = NULL,
    DORC_groups = NULL, cluster_DORCs = TRUE, cluster_TFs = TRUE,
    color_palette = grDevices::hcl.colors(n = 20, 'RdBu', rev = TRUE),
    limits = c(-2, 2), plot_title = NULL, x_title = 'TFs', y_title = 'DORCs',
    highlight_TFs = NULL, highlight_label_color = 'red', highlight_box_color = 'red',
    highlight_box = TRUE, highlight_box_linewidth = 1, group_tag = 'DORC groups')
```

TF × DORC regulatory-score heatmaps. Neither prints or saves.

The first returns a **`ComplexHeatmap::Heatmap`** (needs `ComplexHeatmap` + `BuenColors`; `...` goes
to `Heatmap()`); missing TF–DORC pairs stay `NA`, and its pearson clustering fails on a
single-row/column matrix.

The second returns a plain **ggplot** with **no extra dependencies** — prefer it unless you need
ComplexHeatmap features. It zero-fills missing pairs and **uppercases the x tick labels**.

⚠️ Their filtering semantics differ: in the `_gg` version, supplying both `DORCs` and `TFs` ignores
`score.cut` entirely, and supplying one applies `score.cut` only to the other axis.

---

## Genome-track plots

### `CoveragePlot_custom()`

```r
CoveragePlot_custom(seurat, region, theme = NULL, title = NULL, subtitle = NULL,
    title_theme = NULL, assay_atac = 'ATAC', extend_upstream = 20000,
    extend_downstream = 20000, annotation_font_size = 1.5, add_links = FALSE,
    gene_symbol = FALSE, expression_feature = NULL, assay = 'SCT',
    heights = c(10, 1, 2, 3), widths = c(10, 1))
```

Returns a patchwork from `Signac::CombineTracks()`. **Does not print, does not save** — the caller
does both.

`region` may be `"chr1-1000-2000"`, `"chr1:1000-2000"`, or a **gene symbol** resolved against the
assay's annotation. Needs readable fragment files.

⚠️ **There is no `group.by`** — the coverage track is grouped by whatever `Idents(seurat)` currently
is. Set idents before calling.

⚠️ Track order is **coverage, gene annotation, peak, link**, though the roxygen `heights` text claims
coverage, peak, annotation, link. `heights` is subset to the non-`NULL` tracks, so with
`add_links = FALSE` only the first three apply.

### `create_tracks_patchwork()` / `create_tracks_patchwork2()`

```r
create_tracks_patchwork2(seurat, regions, theme = NULL, nrow = NULL, ncol = NULL,
    title = NULL, subtitle = NULL, title_theme = NULL, assay_atac = 'ATAC',
    extend_upstream = 20000, extend_downstream = 20000, annotation_font_size = 1.5,
    add_links = FALSE, gene_symbol = FALSE, expression_features = NULL, assay = 'SCT',
    heights = c(10, 1, 2, 3), widths = c(10, 1), add_track_titles = FALSE,
    track_title_suffix = NULL, track_titles = NULL, track_title_theme = NULL)
```

Multi-region grid versions. Both return a patchwork and neither prints or saves. With both `nrow` and
`ncol` `NULL` the grid is `ceiling(sqrt(n))` columns.

⚠️ Their track order is **coverage, peak, gene annotation, link** — different from
`CoveragePlot_custom()`.

Version 1 does not expose `annotation_font_size` (fixed at 1.5). **Prefer version 2**, which adds it
plus the per-track title arguments.

### `AnnotationPlot_custom()`

```r
AnnotationPlot_custom(object, region, assay = NULL, mode = 'gene', sep = c('-', '-'),
    extend.upstream = 0, extend.downstream = 0, annotation_font_size = 1.5)
```

The gene-annotation track alone. Returns a ggplot, or **`NULL`** if the object has no annotation.
Stops unless `mode` is `'gene'` or `'transcript'`.

⚠️ Calls **`Signac:::FindRegion()` and `Signac:::reformat_annotations()`** — unexported internals, so a
Signac upgrade can break it. Needs `Signac`, `GenomicRanges`, `GenomeInfoDb`, `IRanges`, `fastmatch`.

Note the dot-separated argument names here (`extend.upstream`) versus the underscored ones in
`CoveragePlot_custom` (`extend_upstream`). `annotation_font_size` is ggplot text size in mm, so 1.5 is
deliberately small.
