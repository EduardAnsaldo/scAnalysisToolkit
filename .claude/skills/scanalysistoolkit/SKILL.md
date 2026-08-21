---
name: scanalysistoolkit
description: >
  API reference for the scAnalysisToolkit R package. Use when writing, fixing or reviewing R code
  that calls it; when about to hand-roll Seurat, Signac, SingleR, DESeq2, clusterProfiler,
  scRepertoire or CellChat code for single-cell work, to check whether a wrapper already exists;
  or when a scAnalysisToolkit call errors, returns early, writes nothing, or prints nothing.
  Covers all 77 exports — signatures, return shapes, required Seurat metadata columns,
  bare-vs-quoted column arguments, the path/figures_path/tables_path disk-write contract,
  Suggests dependencies and version gating, across differential expression, pathway enrichment
  and redundancy collapse, clustering and annotation, TCR repertoire, ATAC differential
  accessibility, motifs, FigR gene-regulatory networks, CellChat, cNMF and pseudotime. Do NOT
  use to scaffold or fill an analysis .qmd or to drive an end-to-end pipeline (use
  scrnaseq-seurat-analysis), to develop the package itself, or for scanpy/Python workflows.
version: 1.0.0
---

# scAnalysisToolkit — API reference

Helper functions and wrappers for single-cell RNA-seq, single-cell VDJ and single-nuclei ATAC-seq.
**Seurat objects only** — no scran, no scanpy.

Most of these are not pure functions. For a majority of the wrappers, **the primary effect is what
they write to disk** — figures, gene-list CSVs, enrichment directories — and the returned value is
secondary. Read the disk-write contract before assuming a return value is the output.

This skill answers *"does a wrapper exist, what does it take, what does it return, why did my call
fail."* It does not drive a pipeline; `scrnaseq-seurat-analysis` does that.

> **Reference files are self-contained.** Each file under `references/` repeats whatever shared
> contract it needs. You should never have to read two of them to make one call.

## Install and version-gate

```r
# renv::install('EduardAnsaldo/scAnalysisToolkit'); renv::snapshot()   # renv projects
# pak::pak("EduardAnsaldo/scAnalysisToolkit")                          # otherwise
library(scAnalysisToolkit)
```

Never `devtools::load_all()`, never pin an old version.

**The installed build often lags the repo.** New arguments are absorbed by `...` with no error, so a
call against an old build runs and silently does less. Check before relying on anything recent:

```r
'collapse_redundancy' %in% names(formals(pseudobulk))   # FALSE on < 0.11.0
```

If `FALSE`: the eight `collapse_pathway_redundancy` / `plot_pathway_*` / `pathway_topics` exports do
not exist, and `collapse_redundancy = TRUE` / `topics =` do nothing at all. Upgrade rather than
working around it.

## Don't hand-roll it

Reach for the wrapper before writing the underlying call yourself:

| about to write | use instead |
|---|---|
| `FindAllMarkers()` + a dotplot | `top_genes_per_cluster()` |
| `SCTransform` → `RunPCA` → `RunUMAP` → `FindNeighbors` → `FindClusters` | `perform_seurat_clustering()` |
| the same for RNA **+** ATAC with a WNN graph | `process_multimodal_seurat()` |
| `AggregateExpression()` + `DESeq2::DESeq()` | `pseudobulk()` |
| `FindMarkers()` between two conditions | `DEG_FindMarkers_SCT_assay()` (or `pseudobulk()`) |
| `SingleR::SingleR()` + a labelled DimPlot | `annotate_seurat_with_SingleR_Eduard()` |
| `EnhancedVolcano`, or a ggplot volcano by hand | `volcano_plot()` |
| a log-log DE scatter by hand | `scatterplot()` |
| `clusterProfiler::enricher()` over MSigDB | `pathway_overrepresentation_analysis()` |
| `clusterProfiler::enrichGO()` | `GO_overrepresentation_analysis()` |
| `clusterProfiler::simplify()` (it does not work across mixed collections) | `collapse_pathway_redundancy()` |
| `Signac::FindMarkers()` on peaks | `run_differential_accessibility()` |
| `Signac::CoveragePlot()` with extra tracks | `CoveragePlot_custom()` |
| `Signac::FindMotifs()` | `find_enriched_motifs()` |
| `scRepertoire::clonalDiversity()` | `calculate_D50()` |
| a clonotype-overlap circos by hand | `overlap_circos_and_tables()` |
| `CellChat::netVisual_bubble()` | `netVisual_bubble_pathways()` |
| a pseudotime-ordered gene heatmap by hand | `plot_pseudotime_heatmap()` |
| `dir.create()` for a fresh output dir | `create_analysis_directory()` — but read the warning below |

## Calling rules

Six things that are not visible in a signature. Details and per-function tables in
`references/contracts.md`.

**⚠️ Bare vs quoted column arguments.** Every column-name argument in the package is a **quoted
string**, with exactly two exceptions that take **bare, unquoted** names:

```r
extract_cell_counts(grouping_var = seurat_clusters, ...)                    # BARE
calculate_D50(cell_grouping_var = major_cell_types, replicate_var = Samples,
              replicate_group_var = Groups, ...)                            # BARE
top_genes_per_cluster(grouping_var = 'seurat_clusters')                     # quoted
```

Passing a quoted string to an embraced argument does not error — `rlang::englue()` mangles it into
the literal `"\"seurat_clusters\""`, which matches no column, and you get a confusing downstream
failure instead.

**⚠️ `create_analysis_directory()` deletes first.** It calls `unlink(path, recursive = TRUE)` before
`dir.create()`. The practical consequence: **running a DE wrapper twice into the same `path` deletes
the first run's enrichment output** rather than merging. Enrichment directory names are built from
one group name plus the cluster, so two contrasts sharing those collide. Only
`pathway_functional_analysis()` and `Metascape_functional_analysis()` accept `file_suffix`, and only
`pseudobulk()` passes one automatically (from `contrast_labels$file_tag`). Vary `cluster` or pass
`file_suffix` explicitly for every other looped comparison.

**Print vs return is argument-dependent.** `pseudobulk()` normally prints its own scatter and volcano.
With `collapse_redundancy = TRUE` the ClusterProfiler path prints **nothing at all** — not even the
plain over-representation bar plots — and every figure comes back on the result for you to place:
`$enrichment$plot_up`, `$plot_down`, `$redundancy_plots_up$clusters` / `$network` / `$topics`.

**Staged early returns.** `perform_seurat_clustering()` and `process_multimodal_seurat()` return
**lists, not Seurat objects**, and sometimes short ones. `seurat <- perform_seurat_clustering(seurat)`
silently replaces your object with a list. Branch on `$stage` before using `$seurat`:

```r
res <- perform_seurat_clustering(seurat, dimensions = 30, selected_resolution = 0.25)
seurat <- res$seurat
```

**Path arguments are two different things.** `path` is a **root** the function builds a subtree
under (`<path>/figures/`, `<path>/gene_lists/`). `figures_path`, `tables_path`, `results_path`,
`local_path` and `local_figures_path` are **leaf** directories written into directly — and most are
**not created for you**. Three incompatible joining conventions coexist; see `contracts.md`.

**Everything is mouse except Metascape.** `org.Mm.eg.db` is a hard Import, the internal gene
`annotations` table is mouse, and the genome defaults are `mm10` / `GRCm39`. The Metascape family
shells out to a local install hard-coded to human tax id 9606.

## Function index

77 exports plus one S3 method. Column two is what you would reach for it *for*.

### Clustering and annotation — `references/clustering-annotation.md`

| function | reach for it when |
|---|---|
| `perform_seurat_clustering()` | SCTransform → PCA → UMAP → clustering at many resolutions; returns a **list** |
| `process_multimodal_seurat()` | the same for RNA + ATAC with a WNN graph; three-stage interactive pipeline |
| `extract_cell_counts()` | cluster/cell-type frequency tables and composition plots (**bare** `grouping_var`) |
| `top_genes_per_cluster()` | `FindAllMarkers` dotplot + ranked marker TSVs |
| `annotate_seurat_with_SingleR_Eduard()` | SingleR labels against ImmGen, per cluster or per cell |
| `annotate_seurat_with_LLM()` | LLM cluster annotation from markers + PubMed context |

### Differential expression — `references/differential-expression.md`

| function | reach for it when |
|---|---|
| `pseudobulk()` | **default DE.** One profile per sample → DESeq2. Needs ≥2 replicates per group |
| `DEG_FindMarkers_SCT_assay()` | per-cell Wilcoxon on the SCT assay — fallback when there are no replicates |
| `DEG_FindMarkers_RNA_assay()` | the same on the RNA assay |
| `bulk_analysis()` | true bulk RNA-seq from a counts table, not a Seurat object |
| `volcano_plot()` | volcano from a DE results table; always prints and saves |
| `scatterplot()` | log-log group-mean scatter; always prints and saves |
| `plot_deg_counts()` | UP/DOWN bar plot across cell types |
| `plot_gene_lists()` | scatter + volcano restricted to named gene lists |
| `validate_cell_counts()` | check a comparison has enough cells before running DE |
| `add_gene_labels_layer()` | the `geom_text_repel` layer the DE plots use, for custom plots |

### Pathway enrichment — `references/pathway-enrichment.md`

Three backends × five analysis kinds. Pick the row, then the column.

| analysis kind | MSigDB (`pathway_*`) | GO (`GO_*`) | Metascape |
|---|---|---|---|
| over-representation, one gene list | `pathway_overrepresentation_analysis` | `GO_overrepresentation_analysis` | `Metascape_overrepresentation_analysis` |
| over-representation, several lists | `pathway_overrepresentation_analysis_multiple_lists` | `GO_overrepresentation_analysis_multiple_lists` | — |
| GSEA over a ranked table | `pathway_GSEA_analysis` | `GO_GSEA_analysis` | — |
| both directions of a DE result | `pathway_functional_analysis` | `GO_functional_analysis` | `Metascape_functional_analysis` |
| cluster-marker identification | `pathway_functional_analysis_cluster_identification` | `GO_functional_analysis_cluster_identification` | `Metascape_functional_analysis_cluster_identification` |

| also here | |
|---|---|
| `run_DEG_functional_analysis()` | dispatcher — `method = 'ClusterProfiler'` routes to the **MSigDB** column, `'Metascape'` to Metascape. The `GO_*` column is **not reachable** from it |
| `run_overrepresentation_analysis()` | dispatcher for a bare gene list |
| `plot_pathways_heatmap()`, `plot_pathways_heatmap2()` | z-scored expression heatmap of a pathway's genes (genes on x / samples on x) |
| `set_enrichment_color_scale()`, `create_no_data_plot()` | the shared enrichment palette and the empty-result placeholder |

### Pathway redundancy — `references/pathway-redundancy.md` — **needs ≥ 0.11.0**

| function | reach for it when |
|---|---|
| `collapse_pathway_redundancy()` | a few hundred near-duplicate terms → gene-overlap clusters (Jaccard + Louvain) |
| `plot_pathway_clusters()` | bar plot of the collapsed clusters |
| `plot_pathway_cluster_map()` | term-overlap network; `interactive = TRUE` gives a plotly widget |
| `plot_pathway_overlap()` | overlap heatmap — returns a **ComplexHeatmap**, not a ggplot |
| `plot_pathway_topic()` | the terms of one topic, faceted by cluster |
| `plot_pathway_topic_genes()` | which genes drive a topic's terms |
| `pathway_topics()` | the seven built-in vocabularies (`lipid`, `immune`, `stress`, …) |
| `pathway_topic_terms()` | the matching rows, for pulling a topic's full gene set |
| `print.pathway_clusters()` | S3 print method — top clusters at a glance |

### ATAC and regulatory — `references/atac-regulatory.md`

| function | reach for it when |
|---|---|
| `run_differential_accessibility()` | DA peaks between two groups + volcanoes + gene-list CSVs |
| `run_differential_accessibility_FindMarkers()` | the stats core alone, no plots |
| `top_peaks_per_cluster()` | per-cluster marker peaks + ranked TSVs |
| `SE_from_Signac()` | ChromatinAssay → SingleCellExperiment, the entry point for FigR |
| `find_enriched_motifs()` | motif enrichment in a cluster's DA peaks |
| `plot_chromvar_heatmap()` | chromVAR motif activity by group |
| `compare_motive_activity_pairwise()` | motif activity between two clusters + logos + violins |
| `analyze_figr_results()` | plots and tables from a finished FigR run |
| `plot_dorcs_heatmap()` | DORC accessibility across groups |
| `plot_dorc_genes()` | DORC vs RNA side by side for chosen genes |
| `compare_DORC_scores_pairwise()` | DORC scores between two clusters |
| `plotfigRHeatmap_EA()`, `plotfigRHeatmap_EA_gg()` | TF × DORC regulatory-score heatmap (ComplexHeatmap / ggplot) |
| `CoveragePlot_custom()` | coverage + annotation + peak + link tracks for one region |
| `create_tracks_patchwork()`, `create_tracks_patchwork2()` | the same across several regions in a grid |
| `AnnotationPlot_custom()` | the gene-annotation track alone |
| `gene_to_jaspar()` / `ID_to_symbol()` | gene symbol ↔ JASPAR motif ID |

### TCR repertoire — `references/repertoire.md`

| function | reach for it when |
|---|---|
| `plot_clonal_repertoire()` | circle-packing of clonotypes by frequency |
| `calculate_D50()` | D50 diversity per group (**bare** column arguments) |
| `create_V_J_chord_diagram()` | V–J gene-pairing chord diagram |
| `overlap_circos_and_tables()` | clonotypes shared across samples: circos + CSV |
| `create_group_palette()` | consistent sample colours derived from group names |

### Other analyses — `references/other-analyses.md`

| function | reach for it when |
|---|---|
| `run_cellchat_analysis()` | cell–cell communication. **Returns `invisible(NULL)`** — pass `output_file` and `readRDS()` |
| `netVisual_bubble_pathways()` | pathway-level communication bubble plot from a CellChat object |
| `run_cnmf_results()` | load a finished cNMF run onto a Seurat object + program plots |
| `plot_pseudotime_heatmap()` | cells × genes heatmap ordered by pseudotime |
| `apply_integer_breaks()` | integer legend breaks on a `FeaturePlot_scCustom()` output |
| `create_analysis_directory()` | a guaranteed-fresh output directory — **destructive**, see above |

## Exported, but call the wrapper instead

`pseudobulk_de()`, `DEG_FindMarkers_SCT_assay_de()`, `DEG_FindMarkers_RNA_assay_de()`,
`bulk_analysis_de()` and `run_differential_accessibility_FindMarkers()` are the statistics cores
behind their wrappers. They return counts and a results table but skip the plots, the CSVs and the
enrichment. Call them only when you want the numbers and nothing else.

## Documented but not exported

`run_figr_analysis()` has a `man/` page but no `@export`, so `scAnalysisToolkit::run_figr_analysis()`
fails. It is also **broken** — `R/GRN_analysis_FigR.r:300` reads `gr` before it is assigned, and it
uses Seurat-v3 slot writes that fail on a v5 `Assay5` object. Do not call it, even via `:::`.

To run FigR, drive the upstream steps yourself and hand the outputs to the exported
`analyze_figr_results()`:

```r
atac_se <- SE_from_Signac(seurat, assayName = "ATAC")
# then FigR::runGenePeakcorr() -> FigR::dorcJPlot() -> FigR::getDORCScores()
#      -> FigR::smoothScoresNN() -> FigR::runFigRGRN()
analyze_figr_results(figR.d, seurat, dorcMat.s, RNAmat.s, umap.d, top_dorcs_genes)
```

`run_standard_enrichment()` is likewise unexported, and is dead code — it calls
`pathways_of_interest_analysis()`, which does not exist anywhere in the package.

## Where to read next

| file | holds |
|---|---|
| `references/contracts.md` | per-function lookup tables: required metadata columns, `Suggests`, path arguments, output layout, version gates |
| `references/troubleshooting.md` | symptom → cause → fix for the failures these functions actually produce |
| the seven domain files | signature, return shape, disk writes, preconditions and gotchas per function |

## Relationship to `scrnaseq-seurat-analysis`

That skill scaffolds and fills a Quarto analysis end to end, and owns the lab's conventions — folder
layout, `here()` paths, palettes, the checkpoint pattern, chunk granularity, pipeline order and the
demultiplexing choices. This skill owns the calls. If the task is "build me an analysis," start there
and come here for any individual function.
