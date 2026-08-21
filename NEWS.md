# scAnalysisToolkit 0.12.0

* Twelve helpers that had documentation but no `@export` are now callable:
  `gene_to_jaspar()`, `ID_to_symbol()`, `plot_dorc_genes()`, `create_tracks_patchwork()`,
  `create_tracks_patchwork2()`, `create_group_palette()`, `create_analysis_directory()`,
  `validate_cell_counts()`, `plot_gene_lists()`, `add_gene_labels_layer()`,
  `create_no_data_plot()` and `set_enrichment_color_scale()`.

* `create_analysis_directory()` is destructive and its help page now says so: it calls
  `unlink()` before creating the directory, so an existing directory of the same name is
  deleted rather than merged into.

* `run_figr_analysis()` stays unexported — it reads an unassigned variable, and its
  Seurat-v3 slot writes fail on v5 objects. Run FigR's own steps and pass the result to
  `analyze_figr_results()`.

* New agent skill under `.claude/skills/scanalysistoolkit/` covering the calling contract for
  all 77 exports: signatures, return shapes, required Seurat metadata columns, and the
  bare-vs-quoted column-argument rule. It ships with the repository, not the installed package.

# scAnalysisToolkit 0.11.0

* New `collapse_pathway_redundancy()` groups an over-representation result's significant
  terms into non-redundant clusters by gene overlap, which works across the mixed MSigDB
  collections where `clusterProfiler::simplify()` does not.

* New plots for the collapsed result: `plot_pathway_clusters()`,
  `plot_pathway_cluster_map()` and `plot_pathway_overlap()`. Pass `interactive = TRUE` to
  `plot_pathway_cluster_map()` for a plotly network with per-node hover.

* `collapse_pathway_redundancy()` gains `label_by`, controlling how each cluster is named in
  `cluster_name`: `"collection"` (default) takes the member from the most curated collection,
  `"richfactor"` the most specific one, and `"genes"` / `"padj"` the largest and most
  significant. All three plots show `cluster_name`; the tables also keep `cluster_label`, always
  the most significant member. Set it from the DE wrappers with `redundancy_label_by`.

* `plot_pathway_topic()` is now always faceted by cluster. It previously dropped the facet when
  every matching term fell in one cluster, which hid which cluster they were.

* New topic helpers, for whether a theme is one redundant block or several distinct signals:
  `pathway_topic_terms()`, `plot_pathway_topic()`, `plot_pathway_topic_genes()` and
  `pathway_topics()`. A topic is a `pathway_topics()` name such as `"lipid"` or any
  case-insensitive regex.

* `pseudobulk()`, `DEG_FindMarkers_RNA_assay()`, `DEG_FindMarkers_SCT_assay()` and
  `bulk_analysis()` gain `collapse_redundancy` and `topics`.

* With `collapse_redundancy = TRUE` the enrichment path prints nothing, including the
  over-representation bar plots, and returns every plot instead — place them yourself from
  the new `redundancy_plots_up` / `redundancy_plots_down` elements. Figures and tables are
  still written to disk.

* `pathway_overrepresentation_analysis()` lowers its `minGSSize` default from 5 to 4, so smaller
  gene sets are tested. `pathway_GSEA_analysis()` and
  `pathway_overrepresentation_analysis_multiple_lists()` stay at 5.

* `bulk_analysis()` no longer errors on its own defaults: `run_pathway_enrichment` defaulted
  to `FALSE`, which the dispatcher rejected after the DE and plotting had already run. It is
  now `NULL`.

* `DEG_FindMarkers_RNA_assay()` now forwards `...` to pathway enrichment, so its options are
  reachable.

* `igraph` and `plotly` added to `Suggests`; install them to use the new functions.

----

Releases before 0.11.0 predate this file; see the git tags and history.
