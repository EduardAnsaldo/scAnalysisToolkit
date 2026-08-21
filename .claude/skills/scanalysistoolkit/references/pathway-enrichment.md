# Pathway enrichment

Getting an enrichment result. For what to do with one you already have — collapsing redundant terms,
topics, the cluster plots — see `pathway-redundancy.md`.

All arguments are **quoted strings**; nothing here uses `rlang` embrace.

---

## Three backends

| | engine | gene sets | reachable from a dispatcher |
|---|---|---|---|
| `pathway_*` | `clusterProfiler::enricher` / `GSEA` | MSigDB **mouse**: MH, M2:CP:BIOCARTA, M2:CP:REACTOME, M2:CP:WIKIPATHWAYS, M5:GO:BP | yes — `method = 'ClusterProfiler'` |
| `GO_*` | `clusterProfiler::enrichGO` / `gseGO` | `org.Mm.eg.db`, `keyType = 'SYMBOL'` | **no** |
| `Metascape_*` | the external `msbio` CLI | Metascape's own | yes — `method = 'Metascape'` |

Two differences that bite:

- **`pathway_*` converts symbols to Entrez internally** (`bitr()`), then back with `setReadable()`.
  **`GO_*` takes gene symbols directly** — do not pre-convert for it.
- **`GO_*` is unreachable from `run_DEG_functional_analysis()` / `run_overrepresentation_analysis()`.**
  `method = 'ClusterProfiler'` routes to the *MSigDB* family. To use GO, call it directly.

⚠️ **Metascape is hard-coded to human.** It shells out with `--source_tax_id 9606 --target_tax_id 9606`
while everything else in the package is mouse. It also needs a local msbio v3-2 install at
`~/msbio_v3-2/`, writes gene lists into `~/msbio_v3-2/data/…`, copies results back, then deletes that
directory. `system()` failures are not checked.

---

## Dispatchers

```r
run_DEG_functional_analysis(results, method = NULL, ...)
run_overrepresentation_analysis(genes, method = NULL, ...)
```

`method = NULL` (the default everywhere) returns `invisible(NULL)` immediately — no enrichment. Valid
values are `'ClusterProfiler'`, `'Metascape'`, or both; `match.arg` partial-matches, and more than two
is an error. **`TRUE` is not valid** and errors in `match.arg`.

One method returns that method's own value unwrapped; two returns a list named by method.

- `run_DEG_functional_analysis` → `pathway_functional_analysis()` or `Metascape_functional_analysis()`.
- `run_overrepresentation_analysis` → `pathway_overrepresentation_analysis()` (with `genes` binding
  positionally to `significant_genes`, so you must supply `all_genes` and `local_path` through `...`)
  or `Metascape_overrepresentation_analysis()`. Only the Metascape branch prints.

---

## Over-representation, one gene list

```r
pathway_overrepresentation_analysis(significant_genes, all_genes, local_path,
    minGSSize = 4, maxGSSize = 400, filename = '', group = '', nterms_to_plot = 12,
    font_size = 10, label_length = 60, wrap_labels = FALSE, collapse_redundancy = FALSE,
    topics = NULL, redundancy_padj = 0.05, redundancy_min_sim = 0.25,
    redundancy_resolution = 1, redundancy_seed = 1, redundancy_label_by = 'collection',
    redundancy_n_clusters = 15, redundancy_max_clusters = 12, ...)

GO_overrepresentation_analysis(significant_genes, all_genes, local_path, ontology = 'ALL',
    minGSSize = 5, maxGSSize = 400, filename = '', group = '', drop_levels = FALSE,
    levels_to_drop = c(), simplify_function = min, simplify_by = 'p.adjust',
    simplify_terms = TRUE, run_network = FALSE, network_n_terms = 100,
    nterms_to_plot = 50, font_size = 8, ...)

Metascape_overrepresentation_analysis(significant_genes_FC_ordered, local_path, group,
    grouping_var, nterms_to_plot_metascape = 12, filename = '', ...)
```

Note `minGSSize` is **4** for the MSigDB entry point and **5** for the others.

**`pathway_overrepresentation_analysis()` returns a bare ggplot carrying three attributes** —
`enrichment_results` (a `setReadable`-converted `enrichResult`, or `NULL`), `pathway_clusters`, and
`pathway_redundancy_plots`. ⚠️ **Read the attributes before adding any ggplot layer**; `p + theme(...)`
drops them.

```r
p <- pathway_overrepresentation_analysis(sig, universe, local_path = out)
res <- attr(p, 'enrichment_results')
```

It returns a placeholder with all three attributes `NULL`, and writes no CSV, when
`length(significant_genes) <= 4` or when ≤ 4 symbols survive Entrez mapping. It never prints.

`GO_overrepresentation_analysis()` returns nothing useful and **always prints** its dotplot.
`Metascape_overrepresentation_analysis()` returns a ggplot and never prints.

**Files.** All into `local_path`, which must already exist.

| function | writes |
|---|---|
| `pathway_overrepresentation_analysis` | `<filename>Pathway_OverRepresentation_analysis_results.csv`; `<filename>_Pathway_enrichment_analysis_metascape.pdf` (6×6) only when `nrow > 1` — the "metascape" in that name is a copy-paste artefact, this is the MSigDB bar plot |
| `GO_overrepresentation_analysis` | `<filename>GO_OverRepresentation_analysis_results_<ontology>.csv`; plus `…_filtered_…csv` and `<filename>GO overrepresentation_analysis_dotplot_<ontology>.pdf` (10×18) when `nrow > 1` — **note the literal space** in that filename |
| `Metascape_overrepresentation_analysis` | `<filename>__Pathway_enrichment_analysis_metascape.pdf` (6×6, double underscore with the default `filename = ''`), plus everything msbio copies back |

With `collapse_redundancy = TRUE` and `nrow > 1`, `pathway_overrepresentation_analysis` also writes
`<filename>Pathway_redundancy_clusters_barplot.pdf` (8×6),
`…_cluster_network.pdf` (9×7), `…_terms.csv`, `…_terms_summary.csv`, and per topic
`…_topic_<slug>_terms.pdf` (8×7) + `…_genes.pdf` (6×6).

**Gotchas.** `qvalueCutoff` is fixed at 0.25 in the MSigDB path. The bar-plot fill scale is fixed at
`limits = c(0, 20)`. Unmapped symbols are dropped with a `message()`. The redundancy plot title is
always `"<filename>Pathways UP in <group>"` regardless of direction.
`GO_overrepresentation_analysis(run_network = TRUE, simplify_terms = FALSE)` references an object that
is never assigned — masked by `try()`, so it degrades to a printed `try-error`.

---

## Over-representation, several lists

```r
pathway_overrepresentation_analysis_multiple_lists(gene_list, all_genes, local_path,
    minGSSize = 5, maxGSSize = 400, filename = '', nterms_to_plot = 5, font_size = 8,
    grouping_var = '', wrap_length = 60, ...)

GO_overrepresentation_analysis_multiple_lists(gene_list, all_genes, local_path,
    ontology = 'ALL', minGSSize = 5, maxGSSize = 400, filename = '', drop_levels = FALSE,
    levels_to_drop = c(), simplify_function = min, simplify_by = 'p.adjust',
    simplify_terms = FALSE, run_network = FALSE, network_n_terms = 100,
    nterms_to_plot = 5, font_size = 8, grouping_var = '', ...)
```

`gene_list` is a **named list** of symbol vectors. Both return `list(plot, results)` and write a
`…_gene_lists_results…csv` plus a dotplot at 10×18 when `nrow > 1`.

⚠️ **Prefer the MSigDB version.** The GO one has two failure modes with no fallback: with ≤ 1
enriched term `p1` is never assigned and it errors with `object 'p1' not found`, and there is no
`NULL` guard on an empty `compareCluster` result. The MSigDB twin returns a no-data placeholder in
both cases (and prints it).

Note `simplify_terms` defaults to `TRUE` in `GO_overrepresentation_analysis` but `FALSE` in its
multiple-lists variant.

---

## GSEA over a ranked table

```r
pathway_GSEA_analysis(results, local_path, group, minGSSize = 5, maxGSSize = 500,
                      pvalueCutoff = 0.05, ...)
GO_GSEA_analysis(results, local_path, ontology = 'ALL', group, minGSSize = 5,
                 maxGSSize = 500, pvalueCutoff = 0.05, ...)
```

`results` needs `genes`, `log2FoldChange`, `padj`. Both return nothing useful and **print** their
dotplot when `nrow > 1`.

⚠️ In `GO_GSEA_analysis()`, `group` is required but sits **after** `ontology`, so a positional third
argument silently binds to `ontology`. Always name it.

⚠️ Neither CSV filename carries the `filename`/`group` prefix, so running both directions into one
`local_path` overwrites.

---

## Both directions of a DE result

```r
pathway_functional_analysis(results, grouping_var, path = './', FC_threshold = 0.3,
    p_value_threshold = 0.05, group1 = '', group2 = '', run_GSEA = FALSE,
    top_n_genes = NULL, file_suffix = '', collapse_redundancy = FALSE, ...)

GO_functional_analysis(results, grouping_var, path = './', FC_threshold = 0.3,
    p_value_threshold = 0.05, group1 = '', group2 = '', run_GSEA = FALSE,
    top_n_genes = NULL, ...)

Metascape_functional_analysis(results, grouping_var, group2, group1, path = './',
    FC_threshold, p_value_threshold = 0.05, top_n_genes = NULL, file_suffix = '', ...)
```

These split `results` into UP and DOWN by `FC_threshold` and run over-representation on each.

⚠️ **`grouping_var` here is not a column name.** It is free text pasted into the output directory
name — normally the cell type or cluster label.

⚠️ In `Metascape_functional_analysis()`, `FC_threshold` is **required despite sitting after an
argument with a default**, as are `group2` and `group1`.

**`pathway_functional_analysis()` returns an 8-element list:** `plot_up`, `plot_down`, `results_up`,
`results_down`, `redundancy_up`, `redundancy_down`, `redundancy_plots_up`, `redundancy_plots_down`.
The four redundancy elements are `NULL` unless `collapse_redundancy = TRUE`. It **prints `plot_up` and
`plot_down` only when `collapse_redundancy = FALSE`** — in redundancy mode nothing is printed at all.

`GO_functional_analysis()` returns `NULL`, always prints, and has **no `file_suffix` and no
`collapse_redundancy`**. `Metascape_functional_analysis()` prints both directions unconditionally and
has no useful return.

⚠️ **Directories are recreated destructively** via `create_analysis_directory()`, which calls
`unlink()` first:

| function | directories |
|---|---|
| `pathway_functional_analysis` | `<path>/Pathway_functional_analysis_UP_in_<group2>_<grouping_var><file_suffix>` + `group1` twin |
| `GO_functional_analysis` | `<path>/GO_functional_analysis_UP_in_<group2>_<grouping_var>` + twin — no suffix available |
| `Metascape_functional_analysis` | `<path>/Functional_analysis_metascape_UP_in_<g2>_<gv><file_suffix>` and `…_DOWN_in<g1>_<gv><file_suffix>` (missing `_` after `DOWN_in`) |

`topics` is not a formal of `pathway_functional_analysis()` — it rides through `...` into
`pathway_overrepresentation_analysis()`. With `run_GSEA = TRUE`, both directions are handed the whole
`results`, so the two GSEA CSVs are identical apart from their directory.

---

## Cluster-marker identification

```r
pathway_functional_analysis_cluster_identification(scRNAseq, results, path = './',
    object_annotations = '', top_gene_number = 50, per_cluster = FALSE, ...)

GO_functional_analysis_cluster_identification(scRNAseq, results, path = './',
    object_annotations = '', top_gene_number = 50, ...)

Metascape_functional_analysis_cluster_identification(seurat, results,
    identities = 'seurat_clusters', path = './', object_annotations = '', ...)
```

`results` is `FindAllMarkers()` output — needs `cluster`, `gene`, `avg_log2FC`, `p_val_adj`.
`scRNAseq[['RNA']]` supplies the universe via `Features()`.

The MSigDB version returns a ggplot, or a patchwork with `per_cluster = TRUE` (only that branch
prints). The GO version returns the dotplot and does not print, and has no `per_cluster`. The
Metascape version returns a patchwork and prints it.

⚠️ The Metascape version coerces `results$cluster` with `forcats::fct_inseq()`, so **cluster labels
must be numeric-like** — `"CD8 T"` errors. Its `seurat` and `identities` arguments are accepted but
never used, and every cluster incurs a full msbio round trip.

---

## Pathway gene-expression heatmaps

```r
plot_pathways_heatmap (genes_to_plot, seurat, title,
    color_palette = grDevices::hcl.colors(n = 20, 'RdBu', rev = TRUE), grouping_var = 'Samples')
plot_pathways_heatmap2(genes_to_plot, seurat, title,
    color_palette = grDevices::hcl.colors(n = 20, 'RdBu', rev = TRUE), grouping_var = 'Samples')
```

Z-scored expression of a gene set across groups. Both return a **ggplot** — despite the roxygen title
saying ComplexHeatmap, neither touches it. No disk output, no printing.

The only differences: v1 puts **genes on x** with italic labels; v2 puts **samples on x**. v1's legend
title contains a literal `<br>` that only renders under `ggtext::element_markdown()` — with the plain
theme used here it prints as `z-score<br>expression`. **Prefer v2** unless you want genes on x.

Needs an `RNA` assay (`AggregateExpression`). Genes absent from the object are silently dropped. Gene
order follows your input order — no clustering. Z-scoring is per gene across groups, so a
zero-variance gene renders as an empty row, and with only one group level everything is `NaN`.

---

## Shared helpers

```r
set_enrichment_color_scale()          # sets options(enrichplot.colours = viridis(4, direction = -1))
create_no_data_plot(label = 'N/A', text_size = 24)   # the placeholder returned on empty results
```

`create_no_data_plot()` is what every enrichment function returns when nothing enriched. Note that a
placeholder is **never saved to disk** anywhere in the package — the one exception is
`scatterplot(scatter_applicable = FALSE)`, which deliberately saves it so a figure set stays complete.
