# Pathway redundancy

**Requires ≥ 0.11.0.** On an older build these eight functions do not exist, and
`collapse_redundancy` / `topics` on the DE wrappers are silently swallowed by `...` with no error.
Check first:

```r
'collapse_redundancy' %in% names(formals(pseudobulk))
```

An over-representation result of a few hundred significant terms is mostly near-duplicates — the same
genes reported under `GOBP_`, `REACTOME_`, `WP_` and `HALLMARK_` names. `clusterProfiler::simplify()`
cannot help, because it only works within a single GO ontology. This family clusters terms by the
**Jaccard similarity of their enriched genes** and runs Louvain community detection over that graph,
which works across mixed collections.

All arguments are quoted strings or plain values.

---

## Getting there

Usually you do not call `collapse_pathway_redundancy()` directly — you pass `collapse_redundancy = TRUE`
to a DE wrapper and take the results off `$enrichment`:

```r
de <- pseudobulk(seurat, comparison = 'Samples', group1 = 'WT', group2 = 'KO',
                 cluster = 'Tcells', path = out,
                 run_pathway_enrichment = 'ClusterProfiler',
                 collapse_redundancy = TRUE, topics = c('lipid', 'immune'))
```

⚠️ **In that mode nothing is printed** — not even the plain over-representation bar plots. Every
figure comes back for you to place:

```r
de$enrichment$plot_up                            # over-representation bar plot
de$enrichment$plot_down
de$enrichment$results_up                         # enrichResult, or NULL
de$enrichment$redundancy_up                      # pathway_clusters object, or NULL
de$enrichment$redundancy_plots_up$clusters       # ggplot
de$enrichment$redundancy_plots_up$network        # ggplot, or plotly widget
de$enrichment$redundancy_plots_up$topics$lipid$terms
de$enrichment$redundancy_plots_up$topics$lipid$genes
```

A topic that matched nothing in a direction is simply absent from `$topics` (with a warning), and no
file is written for it. `collapse_redundancy` and `topics` do nothing under `'Metascape'`.

---

## `collapse_pathway_redundancy()`

```r
collapse_pathway_redundancy(x, padj = 0.05, min_sim = 0.25, resolution = 1, seed = 1,
    name = NULL, label_by = c('collection', 'richfactor', 'genes', 'padj'))
```

`x` may be an `enrichResult` / `gseaResult`, a data frame, or a **length-1 path to an exported CSV**.
It needs `ID`, `p.adjust` and `geneID`; `Description` and `Count` are derived if absent, and leading
`...N` index columns from `write.csv` are stripped automatically.

⚠️ **`geneID` must hold `/`-separated gene symbols** — i.e. the result must have been through
`setReadable()`. Raw Entrez IDs produce meaningless clusters.

`label_by` (`match.arg`, default `'collection'`) chooses how each cluster is named in `cluster_name`:
`'collection'` takes the member from the most curated collection, `'richfactor'` the highest rich
factor, `'genes'` the largest gene count, `'padj'` the most significant. `'collection'` and
`'richfactor'` warn and fall back to `p.adjust` if the needed column is missing.

**Returns a `pathway_clusters` object** (class `c('pathway_clusters','list')`) with four elements:

| element | holds |
|---|---|
| `terms` | one row per surviving term: the original enrichment columns plus `collection`, `label`, `genes` (list-column), `Count`, `cluster`, `cluster_label`, `cluster_name`, `rank` |
| `clusters` | one row per cluster: `cluster`, `n_terms`, `best_padj`, `cluster_label` (always the most significant member), `n_genes`, `collections`, `cluster_name`, `cluster_name_term`, `rank` |
| `sim` | the square Jaccard matrix, in `terms$ID` order |
| `params` | `padj`, `min_sim`, `resolution`, `seed`, `name`, `label_by` |

**Returns `NULL` with a message** when fewer than two terms pass `p.adjust < padj`. Handle that before
plotting — though the plot functions themselves return placeholders rather than failing.

No disk writes, no printing. Needs **`igraph`**.

`print.pathway_clusters()` gives a one-line summary plus the top 5 clusters.

---

## Plots

```r
plot_pathway_clusters(res, n = 15, title = res$params$name)
```
Bar plot of the collapsed clusters. Returns a ggplot.

```r
plot_pathway_cluster_map(res, max_clusters = 12, title = res$params$name, seed = 1,
                         label_length = 40, interactive = FALSE, hover_gene_width = 60)
```
Term-overlap network. `interactive = TRUE` returns a **plotly htmlwidget** whose tooltips carry the
term name, collection, adjusted p-value and comma-separated genes wrapped at `hover_gene_width`.

⚠️ `ggsave()` cannot take the widget, and in Quarto/knitr a widget only renders when it is the
top-level value of its own chunk. Needs `igraph`, plus `plotly` when interactive.

```r
plot_pathway_overlap(res, title = res$params$name, min_label = 3)
```
⚠️ Returns a **`ComplexHeatmap::Heatmap`**, not a ggplot — `ggsave()` will not work. Use
`ComplexHeatmap::draw()` into a device. Needs `ComplexHeatmap`.

```r
plot_pathway_topic(res, pattern, where = 'label', n = 30, title = NULL)
```
Lollipop plot of one topic's terms, faceted by cluster.

```r
plot_pathway_topic_genes(res, pattern, where = 'label', top_genes = 25,
                         lfc = NULL, padj = NULL, title = NULL)
```
How many of a topic's terms each gene appears in. `lfc` and `padj` are **named numeric vectors keyed
by gene symbol**: `padj` colours bars by `-log10(padj)` (sequential), `lfc` by direction (diverging),
and `padj` wins if both are given. A grey bar means that symbol was not found in the vector — worth
checking, since enrichment round-trips symbols through Entrez. `top_genes` truncates *after* sorting,
so ties at the cutoff are arbitrary; raise it or pass `Inf`.

All of these return placeholders rather than erroring on an empty or `NULL` `res`, and none print or
save.

---

## Topics

```r
pathway_topics()
pathway_topic_terms(res, pattern, where = c('label', 'genes', 'both'))
```

`pathway_topics()` takes no arguments and returns the seven built-in vocabularies: **`lipid`,
`mitochondria`, `translation`, `cell_cycle`, `chromatin`, `immune`, `stress`**. Those names are the
valid shorthand anywhere a `pattern` or `topics` argument is accepted.

`pathway_topic_terms()` returns the matching rows of `res$terms` ordered by `p.adjust`, with
attributes `rx` (the resolved regex) and `n_total`. Use it to pull a topic's full gene set. It returns
`NULL` on an empty `res`, and a 0-row tibble plus a warning when nothing matches. Matching is
case-insensitive; gene matching is anchored `^(...)$`.

### The `topics` argument format

Accepted by `pseudobulk()`, the other three DE wrappers, and `pathway_overrepresentation_analysis()`.
**Only has an effect when `collapse_redundancy = TRUE`.**

- A **character vector** — each element is its own topic, and each is either a `pathway_topics()` name
  or a case-insensitive regex: `topics = c('lipid', 'immune')`.
- A **named list** to OR several regexes into one named topic:
  `topics = list(cytokines = c('interleukin', 'chemokine'))`. The name becomes the filename slug.

⚠️ A `pathway_topics()` name is only resolved when `length(pattern) == 1`. Passing
`c('lipid','immune')` to `pathway_topic_terms()` OR-s them into the literal regex `lipid|immune`
rather than unioning the two vocabularies.

---

## Worked example: placing the plots yourself

The one pattern worth spelling out, because `collapse_redundancy = TRUE` inverts the usual
print-as-side-effect behaviour and a mental model built on the other wrappers will be wrong.

```r
de <- pseudobulk(seurat_subset, comparison = 'Samples', group1 = 'WT', group2 = 'KO',
                 cluster = 'Tcells', path = here('results/DEG_KO_vs_WT'),
                 run_pathway_enrichment = 'ClusterProfiler',
                 collapse_redundancy = TRUE, topics = 'lipid')

# nothing has been printed at this point beyond the scatter and volcano

up <- de$enrichment$redundancy_up          # NULL if fewer than 2 terms survived
if (!is.null(up)) {
  print(up)                                # S3 summary: n clusters, top 5
  print(de$enrichment$plot_up)             # plain over-representation bar plot
  print(de$enrichment$redundancy_plots_up$clusters)
  print(de$enrichment$redundancy_plots_up$topics$lipid$genes)

  # the network is a plotly widget when interactive - give it its own chunk
  de$enrichment$redundancy_plots_up$network
}
```

To go back to the underlying terms of a topic:

```r
lipid_terms <- pathway_topic_terms(up, 'lipid', where = 'label')
unique(unlist(lipid_terms$genes))
```
