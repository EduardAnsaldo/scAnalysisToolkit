# TCR repertoire

All of these expect an object already processed with `scRepertoire::combineTCR()` +
`combineExpression()`, so the metadata carries `CTaa`, `CTgene`, `clonalFrequency` and `cloneSize`.

**Shared contract for this file.** `calculate_D50()` is one of only two functions in the package
taking **bare** column names; everything else here takes quoted strings. None of these create their
output directory. Three of the four need `scRepertoire` installed, and two reach into
`scRepertoire:::.expression2List()` — an unexported internal that a version upgrade can break.

---

## `plot_clonal_repertoire()`

```r
plot_clonal_repertoire(seurat, group_name, figures_path,
    viridis_option = 'D', label_size_range = c(3, 10))
```

Circle-packing of clonotypes sized by frequency. Returns a ggplot.

⚠️ **`group_name` is a *value* of the `Groups` column**, not a column name. (Contrast
`create_V_J_chord_diagram()` below, where the same argument name means a column.) Passing a value not
present in `Groups` yields an empty tibble and a downstream error rather than a clean message.

**Requires** metadata `Groups`, `CTaa`, `clonalFrequency`, `cloneSize`. Needs `packcircles`.

Writes `clonal_repertoire_<group_name>.pdf` (8×8) into `figures_path`, which must exist. Does not
print the plot, but does print a stray character vector of clone sizes.

Call once per group:

```r
purrr::walk(levels(seurat$Groups),
            ~ plot_clonal_repertoire(seurat, group_name = .x, figures_path = figures_path))
```

The `cloneSize` colour map is hard-limited to five levels, and the legend labels resolve to `NULL`
because the column they reference is dropped earlier — a cosmetic bug, not a failure.

---

## `calculate_D50()`

```r
calculate_D50(seurat, cell_grouping_var, replicate_var, replicate_group_var = NULL,
              results_path, figures_path, tables_path, colors = NULL)
```

D50 diversity index — the fraction of clonotypes accounting for half the cells — per group.

⚠️ **`cell_grouping_var`, `replicate_var` and `replicate_group_var` are embraced. Pass bare, unquoted
column names.** A quoted string is mangled by `englue()` into a literal `"\"colname\""`.

⚠️ **`replicate_group_var = NULL` does not fall back to `replicate_var`, despite the documentation.**
`englue()` on the `NULL` default yields the literal string `"NULL"`, so the guard passes and the plot
is built with `aes(x = NULL, color = NULL)`. **Always pass it explicitly** — repeat `replicate_var` if
you want per-replicate.

⚠️ **`cell_grouping_var` and `replicate_var` must be factors.** The function reads `levels()`; on
character columns that is `NULL` and the loops silently produce nothing.

```r
calculate_D50(seurat,
              cell_grouping_var  = major_cell_types,
              replicate_var      = Samples,
              replicate_group_var = Groups,
              results_path = results_path, figures_path = figures_path,
              tables_path = tables_path)
```

Returns the ggplot only — the wide results table is written to disk and `head()`-printed, never
returned. Prints both.

**Writes** `<tables_path>/D50_per_<cell_grouping_var>.csv` (glue-interpolated verbatim, no
`here::here()`) and `D50_per_<cell_grouping_var>.pdf` into `figures_path`.

⚠️ `results_path` is a **required positional argument that is never used** — you must still supply it.

**Other gotchas.** Groups with fewer than 20 cells silently become `NA` and drop out. The y scale is
hard-coded to `c(0, 0.5)`, so higher values vanish with a ggplot warning. The `All` row is the *first*
`orig.ident` split, not the pooled repertoire, despite the label. `colors` must be named to match the
`replicate_group_var` levels. Needs `scRepertoire`.

---

## `create_V_J_chord_diagram()`

```r
create_V_J_chord_diagram(trbv_j_table, figures_path, chain_a = 'TRBV', chain_b = 'TRBJ',
    group_name, trbv_palette = 'Teal', trbj_palette = 'Emrl', link_palette = 'YlOrRd',
    cex = 0.7, major_ticks = 10)
```

⚠️ **Not a Seurat function.** `trbv_j_table` is a matrix or data frame whose **rownames are
`"<V>_<J>"`** (split on a single `_`, so gene names containing extra underscores break it) and which
has **one count column per group**.

⚠️ **`group_name` here is a *column name* of `trbv_j_table`** — the opposite of
`plot_clonal_repertoire()`, where it is a metadata value. It names both the column read and the output
file.

Returns a `recordedplot`. Writes `<chain_a>_<chain_b>_chord_diagram_<group_name>.pdf` via
`here::here()` (project-root relative).

⚠️ It **draws twice** — once into the PDF device and once into your active device to build the
`recordPlot()` — so it always paints the current graphics device.

Call once per group:

```r
purrr::walk(levels(seurat$Groups),
            ~ create_V_J_chord_diagram(trbv_j_table, figures_path = figures_path,
                                       group_name = .x))
```

**Gotchas.** `major_ticks` is documented but never used. `check_installed("scRepertoire")` fires even
though the function never uses `scRepertoire` — a false dependency you must satisfy anyway. Link
transparency divides by the range of interaction values, so a table where all values are equal yields
`NaN`. `circos.clear()` runs on exit but not on entry, so leftover state from a previous call can
corrupt the first plot.

---

## `overlap_circos_and_tables()`

```r
overlap_circos_and_tables(seurat, grouping_variable = 'Samples', results_path = './',
    variables_to_color_by = NULL, cell_types_column = 'cell_types', write_table = TRUE,
    figures_path, tables_path, circle_margin = 0.2, major_ticks = 500, cex = 0.6,
    samples, groups, color_palette = NULL, alpha_col = 0.8, color_across_groups = NULL,
    base_colors = c('skyblue','darkorange','wheat','darkolivegreen','orchid',
                    'seagreen','gold','lavenderblush','tomato'))
```

Circos of clonotypes **shared** across `grouping_variable`, plus a clonotype-by-sample CSV.

⚠️ `figures_path`, `tables_path`, `samples` and `groups` have **no defaults but sit after arguments
that do** — always call by name.

Returns **`NULL`**. Everything is on disk.

**Requires** `orig.ident`, `CTaa`, `CTgene`, plus the columns named by `grouping_variable` and
`cell_types_column`. Pass `cell_types_column = 'major_cell_types'` if that is your column — the
default `'cell_types'` rarely exists.

⚠️ **`grouping_variable` must be a factor** — the code reads `levels()`, and a character column makes
the reshape fail.

⚠️ **`figures_path` must end in a trailing `/`.** The circos filename is built with raw `paste0()`, so
without it the file lands beside the directory with a mangled name. Writes
`Circos_clonotypes_per_<grouping_variable>.pdf` **and** `.png`, plus
`<tables_path>/TCR_data_per_group.csv` via `here::here()` — a fixed filename, so repeated runs
overwrite.

⚠️ `results_path` is accepted and never used.

**Palette gotcha.** With `color_palette = NULL`, groups are mapped to samples by **substring
matching** (`samples[grepl(group, samples)]`), so each group label must be a literal substring of its
sample names; a group matching nothing yields `NA` colours. `samples` and `groups` are only consulted
in that case but are required regardless. Limits: 9 groups, 5 samples per group, before it starts
scavenging arbitrary colours.

Prints one line of circos metadata per sector — noisy but harmless. Uses
`scRepertoire:::.expression2List()`, splitting by `orig.ident` and taking the first element, which is
the whole object in the usual single-run case.

---

## `create_group_palette()`

```r
create_group_palette(samples, groups,
    base_colors = c('skyblue','darkorange','wheat','darkolivegreen','orchid',
                    'seagreen','gold','lavenderblush','tomato'))
```

Returns a **named character vector** mapping each sample to a colour, with shades derived from its
group's base colour. This is the palette `overlap_circos_and_tables()` builds internally; call it
directly to reuse the same colours across other plots, or to inspect the mapping before committing to
it.

Same substring-matching caveat: group names must appear literally inside their sample names.
