# scAnalysisToolkit

An R package of helper functions and wrappers for single-cell RNA-seq, single-cell VDJ and
single-nuclei ATAC-seq analysis. It wraps Seurat, Signac, SingleR, DESeq2, clusterProfiler and
CellChat behind consistent interfaces that also write their figures and tables to disk. Operates on
**Seurat objects only** — no scran or scanpy support.

`README.md` (generated from `README.Rmd` — edit the `.Rmd`) has the feature overview. This file
covers the conventions that are not obvious from reading the code.

## Layout

One file per domain under `R/`, 65 exported functions and one S3 method:

| area | files |
|---|---|
| Differential expression | `de_core.r` (stats), `de_wrappers.r` (pipelines), `de_visualization.r` (volcano/scatter/DEG-count plots) |
| Pathway enrichment | `pathway_dispatchers.r` (method dispatch), `ClusterProfiler_MSigDB_functions.r`, `ClusterProfiler_GO_functions.r`, `metascape_functions.r`, `pathway_redundancy.r` (collapses enriched terms into gene-overlap clusters) |
| Clustering & annotation | `seurat_clustering.r`, `cluster_annotation.r`, `llm_annotation.r` |
| ATAC / regulatory | `differential_accessibility.r`, `motif_enrichment.r`, `GRN_analysis_FigR.r` |
| Other analyses | `cellchat_analysis.r`, `cnmf_analysis.r`, `repertoire_analysis.r`, `pseudotime.R` |
| Plotting | `visualization.r`, `heatmaps.r`, `circlize_functions.r` |
| Internal helpers | `shared_utilities.r` |

## Conventions

**`NAMESPACE` and `man/` are generated — never hand-edit either.** Run `devtools::document()` after
touching roxygen comments. Internal helpers are marked `@noRd` and get no `.Rd` file.

**roxygen2 must be 8.0.0**, matching `Config/roxygen2/version` in `DESCRIPTION`. A different version
reformats all 81 `.Rd` files at once and buries the real diff.

**Non-standard evaluation.** The package was migrated off bare column references, and it is easy to
regress:

- `.data$col` inside **data-masking** verbs — `mutate()`, `filter()`, `summarize()`, `arrange()`,
  `count()`, and ggplot `aes()`.
- Quoted strings or `any_of()` / `all_of()` inside **tidy-select** contexts — `select()`, `pull()`,
  `across()`, `relocate()`, `pivot_*()`, `drop_na()`, and the right-hand side of `rename()`.

A `.data$` in a tidy-select context errors at runtime, so the two are not interchangeable.

**Guard `Suggests` dependencies** with `rlang::check_installed()` before first use, so the package
loads without them. The one deliberate exception is `SummarizedExperiment`, which sits in `Imports`
and is used unguarded: it comes in via DESeq2's own `Depends`, so it is always installed.

**`devtools::check()` baseline is 0 errors / 0 warnings / 2 NOTES.** Both notes are known and
accepted, so two notes means passing — don't chase them:

1. `:::` calls into Signac and scRepertoire (no public equivalent), plus `leidenbase` declared in
   `Imports` but only needed as a runtime backend for `FindClusters(algorithm = 4)`.
2. `run_standard_enrichment()` in `shared_utilities.r` calls `pathways_of_interest_analysis()`, which
   does not exist anywhere in the package. That helper is dead code — nothing calls it.

**There is no `tests/` directory.** Verify changes with a scratch script plus `devtools::check()`. For
anything touching output filenames or plot labels, check the generated files, not just the return
values — several functions' primary effect is what they write to disk.

**Output filenames encode the comparison.** Functions in the DE path append a contrast-specific tag
so several comparisons can share one output directory. Note that
`create_analysis_directory()` (`shared_utilities.r`) calls `unlink()` first, so a directory-name
collision deletes the earlier result rather than merging — new writes must stay distinguishable.

## Version bumps and releases

**Never edit `Version:` in `DESCRIPTION` by hand.** A bump is five steps ending in a GitHub release —
a bumped version that is not tagged and released is an unfinished job. Run all of it **on the default
branch**, so merge the feature branch first:

```sh
git switch main && git merge my-feature     # 1. default branch, clean tree
Rscript -e 'usethis::use_version("minor")'  # 2. bump: major | minor | patch
# 3. write the notes under the NEWS.md heading use_version() just added
git commit -am "Bump to 0.12.0"             # 4. commit the bump and the notes together
git tag -a v0.12.0 -m "scAnalysisToolkit 0.12.0"
git push --follow-tags                      # 5. tag, push, then release
Rscript -e 'usethis::use_github_release()'
```

Each of these has bitten us at least once:

- **`use_version()` needs `which`.** With no argument it asks interactively; under `Rscript` it
  silently returns and the version is left unchanged, with no error. Pass `"major"`, `"minor"` or
  `"patch"`.
- **`use_version()` needs a clean tree**, or it aborts with `User input required, but session is not
  interactive`. That is why the bump is step 2, not step 1.
- **Tag with `-a`.** Every tag here is annotated. Push it *before* `use_github_release()`: that
  function pushes the branch but not tags, and if the tag is missing GitHub creates its own
  lightweight one instead.
- **The tag must carry both the bumped `DESCRIPTION` and that version's `NEWS.md` section.** Step 4
  puts them in one commit; a separately committed bump splits them, and tagging that commit
  publishes an empty release. Tag the tip, or the merge commit — never an intermediate bump. Two
  commands to check:

  ```sh
  git show <target>:DESCRIPTION | grep '^Version:'
  git show <target>:NEWS.md | head -1        # heading must be that same version
  ```

- **`use_github_release()` will not run from a detached `HEAD`** (`Detached head; can't continue`),
  so check out the branch, not the SHA. It takes the version from `DESCRIPTION` at `HEAD` to build
  the tag name, pushes the current branch, and publishes immediately with no confirmation prompt. If
  the tag already exists on GitHub, the API ignores the target commit and attaches the release to
  that tag.

**`NEWS.md` is the release body** — `use_github_release()` publishes whatever sits under the topmost
heading, so an empty section ships an empty release.

Keep it **short and strictly user-facing**: what someone calling these functions has to know. A new
or removed function or argument, a changed default, a behaviour change they would notice, a bug that
was biting them, a new dependency they must install. One line each, no subheadings unless the
release is genuinely large. Leave out refactors, internal helpers, NSE cleanups, roxygen and comment
edits, test and tooling changes, and anything only visible reading the diff — that is what the git
history is for. If a bullet needs a paragraph to explain, it probably belongs in the function's
documentation instead.

Before tagging, confirm the notes are complete: diff exported signatures, argument defaults and
`DESCRIPTION` against the previous tag, and account for every change as either noted or new in this
version. Anything new ships under its own entry and needs no separate note — only changes to
already-released behaviour do.

## Code style

- Prefer dplyr/tidyverse over base R.
- Use the native pipe `|>` for chaining.
- Prefer `mutate()`, `filter()`, `select()`, `case_when()`, `arrange()`, `summarize()` over base
  equivalents.
- Use purrr's `map()` family (`map()`, `map_chr()`, `map_dbl()`, `walk()`, …) instead of `for` loops,
  `lapply()`, `sapply()`, or the `apply()` family.
- Follow the tidyverse style guide.
- Prefer `usethis` and `devtools` for package tasks — e.g. `usethis::use_version()`
  (see [Version bumps and releases](#version-bumps-and-releases)), `devtools::document()`,
  `devtools::check()`.
