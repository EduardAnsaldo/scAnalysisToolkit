# scAnalysisToolkit

An R package of helper functions and wrappers for single-cell RNA-seq and ATAC-seq analysis, for
internal use by the DeCarvalho and Constantinides labs at Scripps Research. It wraps Seurat, Signac,
SingleR, DESeq2, clusterProfiler and CellChat behind consistent interfaces that also write their
figures and tables to disk. Operates on **Seurat objects only** — no scran or scanpy support.

`README.md` (generated from `README.Rmd` — edit the `.Rmd`) has the feature overview. This file
covers the conventions that are not obvious from reading the code.

## Layout

One file per domain under `R/`, 57 exported functions:

| area | files |
|---|---|
| Differential expression | `de_core.r` (stats), `de_wrappers.r` (pipelines), `de_visualization.r` (volcano/scatter/DEG-count plots) |
| Pathway enrichment | `pathway_dispatchers.r` (method dispatch), `ClusterProfiler_MSigDB_functions.r`, `ClusterProfiler_GO_functions.r`, `metascape_functions.r` |
| Clustering & annotation | `seurat_clustering.r`, `cluster_annotation.r`, `llm_annotation.r` |
| ATAC / regulatory | `differential_accessibility.r`, `motif_enrichment.r`, `GRN_analysis_FigR.r` |
| Other analyses | `cellchat_analysis.r`, `cnmf_analysis.r`, `repertoire_analysis.r`, `pseudotime.R` |
| Plotting | `visualization.r`, `heatmaps.r`, `circlize_functions.r` |
| Internal helpers | `shared_utilities.r` |

## Conventions

**`NAMESPACE` and `man/` are generated — never hand-edit either.** Run `devtools::document()` after
touching roxygen comments. Internal helpers are marked `@noRd` and get no `.Rd` file.

**roxygen2 must be 8.0.0**, matching `Config/roxygen2/version` in `DESCRIPTION`. A different version
reformats all 72 `.Rd` files at once and buries the real diff.

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

**Never edit `Version:` in `DESCRIPTION` by hand.** A bump is always the same five steps, ending in
a GitHub release — a bumped version that is not tagged and released is an unfinished job. Run the
whole thing **on the default branch**, so merge the feature branch first:

```sh
git switch main && git merge my-feature     # 1. release from the default branch, clean tree
Rscript -e 'usethis::use_version("minor")'  # 2. bump: major | minor | patch
# 3. write the release notes under the new NEWS.md heading use_version() just added
git commit -am "Bump to 0.12.0"             # 4. commit the bump and the notes
git tag -a v0.12.0 -m "scAnalysisToolkit 0.12.0"
git push --follow-tags                      # 5. tag and push, then:
Rscript -e 'usethis::use_github_release()'
```

Tag with `-a`: every tag in this repo is annotated, and a lightweight one would diverge from the
tag GitHub creates for the release. Push the tag *before* `use_github_release()` — it pushes the
branch but not tags, and if the tag is missing GitHub silently creates its own lightweight one.

Two things make this fail silently or noisily under `Rscript`, both verified:

- **The `which` argument is required.** `use_version()` with no argument asks interactively which
  component to bump; non-interactively it just returns and the version is left unchanged, with no
  error. Always pass `"major"`, `"minor"` or `"patch"`.
- **The working tree must be clean.** With uncommitted changes `use_version()` calls
  `challenge_uncommitted_changes()`, which aborts with `User input required, but session is not
  interactive`. Commit first — that is why the bump is step 2, not step 1. Do not work around this
  by editing `DESCRIPTION` instead.

The tag must point at the commit that carries the bumped `DESCRIPTION`, so it cannot be created
before that commit exists.

**`NEWS.md` is the release notes.** `use_version()` adds a heading for the new version to it, and
`usethis::use_github_release()` publishes whatever sits under the topmost heading as the GitHub
release body — so an empty section ships an empty release. Fill it in before tagging.

Keep it **short and strictly user-facing**: what someone calling these functions has to know. A new
or removed function or argument, a changed default, a behaviour change they would notice, a bug that
was biting them, a new dependency they must install. One line each, no subheadings unless the
release is genuinely large. Leave out refactors, internal helpers, NSE cleanups, roxygen and comment
edits, test and tooling changes, and anything only visible reading the diff — that is what the git
history is for. If a bullet needs a paragraph to explain, it probably belongs in the function's
documentation instead.

**`use_github_release()` will not run from a detached `HEAD`** — `challenge_non_default_branch()`
aborts with `Detached head; can't continue`, so check out the branch rather than the SHA. It reads
the version from `DESCRIPTION` for the release name and tag, uses `HEAD` as the target commit, and
publishes immediately with no confirmation prompt. If the tag already exists on GitHub the API
ignores the target commit and attaches the release to that existing tag. Releases before 0.11.0
have tags but no `NEWS.md` section.

## Code style

- Prefer dplyr/tidyverse over base R.
- Use the native pipe `|>` for chaining.
- Prefer `mutate()`, `filter()`, `select()`, `case_when()`, `arrange()`, `summarize()` over base
  equivalents.
- Use purrr's `map()` family (`map()`, `map_chr()`, `map_dbl()`, `walk()`, …) instead of `for` loops,
  `lapply()`, `sapply()`, or the `apply()` family.
- Follow the tidyverse style guide.
- Prefer `usethis` and `devtools` for package tasks — e.g. `usethis::use_version()`
  (see [Version bumps](#version-bumps)), `devtools::document()`, `devtools::check()`.
