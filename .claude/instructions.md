# Code Style Preferences

## R Code Style

- Always use dplyr/tidyverse syntax instead of base R when working with R code
- Use the native pipe operator `|>` for chaining operations
- Prefer `mutate()`, `filter()`, `select()`, `case_when()`, `arrange()`, `summarize()` over base R equivalents
- Use `map()` family functions from purrr (map(), map_chr(), map_dbl(), walk(), etc.) instead of for loops
- When iteration is needed, use `map()` or `walk()` instead of `for`, `lapply()`, `sapply()`, or `apply()` family functions
- Follow tidyverse style guide conventions
