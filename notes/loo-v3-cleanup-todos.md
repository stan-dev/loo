# loo v3.0.0 clean-up: open TODOs

> Last update: 2026-09-29
> 
> Author: Florence Bockting

The PRs and issues that relate to this work are in
[loo-v3-related-prs-issues.md](loo-v3-related-prs-issues.md).
The names of the documents are in
[loo-v3-documentation.md](loo-v3-documentation.md). Use only those names.

## Tutorials (new in loo v3)

### `vignettes/articles-online-only/overview-measures.Rmd`

*Model performance: Overview of predictive measures*. PR #363 "pred_measure".

- [ ] Change the YAML `title:` and `\VignetteIndexEntry{}` to the name above.
  Now: "Overview of predictive measures". Waits for: nothing.
- [ ] Search `(TODO-Vehtari et al., 2026)` in the introduction. Cite
  *Predictive measures: A formula reference*. Waits for: arXiv ID.
- [ ] Search `our supplement` in the introduction. Cite
  *Predictive measures: A formula reference*. Waits for: arXiv ID.
- [ ] Search `[TODO-article-suppl]` in the "Aggregation" bullet below the table.
  Link the formula reference. Waits for: arXiv ID.
- [ ] Search `[TODO: arXiv ID]` in the reference list. Waits for: arXiv ID.
- [ ] Search `[TODO-article]` in the callout "What this article does not cover".
  Link *Model performance: Predictive schemes with `pred_measure()`*.
  Waits for: nothing.
- [ ] Search `[TODO: schemes and pred_measure() article]` in the reference list.
  Link the same article. Waits for: nothing.
- [ ] Search `[model-comparison TODO]` in the callout. Link
  *Model comparison: Case Study with `model_compare()`*. Waits for: issue
  #401.
- [ ] Search `](model-comparison.html)`. Use the link text
  *Model comparison: Explanation of `model_compare()`*. That article exists
  only on `integrate-loo_compare`, so the link is broken on `pred_measure`.
  Waits for: #380 merged.

Use the `articles/articles-online-only/` URL form for all links.

### `vignettes/articles-online-only/pred-measure-workflow.Rmd`

*Model performance: Predictive schemes with `pred_measure()`*. PR #363.

- [ ] Change the YAML `title:` and `\VignetteIndexEntry{}` to the name above.
  Now: "Predictive schemes with pred_measure()". Waits for: nothing.
- [ ] Search `[TODO: arXiv ID]` in the reference list. Waits for: arXiv ID.
- [ ] Search `](model-comparison.html)`. Use the link text
  *Model comparison: Explanation of `model_compare()`*. Waits for: #380
  merged.
- [ ] Search `LOAD-BRMS-GITHUB.txt`. Delete the `child=` chunk that loads brms
  from GitHub. Waits for: brms on CRAN.

### `vignettes/articles-online-only/model-comparison.Rmd`

*Model comparison: Explanation of `model_compare()`*. PR #380.

- [ ] Change the YAML `title:` and `\VignetteIndexEntry{}` to the name above.
  Now: "Model comparison: Using the `model_compare()` function.".
  Waits for: nothing.
- [ ] Search `TODO-ARXIV` and `TODO-Vehtari`. Waits for: arXiv ID.
- [ ] Search `the supplement`. Replace it with the name of the
  formula reference. Waits for: arXiv ID.
- [ ] Search `Computing predictive performance measures` in "See also". Use
  *Model performance: Predictive schemes with `pred_measure()`*.
  Waits for: nothing.
- [ ] Search `TODO-case study` and `[TODO-LINK]`. Link
  *Model comparison: Case Study with `model_compare()`*. Waits for: issue
  #401.

## Tutorials (existing)

`loo_compare()` warns once per session (`.deprecate_once()`). A vignette build
therefore shows the warning.

### `vignettes/loo2-large-data.Rmd`

- [ ] Search `loo_compare(loo_ss_1, loo_ss_2)`. `integrate-loo_compare`
  replaced the other 3 calls. Replace it with `model_compare()`. Check the
  output against #393 first. Waits for: #363 merged, subsampling decision.

### `vignettes/loo2-elpd.Rmd`

- [ ] Search `elpd(log_pd` (two calls) and "using `elpd()`" (two sentences).
  Replace them with `measure_elpd()` or `test_pred_measure()`. Waits for: #363 merged.

### All `vignettes/*.Rmd`

- [ ] Grep for `crps(`, `scrps(`, `loo_crps(`, `loo_scrps(` and
  `loo_predictive_metric(`. Do the two files above first. Waits for: #363
  merged.

## Other

### `R/loo-glossary.R`

*The `loo` glossary*. PR #380.

- [ ] Add the general terms: measure, metric, score, utility, loss.
  Waits for: nothing.

## Project management

### `vignettes/migration-guide.Rmd`

*Migration guide*. PR #379, extended on #363, #380 and #378.

- [ ] Merge the four versions. `loo-v3.0.0` has 182 lines, `pred_measure` 307,
  `integrate-loo_compare` 337 and `parallelization` 216. #363 and #378 both
  extend the file, so the second merge into `loo-v3.0.0` conflicts.
  Waits for: #363 merged, #378 merged.
- [ ] Decide if the sections "Maintainer checklist" and "`_pkgdown.yml`
  reference (this branch)" belong in a user vignette. Waits for: nothing.

### `vignettes/articles-online-only/roadmap.Rmd`

*Development roadmap*. PR #396.

- [ ] Delete `vignettes/articles-online-only/roadmap.html` from the PR. The
  PR's own TODO asks for this. Waits for: nothing.

## R code

### `R/pred_measure.R`

PR #363.

- [ ] Search `overview of scores and metrics`. Use
  *Model performance: Overview of predictive measures*. Waits for: nothing.
- [ ] Search `pred-measure workflow article` (6 places). Use
  *Model performance: Predictive schemes with `pred_measure()`*.
  Waits for: nothing.

### `R/pred_measure-compute.R`

PR #363.

- [ ] Search `See developer notes on computation`. Users cannot see `notes/`.
  Cite the formula reference, or delete the `@note`. Waits for: nothing.
- [ ] Search `not yet implemented`. Decide: keep `group_ids` as a reserved
  argument, or remove it. If you remove it, also delete the test in
  `tests/testthat/test_pred_measure.R` (search `group_ids errors`).
  Waits for: nothing.
- [ ] PR #380. Search `include this correction`. This is the warning for moment
  matching and `reloo()`: only `elpd`, `mlpd` and `ic` are correct. Decide:
  ship v3.0.0 with the warning, or wait for #399. Waits for: #399.

### `R/crps.R`

PR #363.

- [ ] Resolve the conflict with #398. Both PRs edit this file. If #398 merges
  first, merge `master` into `loo-v3.0.0`, then `loo-v3.0.0` into
  `pred_measure`. Waits for: nothing.

## Tests

### `tests/testthat/data-for-tests/test_data_generation.R`

PR #363, #380. In 2026-09 `data-for-tests/` is 5.9 MB. CRAN allows a 5 MB
tarball.

- [ ] Search `saveRDS(`. Add `compress = "xz"` to each call. Then re-save the
  `.Rds` fixtures. The earlier re-save was never committed. It gave
  5.85 → 4.13 MB. Waits for: #363 merged (C2).
- [ ] Search `roaches_compare = 110` in `N_KEEP`. Change it to 53. Then
  re-create `test_data_roaches_compare.Rds`. Waits for: #363 merged (C2).

### `tests/testthat/test_compare.R`

PR #380.

- [ ] Search the assertion that expects `diag_diff` to be `""`. With 53
  observations it becomes `"N < 100"`. Rewrite it. Do this with the `N_KEEP`
  change above.

## Package files

### `NEWS.md`

PR #363.

- [ ] Search `[supported_measures_list()]`. Use backticks, as for every other
  name. Waits for: nothing.
- [ ] Search `vignette("migration-guide")`. Add `package = "loo"`, as in the
  roxygen. Waits for: nothing.

### `.Rbuildignore`

- [ ] Add `^notes$`. `loo-v3.0.0` does not exclude `notes/` yet.
  Waits for: nothing.

### `.github/workflows/pkgdown.yaml`

PR #363.

- [ ] Search `paul-buerkner/brms`. Replace it with `any::brms`.
  Waits for: brms on CRAN (C7).

## Branches and PRs (no document)

- [ ] Merge `pred_measure` into `integrate-loo_compare`. Waits for: nothing.
- [ ] Finish the description of #380. The draft is outside the repo. Waits for: nothing.
- [ ] Update the description of #363. "Known limitations" says
  "`loo_compare` integration still outstanding". #380 does that. Mark #220 as
  fixed, not partial. Waits for: #380 merged.
- [ ] Add "Fixes #281", "Fixes #223", "Fixes #213", "Fixes #201",
  "Fixes #135" and "Fixes #220" to the description of #379. GitHub does not
  close them from #363, because #363 does not merge into `master`.
  Waits for: #363 merged.
- [ ] Measure the tarball of `master` with `R CMD build`. It gives the size
  budget. Waits for: nothing.
- [ ] Run `R CMD build` with the vignettes on this branch. Only this gives the
  real size. Waits for: the two fixture items above.