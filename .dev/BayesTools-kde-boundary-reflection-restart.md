# BayesTools KDE Boundary Reflection Restart Notes

Date: 2026-05-31

## Objective

Switch sample-based KDE calculations to use density reflection when the parameter
has a true finite support boundary. Do not infer reflection boundaries from
plotting ranges, `xlim`, `x_range`, sample minima/maxima, quantile display ranges,
or prior-density grid ranges.

Do not update saved visual snapshots/figures as part of the implementation. Let
visual tests fail so the figure differences can be reviewed manually.

No Git commands were run during the previous work.

## Core Design

Add one internal helper, preferably in `R/priors-density.R` near the density
helpers:

- `.density_kde_boundary(x, n, from, to, bounds, na.rm, ...)`
- `from` / `to` are only the evaluation grid.
- `bounds = c(lower, upper)` is the true support.
- If both bounds are infinite, fall back to ordinary `stats::density()`.
- If one or both bounds are finite, reflect the samples across only those finite
  true support boundaries.
- Preserve the bandwidth chosen by `stats::density()` on the original samples.
- Return a standard `"density"` object and attach
  `attr(out, "boundary_reflection")`.
- Keep zero-height endpoint insertion logic separate. It is a plotting support
  marker, not the boundary correction itself.

Also add a small helper for weightfunction component support:

- `.density.prior_weightfunction_component_bounds(component)`
- `point`: `c(location, location)`
- `beta`: `c(0, 1)`
- independent omega prior on omega scale: prior truncation
- independent omega prior on log-omega scale: exponentiated prior truncation,
  with lower `0` and upper `Inf` when truncation is infinite
- `one_minus_product_beta`: `c(0, 1)`

## Callsites To Wire

Use reflection only where true support is known.

### `R/priors-density.R`

- `.density.prior.simple()`
  - For forced/sample KDE, pass `bounds = c(x$truncation$lower, x$truncation$upper)`.
  - `x_range` remains an evaluation/plotting range.

- `.density.prior.weightfunction(..., individual = TRUE)`
  - Use per-column bounds from `.weightfunction_marginal_components()`.
  - Do not use the prior-wide union truncation for every weight column.

- `.density.prior.orthonormal_or_meandif()`
  - These normally have no finite truncation support. Calling the helper with
    infinite bounds is fine, but do not invent finite bounds from `x_range`.

### `R/JAGS-diagnostics.R`

- `.diagnostics_plot_data_density()`
  - Use `.diagnostics_prior_bounds()` as support bounds.
  - `xlim` and sample range are evaluation ranges only.
  - If custom `transformations` were applied in `.diagnostics_plot_data()`, do
    not reflect unless transformed support is explicitly available. The previous
    attempt attached `attr(model_samples, "custom_transformations")` and used
    infinite bounds when that attribute was true.

### `R/model-averaging-plots.R`

- `.plot_data_samples.simple()`
  - For continuous posterior samples under non-point priors, pass support from
    the non-point prior truncation union.
  - Reflect before any output transformation.
  - Keep point-mass samples excluded from the KDE.

- `.plot_data_samples.weightparameter()`
  - Use `.weightfunction_prior_marginal_components()` to get the component list.
  - Use true component support for reflection.
  - Do not use `.weightfunction_components_range(..., samples = ...)` as support;
    that function is a plotting/tail/sample range.

- `.plot_data_samples.factor()`
  - For simple/treatment/independent bounded factor priors, pass the prior
    truncation bounds.
  - Do not reflect transformed orthonormal/meandif differences unless true
    transformed support is explicitly known.

### Leave Unchanged For Now

- `R/model-averaging-plots.R` `.plot_data_marginal_samples.den()`
  - The available `prior_linear_density` range is a numerical grid, not reliable
    support metadata.

- `R/marginal-distributions.R` `.Savage_Dickey_BF.kd()`
  - No support metadata is available in the function.
  - Boundary-null Savage-Dickey is statistically fragile; do not silently switch
    this inferential calculation to reflected KDE.

## Focused Unit Tests To Add

Use trapezoid integration. Avoid exact curve snapshots.

### `tests/testthat/test-priors-density-numeric.R`

Add tests for:

- Forced sampled KDE for `prior("uniform", list(0, 1))`, `x_range = c(0, 1)`,
  `n_samples = 4000`, `n_points = 512`, `truncate_end = FALSE`.
  - `range(d$x) == c(0, 1)`
  - all finite/nonnegative densities
  - integrated mass near `1`, tolerance around `0.02`
  - endpoint heights clearly nonzero, for example both `> 0.8`
  - `attr(d, "boundary_reflection")` is true

- Forced sampled KDE for `prior("exp", list(rate = 1))`, `x_range = c(0, 5)`,
  `truncate_end = FALSE`.
  - integrated mass near `stats::pexp(5)`, tolerance around `0.02`
  - lower-bound height clearly nonzero, for example `> 0.75`
  - upper plotting endpoint remains small, showing no reflection at artificial
    upper limit

### `tests/testthat/test-JAGS-diagnostic-plot-data.R`

Add a two-chain mock fit with bounded `prior("uniform", list(0, 1))` and samples
inside `(0, 1)`.

Expected:

- chain density grids stay in `[0, 1]`
- finite/nonnegative density
- trapezoid mass near `1`, tolerance around `0.03`
- boundary heights are nonzero; in the previous attempt a robust threshold was
  `> 0.75`
- `attr(chain_density, "boundary_reflection")` is true

### `tests/testthat/test-model-averaging-plots-edge-cases.R`

Add tests for:

- Simple posterior with 20 percent point mass at `0` and 80 percent continuous
  samples in `(0, 1)` under `prior("uniform", list(0, 1))`.
  - point mass is preserved
  - continuous samples exclude exact point samples
  - continuous density integrates to `0.8`, tolerance around `0.03`
  - support is `[0, 1]`

- Omega posterior KDE with exact `omega == 1` samples plus continuous samples in
  `(0, 1)`.
  - point mass equals exact proportion at `1`
  - continuous density integrates to the remaining sample proportion
  - exact `1` samples are not included in the KDE

- Bounded factor posterior under a treatment/independent simple bounded prior.
  - each level density integrates to `1`
  - level metadata is preserved
  - density support is `[0, 1]`
  - boundary reflection attribute is true

## Previous Verification Results

These commands were run before interruption:

- `Rscript -e "devtools::load_all(quiet = TRUE)"`: passed.
- Focused `devtools::test()` runs with `testthat::LlmReporter$new()`:
  - `filter = "priors-density-numeric"`: passed, 44 pass.
  - `filter = "JAGS-diagnostic-plot-data"`: passed after loosening endpoint
    height threshold to `> 0.75`, 90 pass.
  - `filter = "model-averaging-plots-edge-cases"`: passed after matching
    existing factor level-name metadata, 153 pass.
- `Rscript tools/test-profile.R unit`: passed, 3953 pass, 7 expected skips.
- `Rscript tools/test-profile.R visual`: failed only visual snapshots:
  - `test-priors-density.R:38`
  - `test-priors-plot.R:55`
- `Rscript tools/test-profile.R visual-fixture`: failed due to a mix of expected
  snapshot changes and stale/missing fitted-model cache. Snapshot changes were
  reported in JAGS diagnostic/ensemble plot files. Many other failures said:
  "Pre-fitted model cache is missing required artifacts or has stale metadata.
  Run test-00-model-fits.R first."
- `Rscript tools/test-profile.R fit` was started but interrupted by the user
  before completion because branch operations needed to stop.

## Expected Snapshot Failures

Visual failures are expected because reflected KDE changes the density shape near
bounded support. Do not replace snapshot SVGs automatically. Review manually with
`testthat::snapshot_review()` only after the maintainer decides which visual
changes are acceptable.

## Process Note

After the interrupted `fit` profile, multiple `R.exe` workers were still visible
from the shell process list. They were not killed because another agent may have
owned related R work.
