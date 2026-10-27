# Vignettes

Use this guide for files under `vignettes/` and their cached fits under
`models/`. Follow the nearest current vignette when a local pattern differs from
this guide.

## Files and Headers

Inspect the current directory instead of copying a fixed vignette inventory.
Keep sources UTF-8 and use the standard `rmarkdown::html_vignette` YAML header
used by neighboring files. Use only relative project paths, including
`../inst/REFERENCES.bib` and `../inst/apa.csl`.

Use the established setup options and Windows cairo configuration where
applicable. Vignettes may attach packages because they are user-facing
examples. Gate Suggests-only packages with
`requireNamespace(..., quietly = TRUE)` checks.

## Cached Fits

Use `vignettes/vignette-cache.R`.
The cache directory under `models/` matches the vignette filename without
`.Rmd`, and each object is an ordinary `.RDS` file. Only expensive fitted
models and marginal-likelihood results need caching; cheap prior-only examples
do not need caches.

1. Declare `cached_fits <- vignette_cache(...)` in setup, register
   `vignette_cache_setup(cached_fits)`, and use
   `vignette_cache_eval(cached_fits)` as the default chunk evaluation setting.
2. Give each fitting or marginal-likelihood chunk a `cache_fits` option naming
   its output objects, for example `cache_fits = "fit_sleep"` or
   `cache_fits = c("M0", "M1")`. Keep the fitting code itself ordinary R.
3. Knit loads those objects and skips the chunk when its RDS files all exist.
   Otherwise it runs the chunk and saves its objects immediately afterward.
   A chunk with multiple objects regenerates the whole group if any are missing.
   Separate independent fits into separate chunks when individual reuse matters.

Use `knitr::rmarkdown_notangle` as the vignette engine, as in RoBMA, so package
checks do not extract and execute the fitting code. Fitting and result chunks
are skipped during package checks and when declared optional packages are
absent. Ordinary rendering, including the Knit button, generates missing fits
and can therefore take as long as the declared JAGS/Stan fits and marginal
likelihood computations. The `models/` directory remains excluded from the
source package; use builds without vignettes for routine package development.

The R-CMD-check and pkgdown CI jobs set `BAYESTOOLS_VIGNETTE_REFIT=false`.
This disables missing-cache generation in their direct renders and inherited
build/article subprocesses while retaining vignette-rendering checks.
`R CMD build` runs before the package-check gate: ordinary local or release
builds with vignette building enabled can generate missing fits too. Set
`BAYESTOOLS_VIGNETTE_REFIT=false` in the environment for such builds, or use
`R CMD build --no-build-vignettes` when rendering is not needed. An R option
in the calling session does not configure a separate R process.

To regenerate a fit after changing its code, remove its RDS file and Knit again.
Existing files are reused without checking whether the model code changed.
Remove any dependent marginal-likelihood caches too when changing a model.
For rendering without fitting, set `options(BayesTools.vignette_refit = FALSE)`
in the rendering session before calling `rmarkdown::render()`; missing caches
then skip the vignette's evaluated examples. This option overrides the
environment setting when explicitly supplied. Keep the declared seeds and
fitting settings. Do not add source fingerprints, manifests,
migrations, or specialized regeneration runners. Fit models whose objects
retain their formula environment, such as rstanarm fits, in an environment
containing only their data so caches do not store the knitting workspace.

## Code, Citations, and Prose

Make the scientific point prominent and move fitting boilerplate after the
arguments demonstrated by the example. Align related named arguments without
forcing alignment across unrelated nested calls.

- Use Pandoc citations such as `[@key]`, not `\insertCite{}`.
- Verify every citation key in `inst/REFERENCES.bib`; do not invent entries.
- Keep results, mathematical notation, argument names, function names,
  parameter names, and defined abbreviations exact.
- Write concise, active, scientifically precise prose. Avoid marketing language
  and unnecessary section depth.

For setup or prose changes, run focused checks and render the affected
vignettes with `options(BayesTools.vignette_refit = FALSE)`. Missing caches
permit checking the source-only render but do not verify numerical output.
Follow the workspace policy on obtaining approval before regenerating
vignette fits.
