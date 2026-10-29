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
3. Knit loads those objects and skips the chunk when its RDS files all exist
   and their `.key` files match the chunk's current key.
   Otherwise it runs the chunk and saves its objects and keys immediately
   afterward.
   A chunk with multiple objects regenerates the whole group if any are missing.
   Separate independent fits into separate chunks when individual reuse matters.
   Keep printed output, such as `summary()`, out of these chunks: a chunk
   loaded from its cache does not run.

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

The key of a cached chunk hashes its code without comments or formatting, the
deparsed values of the workspace objects the code references, and the keys of
the cached objects it uses. Changing a model's code or its priors therefore
refits that chunk on the next knit, and a refitted model refits the
marginal-likelihood chunks that use it. Keys do not track package versions,
but deparsed objects built by package functions, such as priors, change when
a package changes their structure; refitting an unchanged model with its seed
reproduces the same draws. Removing an RDS file forces a refit of its chunk
only: the key is unchanged, so remove dependent caches too. With refitting disabled, a changed key
loads the existing cache with a warning. Existing caches without a `.key` file
are adopted when knitting with refitting enabled; commit the `.key` files with
their RDS files.
For rendering without fitting, set `options(BayesTools.vignette_refit = FALSE)`
in the rendering session before calling `rmarkdown::render()`; missing caches
then skip the vignette's evaluated examples. This option overrides the
environment setting when explicitly supplied. Keep the declared seeds and
fitting settings. Do not add manifests, migrations, or specialized
regeneration runners beyond these keys. Objects that retain their
formula environment, such as rstanarm fits, serialize the global environment
only by reference, so regenerate caches by knitting into the global
environment (the Knit button and `R CMD build` do); `rmarkdown::render()` with
another `envir` would store the knitting workspace in the cached fits.

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
