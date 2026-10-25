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

Use `vignettes/vignette-cache.R`, following the same pattern as RoBMA.
The cache directory under `models/` matches the vignette filename without
`.Rmd`, and each object is an ordinary `.RDS` file. Only expensive fitted
models and marginal-likelihood results need caching; cheap prior-only examples
do not need caches.

1. Declare `cached_fits <- vignette_cache(...)` in setup and use
   `vignette_cache_eval(cached_fits)` as the default chunk evaluation setting.
2. Load objects with `vignette_cache_load(cached_fits)` in a hidden chunk.
3. Display fitting and marginal-likelihood code with `eval = FALSE`.
4. Save objects with `vignette_cache_save(cached_fits)` in a hidden,
   non-evaluated chunk.

Use `knitr::rmarkdown_notangle` as the vignette engine, as in RoBMA, so package
checks do not extract and execute the fitting code. Ordinary rendering never
runs expensive JAGS/Stan fits or marginal-likelihood computations. Result
chunks are skipped during package checks and when caches or declared optional
packages are absent. The `models/` directory remains
excluded from the source package.

Regeneration is manual: from the vignette directory in a fresh R session with
the current package loaded, run setup, data/prior definitions and fitting
chunks, then the saving chunk; skip the cache-loading chunk. Keep the declared
seeds and fitting settings. Do not add source fingerprints, manifests,
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
vignettes without refitting. Missing caches permit checking the source-only
render but do not verify numerical output. Follow the workspace policy on
obtaining approval before regenerating vignette fits.
