# Vignettes

Use this guide for files under `vignettes/` and their precomputed `.RDS`
objects. Follow the nearest current vignette when a local pattern differs from
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

Every current `.Rmd` has a matching `.RDS` cache in `vignettes/`. Ordinary
vignette rendering must never fit models or compute marginal likelihoods.
Preserve the three-part pattern:

1. validate and load the precomputed cache in setup;
2. display the fitting code with `eval = FALSE`;
3. keep hidden, single-purpose regeneration setup and saving code.

Use `vignettes/precomputed-vignette-cache.R` for the shared cache contract and
the dedicated random-effects cache helper where that vignette requires its
stronger metadata checks. Preserve `purl` controls so regeneration scripts
contain the intended chunks only.

Do not add migrations for stale cached objects. Regenerate them with the current
BayesTools implementation. Use fixed seeds in all regeneration code and record
enough provenance to reject stale or substituted fits. Cheap prior-only
examples should not acquire a fitted-model cache.

Cache validity depends only on source fingerprints: the vignette source, its
cache helper, `DESCRIPTION` without build fields and `Version`, `NAMESPACE`,
package data, and the `R/` and `src/` sources (RandomEffects also fingerprints
its lme4 data sets). Changing any of them makes the cache stale until it is
regenerated. R and package versions are producer provenance: manifests record
them, and validation reports differences in its status and with a message when
the cache is loaded, but they never make a cache stale. Caches therefore build
on newer R versions and after version bumps. Keep integrity checks independent
of the validating session: hash stored serialized bytes, and fingerprint R
objects without the `serialize()` header. Re-serialized output records the
session's R version and native encoding, and fitted models embed the versions
of the namespaces they reference.

The repository tests reject evaluated calls to `JAGS_fit()`,
`JAGS_bridgesampling()`, and `rstanarm::stan_lmer()` during ordinary rendering.
They also enforce BayesTools' marginal-likelihood contract; do not replace it
with direct `bridgesampling::bf()` calls in vignettes.

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

When changing cache contents or fitting code, run `unit`, `fit`, and `fixture`,
then render the affected vignette. Documentation-only prose changes do not
require refitting when the cache contract and displayed results are unchanged.
