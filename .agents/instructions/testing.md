# Testing

Use this guide for changes under `tests/testthat/`, `tools/test-profile.R`, or
the test workflows.

## Development Workflow

Always use the LLM reporter. Prefer the profile runner, for example:

```powershell
Rscript tools/test-profile.R unit
```

For ad hoc `devtools::test()` calls, set `AGENT=1` and pass
`testthat::LlmReporter$new()`. Run the narrowest relevant profile first; do not
use repeated full-suite runs as an iteration loop.

## Test Profiles

- `unit`: fast package-critical tests without live JAGS fitting or visual
  snapshots.
- `fixture`: tests using cached fits, tables, and reference output.
- `visual`: pure `vdiffr` plot tests.
- `visual-fixture`: visual tests that load cached JAGS fits.
- `fit`: centralized live fitting and marginal-likelihood generation. It
  refreshes the cache unless `BAYESTOOLS_TEST_SKIP_REFIT=TRUE` is set
  intentionally for a validated cache.
- `all`: every profile; this is the default CI and coverage profile.

Use the following verification order:

- JAGS fitting, generated syntax, or marginal likelihoods: `unit`, `fit`, then
  `fixture`; add `visual-fixture` when fitted-object plots can change.
- Formula parsing, scaling, design matrices, contrasts, or prediction targets:
  `unit` and `fit`; add `fixture` when cached objects or reference output can
  change.
- Fixture catalogs, cached fits, reference files, or fixture helpers: `unit`,
  `fit`, and `fixture`.
- Pure plotting: `unit` and `visual`.
- Plotting that loads cached fits: `unit`, `fit`, and `visual-fixture`; add
  `fixture` when tables or fixture metadata change.

## Test Organization and Caches

This package uses testthat edition 3. Do not add `context()`.

- Put `skip_if_not_test_profile()` at the top of every new `test-*.R` file.
- Shared helpers live in `tests/testthat/common-functions.R`.
- Profile routing lives in `tests/testthat/helper-test-profiles.R`.
- Only `tests/testthat/test-00-model-fits.R` may fit ordinary models or compute
  marginal likelihoods for cached fixtures. Other tests must load cached fits.
- The expected fit catalog lives in
  `tests/testthat/helper-expected-fit-catalog.R`. Update the fitting file and
  catalog together when a new cached fit is necessary.
- Reuse existing fits whenever possible. Missing or stale required fits are not
  passing release evidence.

The profile runner stores caches below `BAYESTOOLS_TEST_FILES_DIR`, using a
temporary directory by default. Relevant controls are
`BAYESTOOLS_TEST_PROFILE`, `BAYESTOOLS_TEST_FILES_DIR`, and
`BAYESTOOLS_TEST_SKIP_REFIT`.

Run `fit` after fitting, prior/data input, generated JAGS, native distribution,
scaling, registry, or marginal-likelihood changes. The centralized fit file
refreshes the required catalog by default; reuse it only when the existing cache
has been intentionally validated. Do not run `fit` for unrelated plotting,
summary, documentation, or post-fit changes.

Do not modify `GENERATE_REFERENCE_FILES` unless the maintainer explicitly asks.

## Correctness Evidence

- Test behavior, transformations, failure paths, and invariants rather than
  implementation trivia.
- Use analytic identities or independent reference implementations for
  numerical kernels.
- Justify tolerances from numerical or Monte Carlo error. Do not use a broad
  package-wide tolerance merely because it makes a test pass.
- Do not update committed expected results when a test fails. Determine whether
  the implementation or the verified expectation is wrong and involve the
  maintainer before changing a baseline.
- Do not add redundant matrices, samples, fits, or assertions for coverage
  alone.
- Treat Codecov misses as leads, not goals. Reduce reports to missed clusters,
  then add adversarial assertions only for meaningful behavior. Prefer targeted
  `covr` after the relevant profile; report unrelated local coverage failures
  and rely on CI for the final delta.

## Visual Regression

Use the existing `vdiffr::expect_doppelganger()` pattern and the relevant
visual profile. Structural plot-data tests supplement visual snapshots; they do
not replace them.

Never auto-accept or auto-update snapshots. Ask the maintainer to review every
intentional visual change. Keep stochastic plot inputs deterministic so a
snapshot represents rendering behavior rather than random draws.

## Final Verification

Run focused tests after each meaningful change. Before handoff, run the profile
sequence required by the affected subsystem when practical. State exactly what
ran, what did not, and why.
