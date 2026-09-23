# Testing

Use this guide for changes under `tests/testthat/`, `tools/test-profile.R`, or
the test workflows.

In BayesToolsVerse, the shared validation guide defines tests, verification,
and scenarios. Profiles below select existing runner lanes; they do not make
expensive fitting part of the routine development loop. Use the workspace's
configured R and private agent library when available.

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

Start with relevant `unit` tests. Use `fixture` for cached-object behavior and
`visual`/`visual-fixture` for affected plotting. Add `fit` before dependent
profiles when a fitting change invalidates required caches; post-fit or
plotting changes alone do not require refitting. Establish one representative
path before expanding expensive validation.

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
scaling, parameter-map, or marginal-likelihood changes. The centralized fit file
refreshes the required catalog by default; reuse it only when the existing cache
has been intentionally validated. Do not run `fit` for unrelated plotting,
summary, documentation, or post-fit changes.

Do not modify `GENERATE_REFERENCE_FILES` unless the maintainer explicitly asks.

The interactive `test_tests()` runner, loaded by the project `.Rprofile` and
invoked comprehensively by sourcing `.dev/user-tests.R`, caches changed
`test_reference_table()` output beside its baseline as `<name>.new.txt`. After
the test summary, it opens testthat's snapshot reviewer for explicit
Accept/Reject/Skip decisions. Accept replaces the baseline, Reject removes the
candidate, and Skip keeps it for later review. Non-interactive runs retain
candidates and never update baselines.

`test_tests()` runs all five lanes by default. It reuses a validated fit cache;
use `refit = TRUE` to clean and rebuild it. With `filter`, refitting first runs
the centralized `fit` lane and then the selected test files. `regenerate = TRUE`
combines refitting with forced snapshot review. BayesTools has no timing
baselines, so `update_timings = TRUE` fails explicitly. Interactive calls use
the standard progress reporter and leave `AGENT` unset; agent-oriented output
is available explicitly with `reporter = "llm"`.

## Correctness Evidence

- Test behavior, transformations, failure paths, and invariants rather than
  implementation trivia.
- Use analytic identities or independent reference implementations for
  numerical kernels.
- Justify tolerances from numerical or Monte Carlo error. Do not use a broad
  package-wide tolerance merely because it makes a test pass.
- A failing expectation requires diagnosis. Generate a candidate when the
  intended result changed; accept a verified baseline change only after
  maintainer or explicitly delegated review.
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

Retain candidates for review; accept intentional visual changes only after
maintainer or explicitly delegated review. Keep stochastic plot inputs
deterministic so a snapshot represents rendering behavior rather than random
draws.

## Final Verification

Run focused tests after each meaningful change. Before handoff, run the profile
sequence required by the affected subsystem when practical. State exactly what
ran, what did not, and why.
