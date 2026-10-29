# Maintainer scenarios

Scenarios are maintainer-facing analyses of real data with locked text tables,
SVG figures, and timings. They are excluded from routine package tests and from
the built package. The framework mirrors RoBMA's `tests/scenarios/`.

Name each file `test-<scenario>.R` and run it with:

```text
Rscript tools/test-scenario.R <scenario>
```

Use `--refit` to refit wrapped models while comparing against the locked
outputs, `--update` to reuse valid fit caches and replace exercised text and
plot baselines, and `--regenerate` for both. Use `--update-timings` to accept
current timings explicitly, including slower ones. `--list` lists the
scenarios. Review all resulting Git diffs.

The project `.Rprofile` loads `test_scenarios()` for interactive sessions. It
reuses fit caches and compares without replacing text or plot baselines;
`filter` is a regular expression matched against scenario names:

```r
test_scenarios()
test_scenarios(filter = "penicillin")
test_scenarios(filter = "penicillin", refit = TRUE, update_timings = TRUE)
review_scenario_snapshots()
plot_scenario_times("penicillin", unit = "m")
plot_scenario_times(functions = c("JAGS_fit", "JAGS_bridgesampling"))
```

Minimal scenario:

```r
scenario_start("example")

testthat::test_that("Example analysis remains stable", {

  fit <- scenario_fit("fit", {
    JAGS_fit(...)
  })

  scenario_text("estimates", {
    JAGS_estimates_table(fit)
  })

  scenario_plot("posterior", {
    plot_posterior(...)
  })
})
```

Sourcing a scenario file directly in an interactive session also reuses fit
caches, but shows artifact output and immediately creates or replaces its text
and SVG baselines.

- `scenario_fit()` caches a fit under the ignored `cache/` directory, keyed by
  its normalized fitting call; a changed call refits. The key covers only the
  code inside the block, not the values of objects it references, so define
  priors and other model settings inside the block (or refit). Pass a positive
  `cache_version` only when a package-internal change should invalidate an
  unchanged call. A cached BayesTools fit that the current build refuses (for
  example, an older parameter-map version) stops with an error of class
  `BayesTools_scenario_stale_cache`; refit with `--refit` or
  `test_scenarios(filter = "<scenario>", refit = TRUE)`. Bridge sampling run
  inside a fit block is timed as its `fit_marglik` phase.
- `scenario_text()` prints a visible returned value into `results/<scenario>/`.
  Changed text is stored as `<name>.new.txt` for review.
- `scenario_plot()` stores SVG baselines under `_snaps/<scenario>/` and
  requires `vdiffr`. `scenario_agreement_plot()` draws estimate-minus-reference
  differences against the reference with a band of 0.1 reference SD.
- `scenario_time()` times a block without storing its value.

`scenario_text()` and `scenario_plot()` each reset the random seed to 1 before
evaluating an artifact. Text under `results/`, SVG files under `_snaps/`, and
timing baselines under `timings/<scenario>.tsv` are committed; fit caches and
`timings/<scenario>.new.tsv` (the last run) are not.

Timing comparison runs after the complete scenario file finishes. One warning
reports every call more than 25% and at least 0.5 seconds slower than its
baseline, and an unweighted mean regression across calls above 15%; calls
measured below 0.75 seconds are excluded. Missing and improved time and memory
measurements are written automatically; slower ones need `--update-timings`.
`memory_gb` is the peak of R's garbage collector and excludes external
processes such as JAGS.

An interactive `test_scenarios()` run collects all table and figure mismatches
and opens testthat's snapshot reviewer for Accept, Reject, or Skip decisions;
`review_scenario_snapshots()` reopens them. Non-interactive runs retain the
candidates and fail. Runs report locked artifacts that are no longer referenced
by a scenario and never delete them.
