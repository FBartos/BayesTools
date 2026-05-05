# Codex Operating Guide

BayesTools is Codex-first. Treat this file as the canonical agent entry point for repository work. Do not add GitHub Copilot instruction files under `.github/copilot-instructions.md`, `.github/instructions/`, `.github/skills/`, or a Copilot setup workflow; update this file instead. Keep normal GitHub Actions CI under `.github/workflows/`.

## Project Map

BayesTools is an R package for Bayesian analyses, JAGS model automation, and Bayesian model averaging.

- Priors: `R/priors.R` defines the core `prior` S3 class. Distribution implementations live in `R/distributions-*.R` and must provide the relevant `rng`, `quant`, `cdf`, and `pdf` behavior.
- JAGS integration: `R/JAGS-fit.R` wraps fitting through `runjags::run.jags`; `R/JAGS-formula.R` handles formula parsing, data preparation, prior assignment, generated syntax, and initialization.
- Model averaging: `R/model-averaging.R` computes posterior model probabilities, ensemble inference, and Bayes factors, with marginal likelihoods from `bridgesampling`.
- Summary and interpretation tables: `R/summary-tables.R`, `R/model-averaging.R`, and `R/interpret.R` define user-facing summaries and verbal interpretations.
- Validation helpers: `R/tools.R` centralizes input checking.

## Coding Rules

- Read the existing implementation and tests before changing behavior. Match local style and keep edits scoped.
- Use `snake_case` for functions, variables, and files. Use dots only where S3 dispatch or existing API style requires them.
- Prefer explicit assignment with `<-`; do not introduce pipe-heavy code.
- Use the internal validation helpers in `R/tools.R` for standard user-input checks, including `check_bool()`, `check_char()`, `check_int()`, `check_real()`, and `check_list()`. Do not use `stopifnot()` for exported/user-facing validation.
- Use clear `stop()`, `warning()`, and `message()` calls. For user-facing errors, prefer concise messages and `call. = FALSE` where it fits existing style.
- Add or update roxygen2 documentation for exported functions.
- Avoid `setwd()`, absolute local paths, hardcoded credentials, and hidden side effects.
- Qualify non-base package functions in examples and snippets when clarity matters.
- Prefer base data structures and type-stable helpers such as `vapply()` unless the surrounding code establishes a different pattern.

## Testing

Always run tests with LLM-oriented reporting. Prefer `Rscript tools/test-profile.R <profile>`; the runner sets `AGENT=1` and uses `testthat::LlmReporter$new()`. For ad hoc `devtools::test()` calls, set `Sys.setenv(AGENT = "1")` and pass `reporter = testthat::LlmReporter$new()`.

Test profiles:

- `unit`: `Rscript tools/test-profile.R unit` for fast package-critical tests, without real JAGS fitting or visual snapshots.
- `fixture`: `Rscript tools/test-profile.R fixture` for cached model fits, tables, and reference outputs.
- `visual`: `Rscript tools/test-profile.R visual` for pure vdiffr plot tests.
- `visual-fixture`: `Rscript tools/test-profile.R visual-fixture` for visual tests that load cached JAGS fits.
- `fit`: `Rscript tools/test-profile.R fit` for slow fitting and marginal-likelihood tests. It refreshes cached fits by default; set `BAYESTOOLS_TEST_SKIP_REFIT=TRUE` only to reuse a validated cache intentionally.
- `all`: `Rscript tools/test-profile.R all` for full developer verification.

Profile policy:

- Push and pull-request CI should keep using `unit`.
- Use heavier profiles as change-triggered or release verification lanes.
- For JAGS fitting, generated syntax, marginal likelihoods, or `test-00-model-fits.R`, run `unit`, then `fit`, then `fixture`; add `visual-fixture` if fitted-object plots can change.
- For formula parsing, scaling, design matrices, or contrasts, run `unit` and `fit`; add `fixture` if cached objects or reference outputs are affected.
- For fixture registries, cached fits, reference files, or fixture helpers, run `unit`, `fit`, and `fixture`.
- For pure plotting changes, run `unit` and `visual`.
- For plotting changes that load cached JAGS fits, run `unit`, `fit`, and `visual-fixture`; run `fixture` too when tables or fixture metadata are affected.

Test authoring rules:

- This package uses testthat edition 3; do not add `context()`.
- Shared helpers live in `tests/testthat/common-functions.R`.
- Add `skip_if_not_test_profile()` at the top of every new `test-*.R` file.
- Only `tests/testthat/test-00-model-fits.R` may fit models or compute marginal likelihoods for cached fixtures. Other tests should load cached fits.
- Do not modify the `GENERATE_REFERENCE_FILES` flag unless the maintainer explicitly asks.
- Reuse existing cached models whenever possible. Inspect `test-00-model-fits.R` and the model registry before creating new fixtures.

## Common Changes

Adding a prior:

1. Add the distribution path to `prior()` in `R/priors.R`.
2. Implement the distribution methods in the relevant `R/distributions-*.R` file.
3. Add plotting support if required.
4. Add tests using the existing prior-testing helpers in `tests/testthat/test-priors.R`.

Changing JAGS fitting:

1. Edit `R/JAGS-fit.R` for general fitting logic or `R/JAGS-formula.R` for formula handling.
2. Run `Rscript tools/test-profile.R unit`.
3. Run `Rscript tools/test-profile.R fit` to rebuild and verify fit fixtures.
4. Run `Rscript tools/test-profile.R fixture`; add `visual-fixture` when relevant.

## Vignettes

Vignettes live in `vignettes/*.Rmd` and use precomputed model objects to avoid CRAN timeouts.

- Use relative paths such as `../inst/REFERENCES.bib`, `../inst/apa.csl`, and `../models/...`.
- Do not use absolute paths.
- Preserve the three-part pattern: setup/check detection, hidden loading of precomputed models, and a hidden `eval = FALSE` chunk documenting how models are regenerated.
- Do not set model-regeneration chunks to `eval = TRUE` unless intentionally refreshing cached vignette models.
- Use fixed seeds in documented model-generation code.
- Keep prose concise, direct, and scientifically precise. Preserve references, results, mathematical notation, argument names, function names, parameter names, and defined abbreviations.

## Literature-Backed Writing

This project uses a manuscript-local StatsVault workspace. Before writing literature-backed prose, read `.StatsVault/AGENT_INSTRUCTIONS.md` and follow that workflow.

Key workspace files:

- `.StatsVault/AGENT_INSTRUCTIONS.md`: full source-grounding protocol.
- `.StatsVault/project-memory.md`: persistent working memory.
- `.StatsVault/key-papers.md`: core papers for this project.
- `.StatsVault/ACTION_REQUIRED.md`: missing PDFs or unresolved items.
- `.StatsVault/references.bib`: StatsVault-managed bibliography.

Hard rules for literature-backed writing:

- Cite only keys present in `.StatsVault/references.bib`.
- Materialize every paper you rely on before drafting from it.
- Do not invent support. Use `\SV{...}` markers for unresolved evidence gaps.
- Validate manuscript sources before considering a literature-backed draft complete.

## Repository Hygiene

- Agent and workspace coordination files should be committed so the same Codex setup works across development machines.
- Agent-only files must stay out of R source packages through `.Rbuildignore`.
- The small StatsVault coordination files are Git-trackable; bulky local caches, PDFs, query caches, and contribution artifacts remain ignored unless the maintainer explicitly asks to commit them.
- Keep CI workflows in `.github/workflows/`; do not use `.github/` as the place for agent instructions.
