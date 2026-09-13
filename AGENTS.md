# BayesTools

BayesTools owns reusable Bayesian infrastructure: priors and density algebra,
JAGS formula generation and fitting, marginal likelihoods, model averaging,
posterior summaries, and diagnostics. Consumers such as RoBMA own their
domain-specific models and integration.

When this checkout is inside BayesToolsVerse, follow the shared
[workspace guidance](../AGENTS.md), [validation policy](../.agents/instructions/validation.md),
[R environment](../.agents/instructions/r-environment.md), and
[public API contracts](../.agents/instructions/public-api.md).
They own workflow, research tools, agent scratch space, and the distinction
between tests, verification, and scenarios. A standalone checkout retains the
package contracts below; the optional parent workspace is not a dependency.
Check [pending workspace decisions](../.agents/decisions.md) when present.

## Contracts to consult

Read the guide relevant to the change and maintain it with the implementation:

- [Priors](.agents/instructions/priors.md): constructors, distribution methods,
  transformations, provenance, and structural prior-density ordinates.
- [JAGS and formulas](.agents/instructions/jags-formula.md): fitting paths,
  scaling, contrasts, random effects, parameter maps, and marginal likelihoods.
- [Testing](.agents/instructions/testing.md): profiles, fit catalogs, cached
  references, and visual regression.
- [Plotting](.agents/instructions/plotting.md): plot-data and renderer contracts.
- [Vignettes](.agents/instructions/vignettes.md): precomputed caches,
  regeneration, and citation style.

Formula metadata and the stored parameter map are authoritative. Public
selection and summaries use semantic catalog quantities and `parameter_draws()`;
monitored backend coordinates are not automatically public parameters.
Statistical boundaries, prior support, and conditioning targets must retain
their defined meaning. Check the relevant contract before changing behavior.

## Source ownership

- Priors and density algebra: `R/priors*.R`, `R/distributions-*.R`, and
  `R/prior-density-ordinate.R`.
- JAGS runtime, formulas, and bridge sampling: `R/JAGS-*.R`.
- Random effects: `R/random-effects-*.R` and `R/random-group-covariance.R`.
- Marginal inference and model averaging: `R/marginal-*.R` and
  `R/model-averaging*.R`.
- Summaries and validation: `R/summary-tables*.R`, `R/interpret.R`, and
  `R/tools.R`.

Reuse the validators in `R/tools.R`, including `check_bool()`, `check_char()`,
`check_int()`, `check_real()`, and `check_list()`. Follow the established
namespace-import style. Preserve established file families such as `JAGS-*.R`
and S3/public names; use `snake_case` for new ordinary names. Match nearby R
formatting, including the blank line after a function's opening brace.

## Development and backend

Requires R >= 4.3.0 and JAGS >= 4.3.0, through `runjags`/`rjags`. In the
workspace, use its configured R and private agent library for these commands:

```r
devtools::load_all()
devtools::document()
devtools::test(filter = "topic", reporter = "llm")
devtools::check()
```

```text
Rscript tools/test-profile.R unit
Rscript tools/test-profile.R fixture
Rscript tools/test-profile.R fit
```

Start with the affected tests. The `fit` profile refreshes expensive cached
fixtures; run it when fitting inputs or backend behavior changed, following the
testing guide. Visual profiles are `visual` and `visual-fixture`.

Native JAGS distributions live in `src/distributions/`; shared kernels are in
`src/invgamma/`, `src/lkj/`, and `src/nonlocal/`. Registrations are in
`src/BayesTools.cc` and `src/init.c`. Keep registrations and `Makevars*` source
lists consistent when native sources change.

## Documentation and release

- Document exports with roxygen2 and `\insertCite{key}{BayesTools}`;
  vignettes use Pandoc citations.
- For a completed feature, increment the development version in `DESCRIPTION`
  and update `NEWS.md` for the final feature state.
- Preserve reviewed numerical and visual baselines. Generate candidates for
  inspection; accept changes only after maintainer review or explicitly
  delegated review.
- Use `skip_on_cran()` for computationally intensive tests. Keep `AGENTS.md`,
  `.agents/`, and development-only material excluded through `.Rbuildignore`;
  CI workflows belong in `.github/workflows/`.
