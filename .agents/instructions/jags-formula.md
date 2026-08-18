# JAGS, Formula, and Fitted-Metadata Contracts

Use this guide when changing `JAGS_fit()`, generated JAGS syntax, formulas,
scaling, contrasts, random effects, posterior extraction, parameter metadata,
or marginal likelihoods.

## Keep the Fitting Paths Synchronized

`R/JAGS-fit.R` owns the main fitting wrapper. Generated syntax, data,
initialization, monitored nodes, posterior extraction, and bridge-sampling
parameters form one contract. A change to one path must be checked against the
others rather than patched only at the first failing consumer.

- Runtime and settings: `R/JAGS-runtime.R`, `R/JAGS-fit-settings.R`, and
  `R/JAGS-convergence.R`.
- Prior syntax, initialization, and monitors: `R/JAGS-prior-*.R`.
- Formula parsing and prepared designs: `R/JAGS-formula*.R`.
- Bridge context and marginal likelihood: `R/JAGS-bridge-*.R` and
  `R/JAGS-marglik*.R`.
- Posterior extraction: `R/posterior-extraction.R` and
  `R/JAGS-bridge-posterior*.R`.

Do not silently repair malformed covariance matrices, alter prior bounds, drop
formula terms, or substitute a different likelihood target. If a covariance
structure requires intersecting a prior with its mathematical support, make the
normalization explicit and preserve the existing warning or validation
contract. Invalid inputs and unsupported paths must fail with a targeted
condition.

## Formula Coordinates and Scaling

Formula design metadata is authoritative for fixed and random terms. Preserve
the distinction between fitted standardized coordinates, original-scale
display coordinates, unit latent variables, realized group coefficients, and
covariance parameters.

Scaling changes must update design metadata, posterior transformation, prior
transformation, prediction, summaries, and tests together. Do not infer an
unscaling map from sampled values or parameter-name coincidences when the
structural formula metadata is available.

The implementation is split across `R/JAGS-formula-scale.R`,
`R/JAGS-formula-scale-random-sd.R`, `R/JAGS-formula-scale-transform.R`, and
the related design and prediction files. Reuse those maps rather than creating
a second transformation path.

## Parameter Registry

`JAGS_parameter_registry()` is the public, versioned mapping from concrete
posterior columns to formula ownership, semantic role, random-effect block,
fitted scale, monitor status, display label, and internal status. It is the
source of truth for downstream packages.

- Build and attach the registry through `R/JAGS-parameter-registry.R`.
- Use `R/parameter-source.R` and existing accessors rather than parsing names.
- Keep canonical names unique and schema fields type-stable.
- Internal latent, realized, and spike-and-slab implementation coordinates must
  not be presented as original-scale public coefficients.
- Missing, malformed, or unsupported registry versions require refitting with
  the current BayesTools version. Do not add in-memory migrations for stale
  fitted objects without an explicit maintainer decision.

## Random Effects

Random-effect compilation and reconstruction are shared infrastructure. Keep
block names, grouping levels, design columns, correlation structure, scale,
monitor names, and reconstruction metadata consistent across
`R/JAGS-formula-random*.R` and `R/random-effects-*.R`.

Keep the parser's two random-effect families separate:

- `id()`, `diag()`, and `us()` / `un()` are random-coefficient formulas. Plain
  bars default to `us()` and double bars to `diag()`. Intercept controls have
  their formula meaning, while factor bases come from stored contrast metadata;
  `random_block(contrasts = ...)` is the explicit block-level override.
- `cs()` / `hcs()`, `ar1()` / `ar()` / `har()`, and `car()` own an index basis.
  They reject intercept controls and block contrast overrides. Discrete indices
  use persisted factor levels or sorted unique values; `car()` uses actual
  finite numeric distances.

Do not treat `hcs()` as an alias for `us()`: HCS has one common pairwise
correlation, whereas US estimates an unrestricted correlation matrix.

`random_effects_marginal_vcov()` owns posterior `Z G Z'` construction.
`random_effects_marginal_variance_factors()` exposes validated row-aligned
factors for downstream likelihoods that can marginalize only supported blocks.
Do not make a downstream package reconstruct these quantities from raw JAGS
columns.

Preserve explicit new-level policy and the distinction between full covariance
and `diagonal_only` output. Structural unavailability is not permission to
silently substitute zero covariance or a diagonal approximation.

## Verification

Use deterministic syntax and design assertions before live fitting. Follow the
profile order and centralized fixture rules in `testing.md`; do not duplicate
model fitting inside a focused consumer test.
