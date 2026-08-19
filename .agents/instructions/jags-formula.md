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

## Parameter Map

Fitted-parameter metadata has one authoritative, versioned
`BayesTools_parameter_map` with three linked tables:

- `coordinates` maps concrete posterior coordinates, keyed by
  `coordinate_name`, to formula ownership, semantic role, random-effect block,
  fitted scale, monitor status, display metadata, and internal status.
- `quantities` declares selectable semantic quantities, keyed by
  `canonical_name`, together with structured ownership, source provenance, and
  deferred extraction keys.
- `aliases` maps exact accepted selectors to quantity IDs.

Build and attach all three atomically through `R/aaa-parameter-map.R`; the
coordinate and semantic compilers remain pure stages in
`R/JAGS-parameter-registry.R` and `R/JAGS-parameter-catalog.R`.
`parameter_coordinates()` and `parameter_catalog()` are views over the same
stored map. Use these accessors and `R/parameter-source.R` rather than parsing
names.

- Keep coordinate and canonical names unique within their respective schemas,
  and keep every schema field type-stable.
- User-facing summaries, plotting, density estimation, and hypotheses must
  resolve catalog quantities and obtain their draws through
  `parameter_draws()`. Do not promote monitored coordinate rows to public
  aliases.
- Declare source mappings as identity, one-to-one transforms, or composites.
  Record dependencies in extraction keys instead of adding private inputs as
  public catalog rows.
- Internal latent, realized, allocation, LKJ, spike-and-slab, and other
  implementation coordinates remain coordinate-only and must not be presented as
  original-scale public parameters.
- Validate coordinate uniqueness, quantity uniqueness, aliases, extraction
  recipes, and coordinate dependencies atomically. The fit contract stores one
  `parameter_map_version`; there are no separate registry/catalog versions or
  fit attributes.
- Missing, malformed, or unsupported map metadata requires refitting with the
  current BayesTools version. Do not add in-memory migrations for stale fitted
  objects without an explicit maintainer decision.

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

Random-effect catalog names use `(formula) owner: quantity(arguments)`.
Parentheses contain coefficient or parameter names; square brackets contain
factor or index levels. Use public `cor`, while any compact backend `rho`
coordinate remains internal. Total-variance allocations expose `sd_total`,
`var_total`, and `var_prop(...)`; mean-variance allocations expose `sd_common`,
`var_common`, `var_ratio(...)`, and `sd_ratio(...)`.

A bare random formula or unnamed one-entry formula list suppresses a redundant
top-level component prefix. An explicitly named one-entry list retains its
name. In lists with two or more entries, generate missing names as
`component 1`, `component 2`, and so on. Keep allocation `name` as a required
stable backend identifier and use `display_name` and `component_names` for its
public semantic labels.

`random_effects_marginal_vcov()` owns posterior `Z G Z'` construction.
`random_effects_marginal_variance_factors()` exposes validated row-aligned
factors for downstream likelihoods that can marginalize only supported blocks.
`random_effects_marginal_update_plan()` classifies a selected public quantity's
exact scalar covariance update from the parameter map and compiled random
design; downstream optimizations must use this metadata rather than infer an
update form from posterior samples or evaluated covariance matrices.
Do not make a downstream package reconstruct these quantities from raw JAGS
columns.

Preserve explicit new-level policy and the distinction between full covariance
and `diagonal_only` output. Structural unavailability is not permission to
silently substitute zero covariance or a diagonal approximation.

## Verification

Use deterministic syntax and design assertions before live fitting. Follow the
profile order and centralized fixture rules in `testing.md`; do not duplicate
model fitting inside a focused consumer test.
