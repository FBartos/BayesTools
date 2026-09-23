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

Worker connection failures stop fitting retries; new initial values cannot
repair the existing cluster. Preserve the original backend condition. Classify
parallel worker transport by connection/socket wording rather than a closed
phrase list; unmatched wording in that class is fail-closed (no retry), not a
restartable sampler error. If cluster shutdown fails, attempt every remaining worker and close failed
connections before reporting cleanup failure; do not run the runtime finish
callback when worker shutdown failed. Explicit worker-output paths are
call-specific and must not be replayed from a serialized fit.
After a graceful last-valid-fit return, do not recapture `runtime_state` onto
the retained fit from the failed attempt.

Do not silently repair malformed covariance matrices, alter prior bounds, drop
formula terms, or substitute a different likelihood target. If a covariance
structure requires intersecting a prior with its mathematical support, make the
normalization explicit and preserve the existing warning or validation
contract. Invalid inputs and unsupported paths must fail with a targeted
condition.

## Formula Coordinates and Scaling

Fixed-factor contrasts belong to `prior_factor()`, independently of the
intercept. Removing the intercept specifies a structural zero intercept and
preserves the prior-owned factor basis; it does not silently select indicator
or treatment coding as ordinary `stats::model.matrix()` can do.

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
`R/JAGS-parameter-coordinates.R` and `R/JAGS-parameter-catalog.R`.
`parameter_coordinates()` and `parameter_catalog()` are views over the same
stored map. Use these accessors and `R/parameter-source.R` rather than parsing
names.

- Keep coordinate names and quantity IDs unique within their respective
  schemas, and keep every schema field type-stable. `canonical_name` is a
  selector rather than a key: it must be unique per namespace and component,
  because that is what `parameter_catalog_resolve()` narrows on before raising
  a typed ambiguity. A catalog extended by another provider may reuse a
  canonical name for its own view of the same term, and must not be rejected
  at construction for it.
- User-facing summaries, plotting, density estimation, and hypotheses must
  resolve catalog quantities and obtain their draws through
  `parameter_draws()`. Do not promote monitored coordinate rows to public
  aliases.
- Declare source mappings as identity, one-to-one transforms, composites, or
  structural zeros. A factor-level cell with no nonzero contrast weights is
  `structural_zero`; a single weighted cell is `identity`. Do not infer that
  distinction from posterior draws.
  Record dependencies in extraction keys instead of adding private inputs as
  public catalog rows.
- Internal latent, realized, allocation, LKJ, spike-and-slab, and other
  implementation coordinates remain coordinate-only and must not be presented as
  original-scale public parameters.
- Validate coordinate uniqueness, quantity uniqueness, aliases, extraction
  recipes, and coordinate dependencies atomically at map construction. Public
  accessors reuse that result through the map runtime cache rather than
  rebuilding schema checks on every call, as long as the stored coordinate,
  quantity, and alias tables are still the constructed objects. Replacing those
  tables re-runs validation. The fit contract stores one
  `parameter_map_version`; there are no separate registry/catalog versions or
  fit attributes.
- That cache is session-local: the map carries only a `runtime_cache_id`, and
  the entries live in a bounded package registry. Never attach a cache to the
  map itself - it would be serialized into every saved fit and replayed on load,
  long after the code that derived it changed. Consumers store map-derived
  values through `parameter_map_cache()`, supplying a key covering every other
  input they used; BayesTools owns that environment and discards all providers'
  entries whenever the map tables are replaced.
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
  The left side supports continuous slopes, factor slopes, and interactions;
  `1`, `0`, and `-1` control its intercept. `id()` shares one SD across
  independent columns, `diag()` has one SD per independent column, and `us()`
  has one SD per column plus an unstructured correlation matrix. Reuse existing
  fixed-factor contrasts by default when available. For example,
  `us(0 + group | study)` with independent block contrasts has one correlated
  coefficient per group level and no random intercept; `0 + group` alone does
  not require level indicators.
- `cs()` / `hcs()`, `ar1()` / `ar()` / `har()`, and `car()` own an index basis.
  They reject intercept controls and block contrast overrides. Discrete indices
  use persisted factor levels or sorted unique values; `car()` uses actual
  finite numeric distances. CS/HCS accepts one or more discrete columns,
  combining multiple columns by their observed interaction. AR1/HAR accepts
  exactly one discrete column. Discrete indices accept factor, character,
  numeric/integer, or logical values. CAR accepts exactly one finite
  numeric/integer coordinate or an ordered factor with numeric level labels.

Do not treat `hcs()` as an alias for `us()`: HCS has one common pairwise
correlation and level-specific SDs, whereas US estimates an unrestricted
correlation matrix. Persist basis ownership, resolved index levels, design
columns, and public labels for every consumer.

Complete omitted correlation priors after resolving the structure and dimension:
US/UN uses `LKJ(1)`; CS/HCS uses raw `Uniform(-1 / (K - 1), 1)` whose JAGS hull
is closed while evaluation treats the CS singular bound `-1/(K-1)` as open;
AR1/HAR uses raw `Uniform(-1, 1)`; CAR uses raw `Uniform(0, 1)` and includes
`rho = 0`. Do not treat Uniform as an open-interval family. Explicit scalar
priors retain
their Fisher-z default scale. SD magnitude belongs to the outcome model:
BayesTools must not supply a generic scale, and direct `JAGS_formula()` use
requires an SD prior, SD source, or variance allocation.

Random-effect catalog names use `(formula) owner: quantity(arguments)`.
Parentheses contain coefficient or parameter names; square brackets contain
factor or index levels. Use public `cor`, while any compact backend `rho`
coordinate remains internal. Total-variance allocations expose `sd_total`,
`var_total`, and `var_prop(...)`; mean-variance allocations expose `sd_common`,
`var_common`, `var_mult(...)`, and `sd_mult(...)`.
Gate-only roots are slab-plus-inclusion containers and do not publish a root
`sd_total`; the child-split total is the public realized total.
Conditioning a gated total on presence ANDs the binding-factor (parent) gates
that define a positive realized source, not every nested component gate.

Downstream consumers may explicitly request centrally generated simplified
names. Simplification removes only a sole `intercept` argument and permits
omitting an owner only when resolution remains unique; non-intercept arguments
stay explicit. A known group covariance has a fitted `sd`/`var` kernel scale,
not an `sd_mult`/`var_mult`.

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
Allocation-derived component SDs remain affine on `id`/`diag` and
single-column blocks when their public SD is on the fitted coefficient scale.
On correlated multi-column structures they are factor
`column_scale` updates of the selected quantity, not `A + h(sigma) B`.
Formula-scaled allocation-derived SDs, and SDs that depend on several fitted
coordinates, have no exact single-coordinate update and are `unsupported`.
Do not make a downstream package reconstruct these quantities from raw JAGS
columns.

Preserve explicit new-level policy and the distinction between full covariance
and `diagonal_only` output. Structural unavailability is not permission to
silently substitute zero covariance or a diagonal approximation.

## Verification

Use deterministic syntax and design assertions before live fitting. Follow the
profile order and centralized fixture rules in `testing.md`; do not duplicate
model fitting inside a focused consumer test.
