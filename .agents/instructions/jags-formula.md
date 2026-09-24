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

Seeding belongs to the same contract. Initial values are drawn after
`set.seed(seed)`; each chain's `.RNG.seed` comes from `.JAGS_chain_seeds()`
(`sample.int(.Machine$integer.max, chains)` after `set.seed(seed)`, so a
chain's seed does not depend on the number of chains), and automatic restarts
take their seeds from the separate stream of `.JAGS_restart_seeds()`. Every
path that seeds JAGS chains uses these helpers; never derive seeds as
`seed + chain` or `seed + attempt`, which makes adjacent seeds share streams.

Autofit checks convergence of the priors' parameters, the user's
`add_parameters`, and the stochastic nodes that formulas monitor. It leaves
out only the deterministic helper nodes that BayesTools generates for
formulas (correlation matrices, Cholesky factors, bound SDs, derived latent
effects), identified from the persisted formula name maps, unless
`autofit_control$monitor` requests them. A monitored node that the model
syntax defines deterministically and whose draws are identical in every
chain is a structural constant; a stochastic node with constant draws is not
assessable.

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
or treatment coding as ordinary `stats::model.matrix()` can do. An interaction
without one of its lower-order terms (`g:x` in `~ g + g:x`) codes that factor
by level indicators, as `stats::model.matrix()` does (value 2 in the terms
`factors` attribute), and has one coefficient per level. Treatment and
independent priors apply to those level coefficients and name them from the
term design; mean-difference, orthonormal, and ordered priors are defined on
contrast coefficients and are rejected for such a term.

Formula design metadata is authoritative for fixed and random terms. Preserve
the distinction between fitted standardized coordinates, original-scale
display coordinates, unit latent variables, realized group coefficients, and
covariance parameters.

Scaling changes must update design metadata, posterior transformation, prior
transformation, prediction, summaries, and tests together. Do not infer an
unscaling map from sampled values or parameter-name coincidences when the
structural formula metadata is available.

`formula_scale` centers every scaled predictor, also in terms without a free
intercept (maintainer decision): the model is specified on the standardized
scale. `~ 0 + x` with scaled `x` fits `mu = b (x - m) / s`, whose
original-scale intercept is `-b m / s`; original-scale tables report it even
though the standardized intercept is `spike(0)`. A random block without a
random intercept, such as `(0 + x | g)`, likewise implies an original-scale
random intercept `-u m / s`, perfectly correlated with the slope `u / s`; its
summaries report only the block's own coefficients.

Independence is also specified on the centered scale. `(1 + x || g)`, `diag()`,
or `id()` with scaled `x` has independent `u0`, `u1` (SDs `t0`, `t1`); the
original-scale intercept `u0 - u1 m / s` and slope `u1 / s` then have
correlation `-(t1 m / s) / sqrt(t0^2 + (t1 m / s)^2)` (about -0.93 for equal
SDs and `m / s = 2.5`), which such blocks do not report.

Unscaled correlations are `cov / (sd_i * sd_j)` of the transformed
covariance: defined whenever both SDs are positive, including singular draws
(perfect correlation), and missing only when an SD is zero. A rounding excess
of `|r| - 1 <= 1e-8` is set to `sign(r)`; a larger excess is an internal
error. The Cholesky factor and LKJ primitives exist only for positive-definite
correlation matrices, so semantic correlations are read from the unscaled
correlation matrix, not reconstructed from those coordinates. Monitored and
summarized LKJ correlation matrices (the module's `bt_lkj_corr()` and the R
reconstructions through `.bt_lkj_cholesky_L_to_R()`) set their diagonal to
exactly 1 rather than computing it from row products, so `R[k,k]` is a
structural constant in convergence checks.

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
- Square brackets after a factor term always hold a level label, or a cell of
  level labels, never a coordinate position. Every level or cell of a fixed
  factor term (all contrasts, and ordinary factor priors, whose levels are
  1..K) is the quantity `<parameter>[<level>]` with a `factor_level`
  extraction key, and the transformed-summary labels `[dif: <level>]` are its
  aliases. A coordinate is a direct level cell only where the contrast makes
  it so structurally (treatment, independent, first ordered coordinate),
  never by floating-point equality of design rows. Every other coordinate is
  contrast coefficient `<parameter>{j}`, labelled so in tables, diagnostics,
  mixed-posterior columns, and random-slope components (`sd(g{j})`). JAGS
  coordinate names such as `mu_g[2]` remain backend column names, used by
  coordinate-based functions such as `JAGS_materialize_draws()`, but are not
  factor selectors.
- Draws that can be undefined are declared, never inferred from names or
  values: `parameter_draws()` sets the `mcmc.list` attribute `undefined_draws`
  (canonical name to reason, `"correlation"` for original-scale random-effect
  correlations). Summaries such as `ensemble_estimates_table()` accept missing
  draws only for declared columns, summarize the defined draws, and footnote
  their share; any other missing draw is an error.
- Prior draws from `transform_prior_samples()` carry every generated monitor
  that the model defines deterministically from nodes with prior draws
  (allocation-derived SDs, Fisher-z/logit scalar correlations, LKJ factors,
  matrices, and partial correlations), computed with the model's definitions
  and the posterior evaluators, so that catalog quantities built on them
  evaluate on prior and posterior draws alike. Latent group effects, nodes
  derived only from them, and the stochastic auxiliaries of mixture,
  spike-and-slab, and Dirichlet priors have no prior draws; catalog quantities
  that need them (e.g. mixture inclusion) are unavailable from prior draws.
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
  fit attributes. Bump `.bt_parameter_map_version` whenever selectors or map
  semantics change; fits with another version must be refitted.
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
factor or index levels, and curly braces contrast coefficients that are not a
level (`sd(g{j})`). Use public `cor`, while any compact backend `rho`
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

A single grouping factor keeps its declared levels as fitted groups, including
levels without fitting rows, and the fitted term records the observed ones in
`group_observed_levels`. Prediction (conditional and marginal) routes every
unobserved level through `new_levels` exactly like a new level. Blocks with a
known group covariance are the exception: the kernel intends its declared
levels, which keep their fitted coefficients. Multi-variable groupings contain
only observed tuples.

## Verification

Use deterministic syntax and design assertions before live fitting. Follow the
profile order and centralized fixture rules in `testing.md`; do not duplicate
model fitting inside a focused consumer test.
