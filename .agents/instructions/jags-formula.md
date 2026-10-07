# JAGS, Formula, and Fitted-Metadata Contracts

Use this guide when changing `JAGS_fit()`, generated JAGS syntax, formulas,
scaling, contrasts, random effects, posterior extraction, parameter metadata,
or marginal likelihoods.

## Keep the Fitting Paths Synchronized

Unix builds with an explicit JAGS_ROOT or --with-jags-prefix use that prefix's
headers, libraries and version, preserving explicit include/lib/version overrides.
They never inherit unrestricted global pkg-config flags or fall back to unrelated
system headers. No-prefix automatic discovery retains its existing behavior.

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
- Generated deterministic nodes: `R/JAGS-deterministic-nodes*.R`.

Ordinary scalar and vector bridge parameter extraction checks required monitor
names before reading a posterior row and uses `BayesTools_missing_monitored_columns`
(parent `BayesTools_marglik_input`) for absent coordinates. Preserve row-value
shapes and names; a present NA value is distinct from an absent name, and literal
point priors need no monitor.
Bridge sampling checks affine rank of varying finite supplied coordinates before
bound transformations, then checks the transformed coordinates separately with
the existing QR tolerance and diagnostics. Constant columns retain the existing
sampler policy. Before each check, finite nonzero columns are rescaled by their
maximum absolute value in a local diagnostic copy to avoid variance overflow
or underflow; actual draws and sampler bounds remain unchanged. This is not a
general nonlinear dependence detector.

Every deterministic node BayesTools generates (allocation-derived SDs, scalar
and LKJ correlations, publication weights, spike-and-slab and mixture
parameters, formula linear predictors) belongs to one registered family whose
node specification emits its JAGS syntax and carries the R evaluator and the
declared dependencies: `JAGS_formula()` writes the linear predictor from the
node of the fitted design, and `selection_backend_spec()` writes the weights
from the `omega` node specification. `JAGS_deterministic_nodes()` lists the
nodes; `JAGS_deterministic_evaluator()` resolves them once into an evaluator
of draws (repeated or one-draw evaluation), and `JAGS_evaluate_deterministic()`
is that evaluator applied once. A prepared evaluator may locate a family's
columns once per set of draw column names, but shares the family's lookup
rules, arithmetic, and checks (values identical to the family evaluator).
Prior draws, catalog quantities, bridge and marginal-likelihood
parameters, prediction, marginal posteriors of formula parameters,
random-effect unscaling, and convergence-role parents use the family
evaluators. Declared dependencies are the coordinates the R evaluator reads. A
new generated node gets a family and a parity test against the JAGS monitors,
not another evaluator. Intercept priors cannot carry `multiply_by`; it scales
only formula-term coefficients.

Seeding belongs to the same contract. Initial values are drawn after
`set.seed(seed)`; each chain's `.RNG.seed` comes from `.JAGS_chain_seeds()`
(`sample.int(.Machine$integer.max, chains)` after `set.seed(seed)`, so a
chain's seed does not depend on the number of chains), and automatic restarts
take their seeds from the separate stream of `.JAGS_restart_seeds()`. Every
path that seeds JAGS chains uses these helpers; never derive seeds as
`seed + chain` or `seed + attempt`, which makes adjacent seeds share streams.
Seeded functions restore the caller's `.Random.seed` and `RNGkind()`.

Autofit, `JAGS_extend()` and the default `JAGS_check_convergence()` share one
selection read from the fit's `convergence_role` column. The roles are set at
fit time from the prior list, the parsed model syntax, the fully observed
data names and the formula name maps, never from the draws. Generated
deterministic formula monitors are derived and checked only on request. Unit
correlation diagonals and Cholesky constants are structural. A user
`add_parameters` node is structural only when it is fully observed data or its
deterministic definition has no stochastic ancestor. Indicators and inclusion
probabilities come from their priors, not their names. The backend anchor is
auxiliary.

Worker connection failures stop fitting retries; new initial values cannot
repair the existing cluster. Preserve the original backend condition. Classify
parallel worker transport by connection/socket wording rather than a closed
phrase list; unmatched wording in that class is fail-closed (no retry), not a
restartable sampler error. If cluster shutdown fails, attempt every remaining worker and close failed
connections before reporting cleanup failure; do not run the runtime finish
callback when worker shutdown failed. Explicit worker-output paths are
call-specific and must not be replayed from a serialized fit.
After a graceful last-valid-fit return, do not recapture `runtime_state` onto
the retained fit from the failed attempt. A compiled model is continued only
within one call: returned fits carry no runjags `method.options$rjags`, and
`JAGS_extend()` recompiles from the stored chain states, so extending a fit
and a reloaded copy give identical draws.

Do not silently repair malformed covariance matrices, alter prior bounds, drop
formula terms, or substitute a different likelihood target. If a covariance
structure requires intersecting a prior with its mathematical support, make the
normalization explicit and preserve the existing warning or validation
contract. Invalid inputs and unsupported paths must fail with a targeted
condition.

## Formula Coordinates and Scaling

Selected existing terms preserve the fitted terms' variable order, coding,
raw columns, levels, and concrete contrasts on the actual matrix input.
Equivalent interaction orders resolve before prior lookup. Selected replay
requires only its own predictors and checks raw column identities.

Expression-valued scalar point coefficients retain a compact private owner
from the actual full prior declarations and needed model data. Its original
AST and array dimensions are independent of prediction rows. Supported replay
uses the bounded expression grammar in scalar mode, with no loop `i`, through
the registered `point_expression` family; missing parents, syntax, scalar,
index, cycle and owner failures remain distinct typed unavailability.
Valid monitors supply ordinary values first; explicit replay replaces stale
derived values from available parents. A monitor-only coefficient can supply
its dependent formula without claiming the coefficient was recomputed.
Affected missing owners require refitting. These additive owners leave
unaffected formula/map schemas unchanged: expression points were already
derived in map schema 10, and their new declared parents are compiled at
the actual owner. Full public bridge prior-expression refusal remains.

Bridge fixed and random reconstruction use one explicit natural state:
retained scalar constants, decoded ordinary/formula-prior values, genuine
sampled extras, and then available reconstructed formula parents. Natural
owners override transformed row aliases; malformed or contradictory owners
never fall back. Scalar SD and marginal batch/cache routes forward that state.
Prediction data cannot shadow retained multiplier constants.

Fixed-factor contrasts belong to `prior_factor()`, independently of the
intercept. Removing the intercept specifies a structural zero intercept and
preserves the prior-owned factor basis; it does not silently select indicator
or treatment coding as ordinary `stats::model.matrix()` can do. An interaction
without one of its lower-order terms (`g:x` in `~ g + g:x`) codes that factor
by level indicators, as `stats::model.matrix()` does (value 2 in the terms
`factors` attribute), and has one coefficient per level. Such a term does not
use the factor's contrast: it records the independent (identity) coding for
that factor in its `factor_contrasts`, and its prior does not take part in the
agreement of contrast priors across the factor's terms (so `g` mean-difference
with `g:x` independent is valid). When no term codes the factor by its
contrast (no main effect and no contrast-coded interaction, e.g.
`~ g:x + g:z`), the agreement check applies among its indicator-coded terms
instead: they must agree on one contrast. Treatment and independent priors
apply to those level coefficients and name them from the term design;
mean-difference,
orthonormal, and ordered priors are defined on contrast coefficients and are
rejected for such a term. The one exception is a mean-difference or
orthonormal point mass at zero; an ordered prior is rejected even when its
total is zero.

Formula design metadata is authoritative for fixed and random terms. Fixed
and random coefficient formulas may use literal whole interaction-expansion
degrees from 2 through `.Machine$integer.max`, such as `(x + z)^2`. This is the
formula expansion operator, not arithmetic `I(x^2)` or expression arithmetic.
Public degree 0/1 and symbolic powers retain the existing base-R terms refusal;
fractional degrees retain the grammar refusal. Do not reorder formula parsing.
Every factor prior carries complete factor metadata (levels, level names,
contrasts, design, and cell names): formula terms and random-effect SD factor
priors from their design, ordinary factor priors from
`prior_factor_levels()`. Levels and coordinate names are never inferred from
a coefficient count; fitting stops on incomplete factor metadata. Formula
evaluation and prediction (`JAGS_evaluate_formula()`,
`JAGS_predict_formula()`) go only through the stored design: draws without a
fit get it from `JAGS_formula_draws()`, which builds it as `JAGS_fit()` does;
never evaluate from prior-list factor metadata. Formulas stored in designs and
fits keep no caller environment (only the base or empty environment), or a
saved fit serializes the caller's workspace; random-effect term formulas read
their predictors only from data columns, and prediction stops when the data
lack one (`test-fixture-integrity.R` walks every cached fit for this). Preserve
the distinction between fitted standardized coordinates, original-scale
display coordinates, unit latent variables, realized group coefficients, and
covariance parameters.

Replayable data-array expressions pair multiple indices pointwise as JAGS does;
only scalar indices repeat to a common nonzero length. Coordinates must be
finite positive integers within their declared array dimensions. One-index
parameter reconstruction retains its existing vector rules. A coefficient
transform request for an absent formula parameter, including scalar and weighted
prior-density requests, raises `BayesTools_formula_transform_unavailable` with
reason `missing_formula_design` after fit-contract validation.

Scaling changes must update design metadata, posterior transformation, prior
transformation, prediction, summaries, and tests together. Fixed-coefficient
unscaling is derived only from the fitted design that `formula_scale` carries
(`unscale_design`); `formula_scale` without it stops with
`BayesTools_formula_transform_unavailable`, also for random-effect SD columns
alone. Routes that receive the fit take its fitted structure (design, log
intercept, random-effect SD structure) also for a passed `formula_scale`.
Never infer an unscaling map from sampled values or parameter-name
coincidences; tests build `formula_scale` through `JAGS_formula()` (see
`formula_scale_for_test()`).

Formula-design 6, unscale-design 2 and coefficient-transform 3 require refitting
older formats. One private finalizer completes current compiler specs from
actual full priors and owned scalar data for fit/draws-only producers, writing
the same owner to both carriers; fit readers never complete missing ownership.
Explicit prior-only entry points may complete a current compiler copy. Numeric
overrides retain fitted declarations. Original coefficients use
`sum_k A[j,k] * z_k * m_k / m_j` on one immutable original snapshot for every
prefix, with equal canonical multipliers canceled even at zero and declared
zero terms elided. A different needed zero denominator refuses the requested
numeric vector. C has finite static rows and all-NA dynamic rows; A remains
finite. Nonzero ratio arithmetic preserves the ordinary multiplication-first
grouping when its intermediates are normal, then tries ratio-first grouping;
unresolved range loss refuses rather than declaring a zero coefficient.
Weighted original-scale requests refuse nonidentity exp outputs, including
state-dependent maps; a scalar exp target retains its declared transformed law.
Formula values/measures use compiled contributions and fitted rows on
both views. Public marginal `at` values retain fitted/SD units on fitted views
and original units on original views; standardize only the latter to compare
the same physical row. Retain fixed sources, actual multipliers and needed gate/replay
parents, never a sampler/training matrix/random latent closure. Exact posterior
model probabilities and original eligible joint gate frequencies survive output
selection and allocation; unknown certificates remain explicitly unavailable.

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

For a log-intercept formula, `exp(intercept)` aliases the fitted coefficient
quantity; selection alone does not unscale it. Downstream display scale, such
as RoBMA's `standardized_coefficients`, is selected independently.


Fitted-parameter metadata has one authoritative, versioned
`BayesTools_parameter_map` with three linked tables:

- `coordinates` maps concrete posterior coordinates, keyed by
  `coordinate_name`, to formula ownership, semantic role, random-effect block,
  fitted scale, monitor status, display metadata, and internal status.
- `quantities` declares selectable semantic quantities, keyed by
  `canonical_name`, together with structured ownership, source provenance, and
  deferred extraction keys.
- `aliases` maps exact accepted selectors to quantity IDs, each with the label
  parts whose table label is the alias (the alias text itself as one
  component where the alias is no table label of structured parts, e.g. a
  canonical name); `parameter_labels(aliases, "table", vocabulary = )`
  renders random-effect aliases, including the pairwise aliases of `cs()`
  and `hcs()` correlations, under a caller's quantity names. Extension
  aliases may omit their parts.

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
  set by `prior_factor_levels()`: the given names, or 1..K for a count) is
  the quantity `<parameter>[<level>]` with a `factor_level`
  extraction key, and the `[dif: <level>]` labels of transformed
  mean-difference and orthonormal levels are its aliases (transformed
  treatment and ordered levels are labelled by the cell itself). A
  coordinate is a direct level cell only where the contrast makes
  it so structurally (treatment, independent, first ordered coordinate),
  never by floating-point equality of design rows. Every other coordinate is
  contrast coefficient `<parameter>{j}`, labelled so in tables, diagnostics,
  mixed-posterior columns, and random-slope components (`sd(g{j})`). A
  `{j}` selector of a direct level coordinate is refused (in catalog
  resolution and in `hypothesis_parse()` with a catalog) with class
  `BayesTools_selector_unavailable`, naming the level in the same form. JAGS
  coordinate names such as `mu_g[2]` remain backend column names, used by
  coordinate-based functions such as `JAGS_materialize_draws()`, but are not
  factor selectors.
  A continuous coefficient's backend name is not a level label: aliases such
  as `x[mu_x]` and `intercept[mu_intercept]` are unavailable. Parameter-map
  schema 10 requires refitting older maps rather than retaining their aliases.
- Every label is rendered by one renderer (`.bt_label()`, exported as
  `parameter_labels()`) from structured label parts (`R/parameter-labels.R`):
  catalog selectors and aliases, table rows, mixed-posterior column names,
  plot legends, diagnostics titles, and warnings. Catalog quantities store
  their parts in `label_parts`; mixed, transformed-factor, and marginal draws
  store per-column parts, catalog quantity ids, and fitted-coordinate
  dependencies and weights in the `quantities` draw metadata. A column's
  quantity id and parts describe the values it holds: a fitted-scale column
  whose original-scale quantity differs (the SDs of a block with a
  standardized random slope) is its fitted coordinate, and the transform to
  the original scale gives it the catalog quantity it keeps in the internal
  `original_scale_quantities` field. Estimates tables (model, ensemble, and
  marginal) carry the per-row `quantities` attribute (row, quantity id, label
  parts); values transformed after labelling (`posterior_transform()`, the
  `transformations` of model tables) record the transformation in their parts
  and have no quantity id. Consumers map
  columns to coordinates and render labels from these, never by parsing label
  text; original-scale transforms of mixed columns go through the fitted
  design by these coordinates and stop for columns that do not identify them
  or that the design does not contain, and for random-effect SDs that the
  random-effect structure of the passed scale does not cover.
- Catalog quantities declare their exact `support` (from the prior
  provenance of their source coordinates; `NULL` when not derivable, and a
  composite SD's `[0, Inf)` hull is exact only when every scale prior is
  supported on `[0, Inf)`) and `definedness` when the map is built;
  downstream packages read
  them instead of deriving supports.
- Draws that can be undefined are declared, never inferred from names or
  values: `parameter_draws()` sets the `undefined_draws` draw metadata of its
  `mcmc.list` from the catalog's `definedness` (canonical name to reason:
  `"correlation"` for original-scale random-effect correlations,
  `"allocation_active"` for variance shares of gated total-variance
  allocations). Summaries such as `ensemble_estimates_table()` accept missing
  draws only for declared columns, summarize the defined draws, and footnote
  their share; any other missing draw is an error.
- Prior draws from `transform_prior_samples()` carry every generated monitor
  that the model defines deterministically from nodes with prior draws
  (allocation-derived SDs, Fisher-z/logit scalar correlations, LKJ factors,
  matrices, and partial correlations), computed with the model's definitions
  and the posterior evaluators, so that catalog quantities built on them
  evaluate on prior and posterior draws alike. The indicators of mixture and
  spike-and-slab priors and the inclusion probability and slab draws of
  spike-and-slab priors come from the components of their `rng()` draws (the
  stream of the other columns is unchanged); so do ordered-prior totals with
  these nodes of a spike-and-slab or mixture total, whose theta slices share
  one inclusion probability and indicator, as in the fitted model. Latent
  group effects, nodes derived only from them, the component nodes of
  mixtures, and the auxiliaries of Dirichlet priors have no prior draws;
  catalog quantities that need them are unavailable from prior draws.
- Internal latent, realized, allocation, LKJ, spike-and-slab, and other
  implementation coordinates remain coordinate-only and must not be presented as
  original-scale public parameters. Raw estimates tables show them as backend
  coordinates with rendered labels that are not catalog aliases; semantic
  tables omit them, including the gamma allocation nodes of ordered-factor
  priors.
- Ordered-factor terms are summarized like the other contrasts: transformed
  tables show the level effects (without design-fixed zero levels and without
  the total, which the last level repeats); untransformed model tables show
  the sampled parameters, the total and the normalized allocation shares
  (display rows derived from the gamma nodes, not catalog quantities), on the
  fitted scale. Mixed producers retain these primitive sources with their
  selected model and row; raw ensemble tables use their complete declared-model
  union, including fixed allocations and undefined shares in absent models.
  `JAGS_ordered_parameter_spec()` owns source names, labels, tensor projections
  and exact scalar point states. Every available ordered semantic value uses
  its retained primitive projection, including continuous cumulative effects;
  declared view transformations and undefined rows retain their meaning.
  Formula projections receive fitted weights and reuse the aligned fitted
  projection context directly, including ordinary coefficients, without
  reapplying an original-scale coefficient view. Select owners from nonzero
  fitted weights and retain every provided source of each contributing prefix;
  unrelated prefixes are excluded. Scaling belongs to the registered prefix's
  emitted scale record. Contexts agree within a prefix; combined-prefix targets
  union aligned sources by authoritative names and refuse contradictory overlaps.
  Raw ordered producers canonicalize each prefix from its retained primitive
  sources unless that prefix is being unscaled through its own emitted scale
  record. An unscaled prefix therefore canonicalizes raw values and atoms under
  either transform_scaled setting, even when another prefix is scaled. Scaled mixed
  objects without needed fitted ordinary sources are unavailable and must be
  recreated from their source fits; never combine fitted weights with original-scale
  draws or reconstruct fitted draws by inverting the scale map.
  Registered ordered families replay Gamma
  allocations and coefficient chains; `JAGS_ordered_density_kernel()` is the
  separate active-component numerical kernel, not a bridge density. Numeric
  literals and emitted total syntax are persisted as literal/emission evidence
  under the fitted parameter-map contract version; changes to their formatting
  require a deliberate contract-version decision, never an in-memory migration
  or acceptance of source-incomplete fits. Missing provenance requires a real
  refit. Continuous total families remain continuous with expression parameters;
  Bernoulli expression probabilities retain the declared numeric support.
  Expression point locations provide fitted snapshot values only unless a
  certified ancestor recipe exists, and never acquire point atoms from their
  wrapper or observed constancy.
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
  entries whenever the map tables are replaced. The transformed prior densities
  of `as_mixed_posteriors(transform_scaled = TRUE)` live there, for the two
  most recently used input sets of each fit.
- Repeated pure work shares one session memo (`.bt_content_memo()`): validators
  of catalogs, selections, hypothesis ASTs, and label parts
  (`.bt_validate_once()`) and the formula coefficient transform keep recent
  inputs by reference and recognise them only when `identical()`. A modified
  object is validated or rebuilt again; nothing is marked on the object itself.
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

Complete omitted correlation priors after resolving the structure and dimension.
Missing local covariance subcalls inherit the shared covariance; supplied
covariance constructors replace it completely with their constructor defaults,
including removal of its covariance-only SD fallback. Parent prior printing uses
the same resolved block and labels dimension-dependent correlation as
formula-owned. Non-point vector SD laws are unsupported and refuse before
nonnegative compiler mutation; scalar and deterministic point behavior remains.
The scalar `cor_scale` omission/error rule applies to multi-column scalar
correlation structures only. US/UN uses LKJ; ID/DIAG and single-column blocks
have no applicable scalar correlation.

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
Known-group covariance marginal sampling supports noncentered multi-column
ID/DIAG/US blocks with the declared group kernel and each posterior draw's
coefficient covariance; the one-column sampling stream remains unchanged.

A single grouping factor keeps its declared levels as fitted groups, including
levels without fitting rows, and the fitted term records the observed ones in
`group_observed_levels`. Prediction (conditional and marginal) routes every
unobserved level through `new_levels` exactly like a new level. Blocks with a
known group covariance are the exception: the kernel intends its declared
levels, which keep their fitted coefficients. Multi-variable groupings contain
only observed tuples.

Grouping identity uses tagged exact numeric keys (17 significant digits with
a '.' decimal mark and equivalent signed zeros) or tagged UTF8 character keys.
Display labels remain ordinary labels, disambiguated with '[<key>]' only when
needed. Prediction matches exact keys first. Numeric input can display-match
only character-owned training levels; character/factor input can match a unique
fitted label, but ambiguous raw labels refuse before the new-level policy.
The existing structured metric/index replay has its separate unique-label
fallback policy; grouping ownership does not change that policy. The current
formula schemas above require persisted component ownership and refitting
older designs.

## Verification

Use deterministic syntax and design assertions before live fitting. Follow the
profile order and centralized fixture rules in `testing.md`; do not duplicate
model fitting inside a focused consumer test.

`JAGS_with_draws()` creates a nonsampler `BayesTools_draws_view` with only
`original_fit` and `mcmc` slots. It retains the full original fit, copies
nonruntime analysis attributes and the unchanged parameter map, and clears
fit-level draw metadata and cached posterior estimates. Local view overrides
do not modify the original. Descriptive readers use supplied draws and actual
chain geometry; pooled coda draws use chain-major indices, while coda lists
retain chain timing. Extension, model convergence and bridge sampling refuse
views with classed errors before runtime/default arguments are forced. Extend
`fit$original_fit` and regenerate the view; missing view coordinates require
regeneration, while ordinary fitted-metadata errors keep their refit contract.
Uncanceled static ratio terms with subnormal stored weights use the original
basis, source value and declared multiplier ratio for draws and fixed targets.
Uncertified continuous laws from such scales are explicitly unavailable; do
not infer an exact law from draws. Genuine supplied single `mcmc` matrices keep
their own recorded timing through parameter extraction; plain matrices use
the established default timing.
