# JAGS, Formula, and Fitted-Metadata Contracts

Use this guide when changing `JAGS_fit()`, generated JAGS syntax, formulas,
scaling, contrasts, random effects, posterior extraction, parameter metadata,
or marginal likelihoods.

## Keep the Fitting Paths Synchronized

Mixtures of ordered prior containers retain constructor and RNG support, but
JAGS formula binding refuses with the ordered-unavailable condition family.
Supported mixture or spike-and-slab behavior belongs on the ordered total.

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
A bridge target declares its coordinates through the stored or supplied
non-formula priors, the formula priors, and `add_parameters`; a target with none
is refused (`BayesTools_marglik_input`), and only a target whose declared
coordinates are all point-fixed is evaluated exactly. For a fit that stores its
priors, a supplied non-formula prior must have the semantics of the stored prior
of the same name (`BayesTools_refit_required` otherwise, as for formula priors);
supplied names the fit does not store, and stored names omitted from
`prior_list`, stay allowed because callers may integrate omitted parameters in
their callback.
For a block fitted with `random_effects_compile(marginalized = ...)`, the SD
coordinates and priors stay in the bridge target, but under `bridge_context =
FALSE`/`"none"` (the default) the likelihood callback receives neither the SD
values nor the covariance. It reconstructs them from primitive coordinates
(`add_parameters` or its own evaluator), or requests `"marginal"` or `"full"`
(`TRUE`); the bridge does not refuse the default, because callers (RoBMA) supply
their own reconstruction. A refusal with an opt-out is a maintainer decision.
Bridge sampling checks affine rank of varying finite supplied coordinates before
bound transformations, then checks the transformed coordinates separately with
the existing QR tolerance and diagnostics. Constant columns retain the existing
sampler policy. Before each check, finite nonzero columns are rescaled by their
maximum absolute value in a local diagnostic copy to avoid variance overflow
or underflow; actual draws and sampler bounds remain unchanged. This is not a
general nonlinear dependence detector.
Both refusals use `BayesTools_marglik_rank_deficiency` (parent
`BayesTools_marglik_input`), with `rank`, `n_coordinates` after finite/varying
filtering, `n_draws` including duplicate rows, and complete available
`dependent_coordinates` from the QR-pivot tail (otherwise `character(0)`).
`reason = "insufficient_draws"` records `n_draws <= n_coordinates`;
`"dependent_coordinates"` records deficient numerical sample rank otherwise.
These are nonexclusive observations, not certified structural causes or
guaranteed remedies; pivot names are basis-dependent sampled redundancy.
Marginal-likelihood `success` records finite-estimate availability; finite
repetitions remain aggregated even when the known total iteration count
exceeds `maxiter`. Optional schema-1 field `iteration_limit_exceeded` is TRUE
for any known excess, FALSE when all counts are known without excess, and NA
when counts are absent/partly unknown without known excess (also for manual
and exact results). An absent historical field means unavailable telemetry.
Legacy repetition `within_maxiter` remains TRUE for an unknown count, and its
`success` remains finite AND within_maxiter. No field certifies convergence:
upstream totals may include successful retries. Known excess always warns
with `BayesTools_marglik_iteration_limit` (parent
`BayesTools_marglik_diagnostic`; affected repetition rows, maxiter and upstream
warnings) independently of `silent`, which controls progress. Flagged results
print a visible diagnostic; scalar extraction/nonfinite policy stay unchanged.
A `NA` from the bridge likelihood callback is no density: at a posterior draw
the sampler route stops (`BayesTools_marglik_callback_na`, parent
`BayesTools_marglik_input`, with the number of draws), and at proposal draws
it is evaluated as `-Inf` with a warning of the same class that reports the
share of proposal draws. The counts are read from the warnings of
`bridgesampling` (it replaces the `NA` itself), which a regression test pins
against the installed version; the exact route and `NA` priors keep their
refusals. The policy covers `method = "normal"`: `bridgesampling` stops on an
`NA` with `method = "warp3"` with an unclassed upstream error (documented, and
pinned by a test), so a `warp3` callback must not return `NA`.

Every deterministic node BayesTools generates (allocation-derived SDs, scalar
and LKJ correlations, publication weights, the experimental p-hacking nodes,
spike-and-slab and mixture
parameters, the coordinates of ordinary Dirichlet and multivariate-t priors,
formula linear predictors) belongs to one registered family whose
node specification emits its JAGS syntax and carries the R evaluator and the
declared dependencies: `JAGS_formula()` writes the linear predictor from the
node of the fitted design, and `selection_backend_spec()` writes the weights
from the `omega` node specification and the p-hacking constants, null masses,
compositions and `sel_vector_rule` from the `phacking` node specifications (one
evaluator, `.bt_dnode_phacking_quantity()`, serves the registry, bridge sampling
and the marginal-likelihood parameters; a prior without a p-hacking branch
registers no such node, and `selection_backend_spec()` writes its constant
`sel_vector_rule` with the same helper). The `names` and `global_breaks` of
`selection_backend_spec()` change only the syntax it emits: the evaluators
(`JAGS_marglik_parameters()`, `JAGS_marglik_priors()`, bridge sampling) read the
default names and the grid that `JAGS_fit()` emits.
`JAGS_deterministic_nodes()` lists the nodes;
`JAGS_deterministic_evaluator()` resolves them once into an evaluator
of draws (repeated or one-draw evaluation), and `JAGS_evaluate_deterministic()`
is that evaluator applied once. A prepared evaluator may locate a family's
columns once per set of draw column names, but shares the family's lookup
rules, arithmetic, and checks (values identical to the family evaluator).
Prior draws, catalog quantities, bridge and marginal-likelihood
parameters, prediction, marginal posteriors of formula parameters,
random-effect unscaling, and convergence-role parents use the family
evaluators. Declared dependencies are the coordinates the R evaluator reads. A
new generated node gets a family and a parity test against the JAGS monitors,
not another evaluator. Bias-term specs retain their declared branch count;
finite indicators must exactly index that declaration before inactive zeros
are constructed. Missing or inconsistent fitted declarations require refitting. Intercept priors cannot carry `multiply_by`; it scales
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

Coefficient-transform target schemas include the existing `state_dependent`
map type. Fixed exponential outputs must be finite and strictly positive;
underflow remains unavailable rather than a structural zero. Manipulated
continuous marginal grids use `c(-1, 0, 1)` in fitted/SD or original view units.
Public linear-target coefficient/threshold range loss retains the existing
typed hypothesis numerical refusal, without rounding or clamping.

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
a coefficient count; fitting stops on incomplete factor metadata. Refusals of
`JAGS_formula()` inputs that callers may need to tell apart are
`BayesTools_formula_transform_unavailable` errors with a `reason`: a factor with
a single level in the data, which has no contrasts to code
(`single_level_factor`, raised before base R's contrasts error; structured
`cs`/`hcs`/`ar1`/`har` indexes may have one level), a logical predictor that
`model.matrix()` would code by indicators in a factor interaction, fixed or
random (`logical_factor_interaction`), a random-only factor with a contrasts
matrix or an unsupported contrast name (`unsupported_random_factor_contrast`),
and a homogeneous scaled block (`homogeneous_scaled_block`, below).
`transform_scale_samples()` reports a `formula_scale` entry that omits a
standardized predictor of the fit as `incomplete_formula_scale`. Formula
evaluation and prediction (`JAGS_evaluate_formula()`,
`JAGS_predict_formula()`) go only through the stored design: draws without a
fit get it from `JAGS_formula_draws()`, which builds it as `JAGS_fit()` does;
never evaluate from prior-list factor metadata. The `log(intercept)` flag of an
evaluation is the stored flag of the design, also for an explicit formula
without the attribute; an explicit attribute that disagrees with it stops with
`BayesTools_formula_transform_unavailable` (reason `log_intercept_mismatch`; as
`marginal_posterior()` requires them to agree).
Formulas stored in designs and
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

Formula-design 6, unscale-design 2 and coefficient-transform 4 require refitting
older formats. One private finalizer completes current compiler specs from
actual full priors and owned scalar data for fit/draws-only producers, writing
the same owner to both carriers; fit readers never complete missing ownership.
`JAGS_formula_point_finalize()` and `JAGS_formula_scale_finalize()` are the
validated public wrappers of the point and scale finalizers for downstream
packages that rebuild a fitted design; the expression finalization and the
formula-data merge of `JAGS_fit()` stay private.
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

Fixed and random original-scale maps expand each numeric monomial into its
declared lower-order subsets. Exactly zero means never drop a variable;
every nonzero induced coefficient is retained, however small. Unscaled
logicals remain categorical and scaled logicals are numeric. Categorical
functions use the persisted coding on the existing synthetic grid. A complete
positional/name map is accepted only when every induced component has one
identically coded partner; this preserves redundant independent-factor maps.
Otherwise all numeric-subset categorical groups must have full rank before
any missing-term refusal. Solve first in a complete source-variable hierarchy
basis, then a certified 0/1 indicator-partition lift when available, otherwise
the full categorical group. The lift requires complete declared source cells
crossed once with each observed outside cell, without requiring unneeded
outside logical interactions. Excluded basis columns are structurally zero;
computed within-basis nonzeros are never snapped. Normalized categorical
square solves verify all grid rows with the floating backward-error bound;
numeric affine coefficients are applied afterwards. Whole base `prod()` and
safe ratio regrouping preserve representable nonzeros or raise
`nonfinite_transform`, including on the direct random name-paired fallback.
Fixed missing lower-order effects refuse even when their coefficients are tiny.
Eigen-generated contrasts do not imply exact stored-double orthogonality;
tiny within-basis coefficients can have poor relative accuracy in doubles.
Parameter-map 16 and coefficient-transform 4 invalidate older classifications;
refit and recreate detached derived objects. Design 6/unscale-design 2 remain.

Random-effect SDs and correlations of a block with a scaled slope follow the
block's own design, as the fixed coefficients do: when the standardized block
columns equal the original-scale columns times a matrix `A` (verified on the
synthetic grid of `.bt_formula_unscale_design_transform()`), the coefficients
transform as `u_o = A u_s` and the covariance as `A Sigma A'`. A factor that
enters the block as a main effect and in an interaction with the scaled slope
(`diag(1 + f + f:x | g)`) is representable, so the induced shifts of the
intercept and of the factor levels are reported. Every route that has the fit
passes its fitted formula design (`transform_scale_samples()`,
`transform_prior_samples()`, estimates tables, mixed posteriors, catalog
random summaries); the block is found by its block name and SD coordinates.
A block with no such matrix in its own columns (`(0 + x | g)`, `diag(1 + f:x | g)`,
`us(0 + f + x:f | g)`: a term would have to be induced) keeps the documented
name-paired map that drops it only when every structurally nonzero component
is retained once with identical categorical coding or is unrepresentable by
its target numeric group. The existing categorical solver checks omitted
components when that group exists; exact-zero means induce no component.
Mixed-coding or omitted representable contributions refuse the whole
original-scale block with `original_scale_not_representable`, rather than
returning a partial map. Only the design transform's
`original_scale_not_representable` refusal selects the name-paired fallback;
a refusal from its certificate is final, and all other failures propagate.
For a scaled slope that interacts with a factor whose level
intercepts are not columns of the block (`diag(1 + f:x | g)`), the induced
level-specific intercept shifts are not block columns: the original-scale SD
of the intercept is the SD of the group effect at the predictor mean, not the
SD of the intercept of a level, while the slope SDs are exact.
`JAGS_estimates_table()` states this in a footnote of the intercept SD row, and
the `transform_scaled` roxygen of the summary functions documents it. The name-paired map takes each column as its own coordinate by
position: the indicator columns of a logical predictor's interaction with a
scaled slope share one SD leaf and put no weight on each other, and the leaf's
SD is their common per-column value. Hand-built SD-leaf metadata without the
fitted design pair names.
A block whose design is linearly dependent (`diag(1 + x + x:f:f2 | g)`: a
three-way term without its margins) has no unique original-scale map
(`original_scale_not_identified`). The fit and its parameter map are never lost
to it: the block's SDs are composite catalog quantities without an
original-scale law, and every request for original-scale values stops with
`BayesTools_formula_transform_unavailable`, reason
`original_scale_not_identified` (`transform_scaled = TRUE` tables,
`transform_scale_samples()`, `parameter_draws()` and the mixed posteriors of
those SDs, the prior law of an SD on the original scale). The estimates table
without `transform_scaled` leaves the block's original-scale random-effect rows
out and says so in its warnings; `random_effects_summary = "raw"` shows the
fitted SDs. An SD shared by several columns of a block (the indicator columns
of a logical predictor's interaction with a scaled slope) has one original-scale
value only when its columns transform to the same SD, which a correlated block
never guarantees: otherwise original-scale requests stop with
`BayesTools_formula_transform_unavailable`, reason
`original_scale_not_representable`, and the estimates table without
`transform_scaled` leaves the block's original-scale rows out with a warning
in the same way.
For a diagonal block, shared-column equality compares variance weights on
each source SD leaf separately with the existing relative tolerance 1e-12.
Zero/nonzero weights differ exactly, however large another leaf's weight is;
each leaf has its own independent SD scale. Common tiny weights and zero/zero
pairs retain their equality, and correlated blocks retain their refusal.
The ensemble and marginal tables with `transform_scaled = TRUE` take the block
design from the formula state of the mixed posteriors (`as_mixed_posteriors()`;
`mix_posteriors()` attaches a state whenever a requested parameter is a fixed
coefficient or a random-effect SD of a fitted formula, and none for a mixture
of SDs alone that includes a model without a fitted formula design), one design
per model of a mixture. Samples without it unscale a block exactly name-paired
as before, and refuse the SDs of any other block (a scaled column without its
lower-order partners, or two scaled columns that share one term name, as the
indicator columns of a logical predictor's interaction with a scaled slope) with
`BayesTools_formula_transform_unavailable`, reason
`random_effect_design_unavailable`, not with the values of the pairing.

The prior density and prior mass of a random-effect SD of
`as_mixed_posteriors(transform_scaled = TRUE)` are those of the original-scale
SD, never the fitted-scale law. The stored prior-density context carries the
fitted block design (`formula_design`: the `.bt_formula_design_random_unscale_view()`
of each parameter) next to `formula_scale`, and its conditional components
inherit it. An SD whose row of the block's original-scale map has one nonzero
entry, on its own coordinate (`.random_sd_original_scale_laws()`; every column
of a shared leaf must give the same entry), is a pure rescaling: its law is the
fitted SD's law rescaled by `|c|` (exact for every prior kind, applied as the
SD's weight in `.prior_density_context_standardized_weights()`). Any other SD of
a block with a scaled column (the intercept SD with a non-centered scaled slope,
factor-level SDs that absorb induced shifts) combines several fitted SDs, and so
does an SD with any further nonzero entry, however small beside its diagonal:
small predictor units alone can make an induced entry consequential, and SD
priors can also differ in scale by orders of magnitude. The exact zeros of the
map that the draws use
(`.random_sd_unscale_matrix()`) are the one zero rule
(`.random_sd_unscale_row_nonzero()`): the random-effect wrapper returns the
design transform or documented name-paired fallback unchanged. Draws, the prior
law of an SD, the catalog support of a composite SD (which SDs it reads), and the
prior law of a correlation read the same retained rows. The shared design
helper uses the structural affine/categorical representation above; the random
wrapper adds no magnitude threshold. An exactly zero predictor mean keeps a
pure intercept rescaling. Centering by subtraction can leave a numerically
nonzero mean: if its induced entry is retained, the intercept SD is composite,
its transformed draws remain available, and its prior law is refused. The prior
density of a composite is unavailable with cause `random_sd_prior_law_unavailable`, and with cause
`random_effect_design_unavailable` when the context holds no design of its block
(an older object, or `plot_transformed_prior()`, which has no fit). SDs of blocks
without a scaled column and of structured blocks keep their fitted law. The
refusal is raised where the SD's weight is evaluated, never when a context over
fixed columns is built, and `as_mixed_posteriors()` marks the `prior_density`
(and the support the same evaluation refuses) of the SD's element, so
`hypothesis_BF()`, `marginal_posterior(prior_samples = TRUE)` and
`plot_posterior(prior = TRUE)` refuse with
`BayesTools_formula_measure_unavailable`; the draws, the estimate tables and
`plot_posterior(prior = FALSE)` do not read the prior law. The prior draws of
the level comparisons unscale the SDs with the context's design.

Independence is also specified on the centered scale. `(1 + x || g)` or
`diag()` with scaled `x` has independent `u0`, `u1` (SDs `t0`, `t1`); the
original-scale intercept `u0 - u1 m / s` and slope `u1 / s` then have
correlation `-(t1 m / s) / sqrt(t0^2 + (t1 m / s)^2)` (about -0.93 for equal
SDs and `m / s = 2.5`), which such blocks do not report. A homogeneous-SD block
(`id()`, `covariance = "id"`, `diag(hom = TRUE)`) with more than one column and
a scaled-predictor column has no original-scale SDs (one shared SD, different
column SDs) and is refused by `JAGS_formula()`, whether its shared SD is
sampled or fixed by a point prior, with `BayesTools_formula_transform_unavailable`
(reason `homogeneous_scaled_block`; the remedy is a heterogeneous structure or
an unscaled slope, so it holds inside `JAGS_formula()` as well as for the
unscaling of hand-built metadata); single-column homogeneous blocks are
unaffected.

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
  The mapping of a random-effect SD or variance comes from its own persisted
  covariance transformation, not from the coordinates its block evaluator
  reads: in a scaled block an SD whose row of the original-scale map
  (`.random_sd_original_scale_laws()`) has one nonzero entry, on its own
  coordinate, is a one-to-one transform of that coordinate with the entry as
  its scale (also when the block has other SDs), and any other SD is a
  composite.
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
  never by floating-point equality of design rows. In a random interaction, a
  factor that R codes by level indicators (the stored coding of the term: its
  margin is absent from the block, as in `f:g` without `f` or `g`, or `f:x`
  without `x`) counts as independent whatever its declared contrast, so the
  columns of an interaction-only term are level cells (`f[a]:x`), while a
  factor coded by its eigendecomposition contrast keeps the tag gate (labels
  follow the coding, not the declared contrast tag). Every other coordinate is
  contrast coefficient `<parameter>{j}`, labelled so in tables, diagnostics,
  mixed-posterior columns, and random-slope components (`sd(g{j})`). A
  `{j}` selector of a direct level coordinate is refused (in catalog
  resolution and in `hypothesis_parse()` with a catalog) with class
  `BayesTools_selector_unavailable`, naming the level in the same form. JAGS
  coordinate names such as `mu_g[2]` remain backend column names, used by
  coordinate-based functions such as `JAGS_materialize_draws()`, but are not
  factor selectors.
  The SD columns of a random interaction whose factors are all contrast-coded
  (meandif, orthonormal) and whose margins are present in the block
  (`diag(1 + f * x | id)`, `diag(1 + f * g | id)`) are contrast coefficients of
  the term, `<term>{j}` (`sd(f:x{1})`, `sd(f:g{2})`), in the mixed-posterior
  columns (`mu__xREx__id_f__xXx__x{1}`) as in the catalog; an interaction with
  an independent or treatment factor, or with a factor whose margin is absent
  (`diag(1 + f:x | id)`), keeps its level cells.
  The indicator columns of a logical predictor's interaction with a continuous
  one whose margin is absent from the block (`flagFALSE:x`, `flagTRUE:x` of
  `(1 + flag + flag:x | id)`) share one SD, labelled by the term (`sd(flag:x)`),
  and are labelled by their level cells (`flag[FALSE]:x`, `flag[TRUE]:x`) in
  the correlations, so that every correlation name identifies one pair of
  columns.
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
  supported on `[0, Inf)`; the SDs, variances, totals, and variance shares of
  allocations with inclusion gates or a scale prior with a point component list
  the atoms of their own prior density as `points` of a `"mixed"` support,
  derived from the same measure as the density; where such a chain has no
  derivable atoms (an external SD source without a BayesTools prior, a
  composite SD outside its allocation component), the hull is declared inexact
  instead of an exact atomless interval; the original-scale correlation of a
  scaled block whose chain has gates or point-mass SD components has atoms at
  -1 and 1 that the catalog does not place, so its `[-1, 1]` is inexact too, while
  a fitted-scale correlation, undefined where an SD is zero, keeps an exact
  support) and `definedness` when the map is built;
  downstream packages read
  them instead of deriving supports.
- Draws that can be undefined are declared, never inferred from names or
  values: `parameter_draws()` sets the `undefined_draws` draw metadata of its
  `mcmc.list` from the catalog's `definedness` (canonical name to reason:
  `"correlation"` for random-effect correlations (every LKJ correlation,
  fitted- or original-scale, and the original-scale correlations of a block
  with scaled columns, are undefined where an SD is zero; a correlation whose
  block has no readable SD draws, because its external SD source is
  row-indexed or not monitored, keeps its fitted draws, since nothing
  establishes where it is undefined),
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
  these nodes of a spike-and-slab or single-slice mixture total. Multi-slice
  spike-and-slab totals share one inclusion probability and indicator, as in
  the fitted model; general mixture totals require one slice. Joint ordered
  prior draws reuse the first allocation under each declared shared `(id, factor)`
  key, while consuming every original per-prior RNG stream. Standalone prior
  calls remain independent. Latent
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
- The `fitted_scale` of a fixed coordinate is `fitted_standardized` only when its
  design involves a standardized predictor (the intercept of a formula with a
  standardized predictor, and terms with a standardized component); every other
  coefficient keeps `fitted_original`, also in a formula that standardizes
  other predictors or has formula-scale metadata for another reason. The nodes
  that a coefficient prior declares besides the coefficient (the indicator,
  inclusion probability and slab of a spike-and-slab prior, the indicator of a
  mixture) take the fitted scale of their coefficient, found from the prior's
  declared nodes and never from their names.
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
  They reject block contrast overrides. `cs()`, `hcs()`, `ar1()`, and `har()`
  reject intercept controls; `car()` has no intercept, so it rejects an explicit
  `1` and accepts the redundant `0` and `- 1`. Discrete indices
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
the same resolved block and marks a correlation that depends on the formula
block's dimension "(when the formula block has a correlation parameter)".
Non-point vector SD laws are unsupported and refuse before nonnegative compiler
mutation; scalar and deterministic point behavior remains.
The scalar `cor_scale` omission/error rule applies to multi-column scalar
correlation structures only. US/UN uses LKJ; ID/DIAG and single-column blocks
have no applicable scalar correlation.

US/UN uses `LKJ(1)`; CS/HCS uses raw `Uniform(-1 / (K - 1), 1)` whose JAGS hull
is closed while evaluation treats the CS singular bound `-1/(K-1)` as open;
AR1/HAR uses raw `Uniform(-1, 1)`; CAR uses raw `Uniform(0, 1)` and includes
`rho = 0`. Do not treat Uniform as an open-interval family. Explicit scalar
priors retain
their Fisher-z default scale. In double precision a Fisher-z coordinate beyond
about 19 (a logit coordinate beyond about 36.7) in absolute value maps to
`|rho| = 1`: centered fits stop in JAGS ("Invalid parent values"), noncentered
draws are refused when the correlation is reconstructed, and a CAR prior whose
whole transformed support saturates is refused at compile time; narrower
correlation priors avoid it. The compiler applies no cutoff or truncation. SD magnitude belongs to the outcome model:
BayesTools must not supply a generic scale, and direct `JAGS_formula()` use
requires an SD prior, SD source, or variance allocation.

Random-effect catalog names use `(formula) owner: quantity(arguments)`.
Parentheses contain coefficient or parameter names; square brackets contain
factor or index levels, and curly braces contrast coefficients that are not a
level (`sd(g{j})`). Use public `cor`, while any compact backend `rho`
coordinate remains internal. Total-variance allocations expose `sd_total`,
`var_total`, and `var_prop(...)`; mean-variance allocations expose `sd_common`,
`var_common`, `var_mult(...)`, and `sd_mult(...)`.
Inclusion quantities follow the same labels: `inclusion(sd)` for the SD shared
by all columns, `inclusion(sd(f3:h3))` for the SD prior of a term (the term of
the prior's own name, never a backend `__xXx__` name or a list of cells),
`inclusion(sd(x)[component])` for a mixture component, and, for an allocation
component, `inclusion(<component>)` with the allocation's public
`component_names`, as `var_prop(<component>)` has them.
Gate-only roots are slab-plus-inclusion containers and do not publish a root
`sd_total`; the child-split total is the public realized total.
Conditioning a gated total on presence ANDs the binding-factor (parent) gates
that define a positive realized source, not every nested component gate. A
total of an allocation without parent gates whose components all have
inclusion gates is conditioned on the OR of the component gates (at least one
active); one whose components all have inclusion gates and which also has
parent gates combines both kinds of gates and is refused
(`parameter_mixed_posterior()` stops with "combines parent-allocation and
component inclusion gates"); a total without parent gates and with an ungated
component has no inclusion gate and is not conditionable.

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
One conditioning rule holds for every structure (plan field
`held_fixed = "allocation_weights"`): setting an allocation-derived SD to `v`
holds the allocation weights and every non-allocation coordinate fixed and moves
the allocation scale `tau` (`SD_k = tau * m_k(w)`, so `tau' = v / m_k(w)`, with
`m_k(w)` the product of the square roots of the weights on the path from `tau`
to the SD), the conditional of the variance prior declared on the independent
`(tau, w)`: the moved coordinate has the `tau` prior with Jacobian `1 / m_k(w)`
for an SD value and `1 / m_k(w)^2` on the variance parameterization of the scale
for a variance value. The scale reaches every block whose SD binding uses it,
through nested allocations (a split root with one `<block>_components`
allocation per block) or a shared `random_sd_source()`; the plan's `blocks`
lists all of them, the selected block first. One reach computation serves every
quantity that is the allocation scale times weights: the component `sd`/`var`
and the aggregates `sd_common`/`var_common`/`sd_total`/`var_total` of a nested
allocation (or of one with a shared source) move the root scale of their
allocation chain (`allocation$source`) with every weight fixed, so their plans
list the blocks of the sibling allocations as well. A reached block whose
columns are all driven by `tau` has its covariance multiplied by `c^2`,
`c = v / SD_k`, on every covariance structure, so the plan is the affine `scale`
plan (`A + h(v) B`) of the summed contribution of `blocks` whenever the scale
drives every column of every reached block, when the public SD is on the fitted
coefficient scale. If it drives only part of the columns of the one block it
reaches, the plan of a component SD is the factor `column_scale` plan with
`driven_columns`, and the grid carries `column_factor` for them, provided the
selected column is among the driven ones (a column the scale does not drive
would not move: `unsupported`, reason `selected_column_not_driven`); a scale that
reaches several blocks and only part of their columns, or part of the columns of
a block for an aggregate (which has no selected column), has no exact scalar
update (`unsupported`, reason `partial_scale_reach`). Downstream priors for the
update (RoBMA's IWMDE prior side) use the same conditional and move all reached
columns. A direct per-column SD keeps `column_scale` with its siblings fixed and
has no `held_fixed`.
Formula-scaled allocation-derived SDs, and SDs that depend on several fitted
coordinates, have no exact single-coordinate update and are `unsupported`.
Do not make a downstream package reconstruct these quantities from raw JAGS
columns.

Preserve explicit new-level policy and the distinction between full covariance
and `diagonal_only` output. Structural unavailability is not permission to
silently substitute zero covariance or a diagonal approximation.
Posterior-indexed row SD sources require explicit `fitted_rows` for supplied
prediction data. Only the unavailable new-row capability branch uses
`BayesTools_prediction_new_rows_unavailable` (parent
`BayesTools_prediction_unavailable`), with reason
`posterior_indexed_row_source`, block, displayed source label, affected
prediction-row indices and caller context. Evaluable source callbacks return
before this branch. The refusal applies only to new rows whose source is read,
under `new_levels = "sample"`: rows zeroed under `"zero"` read no source and
their `fitted_rows` entries are validated but ignored, so callers pass
`new_row & <policy is "sample">` to the fitted-rows check. Preserve
missing/malformed mapping validation precedence
and existing invalid-source and missing-monitor conditions; those do not mean
new-row capability is unavailable.
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
