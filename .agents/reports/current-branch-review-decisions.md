# Current branch review: maintainer decisions

This file records issues found during the review whose fixes would change public
semantics, serialized metadata, numerical policy, or test/cache policy. Confirmed
implementation defects are fixed in separate commits and are not repeated here.

## D01. No-intercept factor formulas

**Issue.** For contrasts such as treatment, mean-difference, orthonormal, and
ordered coding, `~ factor - 1` currently retains a neutral fixed intercept and
uses only `n - 1` factor columns. This is not R's full-rank cell-means design.

**Impact.** Users familiar with R formula semantics may fit a lower-rank model
than intended, and comparisons with `lm()`, `lme4`, or other formula interfaces
can be misleading.

**Suggested change.** Either implement true no-intercept cell-means semantics,
including the corresponding prior parameterization, or explicitly reject and
document no-intercept factor formulas until that parameterization is supported.

## D02. Formula calls, offsets, dot expansion, and environments

**Issue.** Formula replay is based mainly on literal column names. Calls such as
`I(x^2)`, `offset(o)`, and `y ~ .` are not fully supported, and stored formula
designs do not preserve lexical environments.

**Impact.** Valid R formulas can fail late or replay differently for prediction
and bridge sampling.

**Suggested change.** Prefer one shared `model.frame()`-based construction and
replay path with explicit offset handling and a documented environment policy.
The smaller alternative is to reject unsupported calls and dot expansion early
with a precise error.

## D03. Positive support for log-intercept priors

**Issue.** A formula with `log(intercept)` can accept an explicit prior whose
support includes zero or negative values, such as an untruncated normal. Existing
tests currently rely on this broader behavior.

**Impact.** Generated syntax may evaluate `log()` outside its domain, and bridge
normalization may not describe a valid target.

**Suggested change.** Require strictly positive support recursively across
simple, mixture, spike-and-slab, and point priors. This is mathematically safer
but is an intentional compatibility break.

## D04. Cross-formula random-source dependencies

**Issue.** A row-source callback may read another formula parameter that itself
has sampled random contributions. Public and compiled reconstruction add random
effects sequentially, so changing `formula_list` order can change the likelihood.
Opaque callbacks do not expose a complete dependency graph.

**Impact.** Equivalent named formula lists can yield different marginal
likelihoods solely because of list order.

**Suggested change.** Add explicit dependency metadata and topologically order
reconstruction. A smaller policy is to reject callbacks that depend on formula
outputs with sampled random effects. Documenting order dependence alone is not
recommended.

## D05. Authority of fitted versus supplied row-source data

**Issue.** During bridge reconstruction, supplied formula data and callbacks can
replace fitted row-source data or a fitted callback without a semantic equality
check. A callback absent from the fitted object is currently grafted from the
rebuild input as a compatibility path.

**Impact.** Bridge sampling can evaluate a different likelihood from the model
that produced the posterior.

**Suggested change.** Give fitted source data and fitted callbacks authority,
reject conflicting supplied values, and only graft a supplied callback when the
fitted callback is absent. Decide whether older fits without source snapshots
remain supported through an explicit compatibility mode.

## D06. Prediction row identity beyond simple subsets

**Issue.** Mixed fitted/new prediction rows now preserve their original row
positions, but row-indexed posterior sources fundamentally identify observations
by position. Arbitrary reordering, duplication, or a new dataset with unrelated
row identity has no explicit mapping contract.

**Impact.** A syntactically valid prediction can attach fitted row-specific
scales to the wrong observations.

**Suggested change.** Store and require stable observation keys for row-indexed
sources, or restrict them to unchanged/subsetted fitted rows and reject
ambiguous prediction data.

## D07. Marginal sampling with known group covariance

**Issue.** `marginal_method = "sample"` rejects random blocks with known group
covariance even when all requested levels were fitted. The covariance method
supports the same block.

**Impact.** A documented marginalization method is unavailable for an important
supported covariance specification, and the current error refers only to new
levels.

**Suggested change.** Implement joint group draws using the fitted-level subset
of the known group covariance combined with coefficient covariance. Otherwise,
reject this combination earlier and document the limitation accurately.

## D08. Collision-free random-effect identifiers

**Issue.** Grouping interactions are encoded with `:` and some internal names
replace `:` with `__xXx__`. Real labels containing those strings can collapse
distinct tuples or produce duplicate normalized names.

**Impact.** Distinct groups can silently become one statistical group; metadata
and name-based summaries can also become ambiguous.

**Suggested change.** Separate display labels from internal tuple identities and
use a length-prefixed or otherwise reversible encoding. Rejecting reserved
characters is simpler but unnecessarily restricts valid factor levels and still
requires a migration policy for serialized metadata.

## D09. Dense random-design architecture

**Issue.** Several structured and group-local compilation paths allocate and
retain dense `n x K` model matrices before applying complexity checks. Independent
new-level sampling likewise draws every design column even when each row uses a
small subset.

**Impact.** High-cardinality factors can exhaust memory before the code reaches
the scalable structured representation.

**Suggested change.** Build sparse/indexed term metadata directly and generate
dense matrices only for explicitly requested small outputs. This is an
architectural change and should include performance budgets and compatibility
tests for stored designs.

## D10. Correlation and kernel boundary policy

**Issue.** Several covariance families distinguish mathematical support,
floating-point representability, and stable factorization differently. Examples
include LKJ endpoint matrices and CS/HCS near the lower bound. Kernel symmetry
is currently checked with a fixed `all.equal(..., tolerance = 1e-10)` rule: a
very small matrix with substantial relative asymmetry can be accepted and
symmetrized while a rescaled copy is rejected. Exact duplicate CAR coordinates
are rejected, but distinct values such as `c(0, 1e-320)` can still produce
`rho^gap == 1` and a numerically singular correlation.

**Impact.** The same conceptual boundary value may be accepted in one path,
rejected in another, or fail later during factorization.

**Suggested change.** Define one package-wide policy with separate open
mathematical bounds and representable computational bounds, plus an explicit
scale-aware symmetry rule and coordinate-separation tolerance. Decide whether
near-coincident CAR coordinates should be rejected or handled by a stable
latent-only transform. Apply the policy consistently to constructors, JAGS
syntax, initialization, reconstruction, and marginal covariance.

## D11. Fit retry and extension semantics

**Issue.** Retryable error classes, total time budgets, and whether an extension
failure should discard the last valid fit are not explicit. In addition,
`JAGS_extend(seed = ...)` cannot truly reseed already-running chains in the same
sense as a fresh fit.

**Impact.** Long fits may restart unexpectedly, exceed intended budgets, or lose
a usable result after a failed extension. The seed argument can promise stronger
reproducibility than the backend can provide.

**Suggested change.** Define a retry condition class and a total wall-clock
budget; preserve the last valid fit on extension failure unless strict mode is
requested. Deprecate or precisely document extension seeding according to
backend capabilities.

## D12. Undefined convergence diagnostics

**Issue.** Undefined `Rhat`, ESS, or MCSE values can currently be treated as
passing convergence checks in some paths. Related diagnostic helpers do not
have a consistent contract for one-observation densities, constant chains, or
models whose requested diagnostic monitor set is empty.

**Impact.** Degenerate or too-short chains may be reported as converged because
a diagnostic could not be computed.

**Suggested change.** Default undefined diagnostics to a distinct
`not_assessable` result and fail convergence unless the caller explicitly opts
into ignoring that diagnostic. Define the same result for degenerate density,
autocorrelation, and empty-monitor inputs instead of exposing low-level
bandwidth, range, or subscript errors.

## D13. Failed marginal-likelihood models

**Issue.** Some model-averaging paths convert failed or non-finite marginal
likelihoods to `-Inf`, effectively assigning the model posterior probability
zero.

**Impact.** Numerical failure becomes statistical evidence against a model and
can silently alter model probabilities and Bayes factors.

**Suggested change.** Preserve failure as `NA` with a classed condition and make
the caller choose between dropping failed models, aborting, or explicitly
treating them as zero-evidence models.

## D14. Positivity-preserving posterior sampling

**Issue.** `preserve_positive`-style fallbacks can alter sampled values to keep
them inside a requested support rather than reporting that the approximation
produced invalid draws.

**Impact.** Posterior summaries can reflect an undocumented sampling
transformation rather than the represented distribution.

**Suggested change.** Make any correction strategy explicit and opt-in. The
default should use a support-respecting sampler or fail with a diagnostic that
identifies the approximation and invalid-draw rate.

## D15. Deterministic finite-grid prior tails

**Issue.** Raw Gaussian KDE tails and exact scalar-prior point densities are now
handled without treating the sample/evaluation range as compact support.
However, deterministic `prior_linear_density` grids are finite, renormalized
numerical approximations built from tail quantiles. They retain neither omitted
tail mass nor enough source provenance to recover density or CDF values outside
the grid.

**Impact.** Treating an outside-grid value as exact support zero is wrong, but
extrapolating a deterministic grid is not statistically identified either.
Point or region Bayes factors can therefore depend on an undocumented choice.

**Suggested change.** Carry exact support plus source PDF/CDF provenance and
recompute or extend the grid at the queried value where possible. When that
provenance is unavailable, reject the query with an explicit "outside numerical
approximation range" error; do not silently assign zero or invent KDE tails.

## D16. One-sided weight-function marginal inference

**Issue.** Some exported marginal-inference paths do not implement one-sided
weight-function components even though the priors are supported for fitting.

**Impact.** A fitted supported model can fail only when users request marginal
inference or derived hypotheses.

**Suggested change.** Either implement the missing marginal distribution or
declare the combination unsupported in the exported function documentation and
validate it before computation.

## D17. Formula priors without reconstruction semantics

**Issue.** Unsupported fixed-prior classes, including some uses of
`prior_none()`, can reach formula reconstruction where the current fallback is a
zero effect. It is unclear whether such priors mean an externally defined node
or an omitted coefficient.

**Impact.** The reconstructed likelihood can silently omit a term that generated
syntax still references.

**Suggested change.** Define the meaning of `prior_none()` in formulas. If it is
external, require and validate the external node; otherwise reject it during
`JAGS_formula()` construction. Reconstruction should never silently substitute
zero for an unknown prior class.

## D18. Contrast bases and serialized factor metadata

**Issue.** Supporting a different contrast family for a factor's main effect and
interaction requires term-specific bases. Older serialized fits can also lack
the newer factor metadata, and the current orthonormal basis comes from an
eigendecomposition with platform-dependent orientation in a repeated-eigenvalue
subspace.

**Impact.** Broad contrast support can make coefficient interpretation
ambiguous; old fits may stop replaying; raw coefficients or fixtures may vary by
platform.

**Suggested change.** Decide whether contrast bases are global per factor or
term-specific. Introduce a metadata schema version and explicit legacy migration,
and use a deterministic Helmert/QR-based orthonormal basis for new schemas.

## D19. Random-effect summary ownership and labels

**Issue.** Raw monitor names, semantic random-block names, factor level labels,
and allocation parameters can overlap or use prefix relationships. Several
summary filters infer ownership from strings rather than a canonical parameter
registry.

**Impact.** Columns can be assigned to the wrong block, duplicated, omitted, or
displayed with ambiguous labels.

**Suggested change.** Persist a canonical parameter-to-block registry in fitted
metadata and make summaries/filtering consume it. Keep string heuristics only as
an explicitly tested legacy fallback.

## D20. Marginal-likelihood computation contract

**Issue.** Documentation and implementation do not fully define the reported
scale, multiple bridge repetitions, multiple chains, or the handling of
non-finite per-repetition results.

**Impact.** Callers may combine values on the wrong scale or receive unstable
results whose aggregation depends on incidental vector behavior.

**Suggested change.** Document and encode the scale in the returned object,
define deterministic aggregation and diagnostics across repetitions/chains, and
return all repetition-level results alongside the selected summary.

## D21. Cache and test-layout policy

**Issue.** The RandomEffects vignette cache checks readability and object classes
but has no source/schema/package fingerprint. Some real-fit tests live outside
`test-00-model-fits.R`; visual snapshot paths exceed common source-tar path
limits; and `.StatsVault/` ignore rules conflict with the stated coordination
file policy.

**Impact.** CI can reuse stale fitted objects, local and CI profile behavior can
diverge, source builds emit portability warnings, and coordination state may not
be reproducible across machines.

**Suggested change.** Add a versioned cache manifest derived from all generator
dependencies, decide whether the single-fit-file rule remains mandatory, shorten
snapshot naming/layout, and explicitly list which `.StatsVault` coordination
files are tracked versus ignored.

## D22. Explicit summary schema labels

**Issue.** Explicit summary schemas can override interval labels independently
of the actual probability columns. This may be intentional presentation control
or may allow a label that misstates the computed interval.

**Impact.** A published table can describe a different interval mass/type from
the values it displays.

**Suggested change.** Either derive labels exclusively from interval metadata,
or validate user-supplied labels against that metadata and require an explicit
override flag for deliberately custom wording.

## D23. Replay of legacy scaled formula designs

**Issue.** Older serialized formula designs can lack `source_data`. Prediction
then falls back to `model_frame`, but there is no metadata flag saying whether
that frame is already on the fitted model scale. Applying the stored scaling
again can double-scale predictors.

**Impact.** Predictions from an otherwise readable legacy fit can differ
substantially from predictions made when the fit was created, without an error.

**Suggested change.** Version the formula-design schema and record the scale of
every stored data field. For unversioned designs, either provide an explicit
migration rule that treats `model_frame` as model-scale data or reject replay
with instructions to refit; do not infer scale from the absence of
`source_data`.

## D24. Scale of monitored group-specific coefficients

**Issue.** `transform_scale_samples()` converts fixed-effect coefficients back
to their original predictor scale but leaves monitored
`__xREx__..._xRE_COEFx` group-specific coefficients on the fitted standardized
scale.

**Impact.** A returned posterior object can mix coefficient scales, and tables
or downstream calculations may interpret group-specific slopes as
original-scale effects.

**Suggested change.** Decide whether raw group-specific monitors are covered by
the public transformation promise. If they are, apply the same fixed-effect
unscaling matrix independently within every group and version any affected
table/sample schema. Otherwise, keep them explicitly labeled as standardized
and exclude them from APIs that promise original-scale coefficients.

## D25. Ordered priors with both atoms and continuous mass

**Issue.** For an ordered total prior that mixes point and continuous
components, each nonconstant ordered level is currently represented by one KDE.
This smears the point mass into a continuous bump because the ordered density
output has no per-level mixed-measure representation.

**Impact.** Density values and plots misrepresent the discrete probability
mass, particularly for spike-and-slab totals.

**Suggested change.** Extend ordered density output to carry continuous curves
and atoms separately for every level, and teach plotting and marginal-inference
consumers to preserve both parts. Until that schema exists, reject mixed-measure
ordered density requests rather than smoothing the atoms.

## D26. Random-formula transformation and grouping semantics

**Issue.** Inline transformed random slopes are treated as missing variables,
and the ordering of `g1:g2` grouping levels differs from base R and `lme4`.
Supporting arbitrary calls would also require a stable replay environment.

**Impact.** Familiar R random-formula syntax can fail unexpectedly or produce a
different coefficient/group order, affecting prior alignment and serialized
metadata.

**Suggested change.** Define a supported random-formula grammar and either use a
shared `model.frame()`-based implementation with explicit environments and
base-compatible interaction ordering, or reject transformations and document a
package-specific ordering. Any ordering change needs a metadata migration.

## D27. Hypothesis grammar boundaries

**Issue.** The hypothesis language does not define whether a complete relation
may be wrapped in parentheses or negated with unary `!`. Escaped parameter names
that equal R constants such as `Inf`, `NaN`, `NA`, or `TRUE` also collide with
literal-token diagnostics.

**Impact.** Expressions that are syntactically natural can be accepted,
rejected, or diagnosed differently depending on superficial spelling.

**Suggested change.** Publish a small formal grammar. Decide whether
parenthesized and negated relations are supported, and make escaped identifiers
take precedence over reserved literal tokens while retaining strict rejection
of unescaped constants where they are not meaningful.

## D28. Structured-local subset contract

**Issue.** Structured-local latent-name extraction accepts `n_groups` and
`n_columns` but currently ignores them, unlike dense random layouts. It is
unclear whether zero and prefix requests should return a subset or whether those
arguments are meaningless for this representation.

**Impact.** Internal callers can receive more posterior coordinates than
requested, and future optimizations may rely on a subset contract that the
current implementation does not honor.

**Suggested change.** Either make structured-local extraction filter by its
stored local group/column indices, including zero-size requests, or remove the
dimension arguments and expose a separate explicit full-layout method. Add
schema tests before relying on either behavior.

## D29. Exact references for stochastic fit outputs

**Issue.** A successful cache refresh regenerated all 73 fit objects and 22
marginal-likelihood objects, after which the fixture profile reported 70 exact
text-reference mismatches. Every failure was a sampled numeric value or
diagnostic (posterior summaries, model probabilities, marginal likelihoods,
ESS, or R-hat); table structure and semantic assertions passed. The current
policy forbids agents from enabling reference generation without maintainer
approval.

**Impact.** A valid fit refresh can make the fixture lane fail even when the
implementation and table schema are unchanged. Updating exact files blesses
one stochastic realization, while leaving them unchanged prevents a refreshed
cache from passing locally.

**Suggested change.** Decide which outputs are true golden values. Keep exact
references for deterministic formatting and schema, but compare stochastic
statistics with documented tolerances or invariant assertions. If the current
exact-snapshot policy is intentional, explicitly approve regenerating and
reviewing the 70 affected references from the validated cache as one controlled
update.
