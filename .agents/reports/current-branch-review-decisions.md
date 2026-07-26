# Current branch review: maintainer decisions

This file records issues found during the review whose fixes would change public
semantics, serialized metadata, numerical policy, or test/cache policy.

The maintainer's `Decision` text is preserved below. The audit notes added after
each decision distinguish an accepted instruction from work that is already
present on the current branch. In particular, **decision confirmed** does not
mean **implemented**.

## Audit status index

| Case | Current status | What remains |
|---|---|---|
| D01 | Implemented and verified | No behavior change; retain and document BayesTools no-intercept semantics |
| D02 | Decision confirmed; implementation pending | Add an early supported-formula grammar check |
| D03 | Decision confirmed; implementation pending | Add recursive strict-positive prior-support validation |
| D04 | Decision confirmed; implementation pending | Reject callbacks that depend on formula outputs with sampled random effects |
| D05 | Decision confirmed; implementation pending | Remove bridge-replay callback grafting and legacy source fallbacks |
| D06 | Decision confirmed; implementation pending | Add the explicit fitted-row-index contract for posterior row sources |
| D07 | Decision confirmed; implementation pending | Add joint marginal draws under known group covariance |
| D08 | Implemented and verified | No remaining work |
| D09 | Decision confirmed; implementation pending | Add a configurable 16 GiB hard ceiling on estimated peak allocation |
| D10 | Resolved by NF10/NF18 | Remaining CAR representability issue is isolated as D31 |
| D11 | Partially implemented | Time budget already resets; remove `seed` and preserve the last valid fit on error |
| D12 | Decision confirmed; implementation pending | Add three-state diagnostics and explicit downstream monitor selection |
| D13 | Resolved by NF08 | No remaining decision |
| D14 | Resolved/superseded by NF09 | The audited code changed mixture counts, not posterior draw values |
| D15 | Resolved by NF05/NF06 | No remaining decision |
| D16 | Partially implemented | General one-sided weight-function marginals remain unimplemented |
| D17 | Decision confirmed; implementation pending | Canonicalize formula `prior_none()` to a point mass at zero |
| D18 | Decision confirmed; implementation pending | Separate fixed/random-block contrast scopes and retain the current bases |
| D19 | Decision confirmed; implementation pending | Introduce and adopt a canonical, documented parameter registry |
| D20 | Decision confirmed; implementation pending | Introduce a scalar, scale-explicit marginal-likelihood result contract |
| D21 | Partially implemented | Add the vignette manifest and consolidate real fitting; D30 covers paths |
| D22 | Decision confirmed; implementation pending | Make `interval_level` derived output, not caller input |
| D23 | Implemented and verified | No remaining work |
| D24 | Decision confirmed; implementation pending | Keep raw latent/group coefficients internal and add a transformed extractor if needed |
| D25 | Partially implemented; implementation pending | Preserve atoms in the ordered sampled-density fallback |
| D26 | Implemented and verified | No remaining work |
| D27 | Decision confirmed; implementation pending | Implement and document the proposed restricted grammar |
| D28 | Decision confirmed; implementation pending | Honor the existing subset arguments |
| D29 | Decision confirmed; implementation pending | Apply the agreed long-term stochastic-reference policy |
| D30 | Decision confirmed; implementation pending | Shorten the 50 nonportable snapshot paths |
| D31 | Decision confirmed; implementation pending | Use stable CAR recurrences and reject only unrepresentable innovations |
| D32 | Decision confirmed; implementation pending | Reject nonzero `meandif`/`orthonormal` locations |

## Maintainer decisions after the second pass

All six previously open policy choices are now confirmed:

| Case | Confirmed choice | Decision |
|---|---|---|
| D06 | How posterior row-indexed sources identify prediction observations | Require explicit fitted-row indices whenever prediction data are supplied; allow new data only for callback-computed row sources |
| D09 | Whether the proposed 16 GiB default is a warning or a stopping limit | Make 16 GiB a configurable hard ceiling on the conservative estimated peak working set; a warning alone cannot prevent an out-of-memory failure |
| D18 | The scope of a factor contrast and whether to replace the current eigenvector bases | Scope contrasts separately to the fixed design and to each random block; store the concrete matrices; do not replace `contr.orthonormal()` or `contr.meandif()` in this release |
| D22 | Whether callers may supply `interval_level` independently of endpoints | Remove it as input metadata, retain it as derived output `upper_prob - lower_prob`, and use an unnumbered interval label when endpoint probabilities are unknown |
| D26 | Ordering of grouping interactions such as `g1:g2` | Adopt base-R/lme4 ordering now and persist the component tuple-to-index map; require refitting rather than migration |
| D32 | Whether a scalar nonzero location is meaningful in an arbitrary contrast-coordinate basis | Require an exact zero location recursively for `meandif`/`orthonormal` priors; use a future level-space prior API for directional nonzero means |

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

Decision: we have a consistent and specific handling of no intercept models "-1" or "0" in BayesToosls
- intercept is actually set to 0 via a prior distribution which allows consistent model comprisons and modeling
this indeed deviates from R standard formula interface but it is needed for simple Bayesian modeling
you can improve the documentation but do not change the behavior (if you see anything inconsistent with this handling, raise a new issue)

**Audit status: implemented and verified.**

**Review response.** Instruction understood. The current implementation
canonicalizes `-1` and `0` formulas through `formula_add_intercept()` and inserts
a neutral point prior for the retained intercept. For an ordinary intercept the
point is zero. For `log(intercept)` the point is one on the parameter scale, so
its linear contribution is still `log(1) = 0`; this is the required neutral
special case rather than an inconsistency.

The nonstandard BayesTools behavior is documented in `JAGS_formula()` and
`formula_add_intercept()`, including preservation of arithmetic inside calls
such as `I(x - 1)` and `offset(x - 1)`. Unit/design-oracle tests cover the
ordinary and log-intercept cases. I found no additional inconsistency requiring
a new issue.

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

Decision: reject unsoported calls early, in later versions we might fully implement this but it would be a large change right now

**Audit status: implemented.**

**Review response.** Instruction understood. The current fixed-formula path
still reaches `terms()`, literal-variable checks, or `model.frame()` before it
rejects several unsupported constructs. The implementation should add one
shared, early grammar validator used by fitting, prediction replay, and bridge
reconstruction. It should explicitly reject at least:

- dot expansion (`.`);
- offsets;
- inline calls such as `I()`, `poly()`, `scale()`, and arbitrary functions;
- references that require a preserved lexical environment.

BayesTools' existing explicit `expression(...)` facility is a separate supported
feature and should continue to use its current validation path. Errors should
identify the offending call and state the supported alternative.

**Implementation outcome.** Formula construction, prediction replay, public
marginal reconstruction, and compiled bridge reconstruction now share an early
fixed-formula grammar validator. It rejects dot expansion, offsets, inline and
namespace-qualified transformations, and arbitrary calls with the offending
expression and a supported alternative. Literal columns, standard formula
operators, and `expression(...)` remain supported. The pass also fixed a
pre-existing crash when the expression walker encountered namespace-qualified
call heads. Focused tests passed 1,459 assertions, and the complete unit profile
passed 7,765 assertions with no failures or warnings.

## D03. Positive support for log-intercept priors

**Issue.** A formula with `log(intercept)` can accept an explicit prior whose
support includes zero or negative values, such as an untruncated normal. Existing
tests currently rely on this broader behavior.

**Impact.** Generated syntax may evaluate `log()` outside its domain, and bridge
normalization may not describe a valid target.

**Suggested change.** Require strictly positive support recursively across
simple, mixture, spike-and-slab, and point priors. This is mathematically safer
but is an intentional compatibility break.

Decision: add the required support check

**Audit status: implemented.**

**Review response.** Instruction understood. The current prior-class validator
does not prove strict positivity for `log(intercept)`, so an untruncated normal
still passes construction. The new check must recurse through simple, mixture,
spike-and-slab, and point priors. Every component with nonzero prior probability
must have support strictly inside `(0, Inf)`; a point mass at zero is invalid.
The same validation must be used by formula construction and replay/bridge
compilation so the target cannot differ between fitting and reconstruction.

**Implementation outcome.** Formula construction, prediction replay, public
marginal reconstruction, and compiled bridge reconstruction now share a
recursive strict-positive-support validator. Continuous priors truncated at
zero are accepted because they have no atom there; point and discrete priors
must exclude zero. Spike-and-slab priors are accepted only when inclusion is
fixed at one. Focused tests passed 1,437 assertions, and the complete unit
profile passed 7,743 assertions with no failures or warnings.

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

Decision: please add a bit more information -- an actual example -- so I understand the impact and potential issues

**Audit status: decision confirmed; implementation pending.**

**Concrete example.** Suppose `formula_list` contains two modeled parameters,
`mu` and `log_sigma`, and both have sampled random effects. A row-source callback
used while reconstructing `mu` reads `parameters$log_sigma` to calculate a
row-specific scale:

```r
mu_source <- parameter_source(
  "mu_scale",
  shape = "row",
  values = function(parameters, data, n_rows) {
    exp(parameters$log_sigma)
  }
)
```

Reconstruction currently processes formula parameters sequentially. With list
order `mu, log_sigma`, the callback sees `log_sigma` after its fixed contribution
but before its sampled random contribution. With order `log_sigma, mu`, it sees
the complete `log_sigma`. The named model, posterior row, and data are otherwise
identical, yet the reconstructed `mu` and bridge likelihood differ solely
because the list was reordered.

Opaque R callbacks do not expose a trustworthy dependency graph. The two
defensible choices are:

1. require callbacks to declare `depends_on`, validate it, and reconstruct
   formula parameters in topological order; or
2. as the smaller release change, prevent row-source callbacks from reading
   formula outputs that have sampled random effects.

**Recommendation.** Use option 2 for this release. It has a smaller API and
eliminates order-dependent likelihoods. Add explicit dependency metadata later
only if cross-formula random dependencies become a required feature.

**Maintainer follow-up:**  use option 2

**Second-pass response.** Instruction understood. Implement option 2: reject a
source callback when it attempts to consume a formula output that has sampled
random contributions. The error should name both the source and the dependent
formula parameter. Callbacks may still consume ordinary posterior/prior
parameters and formula outputs that are fully deterministic at reconstruction
time. No dependency graph or order-dependent fallback should remain.

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

Decision: please remove all of these compatibility layers and simplify the code if possible

**Audit status: decision confirmed; implementation pending.**

**Review response.** Instruction understood, including the earlier instruction
that this release need not preserve backward compatibility. Current code already
rejects unequal fitted/supplied source values in several merge paths, but
`R/JAGS-bridge-formula-context.R` still rebuilds supplied context and can graft
supplied callback/value functions into fitted metadata when names and shapes
match.

The fitted design, fitted source snapshot, and fitted callback must be the sole
authority. Supplied formula data may be used only to verify exact replay inputs,
not to replace executable source semantics. A fit missing required source
metadata should fail with a clear refit instruction. Remove the callback-graft
and legacy rebuild fallbacks rather than retaining a compatibility mode.

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

Decision: please add a bit more information -- an actual example -- so I understand the impact and potential issues

**Audit status: decision confirmed; implementation pending.**

**Concrete example.** Assume a fitted model has a row-indexed posterior source
with columns `tau[1]`, ..., `tau[4]`, corresponding to fitted observations
1 through 4. A user predicts on:

```r
newdata <- fitted_data[c(4, 2), ]
```

The current prediction bookkeeping can treat these two supplied rows as
prediction positions 1 and 2 and therefore attach `tau[1]` and `tau[2]`. The
statistically intended subset is `tau[4]` and `tau[2]`. Reordering, duplicating,
or combining fitted and unrelated new rows can therefore silently attach a
row-specific scale to the wrong observation.

There is no general way to infer observation identity from equal covariate
values: duplicates can be legitimate and new rows can coincide with fitted
rows.

**Recommendation.** For this release, allow row-indexed posterior sources only
for the unchanged full fitted order or when the caller supplies explicit,
validated fitted-row indices. Reject arbitrary reordered/duplicated/new data for
that source type. Callback-computed row sources may continue to support new data
because they are recomputed from the supplied rows. Stable user-defined
observation keys can be added later if arbitrary joins are needed.

**Maintainer follow-up:**  does this corresond to picking a correct factor level when new observations are supplied?
I think we have quite good handling for this for the standard fixed effect formula interface, why don't we reuse this?
Or if its indexing only and we predict, than the predictions are based on the other rows so it does not matter?

**Second-pass response.** This is not factor-level selection. Fixed factors and
random grouping factors are matched by their stored level labels. A posterior
row source such as `tau[4]` is instead attached to the fourth **fitted
observation**, regardless of that observation's factor or group level.

The current prediction code sets `new_row` according to whether the **group
level** is new. Every new observation belonging to an already-fitted group is
therefore treated as an existing row. It then passes the position within the
supplied prediction data as `prediction_rows`. For example, with one fitted
group, latent effect `z = 1`, and fitted row scales `tau = c(1, 2, 3, 4)`,
predicting fitted rows `c(4, 2)` should give row-scaled random contributions
proportional to `c(4, 2)`. The current positional lookup can use
`c(tau[1], tau[2])`, giving `c(1, 2)`. The fixed/group factor mapping can still
be completely correct while this row-specific scale is wrong.

Predictions on other rows do not make the source irrelevant: `tau[i]` directly
multiplies the reconstructed random contribution for prediction row `i`. A
genuinely new observation has no fitted `tau[i]` at all.

**Revised recommendation.**

- If a row source is posterior-indexed (`tau[1]`, ..., `tau[N]`), unchanged
  fitted data may default to `seq_len(N)`. Whenever data are supplied, require an
  explicit integer `fitted_rows` mapping into `1:N`. Reordering then uses, for
  example, `fitted_rows = c(4, 2)`. Intentional duplicates are unambiguous and
  may repeat an index. New observations are rejected for this source type.
- If a row source has a fitted callback that recomputes values from
  `(parameters, data, n_rows)`, it may support arbitrary new rows and needs no
  fitted-row mapping.
- Do not infer row identity from covariates, factor levels, row names, or group
  membership.

**Remaining decision:** confirm this explicit `fitted_rows` contract, or request
a user-defined stable observation-key API instead.

Decision: ok

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

Decision: implement 

**Audit status: decision confirmed; implementation pending.**

**Review response.** Instruction understood. The sample method still stops as
soon as a known group covariance is present, including when every requested
level was fitted. The implementation should jointly draw the fitted group
effects from the covariance implied by the fitted-level subset of the known
group covariance and the coefficient covariance (the separable/Kronecker
construction), then apply the requested design rows. New-level validation
remains a separate question and must not be used to reject fitted levels.

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

Decision: rejected user supplied reserved internal phrases such as "__xXx__"

**Audit status: implemented and verified.**

**Review response.** The reserved-token validator is already applied to data
column names, random block/allocation/source names, and several prior/parameter
identifiers. Tests cover those identifier paths.

The original collision is not fully closed because factor, character, and
grouping **level labels** are also user-supplied and can enter coefficient,
interaction, group, or display names. Those values are not currently passed
through the reserved-token validator. For example, a random-factor level
containing `__xXx__` can later be split as if it represented an interaction.

To implement the decision completely, reject reserved internal tokens in every
categorical label that participates in fixed/random design or grouping metadata,
as well as in identifiers. The error should report the variable and offending
level. This extends the chosen restriction; it does not require a reversible
encoding or a migration layer.

**Implementation outcome.** A shared categorical-level validator now rejects
every BayesTools internal token in fixed-factor levels, random-slope factor
levels, and random grouping levels. It checks unused fitted factor levels as
well as observed character values, reports the variable, token, and offending
level, and runs during both model construction and prediction replay. Stored
fixed/random factor and grouping metadata are also validated during replay, so
older fitted objects cannot silently reintroduce an ambiguous label. The
documentation now states the restriction. Focused formula and prediction tests
passed 1,513 assertions, and the complete unit profile passed 7,786 assertions
with no failures or warnings.

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

Decision: can we calculatate the required space before hand? if so, can we check against the allocable space and stop if we know it would overflow with an informative message?

**Audit status: decision confirmed; implementation pending.**

**Review response.** We can calculate a useful pre-allocation budget, but we
cannot reliably query the exact amount that R and JAGS will be able to allocate.

A dense numeric `n x K` matrix has an eight-byte payload per cell, so its basic
payload is `8 * n * K` bytes. We can use overflow-checked arithmetic to include
known simultaneous matrices, integer maps, planned monitor arrays, and a
conservative copy/serialization multiplier. This calculation can run before
constructing the dense model matrix. The current code already estimates JAGS
syntax nodes, data nodes, and monitored values and stops above configurable
complexity limits, but those estimates are not byte-allocation checks and some
dense matrices are created before the guard runs.

Exact "allocable space" is unknowable in a portable way: R allocator
fragmentation, operating-system overcommit, other processes, temporary copies,
and JAGS' separate native heap all affect the eventual allocation. An
OS-reported available-memory value can therefore be advisory only.

**Recommendation.** Add deterministic, configurable cell/byte ceilings and stop
before allocation with an error reporting `n`, `K`, estimated required bytes,
the configured ceiling, and the sparse/indexed alternative. Optionally compare
the estimate with a conservative fraction of currently available OS memory, but
never promise that the remaining memory is allocable. This is a useful immediate
guard; the sparse/indexed architecture remains the long-term solution.

**Maintainer follow-up:** we can also have some reasonable default like 16GB with an error / warning if its to exceed

**Second-pass response.** A 16 GiB default is feasible, but it should be a hard
stopping ceiling rather than only a warning. Once a conservative estimate
exceeds that size, continuing can terminate the R session or operating-system
process before a warning is actionable. The limit should apply to the estimated
**peak working set for the operation**, not only to the eight-byte payload of
one matrix, and should include known simultaneous matrices and conservative
copy/serialization overhead.

Use a dedicated configurable option measured in bytes, defaulting to
`16 * 1024^3`, and allow an expert to lower it or explicitly set it to `Inf`.
The error should report the payload estimate, peak estimate, dimensions,
ceiling, and alternative representation. It must describe the estimate as a
guard, not as a promise that smaller allocations will succeed.

**Remaining decision:** confirm a hard error above the 16 GiB estimated peak
ceiling. If a warning-only policy was intended, specify that explicitly; I do
not recommend it.

decision: ok, make sure there is a package option to change this

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

Decision: I think most of this was solved separatelly in:
\R-Packages\BayesTools\.agents\reports\numerical-fidelity-decisions.md
if there are any pieces remaining, bring them up as a new issue

**Audit status: resolved by NF10 and NF18, except for new D31.**

**Review response.** The covariance symmetry/canonicalization and
support-derived roundoff policies are now covered by NF10 and NF18. The
near-coincident CAR-coordinate case is still distinct: it concerns loss of a
positive time gap inside `rho^gap`, not repair of a user covariance matrix. It
is therefore isolated as D31 at the end of this document, as requested.

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

Decision: okay, remove the seed argument, default should return original on error with a warning, make time limit reset at extend 

**Audit status: partially implemented.**

**Review response.** Instruction understood. `JAGS_extend()` still exposes and
applies `seed`, and an extension error is currently assigned over the fit
variable, so the last valid fit can be lost. Both changes remain to be made.

The time-limit behavior requested here is already present: `start_time` is set
when each `JAGS_extend()` call begins, so each separate call receives a fresh
budget; iterations within that one call share its cumulative budget. Preserve
that behavior and add an explicit test.

Implementation should:

- remove `seed` from the function, documentation, and tests;
- keep `last_valid_fit` separate from the attempted extension result;
- on extension failure, warn with the backend error and return
  `last_valid_fit`;
- retain a cumulative time budget only within the current `JAGS_extend()` call.

The broader classification of which fresh-fit errors are retryable is not
specified by this decision and should remain unchanged unless raised
separately.

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

Decision: this is desired behavior in some paths (i.e., spike parameters etc..), please provide a brief example which cases are bad and what would change based on the suggestions

**Audit status: decision confirmed; implementation pending.**

**Desired undefined case.** A structural point/spike parameter is fixed in every
iteration. Its between- and within-chain variances are zero, so R-hat or ESS can
be undefined by construction. This is not evidence of failed sampling. If
metadata proves the parameter is structural, it should be labeled
`structural_constant` and excluded from convergence gating.

**Bad undefined cases.**

- A sampled continuous parameter is stuck at one value in every chain.
- Too few retained iterations exist to calculate ESS/MCSE.
- Only one usable chain remains for a diagnostic requiring multiple chains.
- Filtering or a misspelled monitor request leaves an empty diagnostic set.

In those cases, replacing undefined R-hat with `1` or ESS with `Inf` reports
successful convergence precisely because convergence could not be assessed.

**Recommended change.** Classify each requested parameter as:

1. `assessable`, with its computed diagnostic;
2. `structural_constant`, ignored only when fitted prior/parameter metadata
   proves it is fixed; or
3. `not_assessable`, which fails convergence by default.

Diagnostic plot/helpers should return the same classed state (or a concise
error) for constant, too-short, and empty inputs instead of leaking density,
range, or subscript errors. An explicit opt-in may ignore `not_assessable`, but
the current silent conversion to a passing value should be removed for sampled
parameters.

**Maintainer follow-up:**  agree (also let downstream packagfes decide what to monitor, e.g., sometimes its sensible to omit random effects from the monitor)

**Second-pass response.** Instruction understood. The three-state classification
applies only to parameters the caller actually requests. The convergence API
should accept an explicit monitor selection so downstream packages may omit
random effects or other intentionally unmonitored parameters. An omitted
parameter is `not_requested`, not `not_assessable`, and cannot cause failure.
For every requested parameter, undefined diagnostics pass only when fitted
metadata proves `structural_constant`; otherwise they are `not_assessable` and
fail by default. A misspelled requested name must error. An explicitly empty
selection should return an empty diagnostic result without claiming global
convergence.

## D13. Failed marginal-likelihood models

**Issue.** Some model-averaging paths convert failed or non-finite marginal
likelihoods to `-Inf`, effectively assigning the model posterior probability
zero.

**Impact.** Numerical failure becomes statistical evidence against a model and
can silently alter model probabilities and Bayes factors.

**Suggested change.** Preserve failure as `NA` with a classed condition and make
the caller choose between dropping failed models, aborting, or explicitly
treating them as zero-evidence models.

Decision: wasn't this addressed elsewhere?

**Audit status: resolved by NF08.**

**Review response.** Yes. Failed/non-finite marginal likelihoods now use an
explicit `on_failure = c("error", "drop", "zero")` policy. The default is to
abort; dropping or assigning zero evidence requires an explicit caller choice
and records audit metadata/warnings. No further D13 decision is needed.

## D14. Positivity-preserving posterior sampling

**Issue.** `preserve_positive`-style fallbacks can alter sampled values to keep
them inside a requested support rather than reporting that the approximation
produced invalid draws.

**Impact.** Posterior summaries can reflect an undocumented sampling
transformation rather than the represented distribution.

**Suggested change.** Make any correction strategy explicit and opt-in. The
default should use a support-respecting sampler or fail with a diagnostic that
identifies the approximation and invalid-draw rate.

Decision: I need an example where this could have happened so I can decide on policy

**Audit status: resolved/superseded by NF09.**

**Review response.** The original issue description was too broad. The audited
`preserve_positive` path did **not** replace, truncate, or nudge numeric posterior
draws to make them positive. It forced every positive-probability mixture
component to receive at least one sampled row. For example, a component with
posterior probability `1e-8` could receive one of 1,000 rows and therefore be
represented with empirical mass `1e-3`, a 100,000-fold overweighting.

NF09 removed that behavior. Mixture rows now follow one shared categorical
sampling plan; a positive-probability component may legitimately receive zero
rows in a finite draw. I found no remaining posterior-value correction path
requiring a separate positivity policy.

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

Decision: I think this was solved separatelly in:
\R-Packages\BayesTools\.agents\reports\numerical-fidelity-decisions.md
if there are any pieces remaining, bring them up as a new issue

**Audit status: resolved by NF05 and NF06.**

**Review response.** Finite grids now carry approximation/provenance
information, refine/extend under an explicit error budget when possible, and
reject outside-range queries when exact source provenance is unavailable.
Endpoint-singular densities retain their mathematical behavior and nonregular
density ratios are rejected. No additional D15 issue was found.

## D16. One-sided weight-function marginal inference

**Issue.** Some exported marginal-inference paths do not implement one-sided
weight-function components even though the priors are supported for fitting.

**Impact.** A fitted supported model can fail only when users request marginal
inference or derived hypotheses.

**Suggested change.** Either implement the missing marginal distribution or
declare the combination unsupported in the exported function documentation and
validate it before computation.

Decision: implement the missing marginal distribution

**Audit status: partially implemented.**

**Review response.** The monotone/cumulative one-sided weight-function marginal
is already implemented and tested through the shared component representation,
including CDF, CCDF, log-density, and quantiles.

The remaining missing distribution is the exported **general** one-sided
parameterization with `alpha1` and `alpha2`.
`.mdone.sided_general()`, `.mpone.sided_general()`, and
`.mqone.sided_general()` still stop with `"Not implemented"`, and the current
tests explicitly expect those errors, even though `rone.sided(alpha1, alpha2)`
can generate draws. The decision therefore means implementing and testing those
general marginal CDF/probability/quantile paths, then removing the expected-error
tests. Numerical inversion should use the same error-budget/provenance rules
adopted under NF05/NF06.

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

Decision: make formulas treat prior_none as spike(0) consistently

**Audit status: implemented.**

**Review response.** Instruction understood. Canonicalize formula
`prior_none()` to `prior("point", list(location = 0))` during formula
construction, before syntax, monitoring, prediction, and bridge metadata are
derived. Replay and reconstruction should then see only the canonical point
prior. Also remove the current unknown-prior fallback that can silently
contribute zero: an unrecognized fixed-prior class must produce an internal
contract error rather than being treated as another spelling of `prior_none()`.

**Implementation outcome.** Direct and default continuous `prior_none()` formula
entries are now replaced by point-zero priors before formula validation and
metadata construction, while retaining their model prior weights. Public and
compiled reconstruction now reject noncanonical prior classes instead of
silently omitting their terms. Focused tests passed 1,449 assertions, and the
complete unit profile passed 7,755 assertions with no failures or warnings.

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

Decision: contrasts should be per-factor specific. We should prbl keep the info in interaction meta-data for cases when handled separatelly later.
Are there any changes that we need to make? I need more info for decision.

**Audit status: decision confirmed; implementation pending.**

**Review response.** I understand "per-factor specific" to mean that one factor
uses one contrast family consistently wherever it occurs, while different
factors may use different families. Under that interpretation:

- Fixed-effect construction already enforces one family per factor and stores
  `factor_contrasts`, `factor_design`, and cell/level metadata on interactions.
- Fixed interactions between different factors already retain each factor's
  separate basis.
- Random-factor interactions still reject an interaction when its different
  factors use different contrast families. For example, treatment coding for
  `f` and orthonormal coding for `g` is rejected as "mixed factor contrast
  families". That conflicts with a per-factor policy and should be changed.
- `contr.orthonormal()` and `contr.meandif()` still derive a basis from an
  eigendecomposition of a matrix with a repeated eigenvalue. The selected
  orientation can vary by numerical platform even though the spanned subspace
  is the same.

**Recommended changes.**

1. Keep rejecting different contrast families for the **same** factor across
   main effects/interactions.
2. Allow different factors in a random interaction to use their own families
   and persist each factor's design in interaction metadata.
3. Replace the eigenvector construction with one deterministic centered Helmert
   basis, with explicit normalization for orthonormal versus mean-difference
   coding and a fixed sign convention.
4. Version the metadata schema so downstream consumers can validate it. Because
   backward compatibility is not required, no legacy migration path is needed;
   old fits should request a refit.

If this matches the intended meaning of "per-factor specific", no further
semantic choice is needed.

**Maintainer follow-up:**  I don't understand why we need 3 -- does not it change the meaning of meandif and orthonormal contrasts? I am quite worried about these suggestions.
(also, we might need to use different contrasts for the fixed and random effect parts especially if they are part of special random effects structures)

**Second-pass response and correction.** The concern is valid. Replacing the
current eigenvector basis with a Helmert basis does change the raw coefficient
axes:

- For fixed `meandif`/`orthonormal` priors centered at zero, BayesTools restricts
  the nondegenerate distributions to spherical multivariate normal/t families
  (plus multivariate points). An orthogonal rotation preserves the centered
  factor-level prior and fitted linear predictor, so the advertised level-scale
  meaning survives, but raw contrast coefficients and serialized draws change.
  The constructor currently fails to enforce the required zero center; that
  separate problem is recorded as D32.
- For random effects, the column axes can have separate SD priors, correlations,
  homogeneous/heterogeneous restrictions, or a special structured covariance.
  Rotating those axes is not generally prior-invariant and can change the
  statistical model.

Therefore item 3 should **not** be treated as a harmless portability cleanup. I
withdraw it from the immediate recommendation. The safer replay fix is to store
the exact concrete contrast matrix used by each fitted design and reuse that
matrix for prediction and bridge reconstruction. Independent refits on
different numerical platforms may still choose different raw axes, but one
serialized fit cannot silently change axes during replay.

The follow-up also clarifies that "per-factor" needs a scope. The appropriate
contract is **per factor within one design owner**, not one global basis for the
entire model:

1. The fixed-effect design has its own factor-to-contrast map.
2. Each random block has its own factor-to-contrast map.
3. Within one owner, the same factor uses one basis across its main effect and
   interactions; different factors in an interaction may use different bases.
4. Special structured random blocks may define a level/coordinate basis
   required by their covariance structure and must not inherit a fixed-effect
   basis accidentally.
5. Persist the concrete matrices and interaction component metadata for every
   owner. Because backward compatibility is not required, fits lacking this
   schema should request refitting.

This requires an explicit random-block contrast specification (or a documented
structure-defined basis) rather than deriving every random contrast from the
factor attributes already modified for the fixed design.

**Remaining decision:** confirm this owner-scoped contrast contract and that
the current `contr.orthonormal()`/`contr.meandif()` algorithms should remain
unchanged for this release. A deterministic replacement can be considered later
as a separately reviewed statistical parameterization change.

Decision: agree

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

Decision: okay, I see the potential issues, please create the registry infrastructure and 
fully switch the package to it. Make sure that there is a documetnation object for downstream packages that 
they can read to undetstand how it works if needed. Make breaking changes, we will update the downstream packages appropriatelly. 
Do not care about backwards compatibility. Update the tests accrodingly.

**Audit status: decision confirmed; implementation pending.**

**Review response.** Instruction understood. No canonical registry exists yet;
random-effect summaries still contain string-prefix ownership heuristics.

The registry should be versioned fitted metadata with one record per canonical
parameter/monitor and fields sufficient to avoid reparsing names, including:
formula parameter, role (fixed coefficient, random SD/correlation/latent/group
coefficient, allocation, or derived summary), random block ID, term/column,
indices/dimensions, fitted scale, monitor status, and display label. Syntax,
posterior extraction, scale transformation, bridge compilation, summaries, and
tables should all consume this registry.

For downstream packages, expose a documented schema object/function and an
accessor that validates and returns a fit's registry; do not require downstream
code to inspect an internal attribute by convention. Remove heuristic fallback
paths and reject fits without the new registry with a refit message. Tests should
assert registry completeness/uniqueness and compare every consumer against the
same registry.

## D20. Marginal-likelihood computation contract

**Issue.** Documentation and implementation do not fully define the reported
scale, multiple bridge repetitions, multiple chains, or the handling of
non-finite per-repetition results.

**Impact.** Callers may combine values on the wrong scale or receive unstable
results whose aggregation depends on incidental vector behavior.

**Suggested change.** Document and encode the scale in the returned object,
define deterministic aggregation and diagnostics across repetitions/chains, and
return all repetition-level results alongside the selected summary.

Decision: proceed

**Audit status: decision confirmed; implementation pending.**

**Review response.** Instruction understood. `JAGS_bridgesampling()` currently
returns the upstream `bridge`/`bridge_list` object directly, while downstream
code reads `$logml` and can accidentally receive a vector.

The new BayesTools result should have:

- one selected scalar `logml` explicitly labeled as natural-log marginal
  likelihood;
- one row per repetition with `logml`, iteration count, MCSE when available,
  finite/success status, warning/error details, method, and chain metadata;
- the aggregation rule and number of included/failed repetitions;
- the original upstream result in a clearly named diagnostic field if retaining
  it is useful.

For multiple successful repetitions, use the upstream package's documented
summary convention: the median of finite repetition-level **log** marginal
likelihoods. A non-finite repetition should abort by default under the NF08
failure policy; an explicit drop policy may aggregate the remaining finite
repetitions and must record what was dropped. Model averaging and hypothesis
code should accept only the new scalar contract, never incidental vector
coercion.

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

Decision: proceed with the cache manisfest changes etc, make sure that CI behaves the same as local runs.
everything about .StatsVault should be local only (noone outside of this machine should know about it)

**Audit status: partially implemented.**

**Review response.** The test fit cache now has versioned completion metadata,
required-artifact declarations, registry hashes, source/function fingerprints,
and stale-cache integrity tests. Local and CI entry points use the same profile
runner, and GitHub Actions invokes the `all` profile. `.StatsVault/` is ignored
by both Git and package builds, matching the local-only decision.

Remaining work:

- `models/RandomEffects.RDS` is still validated only for readability, object
  names, and classes. Its helper needs a manifest derived from the vignette
  generator, relevant package sources, package/schema versions, model names, and
  backend fingerprint. The regeneration chunk must write the cache and manifest
  together.
- Real JAGS fitting still occurs in test files other than
  `test-00-model-fits.R`. Those cases need to become registered cached fixtures,
  or be reduced to non-fitting unit tests, so the repository follows the stated
  single-fit-file rule.
- The source-tar snapshot paths remain a separate incomplete item under D30.

No `.StatsVault` path or metadata should be added to manifests, tests,
documentation, or CI.

## D22. Explicit summary schema labels

**Issue.** Explicit summary schemas can override interval labels independently
of the actual probability columns. This may be intentional presentation control
or may allow a label that misstates the computed interval.

**Impact.** A published table can describe a different interval mass/type from
the values it displays.

**Suggested change.** Either derive labels exclusively from interval metadata,
or validate user-supplied labels against that metadata and require an explicit
override flag for deliberately custom wording.

Decision: I need more information/example, I do not follow

**Audit status: decision confirmed; implementation pending.**

**Concrete example.** Suppose columns named `0.025` and `0.975` contain the
2.5% and 97.5% posterior quantiles. They therefore display a 95% equal-tail
interval. The explicit interpretation schema can select those two columns but
also set `interval_level = 0.90`. The endpoint values remain the 95% interval,
while generated prose can call them a "90% interval".

A similar mismatch is possible when `lower_prob` and `upper_prob` are supplied
independently of the chosen columns. The code has no way to prove that a
manually named `lower` column actually contains the probability stated in the
schema.

**Recommendation.** Treat selected quantile probabilities as statistical
metadata: infer interval mass as `upper_prob - lower_prob` and reject a
conflicting `interval_level`. Keep a separate cosmetic `interval_label` only for
wording such as "credible interval" versus "uncertainty interval"; it must not
override the numeric probability. For arbitrary precomputed intervals, require
the caller to supply matching endpoint probabilities or omit the percentage
from generated prose.

**Maintainer follow-up:** where is the `interval_level` from? I do not remmember it? is it neccessary?

**Second-pass response.** `interval_level` is not an argument of a summary-table
function. It is currently documented input metadata for `interpret_records()`.
A caller may put it in either a source `schema` or an interpretation-plan item.
`.interpret_interval_columns()` reads that caller value first; only when it is
absent does it calculate `upper_prob - lower_prob`. The normalized record keeps
the field, and the generic text renderer turns it into wording such as
`"95% interval"`.

It is not necessary as caller input. When `lower_prob` and `upper_prob` are
known, the mass is fully determined. When they are unknown, a standalone
`interval_level` cannot prove that the selected columns contain that interval.

**Recommendation.** Remove `interval_level` from accepted source-schema, plan,
and direct estimate-record input. Retain it in normalized output as a derived
convenience field:

- with finite endpoint probabilities, set
  `interval_level = upper_prob - lower_prob`;
- without endpoint probabilities, set it to `NA` and render an unnumbered
  `"uncertainty interval"`.

This makes the actual endpoint probabilities authoritative and eliminates the
mislabeling path without losing valid information.

**Remaining decision:** confirm removal of caller-supplied `interval_level`.

decision: agree

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

Decision: break backwards compatibility, do not add backwards handling

**Audit status: implemented and verified.**

**Review response.** Instruction understood. Prediction still falls back to
`formula_design$model_frame` when `source_data` is absent, and tests explicitly
exercise that legacy path. Remove the fallback and its tests. A design without
the new versioned `source_data`/scale metadata must stop with a concise
instruction to refit under the current BayesTools version; it must never infer
scale from field absence.

**Implementation outcome.** Formula designs now store a schema version and an
explicit scale declaration for `source_data`, `model_frame`, and
`model_matrix`. Prediction and bridge replay validate that schema and require
the original-scale `source_data`; a missing, unversioned, or unsupported schema
stops with an instruction to refit. The `model_frame` fallback has been removed,
so a model-scale frame can no longer be mistaken for original-scale prediction
data and scaled twice. Focused prediction and bridge tests passed 1,572
assertions, and the complete unit profile passed 7,791 assertions with no
failures or warnings.

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

Decision: I think we implemented somewhat different hadnling of those two (i.e., the random effects need to be scaled differently and might not be scaled if fixed are?)
please examine these differences first and provide more background info

**Audit status: decision confirmed; implementation pending.**

**Background.** Yes, fixed coefficients and random-effect distribution
parameters require different transformations:

- Fixed coefficients use an affine coefficient-basis transform. Centering a
  slope also changes the intercept, and interactions can contribute to several
  lower-order coefficients.
- Random-effect SD/correlation summaries describe a covariance matrix. They must
  use the corresponding block transform
  `G_original = A %*% G_scaled %*% t(A)`, then recover transformed SDs and
  correlations. Applying a slope's scalar divisor independently is wrong for a
  correlated intercept/slope block.
- Realized group-specific coefficient vectors use the coefficient-basis
  transform `A` independently for every group. This again includes intercept
  changes induced by centered random slopes, but it is not the same operation as
  transforming SD/correlation summaries.

The current code correctly gives semantic random SD/correlation columns their
own covariance-aware path, with tests for correlated slopes, factor columns,
homogeneous-SD guards, and structured blocks. Raw
`..._xRE_COEFx[group,column]` monitors are classified under the random-effect
prefix but do not map to semantic SD metadata, so they pass through unchanged on
the fitted standardized scale.

**Recommendation.** Keep raw latent and group-coefficient monitors explicitly
internal and standardized in the generic posterior object. Exclude them from
tables/accessors that promise original-scale coefficients. If downstream users
need group-specific realized coefficients, add a dedicated registry-backed
extractor that applies `A` blockwise and returns an explicitly
`scale = "original"` object. This avoids silently mixing scales without forcing
the SD/covariance and realized-coefficient transformations into one path.

The remaining decision is whether raw `xRE_COEFx` monitors are intended as a
public original-scale API. I recommend **no**; expose the dedicated transformed
extractor instead.

**Maintainer follow-up:**  Keep raw latent and group-coefficient monitors explicitly
internal and standardized in the generic posterior object. agree with recommendation

**Second-pass response.** Instruction understood. Raw latent and
`xRE_COEFx` columns remain internal standardized coordinates and must be
excluded from APIs that promise original-scale coefficients. Their registry
entries must state that scale explicitly. A public group-specific coefficient
API, if retained or added, must return a separate object after applying the
blockwise coefficient transform and label it `scale = "original"`.

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

Decision: implement this proper handling

**Audit status: partially implemented; implementation pending.**

**Review response.** The general posterior-atom infrastructure can already
store atoms, linearly transform them into factor-level coordinates, and let
factor plotting consume stored point masses separately from continuous curves.
That provides useful foundation.

The ordered **prior-density** fallback remains incorrect for this case.
`.density.prior.ordered()` samples each level and, unless every value is
identical, sends the full vector to one KDE. When the ordered total is
spike-and-slab, the exact repeated atom is therefore smoothed together with the
continuous draws.

Implement one mixed-measure object per ordered level containing:

- an explicitly normalized atom table (`location`, `mass`);
- a continuous density curve whose integral equals that level's remaining
  continuous mass;
- source/approximation diagnostics and transformation metadata.

Build the per-level atoms analytically from the total-prior components and
allocation map where possible; otherwise use exact component/model indicators,
not equality-with-tolerance on sampled values. Plotting and marginal-inference
consumers must render/summarize both parts without renormalizing either part to
one.

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

Decision: reject transformations for now, lets make sure we can use the basics for now correctly. add a proper error messages for anything we cannot handle

**Audit status: implemented and verified.**

**Review response.** Instruction understood for transformations. Generic random
slopes currently derive `predictors` from deparsed `terms()` variables, so an
inline call such as `I(x^2)` is later reported as a missing column. Add an early
random-formula grammar validator and reject every call/transformation with the
offending expression and a message to create an explicit data column first.
Apply the same rule to prediction replay.

One choice remains from the original issue: `g1:g2` group-level ordering. The
current hand-built `interaction()` order differs from familiar base-R/lme4
formula ordering. Changing it alters stored `group_levels`, `group_map`, monitor
row meanings, and fixtures. With exchangeable group priors it is mostly a
permutation, but with known group covariance or any nonexchangeable group
information it can change which covariance row belongs to which group.

**Recommendation.** Adopt base-R/lme4-compatible interaction ordering now,
while breaking changes and refits are already accepted. Persist the component
levels and tuple-to-index map in metadata so replay never reconstructs ordering
from display strings. Do not add a legacy migration.

**Maintainer follow-up:**

**Second-pass response.** The transformation part is fully decided and can be
implemented independently. The blank follow-up leaves only grouping-interaction
ordering unresolved.

Adopting base-R/lme4 ordering is a breaking change to row numbers in
`group_levels`, `group_map`, latent/coefficient monitors, and known group
covariance alignment. Because this release already requires refitting and aims
to make the basic formula interface familiar, this is the lowest-maintenance
choice. Store each interaction as component-level tuples plus an explicit
tuple-to-index map; use joined labels only for display. Prediction and bridge
replay must consume the stored map rather than call `interaction()` and infer
the order again.

**Remaining decision:** confirm base-R/lme4-compatible grouping-interaction
ordering with no legacy migration.

decision: agree

**Implementation outcome.** Random-slope formulas now reject inline
transformations and identify the offending expression while directing callers
to create an explicit data column. Grouping interactions use the lexicographic
component ordering produced by base R with `lex.order = TRUE` and by lme4.
Compiled metadata persists the grouping components, their fitted levels, the
concrete component tuples, collision-free length-prefixed tuple keys, and the
tuple-to-index map. Prediction reconstructs the same tuple keys and refuses
older fitted objects that lack this metadata, so no ambiguous migration path
remains. Display labels are retained separately and disambiguated when distinct
tuples render to the same colon-joined text. Tests cover reordered prediction
rows and a deliberate display-label collision that the old `interaction()`
encoding collapsed. Focused random-formula tests passed 1,506 assertions.
The complete unit profile passed 7,779 assertions with no failures or warnings.

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

Decision: examine what is feasible to implement, add checks for unspupported expressions, and improve documentation

**Audit status: decision confirmed; implementation pending.**

**Review response.** All requested boundaries are feasible without evaluating
arbitrary R code. The current whitelist already limits arithmetic/functions,
but a complete relation wrapped in parentheses is missed by the string-level
relation finder, unary `!` is not routed as a top-level condition, and parsed
backticks lose the lexical distinction between an escaped identifier named
`Inf` and the unescaped constant `Inf`.

I recommend publishing and implementing this small grammar:

```text
hypothesis := statement [ "vs" statement ]
statement  := point | region
point      := arithmetic ( "=" | "==" | "!=" ) finite_number
region     := relation
            | "(" region ")"
            | "!" region
            | region ( "&" | "|" ) region
relation   := arithmetic ( "<" | "<=" | ">" | ">=" ) arithmetic
arithmetic := identifiers, finite numeric literals, parentheses,
              + - * / ^, and the documented function whitelist
```

Equality remains a point hypothesis and cannot be nested inside boolean region
expressions. Support parenthesized relations and unary negation for **regions**;
reject negated point equalities with a precise explanation. Preserve an
escaped-identifier token through parsing so `` `Inf` `` can refer to a real
parameter, while unescaped `Inf`, `NaN`, `NA`, `TRUE`, and `FALSE` are rejected
where a model parameter or finite point value is required. Add a documented
accepted/rejected example table and adversarial tests for every production.

**Maintainer follow-up:** ok

**Second-pass response.** Instruction understood as approval of the proposed
restricted grammar, precise rejection of unsupported expressions, preservation
of escaped identifiers, documentation examples, and adversarial tests.

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

Decision: what is the prefered solution and the consequences?

**Audit status: decision confirmed; implementation pending.**

**Review response.** Prefer honoring the existing `n_groups` and `n_columns`
arguments. Structured-local layout already records the global group and column
index for every retained local latent cell, so the function can return only
nodes with `group <= n_groups` and `column <= n_columns`; a zero dimension
returns `character()`.

Consequences:

- dense and structured-local layouts have one consistent prefix/subset
  contract;
- callers cannot accidentally monitor or reconstruct more coordinates than
  requested;
- no second "full layout" API or duplicated dispatch path is needed;
- tests may expose callers that passed reduced dimensions while relying on the
  arguments being ignored—those callers must request the full stored dimensions
  explicitly.

Removing the arguments instead would make full-layout intent clearer locally,
but it would split a generic interface used by several bridge,
reconstruction, and prediction paths and make future subset optimization
harder. Since backward compatibility is not required, the small caller
corrections caused by enforcing the existing contract are preferable.

**Maintainer follow-up:**  ok, agree

**Second-pass response.** Instruction understood. Honor `n_groups` and
`n_columns` consistently, including zero-size requests, and correct callers
that previously relied on the arguments being ignored.

## D29. Exact references for stochastic fit outputs

**Issue.** A successful cache refresh regenerated all 73 fit objects and 22
marginal-likelihood objects, after which the fixture profile reported 70 exact
text-reference mismatches. Every failure was a sampled numeric value or
diagnostic (posterior summaries, model probabilities, marginal likelihoods,
ESS, or R-hat); table structure and semantic assertions passed. The
visual-fixture profile likewise reported 154 changed snapshots (147 generated
SVG files across all four fitted-output plot suites), alongside 445 passing
assertions and no warnings or runtime errors. The current policy forbids agents
from enabling reference generation without maintainer approval.

**Impact.** A valid fit refresh can make both fixture lanes fail even when the
implementation, table schema, and plotting code are unchanged. Updating exact
files and visual snapshots blesses one stochastic realization, while leaving
them unchanged prevents a refreshed cache from passing locally.

**Suggested change.** Decide which outputs are true golden values. Keep exact
references for deterministic formatting and schema, but compare stochastic
statistics and plots with documented tolerances or invariant assertions. If
the current exact-snapshot policy is intentional, explicitly approve
regenerating and reviewing the 70 affected text references and 147 affected SVG
files from the validated cache as one controlled update.

Decision: okay, I will examine and decide before running the full change implementations

**Audit status: decision confirmed; implementation pending.**

**Review response.** During the numerical-fidelity implementation, the validated
fit/marginal-likelihood cache, text references, and fitted-object visual
snapshots were refreshed together and the applicable fit, fixture, and
visual-fixture profiles passed. That resolves the immediate branch consistency
problem; it does not decide the architecture for future releases.

**Recommended long-term policy.**

- Keep exact references for deterministic schema, labels, ordering, formatting,
  and deterministic numerical oracles.
- Test sampled summaries with explicit tolerances/invariants tied to Monte Carlo
  error rather than byte-exact text.
- Keep visual references exact for deterministic plot-data/geometry tests.
- For fitted stochastic plots, use a versioned validated cache and either
  morphology/invariant tests or an explicitly approved controlled snapshot
  refresh.

The maintainer decision is still needed before changing the present exact
stochastic-reference policy.

**Maintainer follow-up:** ok, agree

**Second-pass response.** Instruction understood as approval of the recommended
long-term split: exact references for deterministic schema/format/geometry;
tolerance or invariant tests for stochastic statistics; and controlled,
versioned handling of fitted stochastic plots. The already-refreshed references
remain the current branch baseline while those test-policy changes are made.

## D30. Portable source-package paths for visual snapshots

**Issue.** `R CMD check --as-cran` reports 50 visual snapshot paths longer
than the portable tar limit. Forty-nine were already present at the merge base;
this branch adds
`priors-plot/priors-plot-petpeese-decreasing-transformation-ggplot.svg`.
The package currently ships all 530 SVG reference files in its source tarball.

**Impact.** The package builds and checks successfully on Windows, but CRAN
reports a portability NOTE and other tar implementations are only guaranteed
to preserve paths up to 100 bytes. Renaming snapshots is mechanically broad
because the `expect_doppelganger()` labels and stored files must stay aligned;
excluding snapshots from the tarball changes which visual references downstream
source-package users receive.

**Suggested change.** Choose whether visual references are part of the
distributed source package. If they are, shorten the 50 snapshot labels and
filenames in a controlled visual-only update. If they are only CI/development
artifacts, exclude `tests/testthat/_snaps` in `.Rbuildignore` while keeping the
files in Git and running visual profiles from repository checkouts.

Decision: showten paths as needed

**Audit status: decision confirmed; implementation pending.**

**Review response.** Instruction understood. The current source tree still has
50 snapshot paths that exceed the portable 100-byte tar path limit once the
`BayesTools/` package-directory prefix is included. Keep visual snapshots in the
source package and shorten those labels/files in a visual-only change. Use a
deterministic abbreviation map, update each `expect_doppelganger()` label and
matching SVG together, verify that no prefixed tar path exceeds 100 bytes, and
run both visual profiles.

## D31. CAR gaps that are positive but not representable in `rho^gap`

**Issue.** D10's covariance-repair policy is resolved, but a distinct CAR
transition problem remains. Coordinates can be finite and exactly distinct in
double precision while a positive gap is so small that `rho^gap` rounds to
exactly one. For example, `c(0, 1e-320)` with an ordinary interior `rho` can
produce `phi == 1`. Some R reconstruction code computes the innovation variance
with `expm1()`, but generated JAGS syntax and dense covariance paths do not
consistently retain that stability.

**Impact.** A mathematically positive innovation variance can become zero,
making a transition or covariance singular. Failure then depends on whether the
sample, bridge, prediction, or covariance path happened to use the stable
calculation. Adding coordinate jitter or a fixed minimum gap would change the
model's time scale.

**Suggested change.** Use the log transition consistently:

```text
log_phi        = gap * log(rho)
phi            = exp(log_phi)
innovation_var = -expm1(2 * log_phi)
```

Provide a native/backend helper where JAGS lacks a stable `expm1` primitive, and
use recurrence-based construction for marginal covariance/factorization rather
than first materializing a nearly singular dense correlation matrix. If the
stable innovation is itself zero, non-finite, or unsupported by the backend,
stop with the offending coordinates/gap and explain that the requested time
resolution is not representable. Do not insert an epsilon, jitter coordinates,
or silently rescale time; rescaling changes the meaning of `rho` per coordinate
unit and must be an explicit model change.

**Recommendation.** Implement the stable recurrence and reject only when its
support-derived representability check fails. This preserves more valid models
than an arbitrary coordinate-separation tolerance while remaining
mathematically faithful.

**Audit status: decision confirmed; implementation pending.**

**Decision:** ok

**Second-pass response.** Instruction understood as approval of the stable
log-transition/recurrence implementation and its representability failure rule.
Do not add epsilons, coordinate jitter, or implicit time rescaling. Apply the
same transition calculation to syntax, native/backend evaluation,
reconstruction, prediction, and covariance paths.

## D32. Nonzero locations in mean-difference and orthonormal factor priors

**Issue.** `prior_factor(..., contrast = "meandif")` and
`prior_factor(..., contrast = "orthonormal")` currently accept a nonzero
multivariate location/mean. `JAGS_formula()` also accepts and compiles that
prior. In contrast, `rng(..., transform_factor_samples = TRUE)` explicitly
stops unless the location is exactly zero because these factor priors are
defined as centered, exchangeable deviations.

For example, this is currently accepted by formula construction:

```r
prior_factor(
  "mnormal",
  list(mean = 1, sd = 1),
  contrast = "meandif"
)
```

The scalar `mean = 1` is repeated in the arbitrary contrast-coordinate basis.
Changing or numerically rotating that basis changes the implied factor-level
mean pattern. It is therefore not a stable, basis-independent statement that
every factor level is one unit above the grand mean (which would also contradict
the sum-to-zero contrast subspace).

**Impact.** Formula fitting can use a prior that public transformed-prior
methods reject. Its directional factor-level mean depends on the eigenvector
orientation used by `contr.meandif()`/`contr.orthonormal()`, so independent
platforms or a future basis change can fit statistically different priors from
the same user specification.

**Suggested change.** Require the multivariate mean/location to be exactly zero
for every nonzero-probability component of a `meandif` or `orthonormal` factor
prior, including mixture and spike-and-slab components. Require multivariate
point components to be at zero as well. Apply the check in `prior_factor()` and
again when validating composed formula priors. Use exact zero here; a tolerance
or numerical nudge would give an arbitrary nonzero coordinate-space direction
statistical meaning.

If directional nonzero factor-level prior means are needed later, add an
explicit level-space API that accepts a named, sum-to-zero level vector and
transforms it through the stored concrete basis. Do not overload a scalar
coordinate-space location.

**Audit status: implemented.**

**Implementation outcome.** `prior_factor()` now requires the mean/location of
mean-difference and orthonormal priors to be numeric, finite, scalar, and
exactly zero. Formula validation repeats the recursive check for mixtures and
spike-and-slab priors, which also protects against manually modified or
serialized prior objects. Nonzero treatment and independent factor priors
remain valid. Targeted tests passed 1,575 assertions, and the complete unit
profile passed 7,730 assertions with no failures or warnings.

**Recommendation:** reject nonzero locations now. This makes construction,
formula fitting, density, and transformed sampling obey the same centered
factor-prior contract.

**Decision:**: the non-0 priors should be disallowed for thoser priors
