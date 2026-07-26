# Numerical and statistical fidelity: maintainer decisions

This document isolates numerical conventions that may be mistaken for statistical
assumptions. It covers both changes introduced during the current branch review
and relevant pre-existing behavior uncovered while reviewing those changes.

The main conclusion is:

- I did **not** find a generic epsilon added to likelihoods, posterior densities,
  variances, probabilities, or Bayes factors merely to keep them finite.
- I did find several tolerances, clipping rules, and fallback values. Some are
  faithful floating-point safeguards, some are approximation policies that need
  to be visible, and some can change the represented quantity or distribution.

For each case, the recommended handling is marked **Recommended**. A decision can
be recorded as, for example, `NF04: accept recommendation`, or replaced with a
different policy. Cases are independent unless their implementation notes say
otherwise.

## Review of the maintainer responses

The maintainer's `Decision` text below is preserved verbatim. Each
`Review response` records either the understood implementation instruction,
additional evidence needed for a decision, or a narrowly scoped remaining
choice.

| Case | Review status |
|---|---|
| NF01 | Implemented and verified |
| NF02 | Implemented and verified |
| NF03 | Implemented and verified |
| NF04 | Implemented and verified |
| NF05 | Implemented and verified |
| NF06 | Implemented and verified |
| NF07 | Implemented and verified |
| NF08 | Implemented and verified |
| NF09 | Implemented and verified |
| NF10 | Implemented and verified |
| NF11 | Implemented and verified |
| NF12 | Implemented and verified |
| NF13 | Implemented and verified |
| NF14 | Implemented and verified |
| NF15 | Implemented and verified |
| NF16 | Implemented and verified |
| NF17 | Implemented and verified |
| NF18 | Implemented and verified |

## Decision index

| Case | Topic | Consequence | Recommended disposition |
|---|---|---|---|
| NF01 | Positive latent initializations | Sampling only, but can affect warmup | Replace the underflow fallback with an explicit interior initialization policy |
| NF02 | Near-diagonal covariance detection | Can discard real covariance in diagonal-only output | Require structural/exact diagonality by default |
| NF03 | Gaussian KDE tail evaluation | Changes the density estimator, not the posterior target | Retain with explicit extrapolation diagnostics and a strict option |
| NF04 | Tiny hypothesis coefficients treated as zero | Can change a tested linear combination and its prior | Use exact structural zeros; never drop a coefficient solely because it is small |
| NF05 | Finite prior grids and omitted tails | Can change prior ordinates and Bayes factors | Make approximation error and omitted mass explicit; extend adaptively |
| NF06 | Endpoint-singular densities | Can replace an infinite boundary density with zero | Preserve singularities numerically and reject nonregular density ratios |
| NF07 | Silent positive-semidefinite projection | Can change predictive distributions | Clamp only roundoff-sized negative eigenvalues; otherwise error |
| NF08 | Failed marginal likelihoods converted to `-Inf` | Converts computation failure into evidence against a model | Abort by default; require an explicit failure policy |
| NF09 | Forced-positive mixture sample counts | Can heavily overweight very small model probabilities | Sample model indicators from the mixture; make stratification opt-in |
| NF10 | Covariance and CAR tolerance/repair rules | Can change a supplied covariance or accept inconsistent inputs | Use one scale-aware validation and canonicalization policy |
| NF11 | Savage-Dickey atom and support identity | Can apply a density ratio where it is invalid | Use exact model metadata and exact null identity |
| NF12 | Weight-function cut coalescing and snapping | Can remove narrow bins or change boundary assignment | Preserve cuts and one explicit interval-closure convention |
| NF13 | Near-one-hot matrices compiled as exactly one-hot | Can change the design matrix used by the model | Require proof of structure or fall back to the dense exact path |
| NF14 | Simplex sum tolerances | Can accept off-simplex values as valid | Canonicalize only roundoff-scale drift and reject larger discrepancies |
| NF15 | Fitted/rebuilt design equivalence tolerances | Can allow marginal-likelihood replay with different data/designs | Treat fitted values as authoritative and compare with scale-aware diagnostics |
| NF16 | Plot-only offsets and spike tolerances | Can visually merge distinct values | Use metadata and representable interior values |
| NF17 | Automatic centered/noncentered thresholds | Affects computation, not the mathematical target | Retain, expose, and record the resolved choice and reason |
| NF18 | Native clamps and faithful safeguards | Usually protects the intended target from roundoff | Retain only support-derived, roundoff-bounded clamps |

## NF01. Underflowed Gamma/Dirichlet latent initializations

**Origin.** The current branch added the fallback to additional cumulative
selection initialization paths in `f398d35`; equivalent behavior already existed
in other Dirichlet initialization paths.

**Current behavior.** In `R/JAGS-prior-inits.R` and
`R/selection-backend-helpers.R`, non-finite or non-positive values returned by
`stats::rgamma()` are replaced by `.Machine$double.xmin`, approximately
`2.225074e-308`.

**Mathematical impact.** This does not alter the prior or posterior density. JAGS
initial values are starting states, not observations or fixed parameters.
Nevertheless, `.Machine$double.xmin` is not a draw from the requested Gamma
distribution and can be a very poor starting state. With extremely small shape
parameters, underflow to zero is common rather than exceptional, and initializing
multiple latent Gamma variables at the same minimum value can create unstable
normalization, slow adaptation, or misleading convergence behavior.

**Recommended.** Replace the fallback with a named, documented interior
initialization policy:

1. Prefer a deterministic interior point whose normalized proportions match the
   Dirichlet prior mean when that point is representable.
2. Rescale the latent Gamma vector to an order-one common scale; its normalized
   proportions, rather than its arbitrary common magnitude, are the important
   starting feature for the simplex representation.
3. If the intended proportions cannot be represented, use a neutral interior
   point and emit a classed initialization diagnostic that identifies the extreme
   hyperparameters. Do not describe this fallback as a prior draw.
4. Preserve user-supplied valid initial values and test extreme-shape behavior
   through adaptation, not only object construction.

Rejecting every mathematically valid but computationally extreme shape would be
unnecessarily restrictive. A strict option may reject such priors when the user
requires initialization to reflect prior proportions.

**Decision:** this seems like an extreme edge case,
instead of custom one-by-one hadnling, should we hjust have a fallback for the case that rng init values are NA, -> throw an immediate warning and fall back to the mean of the histribution (this should cover 99% of cases)

**Review response:** Understood, with one necessary qualification. A shared
fallback is preferable to distribution-path-specific fixes here. The trigger
should be any invalid initialization draw, not only `NA`: `NaN`, infinite,
non-positive, or otherwise outside the latent variable's support must also
trigger it.

For the affected latent Gamma variables, `Gamma(shape = alpha, rate = 1)` has
mean `alpha`, so the mean is finite, positive, and available from already
validated hyperparameters. For a Dirichlet latent vector, use those means and
apply one common rescaling to an order-one magnitude. The rescaling preserves
the mean proportions `alpha / sum(alpha)` while avoiding a second underflow
during normalization. Emit an immediate warning saying that RNG initialization
failed and a deterministic mean-proportion initialization was used. If the
fallback itself is not finite and strictly interior, error rather than applying
another arbitrary constant.

This should be a shared **initialization** helper, not a package-wide rule that
every distribution falls back to its mean: some distributions have no finite
mean, and some have a mean on a support boundary.

**Status:** Implemented and verified.

**Decision update:** ok, maybe change to median since it will be defined even in cases mean is not.

## NF02. Near-diagonal covariance treated as diagonal

**Origin.** Commit `33a6dad` made an existing tolerance scale-invariant.

**Current behavior.** `R/random-effects-marginal-covariance.R` classifies a row
covariance as diagonal when every off-diagonal entry is no larger than
`sqrt(.Machine$double.eps) * sd_i * sd_j`. On the correlation scale, correlations
up to roughly `1.49e-8` are therefore treated as zero.

**Mathematical impact.** Scale invariance is an improvement over an absolute
covariance tolerance, but the test still turns small real correlations into
structural zeros. If this status authorizes a diagonal-only representation, the
returned object omits covariance that exists in the represented model. This is
not merely a floating-point repair unless the off-diagonal entries are known to
arise from roundoff around a structurally diagonal construction.

**Recommended.**

- Make structural diagonality the default requirement. Prefer construction
  metadata proving independence; otherwise require exact numeric off-diagonal
  zeros.
- If approximate diagonality is useful for storage or speed, expose it as an
  explicit option such as `allow_approximate_diagonal`, with a user-configurable
  correlation tolerance.
- Return the maximum omitted absolute covariance and correlation in metadata, and
  warn or error when a requested error budget is exceeded.
- Do not reuse this approximation status when an exact covariance is required for
  inference.

**Decision:** I need a bit more detail, when does this happen? does it alter user input? is it used to identify the type based on other input?

**Review response:** This check occurs only in
`random_effects_marginal_variance_factors()`, whose default
`require_diagonal = TRUE` protects downstream likelihoods that can consume only
one variance multiplier per observation. It is not used to choose the fitted
random-effect type, and it does not mutate the user's covariance input.

The helper first derives an observation-space covariance factor:

- without a known group kernel, rows belonging to the same group have covariance
  `z_i * z_j`;
- with a known group kernel `K`, rows have covariance
  `K[group_i, group_j] * z_i * z_j`.

It then decides whether the full row covariance can be represented by its
diagonal alone. For example, two random-intercept observations from the same
group have factor

```text
[1 1]
[1 1]
```

and must be rejected by a diagonal-only consumer. With the current tolerance, a
real but very small off-diagonal correlation (up to about `1.49e-8`) is instead
classified as zero, so the helper can authorize returning only the diagonal.
That is the potential information loss. The separate
`random_effects_marginal_vcov(diagonal_only = TRUE)` API merely answers an
explicit request for marginal variances and does not claim that the full
covariance is diagonal; it is not the problematic path.

**Suggested decision:** Require structural diagonality when that fact is
available from the group map/kernel, and otherwise require exact numeric
off-diagonal zeros when `require_diagonal = TRUE`. Retain
`require_diagonal = FALSE` for callers that knowingly want only the row
multipliers despite correlated rows. Do not add an approximate default.

**Status:** Implemented and verified.

**Decision update:** ok, accept the strict policy

## NF03. Gaussian KDE tail extrapolation

**Origin.** Commit `4b8b3eb` stopped treating the finite KDE evaluation grid or
sample range as exact compact support.

**Current behavior.** When a finite null value lies outside the grid produced by
`stats::density()`, `.density_kde_gaussian_height()` evaluates the Gaussian KDE
formula directly. Boundary-reflected KDEs use the corresponding reflected
kernels. Values outside known exact support are still assigned density zero.

**Mathematical impact.** This is mathematically faithful to the selected Gaussian
KDE: Gaussian kernels have nonzero tails everywhere. It is not an exact posterior
ordinate, and a far-tail Savage-Dickey estimate can be dominated by bandwidth
choice and floating-point underflow. Replacing it with zero at the plotting-grid
edge would be less faithful, because that edge is arbitrary.

**Recommended.** Retain direct Gaussian-kernel evaluation, but make estimator
uncertainty visible:

- record whether extrapolation occurred, the null's distance from the nearest
  sample and sample range in bandwidth units, the bandwidth rule, and whether
  boundary reflection was used;
- retain the current warning when samples do not span the null;
- add a strict policy that errors instead of extrapolating, while keeping `warn`
  as the compatibility default;
- do not introduce a fixed distance cutoff at which the ordinate is set to zero;
- prefer a valid precomputed ordinate when one is available.

**Decision:** keep behavior, add immediate warning if the extrapollation happens 

**Review response:** Understood. At the moment the warning can be accumulated
and emitted later. The instruction is to emit one
`warning(..., immediate. = TRUE, call. = FALSE)` at the point the code chooses
direct Gaussian-kernel evaluation outside the finite KDE grid. The returned
diagnostic metadata should still record that extrapolation occurred and its
distance in bandwidth units.

**Status:** Implemented and verified.

## NF04. Nonzero hypothesis coefficients classified as exact zero

**Origin.** Pre-existing behavior exposed by the current branch's expanded
hypothesis and prior-density paths.

**Current behavior.** `.prior_linear_density_zero_tol()` returns
`sqrt(.Machine$double.eps)`, about `1.49e-8`. Several functions use this value to:

- remove coefficients from active linear-combination groups;
- treat a weighted scalar prior as a point mass at zero;
- treat a vector combination with a small norm as a point mass at zero;
- decide which prior columns or conditional components participate; and
- classify a very small PET/PEESE scale as exactly zero.

The posterior sample calculation does not consistently make the same
replacement. Consequently, the posterior quantity can include a small
coefficient while the derived prior treats it as zero.

**Mathematical impact.** This directly changes the tested quantity. For example,
`1e-9 * beta` is a rescaled continuous random variable, not a point mass at zero.
Its scale may be small, but its density at a null and its Bayes factor do not in
general equal those of an exact zero. The fixed threshold is also not invariant
to unit changes or algebraically equivalent reparameterizations.

**Recommended.**

1. Use `weight == 0` for structural absence in hypothesis semantics.
2. Keep grid-coordinate coalescing tolerances separate from coefficient
   semantics; they solve different problems and must not share one helper.
3. For a very small nonzero coefficient, compute the scaled distribution using a
   stable change of variables or log scale. If the requested ordinate cannot be
   represented, return an explicit numerical-range error or approximation
   diagnostic rather than a point mass.
4. Add invariance tests showing that multiplying a quantity and its null by a
   nonzero constant gives the corresponding equivalent inference, including
   coefficients below `sqrt(machine epsilon)`.

This is the highest-priority fidelity correction in this report.

**Decision:** agree, it would be incorrect to using tolerance for this task, could/shall we just use information about the structure (i.e., something is multiplied by a coefficient of 0 as set by the user -- spike(0) -- to proceed with those removals rather than using these odd tolerance hacks?)

**Review response:** Yes. Use exact structural information and never numerical
proximity:

- if the user-specified coefficient in a linear expression is exactly `0`, the
  term can be removed algebraically;
- if a design/compiler object records an exact structural zero, that metadata
  can be used;
- a `spike(0)` is slightly different: it is a model component with a point mass
  at zero, not an instruction to remove the parameter from every component. Its
  zero contribution may be simplified within the spike component, but its
  posterior/prior component mass must remain in the mixed-measure metadata;
- every nonzero coefficient, however small, stays in the calculation.

Thus, exact coefficient-zero simplification and point-mass propagation should
have separate helpers. No tolerance should define either semantic condition.

**Status:** Implemented and verified.

## NF05. Finite grids, tail truncation, clipping, and renormalization

**Related existing decision.** This refines D15 in
`current-branch-review-decisions.md`.

**Current behavior.** Deterministic `prior_linear_density` calculations use:

- a default grid of 4096 points;
- a default tail probability of `1e-4`;
- a 1024-point product grid;
- truncation to finite source-quantile ranges;
- clipping negative FFT artifacts to zero;
- replacement of non-finite density evaluations by zero; and
- renormalization of the remaining finite-grid mass to one.

**Mathematical impact.** These choices define a numerical approximation, not the
exact transformed prior. Renormalization hides omitted tail mass, and assigning
non-finite endpoint evaluations zero can suppress integrable singularities.
Derived prior ordinates and region probabilities can therefore depend on grid
size and tail probability. That dependence is especially consequential for
Savage-Dickey and hypothesis Bayes factors.

Clipping a tiny negative FFT artifact is defensible only when its size is bounded
by a numerical error estimate. Clipping every negative result without measuring
the discrepancy can conceal a failed approximation.

**Recommended.**

- Carry exact support and source PDF/CDF or transformation provenance with the
  grid.
- Track intended tail mass, captured numerical mass, clipped negative mass,
  normalization factor, and a grid-refinement error estimate.
- Adaptively extend/refine the grid at the queried point and until the requested
  ordinate or probability meets an explicit tolerance.
- Do not silently reinterpret an outside-grid point as outside support.
- If provenance is insufficient to extend the approximation, return an
  `"outside numerical approximation range"` error.
- Permit renormalization only as an explicitly reported quadrature correction,
  not as a way to erase an unknown amount of omitted mass.

**Decision:** I need a small example to understand where/how this happens and what is the impact

**Review response:** A minimal example is
`theta = beta_1 + beta_2`, with independent standard-normal priors. The exact
prior is `theta ~ Normal(0, sqrt(2))`.

With the current default `tail_prob = 1e-4`, each source grid ends at
approximately `+/-3.719`. Each grid captures probability `0.9998` and is then
renormalized as though that finite interval had probability one. Convolution
therefore uses two truncated-and-renormalized normals and produces a derived
grid ending at approximately `+/-7.438`.

Two consequences illustrate the scale of the issue:

- At `theta = 0`, the exact density is `0.28209479`. Even before finite-grid
  discretization, convolving the two truncated-and-renormalized sources gives
  approximately `0.28220762`, about `0.04%` too high. This is small near the
  center under the defaults, but it is a systematic consequence of silently
  renormalizing omitted tails.
- At `theta = 8`, the exact density is approximately `3.17e-8`, but the current
  finite-grid lookup returns zero because 8 lies beyond the numerical grid.
  The exact normal still assigns approximately `1.44e-7` total probability
  beyond the two grid endpoints. In a point-null Bayes factor, zero versus a
  small positive prior ordinate is a qualitative, not merely cosmetic,
  difference.

The effect can be larger for heavy-tailed, product, or nonlinear transformed
priors. Increasing the default grid alone only moves the arbitrary boundary.

**Policy options:**

1. **Conservative minimum:** store omitted/captured mass and grid provenance,
   report refinement diagnostics, and error with `"outside numerical
   approximation range"` rather than returning zero outside the grid.
2. **Adaptive evaluation (recommended):** do option 1 and, when an ordinate or
   probability is requested, extend/refine the grid until a documented error
   criterion is met; error if it cannot converge.
3. **Defaults only:** use wider/finer default grids and document the
   approximation. This reduces common-case error but does not resolve the
   mathematical ambiguity, so it is not recommended as the sole change.

**Status:** Implemented and verified.

**Decision update:** use option 2, make sure to add tests and verify on a couple of examples when implementing because this might be non-trivial and error prone

## NF06. Endpoint-singular prior densities replaced or avoided

**Origin.** Pre-existing behavior encountered through newly expanded derived
prior and random-effect summary paths.

**Current behavior.**

- One ordered-Dirichlet density integral starts at `.Machine$double.eps` rather
  than zero and uses `rel.tol = 1e-5`.
- Scaled-Beta summary densities may move both endpoints inward by
  `sqrt(.Machine$double.eps) * scale`.
- Infinite density values on the resulting grid are replaced by zero.

**Mathematical impact.** A Beta-family density can be integrably infinite at an
endpoint. Its mathematical value there is not zero. Moving the integration bound
inward may be a sound quadrature technique only if the omitted interval is
handled analytically or its mass is bounded. Replacing infinity with zero changes
the density shape and is particularly dangerous if the attached prior density is
later used in a point-null Bayes factor.

At a singular support boundary, the ordinary Savage-Dickey density ratio may be
nonregular or undefined. Producing a finite result by changing the endpoint
density is not a valid repair.

**Recommended.**

- Integrate on a transformed open-domain scale, such as logit or a
  singularity-removing substitution, so the complete interval is represented.
- Preserve exact support and the fact that a boundary density is infinite.
- Use one-sided limiting calculations only when the statistical method defines
  them.
- Reject ordinary Savage-Dickey evaluation at a nonregular/singular boundary
  with a targeted explanation.
- Report integration error and fail when the requested accuracy is not achieved;
  do not turn non-finite integrand values into zero generically.

**Decision:** are we overengineering this a bit? I would assume that we can get close forms for those priors from the beta functions and we won't require numerical integration?

**Review response:** The original recommendation was too broad. There are two
different paths:

1. For the random-effect summary transformations, there is already a closed
   form. If `W ~ Beta(a, b)` and `Y = sW`, then
   `f_Y(y) = dbeta(y / s, a, b) / s`. If `Y = sqrt(sW)`, then
   `f_Y(y) = dbeta(y^2 / s, a, b) * 2y / s`. No numerical integration is needed.
   The current problem is that these formulas are evaluated only on a finite
   grid and an infinite endpoint value is replaced by zero. The maintainable fix
   is to retain a callable analytic evaluator and exact support; use a finite
   grid only for plotting, and reject an ordinary Savage-Dickey ratio at a
   singular boundary.
2. The ordered-prior product path is genuinely different. It has
   `Y = T * C`, where `C` is Beta but `T` may be an arbitrary supported prior.
   Its density contains
   `integral f_T(y / c) * f_C(c) / c dc` over `0 < c < 1`. There is no universal
   Beta-only closed form because `f_T` is arbitrary. This path currently serves
   direct ordered-prior density/plot evaluation rather than changing the JAGS
   target. Numerical quadrature remains appropriate there, but it should use a
   transformed domain and check its integration error instead of silently
   dropping an epsilon-sized interval.

**Narrowed recommendation:** Use the analytic evaluator for the scaled-Beta
summary path. Keep checked transformed quadrature only for the arbitrary-prior
product path, and treat that lower-impact plotting path as a later priority.

**Status:** Implemented and verified.

**Decision update:** accept narrowed recommendation

## NF07. Silent positive-semidefinite projection in prediction

**Origin.** Pre-existing prediction helper uncovered during random-effect review.

**Current behavior.** `.bt_random_effect_mvn_group_draws()` performs a symmetric
eigendecomposition and applies `pmax(eigenvalues, 0)` before constructing draws.
Every negative eigenvalue, regardless of magnitude, is silently set to zero.

**Mathematical impact.** For a matrix that is materially indefinite, this samples
from a different covariance matrix. The repair can hide upstream errors in
covariance reconstruction and produce overconfident or lower-rank predictions.
For a negative eigenvalue of roundoff magnitude, setting it to zero is a normal
floating-point safeguard.

**Recommended.**

1. Validate finiteness and symmetry before decomposition.
2. Derive an eigenvalue tolerance from matrix dimension, machine precision, and a
   matrix norm.
3. Clamp only eigenvalues in `[-tolerance, 0)` and record the largest correction.
4. Error with the minimum eigenvalue and tolerance when any eigenvalue is below
   `-tolerance`.
5. Where positive definiteness is mathematically required, use a Cholesky
   factorization and fail on genuine singularity rather than silently lowering
   rank.

Do not use a generic nearest-positive-semidefinite projection by default; that is
a modeling change unless explicitly requested.

**Decision:** agree

**Review response:** Understood. Only negative eigenvalues whose magnitude is
bounded by a derived floating-point error estimate may be clamped to zero;
material indefiniteness must error rather than change the covariance.

**Status:** Implemented and verified.

## NF08. Marginal-likelihood failure assigned zero model probability

**Related existing decision.** This is D13 with a concrete default policy.

**Current behavior.** `.model_averaging_margliks()` replaces `NA` log marginal
likelihoods by `-Inf`. Downstream model averaging then gives those models
posterior probability zero. A genuine out-of-support parameter evaluation may
correctly have log density `-Inf`; a failed marginal-likelihood computation is a
different event.

**Mathematical impact.** Computational failure becomes infinitely strong evidence
against the model. Posterior model probabilities and inclusion Bayes factors can
change silently, sometimes substantially.

**Recommended.**

- Keep `NA`/error as a classed computational failure.
- Default to aborting model averaging whenever a positive-prior-probability model
  lacks a valid finite or legitimate `-Inf` marginal likelihood.
- Expose an explicit `on_failure` policy:
  - `"error"`: recommended default;
  - `"drop"`: remove failed models and renormalize both prior and posterior model
    spaces, with a prominent warning and returned audit metadata;
  - `"zero"`: deliberately assign zero evidence, permitted only by explicit
    request and labeled as an assumption.
- Preserve legitimate `-Inf` results only when the evaluator can identify them as
  a valid mathematical result rather than a catch-all failure sentinel.

**Decision:** agree

**Review response:** Understood. Computational failure remains distinct from a
legitimate mathematical `-Inf`, and the default is to abort model averaging.
Any `"drop"` or `"zero"` behavior must be explicitly requested and audited in
the returned object.

**Status:** Implemented and verified.

## NF09. Positive-probability components forced to receive a sample

**Related existing decision.** This corrects the description in D14. The current
behavior changes component **counts**, not the numeric values of already drawn
samples.

**Current behavior.** `.mixture_sample_counts(..., preserve_positive = TRUE)`
starts with deterministic largest-remainder allocation and transfers one draw to
every positive-probability component when possible. For example, with 10,000
draws, a component of probability `1e-12` can receive one draw and therefore an
empirical weight of `1e-4`.

**Mathematical impact.** The resulting unweighted sample is not a random sample
from the stated finite mixture and can dramatically overweight very small
components. Posterior quantiles, tail probabilities, densities, and plots based
on that sample may be affected.

**Recommended.**

- For inferential samples, draw model/component indicators from the exact
  categorical distribution, equivalently via `rmultinom()` or `sample.int()`.
  A rare component may legitimately receive zero observations in a finite sample.
- If deterministic allocation is required, use largest remainder without forced
  positivity and document that the result is a deterministic approximation.
- Offer forced representation only as an explicit stratified diagnostic or plot
  option. Such output should carry observation weights restoring the requested
  mixture probabilities and must not be labeled as an ordinary posterior sample.
- Preserve seed control and test empirical component frequencies over repeated
  samples.

**Decision:** agree, issue that needs to be solved, I think that original reason for smth like this 
was that we use the mixing functions across multiple paramaters from the same model and we were concerned 
of having the same "subset" from the models when mixing them together. Veryfy that the recommencded changes (1) does not create an issue with it

**Review response:** Verified. The joint-alignment requirement is real, but it
does not require forced-positive counts.

The current code gives each parameter the same seed and deterministic component
counts, so corresponding parameters ordinarily receive matching `models_ind`
and `sample_ind`. Replacing the count calculation independently inside each
parameter helper would be fragile. Instead, exact categorical mixing should be
implemented with one shared sampling plan:

1. In `mix_posteriors()`, sample the model/component index for every output row
   once from the exact categorical probabilities.
2. For each selected model, sample the source posterior row indices once.
3. Reuse that same `(models_ind, sample_ind)` plan for every parameter belonging
   to the same mixture.

This preserves both requirements: all parameters in an output row come from the
same selected model, and they use the same within-model posterior row, thereby
retaining within-model joint dependence. A very rare positive-probability model
may still receive zero output rows, as it should in a finite categorical sample.

If a conditional operation gives some parameters a genuinely different model
probability vector, those parameters cannot share the same plan without changing
their distribution. The implementation should create one plan per identical
mixture definition and align only the parameters that mathematically share it.

Tests should assert matching model and source-row indices across multiple
parameters, preservation of a known within-model row relationship, possible
absence of a rare component, and correct empirical frequencies across repeated
seeds.

**Status:** Implemented and verified.

## NF10. Symmetry, diagonal, and covariance repair policies

**Related existing decision.** This supplies a concrete direction for the
symmetry portion of D10.

**Current behavior.**

- Known group covariance kernels are accepted as symmetric under a fixed
  `all.equal(..., tolerance = 1e-10)` rule and then replaced by
  `(K + t(K)) / 2`.
- CAR distance matrices are accepted under separate fixed `1e-12` symmetry and
  diagonal rules, but the accepted matrix is not canonicalized to exact symmetry
  and exact zero diagonal.
- Other covariance paths use still different support and representability rules.

**Mathematical impact.** A fixed absolute/relative hybrid tolerance is sensitive
to scale and dimension. One path changes accepted input, while another accepts a
nearby invalid representation and passes it onward unchanged. Material
asymmetry can therefore be hidden for small-scale matrices, and equivalent
rescaled inputs can receive different treatment.

**Recommended.**

1. Define a package-wide, scale-aware matrix comparison tolerance based on
   machine precision, dimension, and a matrix norm.
2. Distinguish exact model requirements from numerical representation:
   symmetry, zero diagonal, and support remain exact mathematical conditions.
3. When deviations are within the derived roundoff bound, canonicalize to exact
   symmetry/zero diagonal and record the correction.
4. When deviations exceed the bound, error; do not silently project.
5. Apply the same policy to constructors, JAGS data, bridge reconstruction,
   prediction, and marginal covariance.

An optional explicit repair function may handle intentionally noisy empirical
matrices, but it should return the repaired matrix and diagnostics before model
fitting rather than repair silently inside the fit.

**Decision:** I need more infor, when do those repairs happen in practice? what is the source of these issues?

**Review response:** There are two separate situations, and only one normally
changes user-derived input:

1. `random_group_covariance()` accepts a user-supplied group covariance/kernel.
   After ordering/scaling it, the validator accepts near-symmetry under the fixed
   `1e-10` rule and replaces the matrix by `(K + t(K)) / 2`. For example,
   off-diagonal entries `0.3` and `0.300000000001` are replaced by
   `0.3000000000005`. Plausible sources are an externally calculated empirical
   kernel, text/serialization round trips, or manual entry. Standard symmetric
   constructors such as `stats::cov()` and symmetric distance formulas normally
   produce bitwise-matching pairs, and scaling an already symmetric matrix does
   not create asymmetry.
2. CAR distances are built internally from user-supplied time coordinates using
   `abs(outer(time, time, "-"))`. That construction is exactly symmetric with an
   exact zero diagonal. The separate `1e-12` acceptance rule is therefore
   effectively redundant in normal public use; a failure indicates an internal
   reconstruction bug. This path accepts but does not currently repair the
   matrix.

Thus, these tolerances do not identify a model type. The first is a constructor
validation/repair of a user-supplied numerical kernel; the second validates an
internally derived CAR representation.

**Suggested common part:** Assert exact symmetry and exact zero diagonal for the
internally built CAR matrix, because any discrepancy is an internal error. For a
user-supplied group kernel, validate and canonicalize at the public constructor
only, so every downstream path sees one authoritative matrix.

**Remaining policy choice for user kernels:**

1. **Warned roundoff canonicalization (recommended):** derive a scale- and
   dimension-aware roundoff bound, average the two triangles only within that
   bound, emit an immediate warning, and store the maximum correction. Error
   beyond the bound.
2. **Exact rejection:** require the user-supplied matrix to be exactly symmetric
   and ask the user to repair it explicitly. This is simpler and maximally
   transparent, but less convenient for externally computed matrices.

Neither option silently repairs material asymmetry.

**Status:** Implemented and verified.

**Decision update:** option 2 (but if user supplies only lower / only upper triangle matrix, do not average and use the actually speciefied triangle)

## NF11. Savage-Dickey atoms, null equality, and support topology

**Origin.** Mixed current-branch and pre-existing behavior around new posterior
density/support metadata.

**Current behavior.**

- Raw posterior draws are treated as invalidating the density ratio only when
  more than 5% equal the null exactly.
- Declared posterior point masses within
  `sqrt(machine epsilon) * scale` of the null are treated as located at the null.
- Exact support boundaries and point supports use a similar tolerance.
- Support intervals separated by less than that tolerance can be treated as
  connected.
- Posterior ordinate values separated by that tolerance can be rejected as
  duplicates.

**Mathematical impact.** Any genuine point mass at the null invalidates the
ordinary continuous-density Savage-Dickey ratio, even if its mass is below 5%.
Conversely, repeated MCMC values do not by themselves prove a mathematical atom.
A point at `null + 1e-9` is not a point at `null`, and a small but real gap in
support remains a gap. Approximate identity can therefore change whether a Bayes
factor is declared valid and whether exact support excludes the null.

**Recommended.**

- Make declared model/posterior metadata authoritative for atoms and exact
  support.
- Reject the ordinary density ratio for any declared positive atom at the exact
  null; use a mixed-measure Bayes-factor method when atoms are intended.
- Compare semantic null values and point locations exactly after one documented
  canonical serialization step. Do not use proximity to redefine the null.
- Treat repeated raw draws as a diagnostic ambiguity. Without atom metadata,
  warn or require an explicit density policy rather than using a fixed 5%
  definition.
- Preserve support components and gaps exactly. Use tolerance only to assess
  whether floating-point samples are compatible with declared support, never to
  change support topology.

**Decision:** I think the original issue was that we didn't have enough meta-data to know whether 
the created mixed posterior samples have a point at the null value or not. check whether we are 
able to make this fully meta-data driven now. If so, this would allow us to make a proper treatment of these cases.

**Review response:** The metadata infrastructure is now capable of representing
the distinction, but the complete pipeline is **not yet fully metadata-driven**.

What is already available:

- `posterior_density_attribute()` can declare exact atom locations and masses in
  `point_masses`;
- it distinguishes “point masses were explicitly declared and there are none”
  from “point-mass status was not supplied” through
  `point_masses_declared`;
- Savage-Dickey and plotting paths already inspect declared point masses when
  present.

What is still missing:

- ordinary `mix_posteriors()` and `as_mixed_posteriors()` outputs retain
  `prior_list`, `models_ind`, and `sample_ind`, but do not consistently convert
  their known point components into canonical posterior atom metadata;
- transformations currently drop stored posterior density/ordinate metadata
  rather than transforming atom locations and preserving their masses;
- conditional and linear-combination paths do not yet propagate all joint atom
  information;
- `Savage_Dickey_BF()` still applies
  `mean(posterior == null_hypothesis) > 0.05`, even when repeated numeric draws
  are not authoritative evidence of an atom.

There is no conceptual blocker to making it fully metadata-driven. The required
implementation is:

1. At each mixture producer, derive atoms structurally from point-prior
   components and the model/component mapping. For across-model mixtures, use
   posterior model probabilities for atom mass rather than the realized finite
   sample count; this remains correct even when NF09 legitimately draws zero
   rows from a rare point component.
2. Attach an explicit declaration of no atoms when the structure proves the
   posterior is continuous.
3. Propagate atom locations and masses through transformations and conditioning.
   For linear combinations, use joint component metadata; multiplying marginal
   atom probabilities would be incorrect when components are dependent.
4. Compare a declared atom location with the semantic null exactly after one
   canonical numeric representation. Any positive declared atom at that exact
   null rejects the ordinary density ratio, regardless of mass.
5. Remove the 5% rule. Repeated raw values may be reported as a diagnostic, but
   they must not define the measure type.

One compatibility choice remains for legacy or manually constructed posterior
vectors whose atom status is unknown:

1. **Strict (recommended):** point-null Savage-Dickey evaluation errors when
   atom metadata is required but unknown, with instructions to supply an
   explicit atom/no-atom declaration.
2. **Compatibility:** warn that atom status is unknown and continue with the
   continuous KDE at the user's risk. This still removes the arbitrary 5%
   threshold.

**Status:** Implemented and verified.

**Decision update:** make the instractructure fully metadata driven, then make it strict
this is a larger change but well worth the improvements

## NF12. Weight-function cuts coalesced or observations snapped to cuts

**Origin.** Pre-existing model-averaging and selection mapping behavior.

**Current behavior.**

- `.weightfunction_unique_cuts()` merges cuts separated by no more than
  `sqrt(.Machine$double.eps)`.
- `.weightfunction_global_bin_indices()` expands bin containment by the same
  tolerance.
- `.selection_native_step_bin_from_z()` snaps a computed p-value to the first cut
  within `1e-12` before assigning a bin.

**Mathematical impact.** Distinct user-specified cuts define distinct selection
models, even when their bin is very narrow. Merging cuts removes that bin.
Snapping can assign a boundary-near observation to a different weight, and
choosing the first nearby cut is order-dependent when more than one cut is in
range. This is a model-definition change, not merely a plotting approximation.

**Recommended.**

- Preserve every distinct representable user cut.
- Define one interval convention, for example `(lower, upper]`, and use it in R,
  generated JAGS code, native evaluators, summaries, and tests.
- Derive bin membership from original cut indices or evaluate observations and
  cuts on the same stable scale. Do not snap by absolute distance.
- If two requested cuts transform to the same floating-point boundary, report
  that the requested bin is not computationally representable and require the
  user to change it.
- Keep display-label rounding separate from stored numeric cuts.

**Decision:** agree

**Review response:** Understood. Preserve every distinct representable cut and
use one exact interval-closure convention throughout R, JAGS, and native code.
Numerical display rounding remains separate from stored/model cut values.

**Status:** Implemented and verified.

## NF13. Near-one-hot matrices routed through an exact one-hot compiler

**Origin.** The branch's scalable group-local structured random-effect
compilation.

**Current behavior.** `.bt_random_effect_structured_local_layout()` treats model
matrix entries with magnitude at most `sqrt(machine epsilon)` as zero and accepts
the selected entry when it is within that tolerance of one. The compact layout
stores the selected column indices, not the small residual entries or the
near-one multiplier.

**Mathematical impact.** If a matrix entry is genuinely `0.99999999`, or a second
entry is genuinely small but nonzero, the compact model represents a different
design: exactly one selected coefficient with multiplier one. Generated
factor-indicator matrices are normally exact zero/one, so the issue is avoidable
without sacrificing the fast path.

**Recommended.**

- Prefer structural provenance: use the compact compiler only when the design
  builder declares that the matrix is an exact indicator matrix.
- Otherwise require exact finite zero/one entries for the compact path.
- If the condition is not met, fall back to the general dense/sparse design that
  preserves the supplied numeric values; do not quantize the matrix.
- If roundoff canonicalization is needed for a known algebraic source, perform it
  once at design construction, record the maximum correction, and test equality
  of compact and general likelihoods.

**Decision:** how about: use the compact compiler only when the design
  builder declares that the matrix is an exact indicator matrix and never use the compact path otherwise? 
  (if there is ever a sample that is almost 0/1 by accident -- i.e., not because of the structure -- the efficiency gain is more than paid for by complexity and maintenance of the code)

**Review response:** Understood and agreed. The compact path will require trusted
structural provenance from the design builder, for example an
`exact_indicator = TRUE` declaration tied to the builder output. It will not
inspect arbitrary numeric values and opt in merely because they look like
zero/one. If a trusted builder declares the structure but produces a non-finite
or non-exact zero/one value, that is an internal consistency error rather than a
value to quantize. Every object without the declaration uses the general path.

**Status:** Implemented and verified.

## NF14. Fixed tolerances for simplex membership

**Origin.** Pre-existing Dirichlet, ordered-allocation, mixed-posterior, and
random-effect reconstruction behavior.

**Current behavior.** Several paths accept vectors whose sum differs from one by
up to `1e-8` or `sqrt(machine epsilon)`. Examples include Dirichlet log-density
evaluation, fixed ordered allocations, simplex point priors, and reconstructed
Dirichlet allocation samples. Accepted vectors are not consistently normalized.

**Mathematical impact.** A Dirichlet distribution is defined on the simplex.
Using a materially off-simplex vector in its density is not a harmless
approximation. At the same time, decimal user input and floating-point
reconstruction cannot be expected to sum to bitwise exactly one.

**Recommended.**

1. Define one simplex canonicalization helper.
2. Derive its roundoff bound from machine precision, dimension, and
   `sum(abs(x))`; do not use a fixed `1e-8`.
3. For user-specified fixed vectors within that bound, normalize once, retain the
   original sum and correction in metadata, and use the canonical vector
   everywhere.
4. For reconstructed posterior values, prefer a parameterization that satisfies
   the simplex by construction. If a final normalization is needed, permit only
   roundoff-sized correction and record/check its maximum.
5. Reject negative coordinates or sum discrepancies beyond the derived bound.
6. Ensure bridge densities evaluate the same canonical values that the fitted
   model represents.

**Decision:** ok

**Review response:** Understood as acceptance of the recommendation: one
scale-aware canonicalization at the input/reconstruction boundary, with only
roundoff-sized correction and rejection of material off-simplex values.

**Status:** Implemented and verified.

## NF15. Fixed equality tolerances in fitted-model replay

**Origin.** Current-branch bridge-sampling and random-design reconstruction
checks.

**Current behavior.** Several fitted-versus-rebuilt model matrices, scale
metadata, source data, and correlation bounds are considered equal using
`all.equal(..., tolerance = 1e-12)`. If accepted, bridge sampling or
reconstruction proceeds even though the objects need not be identical.

**Mathematical impact.** A tiny difference is often inconsequential, but a fixed
tolerance has no direct relationship to the resulting log likelihood or log
marginal likelihood. It is also scale-dependent. For bridge sampling, proceeding
with data or a design that differs from the fitted model is conceptually more
serious than a display mismatch.

**Recommended.**

- Treat stored fitted design/data values as authoritative whenever they are
  available; reconstruction should recover identity, not choose a nearby
  replacement.
- Compare categorical structure, names, contrasts, indices, and formulas exactly.
- For values intentionally recomputed by deterministic floating-point
  transformations, use a scale-aware forward-error bound specific to that
  transformation and report the maximum discrepancy.
- If exact fitted values are unavailable and the effect of the discrepancy on
  the log joint cannot be bounded, stop rather than proceed under a generic
  `1e-12` rule.
- Consider storing canonical design fingerprints plus the values needed by the
  marginal-likelihood evaluator.

**Decision:** ok, but a bit of a nitpick

**Review response:** Understood as accepted but deliberately lower priority.
Implementation should remain small: prefer authoritative stored fitted
data/design values, compare structural metadata exactly, and replace a generic
fixed tolerance only where replay genuinely recomputes floating-point values.
No broad fingerprinting framework is required unless the existing stored values
prove insufficient.

**Status:** Implemented and verified.

## NF16. Plot-only offsets, spike matching, and numerical integration tolerances

**Origin.** Pre-existing plotting code.

**Current behavior.** Examples include:

- moving transformed plotting bounds inward by a fixed `1e-5`;
- classifying posterior factor samples within `1e-6` of one as the spike at one;
- grouping plotted point locations under scale-based tolerances;
- fixed `1e-8` quantile stopping rules and `1e-7` integration tolerances; and
- clipping computed plotting CDF values to `[0, 1]`.

**Mathematical impact.** These paths do not alter fitting or model probabilities,
but they can visually merge a narrow continuous component with a spike, omit a
meaningful interval on a small scale, or show a slightly inaccurate quantile.
Clipping a CDF overshoot of roundoff size to its known mathematical range is
faithful; using the same clipping to hide a failed integral would not be.

**Recommended.**

- Identify spikes and atoms from component metadata, not numeric proximity.
- Move open plotting bounds to the next representable interior value, or use the
  transformation's mathematical domain, instead of adding `1e-5`.
- Make root/integration tolerances relative to the plotted scale and expose
  convergence/error diagnostics internally.
- Clamp to known mathematical ranges only when the excess is within a derived
  numerical error bound; otherwise warn or fail the plot-data computation.

This is lower priority because the effect is presentational.

**Decision:** agree

**Review response:** Understood. This is presentational work: derive spike
identity from metadata, use representable interior endpoints, and ensure
numerical clamps are bounded by reported computation error.

**Status:** Implemented and verified.

## NF17. Automatic centered/noncentered parameterization thresholds

**Origin.** Current random-effect parameterization infrastructure.

**Current behavior.** The automatic centered-design decision includes
`max_columns = 8`, `min_effective_n = 5`, and `min_rcond = 1e-4`.

**Mathematical impact.** When centered and noncentered implementations are
correct, these thresholds do not change the posterior target. They select a
parameterization intended to improve sampling. They can still affect adaptation,
effective sample size, runtime, and the chance that a finite run meets
convergence criteria.

**Recommended.**

- Retain the heuristic as a computational policy.
- Document that the numbers are tuning defaults rather than statistical
  assumptions.
- Record the resolved parameterization, threshold values, and reason in every
  fitted object and summary diagnostic.
- Preserve an explicit user override.
- Add equivalence tests for the represented log joint and simulation studies for
  the heuristic; tune defaults from those studies rather than presenting them as
  universal constants.

**Decision:** agree

**Review response:** Understood. Retain the heuristic as a computational tuning
policy, preserve the user override, and record both the resolved choice and its
reason without presenting the thresholds as statistical assumptions.

**Status:** Implemented and verified.

## NF18. Native endpoint clamps and safeguards that should remain

**Origin.** Mixed current-branch and pre-existing numerical code.

**Current behavior.** The native LKJ transform clamps a descendant input to
`[0, 1]`, while the stochastic density separately rejects values outside the open
support. Its square-root helper changes values in `(-1e-14, 0)` to zero before
computing `sqrt(1 - x^2)`. Other package safeguards include:

- log-sum-exp normalization of model probabilities;
- survival-tail and `log.p` calculations for extreme probabilities;
- clipping analytically bounded CDF values to `[0, 1]`;
- returning exact truncation bounds for endpoint quantiles;
- selecting representable values immediately inside an open correlation support;
- setting a computed negative variance to zero only when it is within a
  roundoff-scale bound; and
- allowing a machine-precision accumulation error when checking that
  probabilities sum to one.

**Mathematical impact.** These operations generally preserve the intended
mathematical quantity more faithfully than naive arithmetic. The important
condition is that a clamp be derived from a known support or error bound and not
used to convert a materially invalid value into a valid one. The LKJ fixed
`1e-14` rule is the one part that lacks a scale/error derivation.

**Recommended.**

- Retain log-domain calculations, stable tail calculations, exact endpoint
  semantics, and support-derived one-ULP interior bounds.
- Retain probability/CDF validation tolerances that only accept ordinary
  floating-point summation error and do not alter the probabilities.
- Replace the LKJ fixed `-1e-14` allowance with a bound derived from the error in
  `1 - x*x`; error on a materially negative radicand.
- Verify by test that out-of-support stochastic LKJ values always receive
  `-Inf`/rejection before descendant clamping can influence an accepted state.
- Establish the package-wide rule: a silent clamp is allowed only when the
  uncorrected value differs from a known mathematical boundary by a documented
  floating-point error bound. Everything larger must produce a diagnostic.

**Decision:** ok

**Review response:** Understood as acceptance. Faithful log-domain/support
safeguards remain; the fixed LKJ allowance is replaced by a derived arithmetic
error bound, and materially invalid values must not be clamped into support.

**Status:** Implemented and verified.

## Implementation outcome

All maintainer decisions are implemented. No further clarification is needed.
The resulting policies are:

| Case | Implemented handling |
|---|---|
| NF01 | Positive Gamma/Dirichlet initializations use a warned, deterministic median-based interior fallback; unrepresentable cases error. |
| NF02 | Diagonal covariance output requires exact structural diagonality. |
| NF03 | KDE tail extrapolation is explicit in diagnostics and warnings. |
| NF04 | Only exact coefficient zeros are simplified; nonzero coefficients remain in the represented quantity. |
| NF05 | Prior-density grids refine adaptively and report captured mass, tail error, normalization, and FFT clipping diagnostics. |
| NF06 | Scaled-Beta endpoints use analytic semantics; arbitrary products use checked transformed quadrature. |
| NF07 | Covariance eigenvalues are clamped only within a derived roundoff bound; larger violations error. |
| NF08 | Marginal-likelihood computation failures abort by default and require an explicit alternative policy. |
| NF09 | Mixture draws use one shared categorical sampling plan without forced-positive component counts. |
| NF10 | User covariance inputs require exact symmetry, while explicitly lower- or upper-triangular inputs use the supplied triangle. |
| NF11 | Savage-Dickey atom/support identity is strict and metadata-driven, including exact null matching. |
| NF12 | Weight-function cuts and interval assignments use their exact declared values and closure convention. |
| NF13 | Indicator compilation requires exact declared structure; near-one-hot numeric matrices remain dense. |
| NF14 | Simplex canonicalization is limited to roundoff-scale drift; larger discrepancies error. |
| NF15 | Fitted values and structural metadata are authoritative in replay and design comparisons. |
| NF16 | Plot atoms come from component metadata, open endpoints use representable interior values, and mixed measures retain separate probability and density semantics. |
| NF17 | The resolved centered/noncentered parameterization, thresholds, and reason are recorded while preserving user override. |
| NF18 | LKJ/support clamps are bounded by derived floating-point error; materially invalid values are rejected. |

Verification completed with focused adversarial tests, a fresh fitting-cache
run, fixture/reference regeneration, and all package test profiles relevant to
the changed fitting, density, inference, and plotting paths. Detailed final
results:

- unit: 7,720 passed, 0 failed, 0 warnings, 7 profile skips;
- fresh fit/cache rebuild: 14,896 passed, 0 failed, 0 warnings;
- fixture: 9,009 passed, 0 failed, 0 warnings;
- visual: 810 passed, 0 failed, 0 warnings, 11 profile skips;
- visual-fixture: 608 passed, 0 failed, 0 warnings, 26 profile skips;
- targeted `covr` run over seven adversarial test files: 28.91% of the full
  package, with the instrumented tests passing; and
- package check: 0 errors, 0 warnings, 1 already-documented portability note
  for long visual-snapshot paths (D30 in
  `current-branch-review-decisions.md`).
