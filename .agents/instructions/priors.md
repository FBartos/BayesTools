# Priors and Density Semantics

Use this guide when changing prior constructors, distribution methods, density
provenance, transformations, mixtures, or prior-density and Savage-Dickey
ordinates.

## Prior Contract

`prior()` creates the core `prior` S3 family. Distribution aliases and routing
live in `R/priors.R`; constructors and parameter normalization live in
`R/priors-constructors.R`; joint and marginal methods live in
`R/priors-methods-joint.R` and `R/priors-methods-marginal.R`. Specialized
distribution implementations live in `R/distributions-*.R`.

When adding or changing a distribution, keep these paths synchronized:

1. accepted names and aliases in `prior()`;
2. constructor validation, canonical parameters, support, and truncation;
3. applicable `rng`, `quant`, `cdf`/`ccdf`, `pdf`/`lpdf`, and marginal methods;
4. JAGS prior syntax, initialization, monitoring, and bridge-sampling density
   where the distribution is supported for fitting;
5. printing and plotting where applicable;
6. focused unit and visual tests using the existing prior helpers.

Do not claim a joint or marginal method for a prior family unless its
measure and return shape are defined. Point masses, continuous densities,
mixtures, vectors, and factor priors are not interchangeable merely because a
numeric representation can be produced.

Truncation is part of the distribution definition. Preserve exact support and
one-sided boundary behavior; do not replace valid endpoints with nearby values.

`selection_model()` captures fixed source choices and a deferred publication
column reference. Keep `group = NULL` immediately after `weight_rule`; consumers
resolve publication groups only for `"best"` branches. Product selection is
invariant to publication partitions and must not require or bind a group
column. `group` is a deferred column reference rather than a modelling
instruction, and callers parameterize `weight_rule` while forwarding one
`group` for both rules, so a group supplied alongside the product rule is inert
rather than contradictory: keep the combination constructible and let `print()`
report the group as unused. Integration dependencies come from the integrated
sampling and random-effect covariances, independently of publication groups.
Preserve typed reference capture, branch metadata, and source-conditioning
choices.

Generated JAGS syntax and generated initial values are two renderings of one
node layout, so derive both from `.weightfunction_component_node_names()`
rather than rebuilding names on either side. Initial values belong on the
stochastic node; `omega_target` (`omega`, `omega_component_<id>`) is
deterministic and expanded onto the global cut grid, and initializing it either
does nothing or aborts JAGS with a dimension mismatch. Any check that mocks the
model, the sampler, or the backend cannot establish that the two sides still
agree - keep at least one test that compiles the real model.

## Structural Prior-Density Ordinates

`prior_density_ordinate()` classifies the mathematical behavior of a prior at
one exact finite value. Its stable schema and behavior values are public. Keep
`R/prior-density-ordinate.R`, its roxygen documentation, and
`tests/testthat/test-prior-density-ordinate.R` synchronized.

Classification must come from the prior definition and deterministic
provenance. Never infer `regular`, `zero`, `infinite`, `point_mass`, or
`undefined` from:

- prior or posterior samples;
- KDE output;
- a finite grid, interpolation, or FFT clipping;
- probing `value +/- epsilon`;
- binary64 underflow or overflow alone.

A structurally regular ordinate remains regular if its representable density
underflows. Exact point mass takes precedence at the requested value; the
continuous behavior remains diagnostic provenance: every `point_mass` result
(and only such a result) carries `provenance$continuous_behavior`, the
behavior (`regular`, `zero`, `infinite`, `undefined` or `unknown`) of the
measure without its point masses at that value, with `log_density` the log
density of that continuous part (not renormalized). It is documented in the
roxygen `@return`, and consumers (RoBMA's point-test eligibility) read it, so
keep it on every point-mass path, including stored atoms of linear densities.
Unsupported transformations or convolutions return `unknown` rather than a
guessed structural class.

The density-provenance implementation is shared across
`R/priors-density-context.R`, `R/priors-linear-density.R`,
`R/priors-linear-density-combinations.R`, `R/priors-density.R`, and
`R/prior-density-ordinate.R`. Extend this path instead of building a parallel
density algebra. Provenance must remain deterministic and compact: do not store
draws, large grids, fitted objects, environments, or closures capturing them.

Conditional-normal quadratures are independent 1-D integrals: each distinct
design row, model or conditional mixture component, and mixture leaf receives
the full evaluation budget (`n_grid`) and its own convergence check. Each such
integral is split at the multiplier's (other term's) support bounds, at its
declared-prior quantiles (1e-6, 1e-3, .02, .25, .5, .75, .98, and their
complements), and at the location peak of the conditional normal and +-1, 3,
10 local SDs around it when the multiplied SD is at most half of the
absolute value of its mean (always for Gaussian convolutions; beyond that
guard the window is not a peak and is not used). Next to a finite bound with
infinite density, a breakpoint much closer to the bound than the next one
makes QUADPACK count the mass near
the bound twice: there the extreme quantile is skipped, and up to the quartile
on that side a breakpoint is kept only if its distance to the bound is at
least 1e-3 of the next kept breakpoint's distance. Gaussian-peak breakpoints
at least one local SD from the bound are exempt, so a narrow peak next to the
bound keeps its breakpoints (closer peak points, i.e. rounding residues, are
still dropped). Breakpoints closer than
1e-9 of their magnitude (at least 1) are merged, so no piece is only a few
ulps wide, but the merge width never exceeds half the peak's local SD, so the
peak breakpoints are never merged. Every piece gets the full budget, and the
acceptance criterion applies to the summed value and error.
Budgets are never divided among them; the budget only caps how many mixture
leaves are expanded. Known limitations: pure scale mixtures with heavy-tailed
multipliers can stop as non-convergent (mostly with a small multiplied SD); a
piece that QUADPACK flags stops the ordinate even when its value is far below
the absolute tolerance (e.g. a far tail piece reported as "probably
divergent"); a quadrature that evaluates to exactly zero is rejected (it
cannot be told apart from a missed peak), which also stops a mixture with
such a component, e.g. a narrow component far from the value; very narrow
Gaussian peaks can stop as non-convergent where QUADPACK reaches the
floating-point resolution: next to a bound with infinite density away from
zero (e.g. an SD of 1e-10 at the upper bound of a beta), or when only a few
hundred ulps wide; and the scale peak of a value far in the
multiplier's heavy tail (beyond its extreme quantiles) is not a breakpoint,
also with a nonzero multiplied mean when the multiplied SD is not small
against it, so its mass can be missed without a convergence failure. Such
ordinates are small at ordinary scales, but the missed fraction does not
depend on the units of the value, so the absolute tolerance does not bound
it.

Mixture ordinates (model, conditional, and mixture or spike-and-slab terms of
one linear combination) are weighted sums of per-component ordinates, each from
its own exact or regular method, so a numerical grid never spans a density jump
between components. A Gaussian term plus one other continuous scalar term uses
the conditional-normal quadrature wherever it occurs, so the same combination
is never exact in one context and a grid approximation in another. A mixture
height sums its components' exact or regular heights with those of components
that have none; each of the latter uses its own grid (never spanning another
component's jump), all such grids are refined in lockstep, and the documented
refinement criterion applies to the mixture height with the components'
weighted absolute changes (no cancellation between components). A refinement-
change criterion is not an error bound: grid-based heights of components with
singular source densities (e.g. gamma shape < 1) remain approximate within it.

Scalar prior-region probabilities (`R/prior-density-region.R`) follow the
ordinate's classification wherever it has a structural route, for regions
whose relations are linear in the quantity (unions of intervals): point masses
add their exact mass (the exact condition decides strict and inclusive
bounds), scalar priors use their exact distribution function (through a log
source and named monotone output transformations), normal sums their normal
distribution function, and Gaussian convolutions and conditional-normal scale
mixtures the 1-D integral of the other term's density times the Gaussian
probability of the region, with the ordinate's breakpoints (the Gaussian-peak
window, under the same guard and so always for Gaussian convolutions, at every
finite region bound), the full budget per piece, and the same acceptance
criterion on the total; a rejected quadrature stops. Mixture,
spike-and-slab, model and conditional components and distinct design rows are
expanded as for the ordinate and summed with their probabilities; a
combination is evaluated this way only when every component has such a route.
Keep the two paths in step: a combination whose ordinate is structurally
classified has a structural region probability, and one whose ordinate is
`unknown` keeps the grid. A quadrature total may exceed [0, 1] by at most its
absolute error. For a Gaussian convolution G + w T, a piece that lies
entirely outside every window t* +- 10 s / |w| of the finite region bounds is
evaluated exactly as 0 or 1 times its mass under T's declared distribution
function (the Gaussian region probability is constant there to within
2 * Phi(-10) ~= 1.5e-23), so heavy tails and region bounds far beyond T's
quantiles need no quadrature. Known limitation: conditional-normal scale
mixtures integrate every piece, and with a heavy-tailed multiplier (e.g.
half-Cauchy or inverse-gamma(1)) QUADPACK flags the infinite end piece of a
region with an infinite bound ("roundoff error", "probably divergent"), so
such probabilities stop as they did on the grid. A region whose probability
underflows to zero stops like a zero ordinate, and so does a Gaussian
convolution whose region probability comes only from pieces evaluated as 0
(a representable probability below 2 * Phi(-10) of their mass under T, e.g. a
region far below a bounded or nonnegative T).

For the other prior-region probabilities, integrate the continuous grid's
piecewise-linear interpolant up to the exact region boundaries (the comparison
value, or boundaries located by bisection within the cells where a composite
condition changes), normalize by the interpolant's total over the same grid,
and evaluate point masses with the requested strict or inclusive operator. Do
not move a region boundary to nearby grid knots. Keep the existing refinement
gate.

## Savage-Dickey Ordinates

The posterior ordinate of a Savage-Dickey ratio (`R/marginal-savage-dickey.R`)
is a kernel estimate from the continuous draws, but everything around it comes
from declared metadata: posterior atoms, exact support, and for mixtures each
draw's component with the components' exact supports (`posterior_components`:
the model of a `mix_posteriors()` ensemble, or the indicator tuple of the
mixture terms entering the quantity for `as_mixed_posteriors()`). Never infer
atoms, supports, or components from the draws.

- When continuous components have different supports, estimate the ordinate
  per component (boundary-reflected on its own support, zero when that
  support excludes the null, the one-sided limit on a bound) and mix by the
  components' shares of the continuous draws. Shared or unavailable supports
  keep the pooled estimate.
- A null outside the continuous draws gives the finite kernel-tail value with
  a warning that it is not reliable evidence, once per parameter or level;
  never `Inf`.
- A declared posterior point mass at the null leaves the ratio undefined:
  level lists and marginal inference return `NA` with the reason in the
  `"warnings"` attribute and compute the other levels; a scalar call stops.

## Numerical Evidence

Use distribution-library references, analytic identities, or independently
derived results. Cover interior values, exact support boundaries, outside-
support values, invalid parameters, transformations, mixtures, and meaningful
extreme tails. Test identities such as CDF/CCDF complementarity and
quantile/CDF inversion only where they are mathematically defined.

For prior changes, run `unit` first and `visual` when plots change. Add `fit`
and `fixture` when generated JAGS syntax, initialization, monitoring, marginal
likelihoods, or cached fitted objects can change.
