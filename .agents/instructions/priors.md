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
continuous behavior remains diagnostic provenance. Unsupported transformations
or convolutions return `unknown` rather than a guessed structural class.

The density-provenance implementation is shared across
`R/priors-density-context.R`, `R/priors-linear-density.R`,
`R/priors-linear-density-combinations.R`, `R/priors-density.R`, and
`R/prior-density-ordinate.R`. Extend this path instead of building a parallel
density algebra. Provenance must remain deterministic and compact: do not store
draws, large grids, fitted objects, environments, or closures capturing them.

Conditional-normal quadratures are independent 1-D integrals: each distinct
design row, model or conditional mixture component, and mixture leaf receives
the full evaluation budget (`n_grid`) and its own convergence check. Budgets
are never divided among them; the budget only caps how many mixture leaves are
expanded.

Mixture ordinates (model, conditional, and mixture or spike-and-slab terms of
one linear combination) are weighted sums of per-component ordinates, each from
its own exact or regular method, so a numerical grid never spans a density jump
between components. A Gaussian term plus one other continuous scalar term uses
the conditional-normal quadrature wherever it occurs, so the same combination
is never exact in one context and a grid approximation in another. A mixture
height sums its components' exact or regular heights with those of components
that have none; each of the latter uses its own grid (never spanning another
component's jump), all such grids are refined in lockstep, and the documented
refinement criterion applies to the weighted mixture height.

For scalar prior-region probabilities, integrate the continuous grid's
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
