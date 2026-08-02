# Priors and Density Semantics

Use this guide when changing prior constructors, distribution methods, density
provenance, transformations, mixtures, or prior-density ordinates.

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

## Numerical Evidence

Use distribution-library references, analytic identities, or independently
derived results. Cover interior values, exact support boundaries, outside-
support values, invalid parameters, transformations, mixtures, and meaningful
extreme tails. Test identities such as CDF/CCDF complementarity and
quantile/CDF inversion only where they are mathematically defined.

For prior changes, run `unit` first and `visual` when plots change. Add `fit`
and `fixture` when generated JAGS syntax, initialization, monitoring, marginal
likelihoods, or cached fitted objects can change.
