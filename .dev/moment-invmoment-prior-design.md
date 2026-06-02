# Moment and Inverse-Moment Prior Design

Date: 2026-05-31
Restart baseline: `C:/R-Packages/BayesTools-add-random-clean-20260531-205624`
Branch: `0.3.1-add-random`
Verified pushed commit: `3a00e1cd14fbf2c98539dfcef24acea188fda323`

This note records the design decisions for adding Johnson/Rossell nonlocal
moment and inverse-moment priors to BayesTools. It is meant as restart context
for continuing work from the recovered `0.3.1-add-random` worktree. The small
boundary-reflection KDE edits in `C:/R-Packages/BayesTools` are separate and
should not drive the nonlocal-prior implementation.

## Sources and Definitions

The implementation should follow the nonlocal prior parameterizations used by
Johnson and Rossell:

- Moment prior, also pMOM in the scalar Gaussian-base case.
- Inverse-moment prior, also piMOM in the scalar case.

Use `delta = theta - location`, with `location = 0` by default. The prior is
nonlocal at `location`: density is zero at `delta = 0`.

## User-Facing Parameterization

Expose two parameterizations:

1. `mode`, the default and recommended user-facing parameter.
2. `tau`, the canonical scale from the derivation.

The `mode` argument is a positive distance from `location`, not a signed mode.
Internally canonicalize it as:

```r
mode <- abs(mode)
```

That means the user can think in terms of the absolute modal distance from the
null. The density remains symmetric, with modes at:

```text
location - mode
location + mode
```

The parameter object should store the canonical positive distance and the
canonical `tau`, so printing/plotting can show the intuitive `mode` while JAGS
and bridge sampling use a stable canonical scale.

Validation:

- Exactly one of `mode` or `tau` should be supplied.
- `mode` must be finite and nonzero after `abs()`.
- `tau` must be finite and positive.
- `location` must be finite.
- `order` should be a positive integer for both priors.
- Inverse-moment `df` / `nu` must be positive.

Recommended exported distribution names:

```r
prior("moment", list(mode = 0.5, order = 1, location = 0))
prior("moment", list(tau = 0.125, order = 1, location = 0))

prior("invmoment", list(mode = 0.5, order = 1, df = 1, location = 0))
prior("invmoment", list(tau = 0.25, order = 1, df = 1, location = 0))
```

Consider accepting aliases `"inverse_moment"` and `"piMOM"` later, but keep one
canonical internal distribution name. Prefer `"invmoment"` if following the
existing compact distribution naming style.

## Moment Prior

Let `r = order`, `tau > 0`, and `delta = theta - location`.

The scalar moment prior density is:

```text
f(delta | tau, r) =
  delta^(2r) * N(delta | 0, tau) / (tau^r * (2r - 1)!!)
```

Equivalently:

```text
log f =
  2r * log(abs(delta))
  + dnorm(delta, 0, sqrt(tau), log = TRUE)
  - r * log(tau)
  - log((2r - 1)!!)
```

At `delta = 0`, return density `0` and log density `-Inf`.

Modes:

```text
abs(mode) = sqrt(2 * r * tau)
tau = mode^2 / (2 * r)
```

Exact sampling and CDF/quantile representation:

```text
Y ~ chisq(df = 2r + 1)
S ~ Bernoulli(0.5), represented as sign +/- 1
delta = S * sqrt(tau * Y)
```

CDF:

```text
z = delta^2 / tau

if delta < 0:
  F(delta) = 0.5 * P(chisq_{2r+1} >= z)
if delta == 0:
  F(delta) = 0.5
if delta > 0:
  F(delta) = 0.5 + 0.5 * P(chisq_{2r+1} <= z)
```

Quantile:

```text
if p < 0.5:
  delta = -sqrt(tau * qchisq(1 - 2p, df = 2r + 1, lower.tail = FALSE))
if p == 0.5:
  delta = 0
if p > 0.5:
  delta =  sqrt(tau * qchisq(2p - 1, df = 2r + 1))
theta = location + delta
```

Use exact transformations rather than numerical integration for `cdf`, `quant`,
and `rng`.

## Inverse-Moment Prior

Let `k = order`, `nu = df`, `tau > 0`, and `delta = theta - location`.

The scalar inverse-moment prior density is:

```text
f(delta | tau, k, nu) =
  k * tau^(nu / 2) / Gamma(nu / (2k))
  * |delta|^(-(nu + 1))
  * exp(-(tau / delta^2)^k)
```

At `delta = 0`, return density `0` and log density `-Inf`. Direct evaluation
near zero must use log scale first, because the power and exponential terms are
both extreme.

Modes:

```text
abs(mode) = sqrt(tau) * (2k / (nu + 1))^(1 / (2k))
tau = mode^2 * ((nu + 1) / (2k))^(1 / k)
```

For the common piMOM case with `k = 1`:

```text
tau = mode^2 * (nu + 1) / 2
```

Exact sampling and CDF/quantile representation:

```text
a = nu / (2k)
U ~ Gamma(shape = a, rate = 1)
S ~ Bernoulli(0.5), represented as sign +/- 1
delta = S * sqrt(tau) * U^(-1 / (2k))
```

CDF:

```text
u = (tau / delta^2)^k

if delta < 0:
  F(delta) = 0.5 * P(Gamma(a, 1) <= u)
if delta == 0:
  F(delta) = 0.5
if delta > 0:
  F(delta) = 0.5 + 0.5 * P(Gamma(a, 1) >= u)
```

Quantile:

```text
if p < 0.5:
  u = qgamma(2p, shape = a, rate = 1)
  delta = -sqrt(tau) * u^(-1 / (2k))
if p == 0.5:
  delta = 0
if p > 0.5:
  u = qgamma(2 * (1 - p), shape = a, rate = 1)
  delta =  sqrt(tau) * u^(-1 / (2k))
theta = location + delta
```

Handle `p = 0` and `p = 1` as `-Inf` and `Inf`. The inverse-moment prior has
heavy tails, so tests should avoid overly strict finite-sample tail assertions.

## Native/JAGS Strategy

Use the existing compiled BayesTools JAGS module pattern restored on
`0.3.1-add-random`.

This is the preferred implementation:

- Shared C++ math core is the source of truth.
- R `pdf`, `lpdf`, `cdf`, `quant`, and `rng` call native routines or exact R
  wrappers around the same formulas.
- JAGS distributions call the same C++ core.
- Bridge sampling uses the R log-density path, so it stays aligned with JAGS.

This is better than a latent-variable JAGS-only implementation because:

- The density used by JAGS and `bridgesampling` is identical.
- Truncation and direct density evaluation are easier to reason about.
- Future priors can reuse the same native module registration pattern.
- The syntax exposed to users stays compact:

```jags
theta ~ dbt_moment(location, tau, order)
theta ~ dbt_invmoment(location, tau, order, df)
```

The package already ships the JAGS module for LKJ/Cholesky support. Reuse that
module rather than creating a separate module.

Native files to extend:

- `src/BayesTools.cc`: insert new JAGS distribution classes.
- `src/init.c`: register new `.Call` symbols for R.
- `src/Makevars.in`, `src/Makevars.win`, `src/Makevars.ucrt`: add new object
  files.
- `R/zzz.R`: extend native symbol checks if needed.

Recommended new C++ files:

- `src/nonlocal/BTNonlocalCore.h`
- `src/nonlocal/BTNonlocalCore.cc`
- `src/distributions/DBTMoment.h`
- `src/distributions/DBTMoment.cc`
- `src/distributions/DBTInvMoment.h`
- `src/distributions/DBTInvMoment.cc`
- `src/r-nonlocal.cc`

Follow the LKJ structure:

- `src/lkj/BTLKJCore.*` is shared math.
- `src/distributions/DBTLKJCPC.*` is a JAGS distribution.
- `src/r-lkj.cc` exposes R `.Call` wrappers.

## R Integration

Add a new R distribution implementation file, probably:

```text
R/distributions-nonlocal.R
```

Integrate with:

- `R/priors.R`: add `moment` and `invmoment` to `prior()`.
- Distribution methods: `rng`, `cdf`, `ccdf`, `quant`, `pdf`, `lpdf`.
- Prior-to-JAGS syntax helpers.
- Prior initialization helpers.
- Bridge-sampling prior log-density helpers.
- Plotting via existing prior plotting methods.

The generic prior plotting should mostly work once `pdf` is available, but
choose sane default plotting ranges. For `mode`-parameterized priors, default
plot windows should include both modes and enough tail mass, for example:

```text
location +/- max(4 * mode, quantile-based bounds)
```

Inverse-moment plots need care around `location` because the density is zero at
the null but can be numerically sharp near the modes.

## Bridge Sampling

Bridge sampling should evaluate these priors through the ordinary continuous
prior log-density path:

```r
.JAGS_marglik_priors.simple(samples, prior, parameter_name)
```

No discrete indicators are needed. Do not introduce latent sign/gamma
coordinates for the bridge unless the direct JAGS module path proves impossible.

Bounds:

- Moment and inverse-moment priors have full real support for `theta`, except
  the density is zero at `location`.
- For bridge bounds, use `-Inf` / `Inf`.
- Treat exact `theta == location` as log-density `-Inf`.

## Tests

Unit tests should cover:

- Parameter validation and canonicalization.
- `mode <-> tau` conversions.
- Density values at `location`, modes, and tails.
- Numeric integration close to 1.
- CDF/quantile inverses.
- RNG distribution checks using exact transformed chi-square/gamma identities.
- Symmetry around `location`.
- JAGS syntax strings and metadata.
- C++ R wrappers versus pure formula or exact R references.
- Prior plotting data is finite and nonnegative.
- Bridge log-density accepts the priors as continuous priors.

Fit-profile tests should cover:

- Loading the BayesTools JAGS module.
- Sampling simple models with `dbt_moment` and `dbt_invmoment`.
- Truncated use if supported by JAGS syntax, e.g. `T(0,)`.
- Marginal likelihood on a small direct model.

Use the existing profile policy:

```text
Rscript tools/test-profile.R unit
Rscript tools/test-profile.R fit
Rscript tools/test-profile.R fixture
```

Add `visual` only if plot snapshots are introduced or changed.

## Implementation Order

1. Restore/continue from `0.3.1-add-random`, not from current `master`.
2. Add C++ nonlocal core and R `.Call` wrappers.
3. Register new JAGS distributions in the existing BayesTools module.
4. Add R distribution methods and `prior()` integration.
5. Add JAGS syntax generation and initialization support.
6. Add bridge-sampling support through continuous prior log density.
7. Add plotting defaults/range handling if generic plotting is insufficient.
8. Add unit tests, then fit-profile JAGS tests.
9. Run unit tests before touching fixture or visual outputs.

## Open Decisions

- Final canonical user names: likely `moment` and `invmoment`.
- Whether to accept aliases such as `pmom`, `pMOM`, `pimom`, `piMOM`,
  `inverse_moment`.
- Whether `mode` should silently accept negative inputs by `abs()` or emit a
  message/warning. Current preference: silently canonicalize to the positive
  distance and print the positive value.
- Whether `order` must be integer. Current preference: require positive integer
  to match the standard pMOM/piMOM derivations and simplify JAGS validation.
- Whether to expose `df` or `nu` in the R API. Current preference: expose `df`
  and store `nu` internally only if needed.

