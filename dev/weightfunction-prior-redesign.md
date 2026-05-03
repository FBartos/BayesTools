# Weightfunction Prior Redesign

This note records the proposed redesign of `prior_weightfunction()` and the
downstream weightfunction machinery for a breaking major release. The goal is to
replace the current compound distribution names, such as `"one.sided.fixed"`,
with a unified API that separates weightfunction geometry from the prior placed
on the publication weights.

## Motivation

The current weightfunction API combines several independent concepts in one
`distribution` string:

- side geometry: one-sided vs two-sided;
- weight prior family: cumulative Dirichlet, fixed weights, future independent
  weights;
- latent scale: direct omega scale vs log-omega scale.

This makes extension awkward. Adding names such as
`"one.sided.independent.log"` would work locally, but it would keep expanding a
string-dispatch API that downstream code already handles through brittle checks
on `distribution` and `names(parameters)`.

The major release should instead introduce an explicit internal representation
and make all downstream code dispatch on structured fields.

## Proposed Public API

Use `prior_weightfunction()` for the p-value geometry and a separate `weights`
object for the prior on bin weights:

```r
prior_weightfunction(
  side = "one-sided",
  steps = c(.025, .05),
  weights = wf_cumulative(alpha = c(1, 1, 1)),
  prior_weights = 1
)
```

Examples:

```r
prior_weightfunction(
  side = "two-sided",
  steps = c(.05),
  weights = wf_cumulative(alpha = c(1, 1))
)

prior_weightfunction(
  side = "one-sided",
  steps = c(.05, .10),
  weights = wf_independent(prior("beta", list(1, 1)), scale = "omega")
)

prior_weightfunction(
  side = "one-sided",
  steps = c(.05, .10),
  weights = wf_independent(
    prior("normal", list(0, 1), truncation = list(upper = 0)),
    scale = "log_omega"
  )
)

prior_weightfunction(
  side = "one-sided",
  steps = c(.05),
  weights = wf_fixed(omega = c(1, .5))
)
```

The old compound distribution strings can be removed or supported only through a
compatibility shim that constructs the new objects.

## Weight Prior Constructors

### `wf_cumulative(alpha)`

Current cumulative Dirichlet model. The unconstrained latent variables are
gamma draws, normalized to a Dirichlet vector, then accumulated into monotone
weights. One bin is fixed to `omega = 1`.

This replaces current `"one.sided"` and `"two.sided"` specifications with
`alpha`.

### `wf_fixed(omega)`

Fixed publication weights. This is equivalent to point-mass priors on each bin
but should remain a first-class convenience constructor.

Validation should require the vector length to match the number of bins and the
declared reference bin to be exactly `1`.

### `wf_independent(prior, scale = "omega")`

Independent and exchangeable priors on all non-reference bin weights.

For `scale = "omega"`, the prior lives directly on the non-negative relative
weight scale, for example:

```r
wf_independent(prior("beta", list(1, 1)), scale = "omega")
```

For `scale = "log_omega"`, the latent parameter is:

```r
log_omega[j] ~ prior
omega[j] = exp(log_omega[j])
```

The prior must be constrained only as needed by the modeling choice. The
implementation maps `omega = exp(log_omega)`, so unrestricted log-scale priors
can represent relative weights above `1`.

### Future: `wf_branching(...)`

The current nonmonotone one-sided `alpha1` / `alpha2` interface should probably
be replaced by an explicit branching prior family. For example:

```r
wf_branching(
  expected = wf_cumulative(alpha = c(...)),
  unexpected = wf_cumulative(alpha = c(...))
)
```

This should be designed after the main component-local architecture is in place.

## Reference Bin

Make the reference bin explicit:

```r
reference = "most_significant"
```

Default behavior:

- the most significant bin is fixed to `omega = 1`;
- all stochastic priors apply only to non-reference bins;
- `wf_fixed()` may accept a full `omega` vector, but validation should require
  the reference bin to equal `1`.

For one-sided `steps = c(.025, .05)`, bins are:

```text
[0, .025]     omega = 1
[.025, .05]   omega ~ prior
[.05, 1]      omega ~ prior
```

This makes independent priors identifiable and avoids hidden rescaling.

## Canonical Internal Representation

Every weightfunction prior should resolve to a structured object similar to:

```r
list(
  side = "one-sided",
  steps = c(.025, .05),
  bins = data.frame(
    lower = c(0, .025, .05),
    upper = c(.025, .05, 1),
    reference = c(TRUE, FALSE, FALSE)
  ),
  weights = list(
    type = "independent",
    scale = "log_omega",
    prior = prior("normal", list(0, 1), truncation = list(upper = 0))
  ),
  prior_weights = 1
)
```

Downstream code should dispatch on structured fields:

```r
prior$side
prior$weights$type
prior$weights$scale
```

It should not infer behavior from strings such as `"one.sided.fixed"` or from
parameter names such as `alpha`, `alpha1`, or `omega`.

## Two-Sided Priors and Global Mapping

For two-sided priors, `steps` remain two-sided p-value thresholds. For example:

```r
prior_weightfunction(
  side = "two-sided",
  steps = c(.05, .10),
  weights = wf_independent(prior("beta", list(1, 1)))
)
```

This defines bins:

```text
[0, .05]     omega = 1
[.05, .10]   omega ~ Beta(1, 1)
[.10, 1]     omega ~ Beta(1, 1)
```

When two-sided and one-sided priors are mixed, two-sided priors must be expanded
onto the one-sided global cut grid, as the current `weightfunctions_mapping()`
does. The redesign should retain this mapping behavior but implement it using
`prior$side` and bin metadata instead of string matching on `prior$distribution`.

Example:

```text
global cuts: 0, .025, .05, 1

component A steps: .05
component A local bins: [0, .05], [.05, 1]
component A expanded:   [0, .05], [0, .05], [.05, 1]

component B steps: .025, .05
component B expanded: [0, .025], [.025, .05], [.05, 1]
```

## JAGS Architecture

The current bias-mixture code uses one shared latent representation:

```jags
bias_indicator ~ dcat(prior_weights[])
eta[i] ~ dgamma(eta_shape[i, bias_indicator], 1)
omega = eta2omega(eta, ...)
```

This assumes that all stochastic weightfunction components can be expressed as a
cumulative-Dirichlet `eta -> omega` transform. That does not generalize well to
fixed, independent omega-scale, independent log-scale, or branching priors.

The preferred design is component-local omega construction:

```jags
bias_indicator ~ dcat(prior_weights[])

# component 1: no selection bias
for (j in 1:J_global) {
  omega_component_1[j] <- 1
}

# component 2: cumulative Dirichlet
eta_component_2[1] ~ dgamma(...)
eta_component_2[2] ~ dgamma(...)
omega_component_2[...] <- cumulative_dirichlet_transform(...)

# component 3: fixed
omega_component_3[1] <- 1
omega_component_3[2] <- .5

# component 4: independent beta
omega_component_4[1] <- 1
omega_component_4[2] ~ dbeta(a, b)
omega_component_4[3] ~ dbeta(a, b)

# active global omega used by the likelihood
for (j in 1:J_global) {
  omega[j] <-
    omega_component_1[j] * equals(bias_indicator, 1) +
    omega_component_2[j] * equals(bias_indicator, 2) +
    omega_component_3[j] * equals(bias_indicator, 3) +
    omega_component_4[j] * equals(bias_indicator, 4)
}
```

Single-model weightfunctions can use the same component builder without the
indicator switch.

## Suggested Internal Builder

Introduce a helper with a narrow contract:

```r
.build_weightfunction_component_jags(
  prior,
  component_id,
  global_cuts,
  parameter_prefix = "omega_component"
)
```

It should return:

- JAGS syntax for latent priors;
- JAGS syntax for deterministic component-local omega;
- names of sampled latent parameters for monitoring and bridge sampling;
- lower and upper bounds for bridge sampling;
- a log-prior evaluator for marginal likelihood;
- local-to-global bin mapping.

This builder should be used by both single-prior JAGS generation and bias
mixture JAGS generation.

## Marginal Likelihood

Bridge sampling should use the natural latent parameters for each weight prior
family:

- `wf_cumulative()`: bridge on gamma variables, e.g. `eta_component_*`;
- `wf_independent(scale = "omega")`: bridge on sampled omega variables;
- `wf_independent(scale = "log_omega")`: bridge on sampled log-omega variables;
- `wf_fixed()`: no sampled parameters and log prior contribution `0`.

For log-scale independent priors, no Jacobian is needed if bridge sampling uses
`log_omega` as the sampled parameter and the model defines
`omega = exp(log_omega)` deterministically. A Jacobian would only be needed if
the prior density were evaluated on the omega scale.

The marginal likelihood code should dispatch through the same structured
weight-prior object instead of checking `names(prior$parameters)`.

## Prior Plotting and Analytical Densities

The recent analytical prior-plot machinery can extend naturally:

- cumulative Dirichlet marginals are beta distributions or point masses;
- fixed weights are point masses;
- independent omega-scale priors use the nested prior density, CDF, and
  quantile functions directly;
- independent log-scale priors use the transformed density and CDF:

```text
F_omega(x) = F_log(log(x))
f_omega(x) = f_log(log(x)) / x
```

For prior mixtures and model averages, each global omega interval should be
represented as a mixture of marginal components, using the global cut mapping.

## Public Posterior Output

For the major release, make `omega` the stable public posterior output.
Parameters such as `eta`, `eta1`, `eta2`, `log_omega`, and any component-local
latent variables should be treated as private support parameters.

Summaries, plots, diagnostics, and RoBMA extraction should consume aligned
`omega[...]` columns. Marginal likelihood code can still know about the latent
parameters through the component builder.

## Implementation Phases

1. Introduce `wf_cumulative()`, `wf_fixed()`, and `wf_independent()` constructors.
2. Refactor `prior_weightfunction()` to produce the canonical structured object.
3. Rewrite bin construction and `weightfunctions_mapping()` around explicit bin
   metadata.
4. Add component-local JAGS builders for single weightfunction priors.
5. Refactor bias-mixture JAGS generation to switch between component-local
   omega vectors.
6. Refactor bridge-sampling parameter extraction, bounds, log-prior evaluation,
   and transformed parameter reconstruction.
7. Extend analytical prior plotting to independent omega and log-omega priors.
8. Update summaries, diagnostics, and posterior extraction to treat only aligned
   `omega[...]` as public output.
9. Add compatibility shims or clear deprecation errors for old compound
   distribution names.

## Open Questions

- Should old `prior_weightfunction("one.sided", list(...))` calls be supported
  temporarily, or should the major release fail fast with a migration message?
- Should `wf_independent(scale = "log_omega")` accept only priors already
  truncated to `upper = 0`, or should the constructor automatically apply that
  truncation?
- Should fixed priors require the reference bin to be the only bin equal to `1`,
  or only require at least the reference bin to be `1`?
- Should the nonmonotone `alpha1` / `alpha2` prior be represented as
  `wf_branching()` immediately, or temporarily translated into an internal
  compatibility family?
