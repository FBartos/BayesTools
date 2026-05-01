# Mathematical Faithfulness Audit

Date: 2026-05-01

Scope: feature-by-feature audit of implemented mathematical behavior in BayesTools. This is a triage list of potential issues to decide on later; no package source code was changed as part of this audit.

Severity guide:

- High: likely incorrect mathematical result for ordinary supported use.
- Medium: mathematically inconsistent edge case, unsupported option accepted as if supported, or behavior that can materially bias results in some workflows.
- Low: documentation/API mismatch, validation gap, or narrow edge behavior.

## Prior Distributions And Generic Prior Methods

1. High - inverse-gamma variance formula is wrong.
   - Location: `R/priors.R:2742`
   - Current formula uses `scale^2 / (shape - 1)^2 * (shape - 2)`.
   - Standard inverse-gamma variance is `scale^2 / ((shape - 1)^2 * (shape - 2))` for `shape > 2`.
   - Probe: `shape = 4, scale = 2` gives package value `0.8888889`; expected value is `0.2222222`.

2. High - bounded truncation is ignored when deciding whether moments exist.
   - Locations: `R/priors.R:2580-2589`, `R/priors.R:2752-2761`
   - Bounded Cauchy/t or inverse-gamma priors have finite truncated moments even when their untruncated moments do not exist.
   - Current code returns `NaN` before using numerical integration in these cases.
   - Probe: bounded Cauchy on `[-1, 1]` and bounded inverse-gamma on `[0.1, 10]` both returned `NaN` mean/variance.

3. High - truncated Bernoulli probabilities and moments are not conditioned correctly.
   - Locations: `R/priors.R:2056`, `R/priors.R:2078`, `R/priors.R:2592-2597`, `R/priors.R:2764-2771`
   - A Bernoulli truncated to `[0.5, 1]` should condition on `X = 1`: density/mass at 1 should be 1, CDF below 1 should be 0, mean should be 1.
   - Probe with `p = 0.25`, truncation `[0.5, 1]`: package returned `pdf(1)=0.25`, `cdf(0.75)=0.75`, `mean=0`.

4. Medium - uniform priors accept but ignore `truncation`.
   - Locations: `R/priors.R:877-895`, `R/priors.R:2111`
   - General prior docs say `truncation` defines bounds, but uniform treats all uniforms as already truncated and leaves support on `[a, b]`.
   - Probe: `prior("uniform", list(0, 10), truncation = list(2, 5))` still gives `pdf(1)=0.1`, `cdf(5)=0.5`, `mean=5`.

5. Medium - density transformation Jacobians are wrong for some transformations.
   - Locations: `R/priors-density.R:608-616`, `R/priors-density.R:623-635`
   - `lin` uses `1 / b` instead of `1 / abs(b)`, producing negative densities for decreasing transforms.
   - `exp_lin` uses `1 / (b * x)`, which is not the inverse-Jacobian for `exp(a + b * log(x))`; for `a = 0, b = 1` the transform is identity but the Jacobian becomes `1 / x`.
   - Probe: normal prior with `lin(a=0,b=-2)` produced density range `[-0.199470477, -0.007229872]`.

6. Low/Medium - point and multipoint upper-tail documentation does not match implementation at atoms.
   - Locations: `R/distributions-point.R:12-13`, `R/distributions-point.R:94-98`, `R/distributions-mpoint.R:14-15`, `R/distributions-mpoint.R:125-128`
   - Docs describe `lower.tail = FALSE` as `P[X >= x]`, while implementation returns `1 - P[X <= x]`, i.e. `P[X > x]`.
   - Probe: `ppoint(0, location = 0, lower.tail = FALSE)` returns `0`; documented inclusive upper tail would be `1`.

7. Medium - sampled factor-prior density can reference an undefined symbol.
   - Location: `R/priors-density.R:403-405`
   - In the `force_samples = TRUE` path for orthonormal/meandif priors, the code fills `K` using `.get_prior_factor_levels(prior)` but the local object is named `x`.

## Spike-And-Slab And Prior Mixtures

8. High - spike-and-slab variance is wrong for point inclusion probabilities and likely for non-beta random inclusion priors.
   - Locations: `R/priors.R:1427`, `R/priors.R:2712-2726`
   - RNG uses a Bernoulli indicator, so `E[I^2] = E[I]`. The variance calculation uses zero variance for point inclusion and effectively uses `E[p^2]` for general inclusion priors.
   - Probe: normal slab with fixed inclusion probability `0.5` should have variance about `0.5`; `var()` returned `0.25`, while simulation returned about `0.5015`.

9. Medium - `prior_spike_and_slab(..., prior_weights = ...)` appears ignored.
   - Locations: `R/priors.R:439-441`, `R/priors.R:471-478`
   - The argument is documented as prior odds for the spike-and-slab prior, but the returned mixture stores component weights `c(1, 1)` instead.
   - Probe: `prior_spike_and_slab(..., prior_weights = 9)` returned `attr(prior, "prior_weights") == c(1, 1)`.

10. Medium - converting `prior_none()` inside mixtures drops its original `prior_weights`.
    - Locations: `R/priors.R:585-603`, `R/priors.R:633-650`
    - `prior_none(prior_weights = 7)` is converted to a point prior without preserving weight.
    - Probe: `prior_mixture(list(prior(..., prior_weights=2), prior_none(7)))` stored weights `c(2, 1)`, not `c(2, 7)`.

## JAGS Fitting And Bridge Sampling

11. High - formula marginal-likelihood parameter reconstruction does not handle inverse-gamma formula priors.
    - Locations: `R/JAGS-fit.R:687`, `R/JAGS-fit.R:716`, `R/JAGS-marglik.R:408`, `R/JAGS-marglik.R:949`, `R/JAGS-marglik.R:1160`, `R/JAGS-marglik.R:1213`
    - JAGS samples inverse-gamma priors through `inv_<param>` and reconstructs ordinary parameters for non-formula priors.
    - Formula reconstruction reads `samples[[term]]` or `samples[[parameter_intercept]]` directly, so inverse-gamma formula terms are missing or on the wrong coordinate.

12. Medium - bridge sampling with scaled formulas can use the wrong predictor scale unless scaling is passed again.
    - Locations: `R/JAGS-fit.R:292`, `R/JAGS-marglik.R:93`
    - `JAGS_fit()` stores formula scale metadata on the fit, but `JAGS_bridgesampling()` rebuilds formulas from its own `formula_scale_list` argument, defaulting to `NULL`.
    - If the fit used standardized predictors and the bridge call omits the same scale list, the bridge log posterior can evaluate a different model.

13. Medium - vector-prior truncation is accepted but not consistently implemented.
    - Locations: `R/priors.R:934-989`, `R/JAGS-fit.R:734`, `R/JAGS-marglik.R:441`, `R/priors.R:1733`
    - Vector priors store truncation and bridge bounds use it, but JAGS vector syntax does not emit truncation and vector `lpdf()` ignores it.
    - If vector truncation is intended, JAGS sampling, bridge bounds, and target density disagree. If not intended, it should be rejected.

14. Low - bridge `add_bounds` documentation says `"up"` while validation requires `"ub"`.
    - Locations: `man/JAGS_bridgesampling.Rd:56`, `R/JAGS-marglik.R:280`

## Model Averaging, Inference, And Summary Tables

15. High - transformed Stan summaries are mathematically wrong for nonlinear transformations.
    - Locations: `R/summary-tables.R:1413-1426`
    - Samples are transformed and mean/SD are recomputed, but median is transformed from the old median, MC error is transformed like a parameter value, and lower/upper intervals are not recomputed from transformed samples before column renaming.
    - For nonlinear transforms such as `exp`, intervals and MC error are not faithful summaries of the transformed posterior.

16. Medium - deterministic posterior/prior mixture allocation can distort model weights and sample counts.
    - Locations: `R/model-averaging.R:352`, `R/model-averaging.R:425`, `R/model-averaging.R:652`, `R/marginal-distributions.R:704-719`, `R/marginal-distributions.R:759`
    - Posterior mixing uses `round(post_probs * n_samples)` and skips components with rounded count `<= 1`; total draws need not equal `n_samples`.
    - Prior mixing uses `ceiling()` plus truncation; small components can be over- or under-represented.
    - A multinomial allocation or direct model-index sampling would preserve mixture probabilities more faithfully.

17. Medium - all-failed marginal likelihood comparisons can return `NaN` probabilities/BFs.
    - Locations: `R/model-averaging.R:62-64`, `R/model-averaging.R:98-106`, `R/model-averaging.R:133-136`, `R/model-averaging.R:1499-1504`
    - `models_inference()` maps `NA` marginal likelihoods to `-Inf`, but all `-Inf` values make `max(margliks)` also `-Inf`, producing `-Inf - -Inf = NaN`.
    - Probe: `compute_inference(c(1,1), c(-Inf,-Inf), is_null=c(TRUE,FALSE))` returned `post_probs = NaN, NaN` and `BF = NaN`.

18. Low - zero total prior weights are not explicitly rejected before normalization.
    - Locations: `R/model-averaging.R:46`, `R/model-averaging.R:62`, `R/model-averaging.R:94`, `R/model-averaging.R:136`, `R/model-averaging.R:192`
    - Weights are checked as nonnegative, then normalized by `sum(prior_weights)`.
    - Downstream code may error, but the mathematical condition `sum(prior_weights) > 0` should be validated at the boundary.

19. Low - direct probability BF helper has questionable boundary behavior for incoherent manual inputs.
    - Locations: `R/model-averaging.R:1483-1485`
    - If user-supplied prior/posterior probabilities put zero mass on a side, the odds ratio can be undefined; the helper returns hard-coded `Inf` or `0` in some such cases.

## Weightfunction, Selection, P-Hacking, And Bias Features

20. Medium - fixed-weight quantile clamps infinite point quantiles to `1` even when fixed weights may exceed `1`.
    - Locations: `R/distributions-weightfunctions.R:18`, `R/distributions-weightfunctions.R:621-625`
    - Current code permits relative weights above 1, but the fixed-weight quantile path maps `Inf` to `1`.
    - Probe: `mqone.sided_fixed(0, omega = c(1, .5, 1.5), lower.tail = FALSE)` returned `1` for all components, including the point mass at `1.5`.

21. Low - documentation still describes some relative weights as probabilities or constrained to `[0, 1]`.
    - Locations: `R/priors.R:237`, `man/weightfunctions.Rd:87`, `dev/weightfunction-prior-redesign.md:98`, `dev/weightfunction-prior-redesign.md:111`
    - Source allows log-omega priors and weights above 1, but docs/dev notes retain probability language.

## Informed Priors And API Faithfulness

22. Low - medical `prior_informed()` without `parameter` reaches a raw `switch()` error.
    - Location: `R/priors-informed.R:74`
    - This is not a formula error, but the API does not clearly reject a mathematically required argument.

## Checked Areas With No Obvious Mathematical Issue

- Simple JAGS priors use JAGS precision/rate conventions correctly for normal, lognormal, t, gamma, exponential, beta, uniform, and inverse-gamma priors.
- Non-formula inverse-gamma bridge sampling uses the sampled inverse coordinate consistently and reconstructs the original parameter before user likelihood evaluation.
- Continuous truncation normalization via `.prior_C1()`, `.prior_C2()`, and `.prior_C()` looked coherent for continuous distributions when moments are allowed to reach the integration path.
- Gamma `scale` to `rate` and exponential `scale` to `rate` conversions looked correct.
- Cumulative Dirichlet weightfunction beta marginals and RNG looked internally consistent.
- P-hacking null-mass calibration and destination/source mass transfer looked internally consistent with the package tests.
- Inclusion Bayes factors use posterior odds over prior odds in the expected direction for ordinary package-generated inputs.
- `as_mixed_posteriors()` vector K extraction uses `$parameter`, which works through R partial matching against `parameters` by default. This is fragile style, but I did not count it as a mathematical bug.
