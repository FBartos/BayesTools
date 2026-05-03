# Mathematical Faithfulness Audit

Date: 2026-05-01

Scope: feature-by-feature audit of implemented mathematical behavior in BayesTools. This is a triage list of potential issues to decide on later. No package source code was changed as part of this audit.

Coverage: prior constructors and generic methods, prior density rendering, linear-density transforms, JAGS formula construction/evaluation, bridge-sampling parameter reconstruction, posterior extraction, model averaging, marginal distributions, summary tables, weightfunction/selection-kernel/p-hacking/bias helpers, and informed priors.

Severity guide:

- High: likely incorrect mathematical result or failed computation for ordinary supported use.
- Medium: mathematically inconsistent edge case, unsupported option accepted as if supported, or behavior that can materially bias results in some workflows.
- Low: documentation/API mismatch, validation gap, narrow boundary behavior, or missing regression coverage.

## Prior Distributions And Generic Methods

1. High - raw vector-prior generics fall through or return functions/errors.
   - Locations: `R/priors.R:1819`, `R/priors.R:2190`, `R/priors.R:2338`, `R/priors.R:2431`, `R/priors.R:2615`, `R/priors.R:2775`
   - `rng()` and joint `lpdf()`/`pdf()` handle raw `mnormal`, `mt`, and `mpoint`, but `quant()`, `mcdf()`, `mlpdf()`, `mquant()`, `mean()`, `var()`, and `sd()` do not handle non-factor vector priors consistently.
   - Probe: `quant(prior("mnormal", list(0, 1, 2)), .5)` returns a function; `var(...)` returns the `var` generic; `mean(...)` errors with `object 'm' not found`.
   - Decision: Implement

2. Medium - vector dimension `K` validation is bypassed when `allow_NA = TRUE`.
   - Locations: `R/priors-tools.R:117`, `R/priors-tools.R:123`, `R/priors.R:945`, `R/priors.R:969`, `R/priors.R:997`, `R/priors.R:1018`
   - The vector constructors pass `allow_NA = TRUE`, but the helper then skips validation for non-`NA` values too.
   - Probes: `prior("mnormal", list(0, 1, 2.5))` and `prior("mt", list(0, 1, 3, -2))` are accepted. JAGS code can then contain dimensions such as `rep(0, 2.5)` or loops like `1:2.5`.
   - Decision: Desired behavior -- we need to be able to create vector priors based on data length -- so it is fully created only later in code when data length is known

3. High - mixture weights are lost when `none`/`point` priors are coerced inside mixtures.
   - Locations: `R/priors.R:585-603`, `R/priors.R:633-652`, `R/priors.R:637-638`, `R/priors.R:1459`
   - Simple mixtures convert `prior_none()` to a point prior without preserving its `prior_weights`.
   - Factor mixtures convert scalar `point` or `none` priors to factor-point priors without preserving the original model weight.
   - Probes: a simple mixture with weights `2` and `7` stores `c(2, 1)`; a factor mixture with point weight `9` and slab weight `1` stores `c(1, 1)`.
   - Impact: intended model prior odds can silently change.
   - Decision: Fix

4. Low - point and mpoint constructors ignore explicit truncation.
   - Locations: `R/priors.R:913`, `R/priors.R:926`, `R/priors.R:1009`, `R/priors.R:1023`
   - `prior("point", list(1), truncation = list(2, 3))` silently returns a point at `1` with truncation `[1, 1]`.
   - If truncation is documented as conditioning on the supplied interval, excluding the only atom should be rejected or clearly documented as ignored for point priors.
   - Decision: What would you suggest, the truncation for points and uniforms is defined only for consistency with later functions depending downstream? Maybe plotting range, but in fact these should always relly on the actual point or interval range

5. Medium - point-mass upper-tail documentation is inclusive, while implementation is exclusive at atoms.
   - Locations: `R/distributions-point.R:12`, `R/distributions-point.R:94`, `R/distributions-mpoint.R:14`, `R/distributions-mpoint.R:125`, `R/distributions-weightfunctions.R:23`, `R/distributions-weightfunctions.R:500`, `R/distributions-weightfunctions.R:527`
   - Docs describe `lower.tail = FALSE` as `P[X >= x]`, but the implementation computes `1 - P[X <= x]`, i.e. `P[X > x]`.
   - Probe: `ppoint(1, location = 1, lower.tail = FALSE)` returns `0`; the documented inclusive upper tail would be `1`.
   - This also affects fixed/reference weightfunction components that call `ppoint()`.
   - Decision: suggest a solution, think like a senior engineer

6. Low - `JAGS_formula()` validates `parameter` with reversed arguments.
   - Location: `R/JAGS-formula.R:97`
   - The call is `check_char("parameter", parameter)` instead of `check_char(parameter, "parameter")`.
   - Invalid or multi-length parameter names can pass the wrong validation path and produce malformed downstream names.
   - Decision: Fix

7. Medium - `log(intercept)` formula mode does not enforce positive-support intercept priors.
   - Locations: `R/JAGS-formula.R:111-113`, `R/JAGS-formula.R:204-220`, `R/JAGS-formula.R:291-293`, `R/JAGS-formula.R:1030-1044`
   - The formula syntax can generate `log(<intercept>)`, but continuous-prior validation only rejects incompatible prior classes; it does not require the intercept prior to have positive support.
   - With an unconstrained normal/t/etc. intercept prior, JAGS can evaluate `log()` at nonpositive values and marginal/posterior transformations assume a domain that the prior does not guarantee.
   - Decision: do not fix, handled by JAGS

## Prior Density And Linear-Density Features

8. High - log source transforms are ignored in linear-density construction.
   - Locations: `R/priors-linear-density.R:760-762`, `R/priors-linear-density.R:775`, `R/priors-linear-density.R:812`, `R/priors-linear-density.R:884-890`, `R/priors-linear-density.R:1786-1801`
   - The range path passes unnamed scalar source transforms via `Map()`, but the distribution path passes a named scalar such as `c(a = "log")`. Downstream `identical(source_transform, "log")` checks then fail.
   - Probe: for a lognormal source with `source_transforms = c(a = "log")` and `output_transformation = "exp"`, the intended transform is identity on the original lognormal variable. Current construction returns the range of `exp(X)` instead.
   - Impact: log-intercept transformed prior densities can be on the wrong scale.
   - Decision: Fix

9. Medium - forced sampled density paths can return mismatched `x`/`y`.
   - Locations: `R/priors-density.R:134-137`, `R/priors-density.R:235-237`, `R/priors-density.R:410-412`
   - In `force_samples = TRUE`, the code keeps caller-provided `x_seq` but uses only `stats::density(..., n = n_points)$y`.
   - Probe: normal prior with `x_seq` length `5` and `n_points = 100` returns `length(x) = 5`, `length(y) = 100`.
   - For discrete forced sampling, this also converts point probabilities into KDE heights.
   - Decision: I think this is intended behavior. the function is used to provide "density estimate based on samples?" or am I missing something?

10. Medium - direct `density()` on discrete priors fails without an explicit range.
    - Locations: `R/priors-density.R:87-88`
    - The default range uses single-bracket list extraction, so `seq()` receives list objects.
    - Probe: `density(prior("bernoulli", list(.3)))` errors with `default method not implemented for type 'list'`; specifying `x_range = c(0, 1)` works.
   - Decision: Fix

11. Medium/Low - custom density-transformation Jacobian contract is ambiguous.
    - Locations: `R/priors-density.R:35-37`, `R/priors-density.R:608-616`, `man/density.prior.Rd:61-63`
    - Docs say `jac` is the Jacobian of the transformation, but the implementation multiplies by `jac` after replacing `x` with transformed coordinates.
    - Therefore `jac` must be the absolute derivative of the inverse transform on transformed support. A natural custom `exp` transform with `jac = exp` would be mathematically wrong; the implementation expects `1 / y`.
   - Decision: Fix

12. Low - transformed spike-and-slab density stores an untransformed top-level range.
    - Locations: `R/priors-density.R:464-465`, `R/priors-density.R:482-485`
    - Subdensities are transformed, but `attr(out, "x_range")` is set to the original range.
    - Probe with `transformation = "exp"` gave top-level range near `[-2.576, 2.576]` while component ranges were on `[0.076, 13.142]`.
   - Decision: Fix

13. Low - transformed density plot labels can show original prior math.
    - Locations: `R/priors-plot.R:249`, `R/priors-plot.R:436`, `R/priors-plot.R:479`, `R/priors-density.R:33`
    - Plot labels use `print(x, plot = TRUE)` without indicating that a transformed random variable is displayed.
    - The `tanh` transform is described as "Fisher's z transformation"; mathematically, Fisher's z is usually `atanh(r)`, while this transform applies `tanh`.
   - Decision: Fix

## JAGS Formula Construction, Evaluation, And Bridge Sampling

14. High - `JAGS_evaluate_formula()` does not mirror `JAGS_formula()` for `-1` factor formulas.
    - Locations: `R/JAGS-formula.R:122`, `R/JAGS-formula.R:928`, `R/JAGS-formula.R:1025`
    - Construction rewrites `-1` formulas to include a spike-zero intercept and builds the design matrix from that rewritten formula.
    - Evaluation uses the raw user formula, so a formula such as `~ f - 1` can build a two-column contrast basis during fitting but a three-dummy-column basis during evaluation.
    - Probe from the formula slice: `JAGS_formula(~ f - 1, ..., prior_factor(..., contrast = "meandif"))` builds a two-column contrast basis; `JAGS_evaluate_formula(..., ~ f - 1, ...)` errors with `non-conformable arguments`.
   - Decision: Fix

15. High - top-level factor inverse-gamma bridge reconstruction uses the wrong coordinate for `K > 1`.
    - Locations: `R/JAGS-fit.R:691`, `R/JAGS-fit.R:720`, `R/JAGS-marglik.R:450`, `R/JAGS-marglik.R:513`, `R/JAGS-marglik.R:1064`
    - JAGS samples indexed inverse-gamma auxiliaries such as `inv_mu_x[1]`, `inv_mu_x[2]`; bridge bounds correctly use those names.
    - `.JAGS_marglik_parameters.factor()` reconstructs top-level factor parameters by reading `mu_x[i]` directly, which is absent for inverse-gamma factor priors.
    - Probe returned bridge parameters `inv_mu_x[1]`, `inv_mu_x[2]`, but `JAGS_marglik_parameters.factor()` returned `NA, NA` for `mu_x`.
    - Related point-prior note: bridge excludes point factor priors, but this branch also does not reconstruct constants coherently if reached.
   - Decision: Fix

16. Medium - indexed inverse-gamma auxiliaries for factor priors leak into posterior extraction.
    - Locations: `R/posterior-extraction.R:42`, `R/JAGS-fit.R:1527`, `R/JAGS-fit.R:1541`
    - Auxiliary removal drops exact scalar names such as `inv_<par>`, but factor inverse-gamma monitoring produces indexed names such as `inv_mu_x[1]`.
    - Probe with indexed inverse-gamma columns showed they remain after `.remove_auxiliary_parameters()`, so summaries/diagnostics can expose inverse-scale auxiliaries as user parameters.
   - Decision: Fix

17. Medium - formula marginal-likelihood reconstruction can break term names containing `_data`.
    - Locations: `R/JAGS-formula.R:317`, `R/JAGS-formula.R:376`, `R/JAGS-marglik.R:1201`
    - Formula data names are created as `<parameter>_data_<term>`, but bridge reconstruction strips `_data` globally with `gsub("_data", "", ...)`.
    - A predictor named `x_data` maps `mu_data_x_data` to `mu_x`, not `mu_x_data`, so its formula contribution can be dropped or malformed.
   - Decision: Fix

18. Medium - scaled random-slope terms are not back-transformed to the original scale.
    - Locations: `R/JAGS-formula.R:229`, `R/JAGS-formula.R:553`, `R/JAGS-formula.R:1604`, `R/JAGS-formula.R:1762`
    - Continuous predictors are standardized before random-effect formula construction, and random-effect SD parameters are named through `_xREx__`.
    - The unscale parser handles ordinary fixed effects but leaves `_xREx__` SDs unchanged.
    - Derivation: if `z = (x - m) / s`, a random slope `b_g z` implies original slope `b_g / s` and a random-intercept contribution `-b_g m / s`. Leaving the random-slope SD unchanged is not original-scale faithful.
   - Decision: Fix

## Model Averaging, Marginal Inference, And Tables

19. High - zero-prior models can corrupt marginal-likelihood inclusion Bayes factors.
    - Locations: `R/model-averaging.R:46`, `R/model-averaging.R:1595`, `R/model-averaging.R:1598`
    - Zero prior weights are allowed, but `.inclusion_BF.margliks()` subtracts `max(margliks)` over all models, including zero-prior models, and then divides group weighted likelihood sums.
    - Probes: `inclusion_BF(c(0, 1), margliks = c(0, 0), is_null = c(TRUE, FALSE))` returns `NaN`; a zero-prior high-likelihood model can also force `NaN`.
    - Mathematically, zero-prior models should not affect the BF. The stability max should be taken over finite models with positive prior mass, and zero mass for an entire hypothesis should get explicit boundary handling.
   - Decision: Fix

20. Medium - spike-and-slab inference Bayes factors have boundary `NaN` cases.
    - Locations: `R/summary-tables.R:1160-1162`, `R/summary-tables.R:1197`
    - Spike-and-slab inference computes prior and posterior odds directly.
    - If `prior_prob = post_prob = 1`, the calculation is `Inf / Inf`; if both are `0`, it is `0 / 0`.
    - Mixture priors route through `inclusion_BF()`, but spike-and-slab priors do not get the same boundary handling.
   - Decision: Fix

21. Medium - marginal prior densities use averaged linear weights where posterior summaries use a row mixture.
    - Locations: `R/marginal-distributions.R:310`, `R/marginal-distributions.R:397`
    - Formula marginal posteriors split row-specific predictions, but prior densities use `colMeans(linear_weights[...])`.
    - The density of `wbar' beta` is generally not the equal mixture of row-specific densities of `w_r' beta`.
    - Impact: marginal prior heights and therefore `Savage_Dickey_BF()` can be biased for formula marginal inference.
    - Decision: Please explain in detail, this might be a critical problem needing more examination

22. Low - BF equal to `1` is interpreted directionally as evidence against.
    - Locations: `R/interpret.R:115`, `R/interpret.R:123`
    - Probe: `.interpret.BF(1, "effect", NULL)` returns `"weak evidence against the effect, BF = 1.00"`.
    - BF `1` means posterior odds equal prior odds, not evidence against either side.
    - Decision: Fix exact eq to 1 should be absence of evidence

23. Low/Medium - interval wording and ensemble quantile defaults are potentially misleading.
    - Locations: `R/interpret.R:141`, `R/interpret.R:156`, `R/summary-tables.R:76`
    - Interpretive text reports posterior quantiles as `"95% CI"`; `"credible interval"` would be mathematically clearer.
    - `ensemble_estimates_table()` defaults to `probs = c(0.025, 0.95)`, unlike the usual equal-tailed `c(0.025, 0.975)` used elsewhere.
    - Decision: fix, keep CI, fix the error in `probs = c(0.025, 0.95)` to the correct  `c(0.025, 0.975)`

24. Low - `update.BayesTools_table(remove_parameters = ...)` keeps rather than removes.
    - Locations: `R/summary-tables.R:1765`, `R/summary-tables.R:1790`, `man/update.BayesTools_table.Rd:27`
    - Docs say the argument removes parameters from the table, but the implementation subsets with `%in%`, retaining only the named rows.
    - This is mostly an API/table issue, but it can invert what is reported to users.
    - Decision: Fix

## Weightfunction, Selection, P-Hacking, And Bias Features

25. Medium - fixed-weight quantile endpoint clamps to `1` even when fixed weights exceed `1`.
    - Locations: `R/distributions-weightfunctions.R:621-625`, `R/distributions-weightfunctions.R:647-652`
    - Fixed `omega` values are validated as nonnegative relative weights with reference weight `1`; they are not bounded above by `1`.
    - Probe: `mqone.sided_fixed(0, omega = c(1, .5, 1.5), lower.tail = FALSE)` returns `1` for every component, including the point mass at `1.5`.
    - Decision: Fix

26. Low/Medium - omega diagnostic bounds still assume `[0, 1]`.
    - Locations: `R/JAGS-diagnostics.R:427-428`, `R/priors.R:1118-1120`
    - `.diagnostics_prior_bounds()` returns upper bound `1` for weightfunction `omega`, but `prior_weightfunction()` supports fixed and independent weights above `1`.
    - Probe: `range(prior_weightfunction("one-sided", .05, wf_fixed(c(1, 1.5))))` reaches `1.5`, while the diagnostics bound remains `1`.
    - Decision: Fix

27. Low - `report_scale` for p-hacking priors appears unused.
    - Locations: `R/selection-kernels.R:17`, `R/selection-kernels.R:29`, `R/selection-kernels.R:37`, `R/selection-kernels.R:61`, `R/selection-kernels.R:365-367`, `R/priors-print.R:317-330`, `R/summary-tables.R:946-948`
    - The argument is documented, validated, and stored, but backend monitors and summaries look identical for `report_scale = "alpha"` and `"pi_null"`.
    - If it is intended as a reporting-scale switch, downstream summaries and labels should use it.
    - Decision: Fix

28. Low - weightfunction documentation still uses probability or old bounded-weight language.
    - Locations: `man/weightfunctions.Rd:87`, `dev/weightfunction-prior-redesign.md:98`, `dev/weightfunction-prior-redesign.md:111`, `NEWS.md:27`
    - Some docs/dev notes describe fixed `omega` as probabilities or imply `[0, 1]` support, while current implementation supports nonnegative relative weights above `1`.
    - Decision: Fix

## Informed Priors

29. Low - `prior_informed_medicine_names` is not a type-specific validity map.
    - Locations: `R/priors-informed.R:69`, `R/priors-informed.R:1056`, `man/prior_informed_medicine_names.Rd:6`
    - The object has length `57`, while docs still say "all 46 subfields" and imply listed elements are valid names.
    - Validity depends on `type`: for example, listed names such as `"Breast Cancer"` error for `type = "smd"`, and all non-`"Cochrane"` names error for `type = "logHR"`.
    - Decision: Fix

30. Low - `prior_informed()` normalizes `name` and `type`, but not `parameter`.
    - Locations: `R/priors-informed.R:63`, `R/priors-informed.R:67`, `R/priors-informed.R:71`, `R/priors-informed.R:74`
    - `parameter = "Effect"` fails although `name` and `type` get case/spacing normalization.
    - `prior_informed("Cochrane")` also reaches a rough `switch()` error rather than a clear "parameter required" message.
    - Decision: Fix

31. Low - logHR documentation has a scale wording ambiguity.
    - Locations: `R/priors-informed.R:21`, `R/priors-informed.R:36`, `man/prior_informed.Rd:42`
    - The implementation and cited prior are on the log hazard-ratio scale, but the argument item says `logHR` is "for hazard ratios".
    - That can mislead users into applying the prior on HR rather than log(HR).
    - Decision: the scale is logHR, fix

32. Low - logHR informed-prior branch lacks regression-test coverage.
    - Locations: `R/priors-informed.R:550`, `R/priors-informed.R:1016`, `tests/testthat/test-priors-informed.R:107`
    - The test comment says the table covers "dichotomous and time to event", but the loop covers `logOR`, `logRR`, and `RD` only.
    - Manual probe matched the embedded values: `Cochrane` logHR effect is `Student-t(0, 0.13, 2)` and heterogeneity is `InvGamma(2.42, 0.30)`.
    - Decision: Fix

33. Low/watchlist - RD effect priors intentionally use unbounded Student-t support.
    - Locations: `R/priors-informed.R:433`, `tests/testthat/test-priors-informed.R:323`
    - This appears to follow the cited normal-normal approximation/table, so it is not necessarily a bug.
    - If RD is treated as a literal risk difference, support is naturally `[-1, 1]`; for example, the Cochrane RD Cauchy prior with scale `0.034` puts nonzero mass outside that interval.
    - Decision: do not fix, intended

## Checked Areas With No Obvious Mathematical Issue

- Simple distribution parameterizations looked coherent for normal, lognormal, Student-t/Cauchy, gamma, exponential, beta, uniform, Bernoulli, and inverse-gamma priors. Gamma/exponential scale-to-rate conversion and inverse-gamma shape/scale formulas looked correct.
- Continuous truncation normalization and bounded-truncation moment handling now appear coherent in the current code. Bernoulli truncation is conditioned on retained atoms correctly.
- Vector prior truncation is now rejected rather than silently mixing incompatible JAGS and bridge semantics.
- Built-in density Jacobians for `lin`, `exp`, `tanh`, and `exp_lin` looked correct for monotone transforms; the issue above is the custom contract and named source-transform path.
- Factor sampled-density construction no longer has the old undefined-`K` issue.
- Spike-and-slab `rng()`, `mean()`, and `var()` matched the intended Bernoulli-indicator mixture behavior in targeted probes. The remaining spike-and-slab concern is boundary BF handling in summary/inference tables.
- Scalar inverse-gamma handling in JAGS and bridge sampling looked internally consistent: JAGS samples the inverse coordinate, bridge bounds use that coordinate, and scalar reconstruction returns the original parameter.
- Ordinary non-random fixed-effect formula scaling/unscaling matched the expansion of products of standardized predictors, including lower-order interaction contributions.
- Bridge sampling now infers formula-scale metadata from the fit when formula lists are supplied, avoiding the old scaled-formula mismatch.
- Positive-prior model-averaging algebra looked correct: posterior odds over prior odds matches inclusion BF. Conditional model-list inference renormalizes non-null prior/posterior probabilities while keeping the unconditional BF.
- Model/posterior mixture sample allocation has been moved to explicit sample-count helpers and no longer shows the old deterministic rounding/skipping issue in the audited paths.
- Cumulative Dirichlet weightfunction marginals and RNG looked faithful: if `theta ~ Dirichlet(alpha)`, then `omega_j = sum(theta[j:J])` has the expected beta marginal and `omega_1 = 1`.
- P-hacking null-mass calibration and source/destination mass transfer looked internally coherent in focused probes.
- Informed-prior numeric values for psychology, SMD medicine, logOR, logRR, and RD branches matched the embedded tests/sources inspected by the informed-prior reviewer. `tests/testthat/test-priors-informed.R` passed in the delegated run.
- Relevant density tests reported by the density reviewer passed: `tests/testthat/test-priors-density-numeric.R` and `tests/testthat/test-priors-linear-density.R`.
