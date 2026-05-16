# version 0.3.1
### Features
- adds the `prior_random()` interface for formula random effects, including `random_block()` / `random_term()`, `random_covariance()`, `random_monitor()`, `random_new_levels()`, `random_variance_allocation()`, and `allocation_ref()` helpers for specifying random-effect standard deviation priors, covariance structures, monitoring policy, and total-variance allocation priors
- adds lme4-like formula random-effect parsing through `reformulas`, including ordinary and independent random effects, named random-effect blocks, nested grouping expressions, factor random slopes, and structured covariance shortcuts for diagonal, shared-SD independent, unstructured, compound-symmetry, heterogeneous compound-symmetry, discrete AR(1), heterogeneous AR(1), and continuous-time AR(1) random effects
- adds LKJ correlation priors for unstructured random-effect covariance matrices via `prior_lkj()` and `JAGS_lkj_corr_cholesky()`, with a package-shipped compiled JAGS backend and a pure-JAGS syntax fallback
- adds `BayesTools_load_JAGS_module()` and package compilation support for the BayesTools JAGS module used by generated LKJ-Cholesky syntax
- adds `formula_random_prior_list` to `JAGS_fit()` and `JAGS_bridgesampling()` so formula random effects can be fitted and bridge sampled through the explicit `prior_random()` interface
- adds bridge-sampling support for formula random effects by using standardized latent random effects, scalar correlation coordinates, and LKJ primitive coordinates as bridge parameters
- adds semantic random-effect summaries to `runjags_estimates_table()` / `JAGS_estimates_table()` through `random_effects_summary`, `random_effects_metadata`, `remove_random_effects`, `keep_random_effects`, `remove_random_structures`, and `keep_random_structures`
- adds random-effect parameter filter aliases such as `"random"`, `"random_sd"`, `"random_rho"`, `"random_correlation"`, `"random_variance_fraction"`, `"random_allocation"`, and `"random_sd_multiplier"` for estimates tables
- adds Dirichlet simplex priors via `prior("dirichlet", ...)` / `prior("simplex", ...)`, including random generation, log-density, marginal distribution helpers, JAGS syntax, initialization, posterior extraction, and bridge-sampling support
- adds a `RandomEffects` vignette comparing BayesTools formula random effects with lme4 and rstanarm examples

### Changes
- formula random effects now require `prior_random()`; legacy random-effect priors named with `"term|group"` inside formula `prior_list` are no longer supported
- internal random-effect term metadata now uses canonical `structure` fields only; covariance-only random-effect metadata is rejected as malformed rather than repaired downstream
- internal fitted random-effect metadata now requires canonical `homogeneous_sd`, correlation, scalar-correlation, and variance-allocation fields where applicable; malformed metadata is rejected instead of silently defaulted during prediction, summaries, reconstruction, and bridge sampling
- `JAGS_fit()` and `JAGS_extend()` now carry generated `add_parameters`, `required_packages`, and `jags_modules` metadata so formula-generated monitors and JAGS modules remain available during fitting, extension, convergence checks, and parallel execution
- `JAGS_estimates_table(transform_scaled = TRUE)` now derives formula random-effect SD and correlation summaries on the transformed original scale when fitted formula-scale metadata are available
- `JAGS_evaluate_formula()` can evaluate fitted random-effect formulas for existing grouping levels when latent random effects or group-level coefficients were monitored
- `JAGS_check_convergence()` ignores model indicator variables by default and excludes generated auxiliary monitor parameters from convergence checks unless explicitly requested

### Fixes
- improves validation and error messages for malformed random-effect formulas, unsupported covariance structures, conflicting variance-allocation specifications, missing random-effect priors, and unsupported new grouping levels in prediction
- rejects non-diagonal random-effect covariance wrappers combined with `||` syntax consistently, including mixed formulas with multiple random-effect blocks
- validates rebuilt bridge-sampling random-effect variance-allocation metadata against fitted designs and treats empty prior lists as zero log-prior contributions
- expands fixture-cache source hashing to cover random-effect and LKJ R/native sources so stale fitted fixtures are invalidated after random-feature changes
- fixes bridge-sampling bound validation by requiring named lower and upper bounds that match `add_parameters`
- fixes bridge sampling for structured formula random-effect correlations by applying canonical scalar-rho support bounds, treating covariance-boundary proposals as zero density before reconstruction, and preserving fixed scalar-rho point priors in metadata
- fixes bridge-sampling support checks for positive gamma auxiliary coordinates used by inverse-gamma priors, Dirichlet simplex priors, random-effect variance allocations, and cumulative weightfunction priors before derived parameters are reconstructed
- makes semantic random-effect SD summaries fail on missing canonical SD, point-prior, or allocation coordinates instead of silently dropping unreconstructable summaries
- enforces open LKJ primitive-coordinate support before random-effect Cholesky or correlation reconstruction and treats out-of-support LKJ bridge densities as zero density
- makes transformed-scale random-effect SD summaries and direct `transform_scale_samples()` calls require complete correlation coordinates for correlated random-effect blocks instead of silently using diagonal covariance
- makes random-effect allocation summaries fail on missing Dirichlet allocation coordinates instead of silently omitting variance-fraction summaries
- removes auxiliary Dirichlet gamma coordinates from ordinary posterior summaries while keeping them available for marginal-likelihood calculations

# version 0.3.0
### Features
- major refactoring and speed-up of unit tests
- adds support for `__default_factor` and `__default_continuous` priors in `JAGS_formula()` - when specified in the `prior_list`, these are used as default priors for factor and continuous predictors that are not explicitly specified
- adds automatic standardization of continuous predictors via `formula_scale` parameter in `JAGS_formula()` and `JAGS_fit()` - improves MCMC sampling efficiency and numerical stability
- adds `transform_scale_samples()` function to transform posterior samples back to original scale after standardization
- adds `transform_prior_samples()` function to generate and transform prior samples using the same matrix transformation as posterior samples - enables correct visualization of priors on the original (unscaled) predictor scale, including proper handling of the intercept which depends on multiple coefficient priors
- adds `transform_scaled` argument to `plot_posterior()` for visualizing prior and posterior distributions on the original (unscaled) scale when using formula-based models with auto-scaling
- adds `exp_lin` transformation type for log-intercept unscaling in density/plotting functions: `exp(a + b * log(x))`
- adds `log(intercept)` formula attribute for specifying models of the form `log(intercept) + sum(beta_i * x_i)` - useful for parameters that must be positive (e.g., standard deviation) while keeping the intercept on the original scale. Set via `attr(formula, "log(intercept)") <- TRUE`. Supported in `JAGS_formula()`, `JAGS_evaluate_formula()`, and marginal likelihood computation
- adds advanced parameter filtering options to `runjags_estimates_table()`:
  - `remove_parameters = TRUE` to remove all non-formula parameters
  - `remove_formulas` to remove all parameters from specific formulas
  - `keep_parameters` to keep only specified parameters
  - `keep_formulas` to keep only parameters from specified formulas
  - when `bias` is specified in `remove_parameters` or `keep_parameters`, the corresponding bias-related parameters (`PET`, `PEESE`, `omega`, `alpha`, `pi_null`, and `phack_kind`) are automatically included based on the bias prior type
- adds `probs` argument to `runjags_estimates_table()` and `runjags_estimates_empty_table()` for custom quantiles (default: `c(0.025, 0.5, 0.975)`)
- adds `effect_direction` argument to `plot_posterior()`, `plot_prior_list()`, `lines_prior_list()`, and `geom_prior_list()` for PET-PEESE regression plots - use `"positive"` (default) for `mu + PET*se + PEESE*se^2` or `"negative"` for `mu - PET*se - PEESE*se^2`
- redesigns `prior_weightfunction()` around a unified `side`, `steps`, and `weights` specification, with `wf_cumulative()`, `wf_fixed()`, and `wf_independent()` constructors for cumulative Dirichlet, fixed, independent, and log-independent weightfunction priors
- adds p-hacking and composed selection-bias priors via `prior_phacking()`, `prior_bias()`, calibration helpers, and `selection_backend_spec()` for compiling active step/p-hacking backend parameters
- adds error % for inclusion BF calculation

### Changes
- changes quantile column names in `runjags_estimates_table()` and `stan_estimates_table()` from `lCI`/`Median`/`uCI` to numeric values (e.g., `0.025`/`0.5`/`0.975`) for consistency with ensemble summary tables
- implied prior distributions for estimated marginal means, unstandardized coefficients, and PET-PEESE no longer require prior samples 
- implied prior distributions for weightfunction weights now use analytical forms for cumulative Dirichlet, fixed, independent, and log-independent priors, including mixture and model-averaged weightfunctions where possible
- independent weightfunction priors now allow non-reference weights above one via non-negative omega-scale priors or unrestricted log-omega priors
- replaces the legacy dot-named weightfunction prior specifications with the unified weightfunction prior API and updates JAGS generation, marginal likelihood computation, posterior extraction, diagnostics, and summary tables to use the new component-local `omega` representation
- composed selection-bias priors and publication-bias mixtures now support prior sampling and explicit unsupported-operation errors for ambiguous scalar prior generics

### Fixes
- reports inclusion Bayes factors as `NA` when the prior assigns probability 0 or 1 to inclusion, while keeping finite-sample bounds for posterior inclusion probabilities of 0 or 1
- fixes incorrect ordering the printed mixture priors
- fixes formula with no intercepts coded as `0` (instead of only `-1`)
- fixes bug in `.is.wholenumber` with NAs and `na.rm = TRUE`
- fixes ggplot prior spike layers for marginal factor plots with density and point components

# version 0.2.23
### Fixes
- `JAGS_diagnostics` functions now correctly handle factor parameters nested within mixture priors

# version 0.2.22
### Fixes
- `plot_posterior()` function with spike and slab priors 

### Changes
- unifies back-end of `prior_mixture()` and `prior_spike_and_slab()` 

# version 0.2.21
### Fixes
- `JAGS_formula()` function now replaces removed missing intercept with 0 (so the model matrix remains unchanged)
-  resetting `silent = FALSE` argument in the `JAGS_fit()` function now fits the model non-silently again 

# version 0.2.20
### Features
- extending prior functions to accept `expression()` instead of a parameter, such objects can be use to create prior distributions that depend on other parameters in JAGS
- extending the formula interface of `JAGS_fit()` function to accept expressions that are appended as literal text to the generated JAGS formula 
- extending the formula interface of `JAGS_fit()` function to handle uncorrelated random effects via `(x||y)` (lme4-like) notation

### Fixes
- `JAGS_estimates_table` not printing formula prefix when only spike and slab priors are supplied 

# version 0.2.19
### Features
- adds `max_extend` option to `autofit_control` argument in `JAGS_fit()` to limit the number of iterations for the model extension
- adds JASP progress bar integration

### Fixes
- `JAGS_diagnostics_density()` plots for mixture distributions
- prior and posterior `plot_posterior()` for simple `as_mixed_posteriors` objects
- `JAGS_evaluate_formula()` for mixture and spike and slab priors
- set Bayes factors based on alternative only prior distributions to NA
- better handling of posterior samples in `.fit_to_posterior()`

# version 0.2.18
### Features
- adding `prior_mixture()` function for creating a mixture of prior distributions
- adding `as_mixed_posteriors()` and `as_marginal_inference()` functions for a single JAGS models (with spike and slab or mixture priors) to enabling tables and figures based on the corresponding output
- adding `interpret2()` function for another way of creating textual summaries without the need of inference and samples objects
- speedup and improvements to the `runjags_estimates_table()` function

### Fixes
- small fixes for expansion of the RoBMA functionality

## version 0.2.17
### Features
- adding informed prior distributions for dichotomous and time to event outcomes based on Cochrane Database of Systematic Reviews to `prior_informed()` function
- adding bridge object convenience function `bridge_object()` (fixes: https://github.com/FBartos/BayesTools/issues/28)
- adding `Na/NaN` tests for `check_` functions (fixes: https://github.com/FBartos/BayesTools/issues/26)

### Fixes
- ability to run more than 4 chains (fixes: https://github.com/FBartos/BayesTools/issues/20)

## version 0.2.16
### Features
- update an existing JAGS fit with `JAGS_extend()` function
- new element of the `autofit_control` argument in `JAGS_fit()`: `"restarts"` allows to restart model initialization up to `restarts` times in case of failure

## version 0.2.15
### Fixes
- fixing repeated print of previous prior distribution in `model_summary_table()` in case of `prior_none()`

## version 0.2.14
### Features
- adding `contrast = "meandif"` to the `prior_factor` function which generates identical prior distributions for difference between the grand mean and each factor level
- adding `contrast = "independent"` to the `prior_factor` function which generates independent identical prior distributions for each factor level
- `remove_column` function for removing columns from `BayesTools_table` objects without breaking the attributes etc...
- adding empty table functions (https://github.com/FBartos/BayesTools/issues/10)
- adding `remove_parameters` argument to `model_summary_table()`
- adding multivariate point distribution functions
- adding `point` prior distribution as option to `prior_factor` with `"meandif"` and `"orthonormal"` contrasts
- adding `marginal_posterior()` function which creates marginal prior and posterior distributions (according to a model formula specification)
- adding `Savage_Dickey_BF()` function to compute density ratio Bayes factors based on `marginal_posterior` objects
- adding `marginal_inference()` function to combine information from `marginal_posterior()` and `Savage_Dickey_BF()`
- adding `marginal_estimates_table()` function to summarize `marginal_inference()` objects
- adding `plot_marginal()` function to visualize `marginal_inference()` objects

### Changes
- `contrast = "meandif"` is now the default setting for `prior_factor` function
- depreciating `transform_orthonormal` argument in favor of more general `transform_factors` argument 
- switching `dummy` contrast/factor attributes to `treatment` for consistency (https://github.com/FBartos/BayesTools/issues/23)

### Fixes
- zero length inputs to `check_bool()`, `check_char()`, `check_real()`, `check_int()`, and `check_list()` do not throw error if `allow_NULL = TRUE`
- properly aggregating identical priors in the plotting function (previously overlying multiple spikes on top of each other when attributes did not match)
- `student-t` allowed as a prior distribution `name`
- fixing factor contrast settings in `JAGS_evaluate_formula`
- fixing spike prior transformations

## version 0.2.13
### Features
- `runjags_estimates_table()` function can now handle factor transformations 
- `plot_posterior` function can now handle factor transformations 
- ability to remove parameters from the `runjags_estimates_table()` function via the `remove_parameters` argument

### Fixes
- inability to deal with constant intercept in marglik formula calculation
- `runjags_estimates_table()` function can now remove factor spike prior distributions
- marginal likelihood calculation for factor prior distributions with spike 
- mixing samples from vector priors of length 1
- same prior distributions not always combined together properly when part of them was generated via the formula interface

## version 0.2.12
### Features
- `stan_estimates_summary()` function
- reducing dependency on runjags/rjags

### Fixes
- dealing with posterior samples from rstan
- dealing with vector posterior samples
- fixing MCMC error of SD calculation for transformed samples (previously reported 100 times lower)

## version 0.2.11
### Features
- adding Bernoulli prior distribution
- adding spike and slab type of prior distributions (without marginal likelihood computations/model-averaging capabilities)
- new vignette comparing Bayes factor computation via marginal likelihood and spike and slab priors

### Fixes
- when a transformation is applied, JAGS summary tables now produce the mean of the transformed variable (previous versions incorrectly returned transformation of the mean) 

### Changes
- runjags_XXX_table functions are now also exported as JAGS_XXX_functions for consistency with the rest of the code

## version 0.2.10
### Features
- trace, density, and autocorrelation diagnostic plots for JAGS models

## version 0.2.9
### Fixes
- dealing with NaNs in inclusion Bayes factors due to overflow with very large marginal likelihoods

## version 0.2.8
### Fixes
- dealing with point prior distributions in `JAGS_marglik_parameters_formula` function
- posterior samples dropping name in `runjags_estimates_table` function
- `ensemble_summary_table` and `ensemble_diagnostics_table` function can create table without model components

## version 0.2.7
### Features
- `JAGS_evaluate_formula` for evaluating formulas based on data and posterior samples (for creating predictions etc)  
- `JAGS_parameter_names` for transforming formula names into the JAGS syntax

## version 0.2.6
### Features
- `plot_models` implementation for factor predictors
- `format_parameter_names` for cleaning parameter names from JAGS
- `mean`, `sd`, and `var` functions now return the corresponding values for differences from the mean for the orthonormal prior distributions

### Fixes
- proper splitting of transformed posterior samples based on orthonormal contrasts in `runjags_summary_table` function (previous version crashed under other than default `fit_JAGS` settings)
- always showing name of the comparison group for treatment contrasts  in `runjags_summary_table` function
- better handling of transformed parameter names in `plot_models` function

## version 0.2.5
### Features
- `add_column` function for extending `BayesTools_table` objects without breaking the attributes etc...
- ability to suppress the formula parameter prefix in `BayesTools_table` functions with with `formula_prefix` argument

### Fixes
- allowing to pass point prior distributions for factor type predictors

## version 0.2.4
### Features
- adding possibility to multiply a (formula) prior parameter by another term (via `multiply_by` attribute passed with the prior)
- t-test example vignette

## version 0.2.3
### Fixes
- fixing error from trying to rename formula parameters in BayesTools tables when multiple parameters were nested within a component

## version 0.2.2
### Fixes
- fixing layering of prior and posterior plots in `plot_posterior` (posterior is now plotted over the prior)

## version 0.2.1
### Fixes
- fixing JAGS code for multivariate-t prior distribution

## version 0.2.0
### Changes
- ensemble inference, summary, and plot functions now extract the prior list from attribute of the fit objects (previously, the prior_list needed to be passed for each model within the model_list as the priors argument

### Features
- adding formula interface for fitting and computing marginal likelihood of JAGS models
- adding factor prior distributions (with treatment and orthonormal contrasts)

## version 0.1.4
### Fixes
- fixing DOIs in the references file
- adds marglik argument `inclusion_BF` to deal with over/underflow (Issue #9)
- better passing of BF names through the `ensemble_inference_table()` (Issue #11)

### Features
- adding logBF and BF01 options to `ensemble_summary_table` (Issue #7)

## version 0.1.3
### Features
- `prior_informed` function for creating informed prior distributions based on the past psychological and medical research

## version 0.1.2
### Fixes
- `prior.plot` can't plot "spike" with `plot_type == "ggplot"` (Issue #6)
- `MCMC error/SD` print names in BayesTools tables (Issue #8)
- `JAGS_bridgesampling_posterior` unable to add a parameter via `add_parameters`

### Features
- `interpret` function for creating textual summaries based on inference and samples objects

## version 0.1.1
### Fixes
- `plot_posterior` fails with only mu & PET samples (Issue #5)
- ordering by "probabilities" does not work in 'plot_models' (Issue #3)
- BF goes to NaN when only a single model is present in 'models_inference' (Issue #2)
- summary tables unit tests unable to deal with numerical precision
- problems with aggregating samples across multiple spikes in `plot_posterior'

### Features
- allow density.prior with range lower == upper  (Issue #4)
- moving rstan towards suggested packages

## version 0.1.0
- published on CRAN

## version 0.0.0.9010
- plotting functions for models

## version 0.0.0.9009
- plotting functions for posterior samples

## version 0.0.0.9008
- plotting functions for mixture of priors

## version 0.0.0.9007
- improvements to prior plotting functions

## version 0.0.0.9006
- ensemble and model summary tables functions

## version 0.0.0.9005
- posterior mixing functions

## version 0.0.0.9004
- model-averaging functions

## version 0.0.0.9003
- JAGS fitting related functions

## version 0.0.0.9002
- JAGS bridgesampling related functions

## version 0.0.0.9001
- JAGS model building related functions

## version 0.0.0.9000
- priors and related methods
