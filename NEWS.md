# version 0.4.0
### Breaking changes
- requires R >= 4.3.0 and the `reformulas` package; BayesTools now compiles a JAGS module (inverse-gamma, moment, inverse-moment, and LKJ distributions) and requires JAGS 4.x (>= 4.3.0), including its headers and libraries for installation from source
- fits, mixed and marginal posteriors, marginal likelihoods, and weight-function priors created with BayesTools 0.3.0 or earlier (including RoBMA 4.0.0 models) must be refitted or recreated; post-fit functions stop on them with a `BayesTools_refit_required` error
- marginal likelihoods are `BayesTools_marglik` objects (still inheriting from `"bridge"`), and the model-averaging functions accept only these; use `JAGS_bridgesampling()` or `bridgesampling_object(logml)`
- inverse-gamma priors use the module distribution `dbt_invgamma` on the parameter itself instead of a gamma prior on `inv_<parameter>`
- removes the deprecated `transform_orthonormal` argument (use `transform_factors`) and the ineffective `seed` argument of `JAGS_extend()`; pass the arguments that followed them by name
- `check_bool()` rejects `NA` by default and `check_int()` rejects infinite values
- fixed weight-function priors require a reference weight of exactly 1, and cumulative weight-function and Dirichlet priors require concentrations of at least 0.01
- stricter inputs: finite autofit convergence thresholds, named `add_bounds`, factor priors outside formulas need `prior_factor_levels()`, formula random effects need `prior_random()`, and mean-difference or orthonormal priors on indicator-coded interaction terms (e.g., `g:x` in `~ g + g:x`) and `multiply_by` on the intercept prior are rejected
- `JAGS_evaluate_formula()`, `transform_scale_samples()`, and original-scale transformations require the formula design stored by `JAGS_fit()`; attach it to draws without a fit with `JAGS_formula_draws()`
- chain seeds are drawn from `seed` with R's random-number generator instead of `seed + chain`, and `mix_posteriors()` allocates draws with multinomial counts, so seeded results differ from 0.3.0; seeded functions restore the caller's random-number state instead of resetting it
- conditioning on an unknown parameter name is an error instead of a warning
- `Savage_Dickey_BF()` requires an exact, finite, and nonzero prior ordinate at the null hypothesis (classed errors otherwise), returns a finite Bayes factor when the null lies outside the posterior draws, and returns `NA` with the reason for levels with a posterior point mass at the null
- `JAGS_check_convergence()` and autofit classify parameters by their declared roles: constant parameters are structural, sampled parameters with constant draws fail unless `allow_not_assessable = TRUE`, and indicators are identified from their priors rather than their names
- parameter labels are rendered from structured label parts, which changes some labels: contrast coefficients that are not levels are labelled `term{j}`, interaction cells and formula prefixes are named consistently across tables, Bayes factor MC-error warnings name parameters by their table labels, and individual weight-function omegas in plots follow the order of summary tables
- `JAGS_fit()` monitors point priors
- posterior draw metadata are stored in one validated attribute (`posterior_metadata()`); arithmetic drops the metadata, values modified in place stop with `BayesTools_stale_metadata`, and posterior plots take point masses only from declared posterior atoms

### Features
- adds formula random effects via `prior_random()`, `random_block()`, and related helpers, with lme4-like syntax (independent, nested, and interacted grouping, factor random slopes, where a mean-difference or orthonormal factor term is one exchangeable unit with one SD) and diagonal, unstructured (LKJ prior), compound-symmetry, AR(1), heterogeneous, and continuous-time AR(1) covariance structures
- adds random-effect variance allocations with Dirichlet weights and optional inclusion gates (`random_variance_allocation()`), known group covariance kernels (`random_group_covariance()`), external SD sources (`random_sd_source()`), spike-and-slab and mixture priors on random-effect SDs and correlations, and centered or noncentered parameterizations
- adds predictions with random effects (`JAGS_predict_formula()`; `formula_target`, `blocks`, and `new_levels` in `JAGS_evaluate_formula()`), marginal covariance helpers, and exact marginalization of Gaussian random-effect blocks in `JAGS_bridgesampling()`
- adds random-effect summaries to estimates tables (`random_effects_summary` and related filters) on the fitted or original scale, with level-cell correlations and multiple correlations (`cor_total`) of factor terms, and `random_effects_summary_posterior()`
- adds the `RandomEffects` vignette comparing BayesTools formula random effects with lme4 and brms
- adds `prior_ordered()` for ordered-factor priors that split a total effect across level increments with fixed or Dirichlet allocations
- adds Dirichlet simplex priors and moment and inverse-moment nonlocal priors
- adds `prior_factor_levels()` for attaching factor metadata to factor priors
- adds `selection_model()` for specifying how random effects and sampling errors enter selection models, and selection helpers for downstream packages
- adds a versioned parameter map stored with each fit, with `parameter_catalog()`, `parameter_catalog_resolve()`, `parameter_draws()`, `parameter_labels()`, and `parameter_mixed_posterior()` for selecting, extracting, and labelling model quantities
- adds `JAGS_deterministic_nodes()` and `JAGS_evaluate_deterministic()` for the deterministic quantities generated in JAGS models, and `parameter_source()` for referencing existing JAGS nodes
- adds `posterior_metadata()`, metadata constructors, and `posterior_transform()` for transforming posterior draws together with their supports, atoms, and densities
- adds `parameter_prior_density()`, `prior_density_ordinate()`, and `prior_ordinate_status()` for exact prior densities and ordinates of model quantities, including random-effect SDs, variance allocations, and LKJ correlations
- adds `hypothesis_BF()` for point and region hypotheses specified as expressions (including comparisons of factor levels), with the hypothesis parser `hypothesis_parse()` and `hypothesis_linear_target()`
- `hypothesis_BF()` also tests the quantities of a fitted model by the names of its parameter catalog
- adds `hypothesis_BF_model()`, `hypothesis_provider()`, and `register_hypothesis_density_backend()` for downstream model classes and density methods
- adds the `density_method` argument to `Savage_Dickey_BF()`, `plot_posterior()`, and `plot_marginal()` for posterior densities precomputed by downstream packages
- adds `JAGS_formula_coefficient_transform()` and `JAGS_formula_prior_density()` for original-scale prior densities of scaled formula coefficients
- adds `JAGS_formula_draws()` for evaluating formulas on draws without a fit
- adds `runtime_setup`, `runtime_cache`, and `worker_output` to `JAGS_fit()` and `JAGS_extend()`, `jags_modules` to `JAGS_fit()`, and the parallel cluster helpers `JAGS_runtime_cluster()` and `JAGS_runtime_cluster_stop()`
- adds `monitor` and `allow_not_assessable` to `JAGS_check_convergence()` and `autofit_control`, with per-parameter convergence diagnostics
- adds `seed`, `repetitions`, `method`, `cores`, and `nonfinite` to `JAGS_bridgesampling()`, and `JAGS_marglik_priors_rows()` for row-wise prior densities
- adds `on_failure` to the model-averaging functions for models whose marginal likelihood failed
- adds `harrell_davis_quantile()`; posterior PET-PEESE bands use Harrell-Davis quantiles
- adds `as_marginal_inference(compute_BF = FALSE)`
- exports helpers used by downstream packages (Bayes factor formatting, formula scaling, random-effect grouping, and posterior metadata helpers)

### Changes
- `JAGS_estimates_table()` uses `transform_scaled` for both fixed and random effects; original-scale rows that cannot be computed are reported as `NA` with a note
- marginal likelihoods exceeding the iteration limit keep their estimate and warn
- `JAGS_extend()` recompiles the model from the stored chain states, so repeated extensions and extensions of reloaded fits are reproducible
- estimates tables identify inclusion rows from the priors instead of the row labels
- small Bayes factor bounds print as reciprocals (e.g., `<1/10000`)

### Fixes
- prior densities of linear combinations (formula levels, marginal means, original-scale coefficients, and Savage-Dickey ordinates) use exact routes where available instead of a numerical grid, which was off by up to 20% for `multiply_by` products
- `mean()` and `var()` of truncated normal priors are computed analytically; numerical integration was wrong far in the tails and for narrow truncations
- far-tail truncated normal, t, and Cauchy priors sample and invert exactly
- the Savage-Dickey posterior ordinate is the exact kernel density at the null, reflected at exact support bounds; the previous grid estimate overestimated the density of long-tailed posteriors
- original-scale transformations of formulas such as `~ g + g:x` and `~ f/x` are derived from the fitted design; formulas without an original-scale representation stop
- prior densities of monitored formula coefficients no longer apply `multiply_by` scaling
- `expression()` formula terms are included in marginal-likelihood reconstruction
- `log(intercept)` formulas in `marginal_posterior()` and `as_mixed_posteriors(transform_scaled = TRUE)`
- `marginal_posterior()` for intercept-only formulas, predictors entering only through interactions, and model averages in which a model omits a factor term
- naming of interaction coefficients without main effects (e.g., `g:x` in `~ g + g:x`)
- formulas generating duplicate JAGS node names stop with an informative error
- `plot_posterior()` bias prior overlays follow the full conditioning of the samples
- prior overlays weight spike-and-slab priors by the inclusion probability instead of 50/50
- individual omega and weight-function posterior plots with duplicate or zero-weight bias priors
- factor plot labels of interaction cells and secondary axes of mixed density plots
- `marginal_estimates_table()` reports Bayes factor warnings of scalar parameters
- `stan_estimates_table()` handles matrix-valued parameters
- p-hacking null masses no longer become negative for far-tail cut points
- `JAGS_bridgesampling()` stops on rank-deficient bridge coordinates instead of returning a meaningless estimate
- package contrasts resolve without attaching BayesTools

### Performance
- `JAGS_fit()` resolves the JAGS executable once per session (about 0.8 s faster per fit)
- prior density curves and PET-PEESE prior bands are evaluated in batches (prior band plots take seconds instead of minutes)
- `JAGS_bridgesampling()` compiles formula and prior evaluators once per call
- two-bin cumulative weight-function priors are sampled through their Beta marginal
- saved fits no longer contain the compiled rjags model
- `as_mixed_posteriors(transform_scaled = TRUE)` caches the transformed prior densities

### Development
- vignettes generate missing or outdated fitted-model caches when knitted; set `BAYESTOOLS_VIGNETTE_REFIT=false` to render without fitting
- declares `waldo` as a test dependency

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
