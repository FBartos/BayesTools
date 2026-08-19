# version 0.3.1
### Features
- consolidates unreleased formula-random public parameter names as
  `(formula) owner: quantity(parameter[level], ...)`, omitting `owner: ` for
  exactly one random block and retaining it for multiple blocks. Formula lists
  with two or more unnamed
  components use `component 1`, `component 2`, and so on. Each fit stores one
  authoritative, versioned `parameter_map` containing linked backend-coordinate,
  public-quantity, and alias tables; `parameter_coordinates()` and
  `parameter_catalog()` are validated views rather than separately persisted
  metadata. The coordinate field is `coordinate_name`, while `canonical_name`
  is reserved for semantic quantities. Quantity rows declare identity,
  one-to-one, or composite source provenance while keeping backend LKJ and
  allocation coordinates internal. Public correlations use `cor`; aggregate
  allocations distinguish `sd_total` / `var_total` from `sd_common` /
  `var_common`, with components exposed as `var_prop`, `var_ratio`, and
  `sd_ratio`. Variance allocations retain a required stable internal name while
  recording their public owner and component names separately. Linked formula-
  coefficient transforms require the exact current formula-design and
  parameter-map schemas rather than accepting stale versioned metadata.
- adds `JAGS_with_draws()` for replacing fitted backend draws while preserving
  and refreshing the fit's draw geometry, allowing map-defined semantic
  quantities to be evaluated on posterior or simulated-prior coordinates
- routes public random-effect posterior and estimates-table summaries through
  `parameter_catalog()` resolution and `parameter_draws()`, and adds
  `parameter_transform()` plus authoritative forward, inverse, and Jacobian
  evaluators for one-to-one semantic coordinate maps. Standard random-effect
  tables now report only prior-facing quantities, while full tables retain all
  deterministic representations; directly specified SD ratios remain standard.
  `parameter_draws()` can also evaluate a selection on an already materialized
  posterior matrix for downstream summaries.
- adds `random_effects_marginal_factor_states()` so downstream likelihoods can
  reuse the bridge-sampling random-covariance compiler without constructing
  dense draw-by-row-by-row arrays. Supported non-row-indexed blocks now compile
  their SD and correlation metadata once and reconstruct exact factor states
  across all posterior draws in one batch.
- adds `JAGS_marglik_priors_rows()` and
  `JAGS_marglik_priors_rows_evaluator()` for exact row-preserving prior-density
  evaluation and reusable compiled evaluation. Supported scalar, independent
  or treatment factor, and Dirichlet prior lists use vectorized evaluators,
  while other prior families retain the compiled scalar route.
- adds `hypothesis_level_contrast()` for certifying atom-free pairwise level
  contrasts with an exact joint-prior ordinate, and normalizes symbolic
  point equalities such as `theta = phi` to a difference from zero; constant-
  left relations such as `0 > theta` and `0 = theta` are canonicalized to the
  equivalent parameter-left forms
- adds `JAGS_formula_coefficient_transform()` and `JAGS_formula_prior_density()` for versioned fitted-to-original coefficient maps and exact induced prior measures, including structural point values, interactions, log-intercept Jacobians, and model-mixture atoms
- adds `JAGS_formula_predictor_basis()` for exact observation-level affine
  update bases derived from the fitted parameter map and persisted formula
  design, including coefficient ordering, contrasts, scaling, and term
  multipliers. Metadata-declared nonlinear or coupled coordinates are reported
  explicitly so downstream evaluators can retain their complete formula path.
- adds a versioned `hypothesis_parse()` syntax tree with stable rendering, exact symbol discovery and rewriting, parameter-catalog resolution including unquoted non-syntactic public aliases, and direct `hypothesis_BF()` consumption without reparsing expression text
- adds a metadata-only `parameter_catalog()` view over the fitted parameter map, with classed exact resolution, validated provider extensions, and deferred `parameter_draws()` extraction for declared coordinates and random-effect summaries
- adds versioned `JAGS_draw_geometry()` metadata and parameter-map-based `JAGS_materialize_draws()` reconstruction, including exact structural point-prior values, preserved chain timing, valid zero-column public draws, and a private deterministic backend anchor for models with no ordinary monitor
- adds an injective UTF-8 `JAGS_parameter_encode()` / `JAGS_parameter_decode()` semantic identifier, a persisted formula name map, and a single strict `JAGS_fit_contract()` schema so downstream packages can consume formula and fitted metadata without parsing established JAGS column names
- adds `prior_density_ordinate()` for exact-value structural classification of scalar prior and induced `prior_linear_density` ordinates, including continuous limits, point masses, deterministic mixtures, analytic normal combinations, and supported named transformations
- adds the `prior_random()` interface for formula random effects, including `random_block()`, `random_covariance()`, `random_monitor()`, `random_new_levels()`, `random_variance_allocation()`, and `allocation_ref()` helpers for specifying random-effect standard deviation priors, covariance structures, monitoring policy, and total-variance allocation priors
- adds lme4-like formula random-effect parsing through `reformulas`, including ordinary and independent random effects, named random-effect blocks, nested grouping expressions, factor random slopes, and structured covariance shortcuts for diagonal, shared-SD independent, unstructured, compound-symmetry, heterogeneous compound-symmetry, discrete AR(1), heterogeneous AR(1), and continuous-time AR(1) random effects
- adds LKJ correlation priors for unstructured random-effect covariance matrices via `prior_lkj()` and `JAGS_lkj_corr_cholesky()`, using the package-shipped compiled JAGS backend
- adds `BayesTools_load_JAGS_module()` and package compilation support for the BayesTools JAGS module used by generated LKJ-Cholesky syntax
- adds `formula_random_prior_list` to `JAGS_fit()` and `JAGS_bridgesampling()` so formula random effects can be fitted and bridge sampled through the explicit `prior_random()` interface
- adds `JAGS_predict_formula()` for fixed, conditional, and marginalized formula prediction, with selectable random-effect blocks, explicit new-level policies, and covariance- or sampling-based marginal output
- adds block-level `"noncentered"`, `"centered"`, and `"auto"` random-effect parameterizations without changing the prior or semantic output contracts; structured scalar-correlation blocks use scalable direct recurrences and expose dense correlation draws only on explicit reconstruction
- adds `random_group_covariance()` for known covariance or correlation kernels across random-effect grouping levels
- adds `random_effects_compile()` and `random_effects_marginal_variance_factors()` as the public contract for structurally marginalized random-effect blocks; unsupported factor representations signal the classed `BayesTools_random_effects_marginal_variance_unavailable` condition
- adds `random_effects_marginal_vcov()` for posterior observation-level covariance or variance draws implied by formula random effects; `diagonal_only = TRUE` provides scalable marginal-variance extraction without allocating dense draw-by-row-by-row covariance arrays
- adds `fit_backend_fingerprint()` for stable cache invalidation when fitted-model backend code or native binaries change
- adds bridge-sampling support for formula random effects by using standardized latent random effects, scalar correlation coordinates, and LKJ primitive coordinates as bridge parameters
- adds semantic random-effect summaries to `runjags_estimates_table()` / `JAGS_estimates_table()` through `random_effects_summary`, `random_effects_metadata`, `remove_random_effects`, `keep_random_effects`, `remove_random_structures`, and `keep_random_structures`
- adds random-effect parameter filters such as `"random"`, `"random_sd"`,
  `"random_cor"`, `"random_var_prop"`, `"random_var_ratio"`,
  `"random_allocation"`, and `"random_sd_ratio"` for estimates tables
- adds compact print methods for the public random-effect specification helpers, LKJ priors, parameter sources, random-SD sources, and variance-allocation references
- adds `prior_ordered()` for ordered-factor priors that separate a scalar total effect from fixed or Dirichlet allocations across cumulative level increments
- adds Dirichlet simplex priors via `prior("dirichlet", ...)` / `prior("simplex", ...)`, including random generation, log-density, marginal distribution helpers, JAGS syntax, initialization, posterior extraction, and bridge-sampling support
- adds moment and inverse-moment nonlocal priors, including R density/distribution/quantile/RNG helpers and compiled JAGS-module support for prior-only and formula-model fitting
- adds `hypothesis_BF()` for expression-based point and region hypotheses on posterior/prior draws or marginal posterior objects, including comparisons between named factor levels
- adds a `RandomEffects` vignette comparing BayesTools formula random effects with lme4 and rstanarm examples

### Changes
- Diagonal and unstructured random-factor blocks now retain one SD per
  generated coefficient for mean-difference, orthonormal, and ordered bases.
  Scalar SD templates expand to indexed independent priors, while explicitly
  multivariate factor priors retain their joint prior and map to the same
  coefficient coordinates.
- Random-coefficient blocks without an explicit
  `random_block(contrasts = ...)` override now reuse the concrete contrast
  metadata resolved for the same fixed factor.
- Formula random-effect covariance priors are now completed after the parsed
  structure and dimension are known. Omitted US/UN correlations use `LKJ(1)`;
  omitted CS/HCS, AR1/HAR, and CAR correlations use uniform priors over their
  complete admissible raw-correlation intervals. Explicit scalar correlation
  priors retain the existing Fisher-z default scale.
- `JAGS_bridgesampling()` bypasses formula reconstruction and bridge-context
  replay for ordinary non-formula models while retaining their complete prior
  and likelihood target.
- `JAGS_bridgesampling()` can pass an exact nodes-only bridge context, compiles
  invariant random-effect prior, replay, and SD-binding plans once, avoids
  reconstructing allocation nodes already supplied by formula priors, retains
  the complete context option, can select an exact named node subset without
  flattening unrelated state, and forwards an explicit `cores` setting to
  `bridgesampling::bridge_sampler()`. Marginal random-effect covariance replay
  now consumes the already reconstructed natural formula-prior parameters, and
  compiled scalar priors call the same validated primitive density calculation
  without repeating public-interface validation for every bridge row.
- `JAGS_bridgesampling()` can exactly integrate selected fitted sampled
  Gaussian formula random-effect blocks during bridge evaluation. It removes
  only their standardized latent coordinates, retains all covariance
  parameters and priors, and supplies either the full draw-specific `ZGZ'`
  covariance or a validated exact block-factor representation to the
  likelihood callback. The factor contract includes full coefficient
  covariance for known group kernels and row-specific external SD scales. Its
  nodes-only bridge context can omit the coefficient covariance already
  represented exactly by a supplied factor, while the complete generic context
  and covariance-valued evaluator remain available. Sampled SD positions,
  covariance structure, and known-group kernels are compiled once and reused
  without changing source-parameter precedence or covariance reconstruction;
  an opt-in compact factor-state contract separates this invariant geometry
  from draw-varying factors while retaining the full context option.
- Bridge evaluation of structured random-effect factors now reconstructs and
  validates each exact coefficient Cholesky factor once instead of rebuilding
  the same factor a second time for correlation validation, and compiles
  invariant scalar-correlation metadata, coordinates, and structural support
  before repeated bridge states are evaluated. Exact compound-symmetry and
  Markov Cholesky recurrences in this bridge hotspot are evaluated by a
  package-native kernel after the same R-level parameter and support checks.
  The same compiled evaluator now serves posterior reconstruction, prediction,
  covariance summaries, and plotting consumers; the pure-R subset recurrence
  remains the exact reference and diagnostic fallback.
- Batched marginal random-effect factors now handle direct posterior and
  row-indexed SD sources for every supported random structure, cache shared
  allocation replay within a batch, and reconstruct simplex weights from their
  authoritative gamma auxiliaries. Valid zero auxiliary components retain
  their exact boundary weights, while invalid or non-normalized coordinates
  fail without renormalization. Random-effect prediction samples directly from
  persisted Cholesky factors instead of decomposing and repairing reconstructed
  covariance matrices. Fisher-z and bounded-logit scalar correlations use
  their declared transforms exactly; numerically saturated boundary values are
  rejected rather than clamped into the admissible interval.
- The compact bridge factor-state contract now labels exact diagonal and
  Markov coefficient structures. AR1, HAR, and CAR states expose their complete
  coefficient scales, adjacent transitions, and innovation variances alongside
  the unchanged full coefficient factor, allowing downstream exact linear-time
  likelihood evaluation without removing the generic covariance contract.
- `JAGS_bridgesampling()` exposes `repetitions` and `method` alongside its
  existing bridge controls and always derives the effective sample size from
  the fitted chains.
- hypothesis parsing now accepts unquoted colon-separated formula interaction
  level references such as `factor:moderator[level]`
- `as_marginal_inference(compute_BF = FALSE)` now returns averaged and conditional marginal posteriors without computing inclusion Bayes factors
- deterministic scalar, multivariate, and factor point-prior parameters are now monitored and retained as structural posterior columns
- `JAGS_bridgesampling()` now evaluates zero-dimensional fixed-parameter models exactly instead of requiring an artificial sampled parameter
- `JAGS_fit()` and `JAGS_extend()` now carry generated `add_parameters`, `required_packages`, and `jags_modules` metadata so formula-generated monitors and JAGS modules remain available during fitting, extension, convergence checks, and parallel execution
- `JAGS_estimates_table(transform_scaled = TRUE)` now derives formula random-effect SD and correlation summaries on the transformed original scale when fitted formula-scale metadata are available
- `JAGS_evaluate_formula()` can evaluate fitted random-effect formulas for existing grouping levels when latent random effects or group-level coefficients were monitored
- `JAGS_check_convergence()` ignores model indicator variables by default and excludes generated auxiliary monitor parameters from convergence checks unless explicitly requested
- `Savage_Dickey_BF()` and marginal-posterior `hypothesis_BF()` point-null tests now use boundary-reflected KDE ordinates when exact posterior-support metadata is available and validated against the samples; boundary-null Bayes factors for bounded parameters can therefore differ from version 0.3.0 standard-KDE results
- `as_marginal_inference()` conditional marginal summaries use active-subset conditioning: each marginal level conditions only on requested parameters with nonzero weight in that level's linear combination, and levels with no active requested conditionals use the fully averaged context

### Fixes
- preserves exact numeric values when hypothesis accessors consume an existing
  syntax tree, and retains authoritative attached prior-density metadata when
  constructing marginal posteriors. Scalar random blocks no longer require
  unused correlation state, including when their SD source is row-indexed,
  and Dirichlet auxiliary coordinates must be
  strictly positive before simplex reconstruction.
- generates stored LKJ primitive coordinates when drawing formula priors so
  public semantic `cor(...)` summaries can be reconstructed exactly.
- preserves matrix dimensions while validating one-coefficient random-effect
  correlation Cholesky draws.
- evaluates transformed prior density grids on the displayed plotting range,
  including valid zero boundaries, and propagates that range to marginal-prior
  plots instead of exponentiating user-facing axis limits as fitted-scale
  coordinates
- resolves one density-to-probability mapping from `ylim` and `ylim2` for the
  secondary axis and all point masses, and reuses it for base posterior and
  prior line overlays across fitted objects; out-of-range overlays warn instead
  of silently rescaling the active plot
- exposes every fitted factor level and interaction cell as a named coefficient-level quantity, using injective named-cell components, semantic display labels, direct coordinates, structural zeroes, or the persisted term-only contrast transformation for treatment, independent, mean-difference, orthonormal, and ordered encodings; incomplete coordinate maps now fail closed and boundary whitespace remains resolvable
- validates the parameter map atomically, reusing identity random-effect coordinates with their sampled or structural provenance, hiding transformed fitted-scale and private implementation rows, checking every native extraction dependency, and reserving the `BayesTools` provider namespace
- preserves exact finite transformed-support endpoints through structural
  transformation provenance, avoiding round-trip misclassification for
  transformations such as negative powers of truncated positive priors
- corrects the zero-boundary exponent for negative powers of inverse-gamma variables, distinguishing infinite, finite nonzero, and zero transformed densities according to the inverse-gamma shape and power
- ensures targeted `JAGS_check_convergence()` checks still include eligible product-space indicators when `check_indicators = TRUE`, while auxiliary inclusion-probability coordinates remain opt-in
- makes optional model-indicator convergence checks label invariant by diagnosing binary state occupancy and each observed categorical state separately, while recognizing one-state indicators as structural only when the prior fixes their support
- derives formula coefficient transforms from parameter-coordinate rows with the `fixed_coefficient` role, so sampled and marginalized random-effect priors no longer make mixed-formula coefficient transforms and prior densities fail with a source mismatch
- makes `JAGS_extend()` validate the declared fit contract and every preserved metadata object before backend work, preventing stale formula designs from being relabeled as current after an early exit or extension
- restores literal `expression()` formula terms that reference sampled scalar or one-dimensional indexed parameters, including individually monitored indexed coordinates, persisting their data and parameter dependencies for draw-aware fitted/new-data prediction and marginal-likelihood reconstruction while rejecting opaque JAGS-derived nodes that cannot be replayed
- resolves level-qualified hypothesis symbols through semantic catalog components, including multiple factor levels in the same hypothesis, instead of dropping the level and reporting the factor term as ambiguous
- fixes bridge-sampling bound validation by requiring named lower and upper bounds that match `add_parameters`
- hardens JAGS build discovery by adding non-default Unix rpath flags, enforcing JAGS >= 4.3.0 when the version is discoverable, removing unused JAGS major-version compile defines, and making Windows JAGS root detection robust to spaces and semantic version ordering
- hardens fixture test infrastructure by routing hypothesis bridge comparisons through shared cache-currency validation, registering indexed JAGS parameter tests in the unit profile, and expanding random-effect and fitted-metadata source hashes that invalidate fitted-model caches
- resolves package-defined factor contrasts inside the namespace across fixed, random-effect, prediction, and marginal-posterior design matrices, so `BayesTools::` calls do not require attaching the package

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

