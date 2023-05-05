## version 0.2.14
### Features
- adding `contrast = "meandif"` to the `prior_factor` function which generates identical prior distributions for difference between the grand mean and each factor level
- adding `contrast = "independent"` to the `prior_factor` function which generates independent identical prior distributions for each factor level
- `remove_column` function for removing columns from `BayesTools_table` objects without breaking the attributes etc...
- adding empty table functions (https://github.com/FBartos/BayesTools/issues/10)
- adding `remove_parameters` argument to `model_summary_table()`

### Changes
- `contrast = "meandif"` is now the default setting for `prior_factor` function
- depreciating `transform_orthonormal` argument in favor of more general `transform_factors` argument 
- switching `dummy` contrast/factor attributes to `treatment` for consistency (https://github.com/FBartos/BayesTools/issues/23)

### Fixes
- zero length inputs to `check_bool()`, `check_char()`, `check_real()`, `check_int()`, and `check_list()` do not throw error if `allow_NULL = TRUE`
- properly aggregating identical priors in the plotting function (previously overlying multiple spikes on top of each other when attributes did not match)

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
