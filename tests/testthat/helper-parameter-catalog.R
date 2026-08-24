.parameter_catalog_test_fit <- function(chains, prior_list,
                                        formula_design = NULL,
                                        formula_scale = NULL){

  fit <- chains
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- prior_list
  if(!is.null(formula_design)){
    attr(fit, "formula_design") <- formula_design
  }
  if(!is.null(formula_scale)){
    attr(fit, "formula_scale") <- formula_scale
  }
  attr(fit, "parameter_map") <- .bt_build_parameter_map(
    columns = colnames(chains[[1L]]),
    prior_list = prior_list,
    formula_design = formula_design,
    formula_scale = formula_scale
  )
  fit <- .bt_attach_draw_geometry(fit)
  .bt_attach_fit_contract(fit)
}

.parameter_catalog_random_summary_samples <- function(
    model_samples, prior_list, formula_design = NULL,
    mode = c("standard", "full", "raw", "none"),
    formula_scale = NULL){

  model_samples <- as.matrix(model_samples)
  chains <- coda::mcmc.list(coda::mcmc(model_samples))
  fit <- .parameter_catalog_test_fit(
    chains = chains,
    prior_list = prior_list,
    formula_design = formula_design,
    formula_scale = formula_scale
  )
  .bt_parameter_catalog_random_summary_samples(
    fit = fit,
    model_samples = model_samples,
    prior_list = prior_list,
    coordinates = parameter_coordinates(fit),
    mode = mode
  )
}
