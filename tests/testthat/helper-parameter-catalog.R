.parameter_catalog_test_fit <- function(chains, prior_list,
                                        formula_design = NULL,
                                        formula_scale = NULL){

  fit <- chains
  class(fit) <- c("BayesTools_fit", class(fit))
  if(is.list(formula_design)){
    if(is.null(formula_scale)) formula_scale <- list()
    for(parameter in names(formula_design)){
      design <- formula_design[[parameter]]
      scale <- formula_scale[[parameter]]
      if(is.null(scale)) scale <- design$formula_scale
      completed_scale <- .bt_formula_scale_finalize(scale, design,
        prior_list = prior_list, owner_scope = "fit")
      formula_design[[parameter]]$formula_scale <- completed_scale
      formula_scale[[parameter]] <- completed_scale
    }
  }
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

# Synthetic fitted object whose draws have the columns that a fit of
# 'formula_result' monitors; the draws themselves are placeholders.
.prior_monitor_test_fit <- function(formula_result, columns){

  chains <- coda::mcmc.list(coda::mcmc(matrix(
    0.5,
    nrow = 2L,
    ncol = length(columns),
    dimnames = list(NULL, columns)
  )))
  formula_scale <- formula_result$formula_scale
  .parameter_catalog_test_fit(
    chains         = chains,
    prior_list     = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design),
    formula_scale  = if(length(formula_scale) > 0L) list(mu = formula_scale)
  )
}

.prior_monitor_matrix_names <- function(name, K){

  as.vector(outer(seq_len(K), seq_len(K), function(row, column){
    paste0(name, "[", row, ",", column, "]")
  }))
}

# Draws of every public catalog quantity, evaluated on 'samples'.
.prior_monitor_catalog_draws <- function(fit, samples){

  catalog <- parameter_catalog(fit)
  public  <- catalog$quantities[!catalog$quantities$internal, , drop = FALSE]
  draws <- lapply(seq_len(nrow(public)), function(i){
    selection <- parameter_catalog_resolve(
      catalog,
      alias     = public$canonical_name[[i]],
      namespace = public$namespace[[i]],
      component = public$component[[i]]
    )
    parameter_draws(fit, selection, model_samples = samples)
  })
  names(draws) <- public$canonical_name
  draws
}
