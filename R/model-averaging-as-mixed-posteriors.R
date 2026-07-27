#' @title Export BayesTools JAGS model posterior distribution as model-average posterior distributions via \code{mix_posteriors}
#'
#' @description Creates a model-averages posterior distributions on a single
#' model that allows mimicking the [mix_posteriors] functionality. This function
#' is useful when the model-averaged ensemble is based on [prior_spike_and_slab]
#' or [prior_mixture] priors - the model-averaging is done within the model.
#'
#' @param model model fit via the [JAGS_fit] function
#' @param conditional a character vector of parameters to be conditioned on
#' @param conditional_rule a character string specifying the rule for conditioning.
#' Either "AND" or "OR". Defaults to "AND".
#' @param force_plots temporal argument allowing to generate conditional posterior samples
#' suitable for prior and posterior plots. Only available when conditioning on a
#' single parameter.
#' @param transform_scaled whether to transform samples from standardized (scaled) to
#' original (unscaled) scale. When \code{TRUE}, posterior samples are
#' transformed, and the result can be directly passed to [plot_posterior] which will
#' automatically detect the transformation and use transformed deterministic prior densities.
#' Requires a model fitted with \code{formula_scale_list}. Defaults to \code{FALSE}.
#' @param n_prior_samples controls the numerical grid used for transformed
#' prior densities when \code{transform_scaled = TRUE}. Defaults to 10000.
#' @inheritParams ensemble_inference
#'
#' @return \code{as_mix_posteriors} returns a named list of mixed posterior
#' distributions (either a vector of matrix).
#'
#' @seealso [mix_posteriors]
#'
#' @name as_mixed_posteriors
#' @export
as_mixed_posteriors <- function(model, parameters, conditional = NULL, conditional_rule = "AND", force_plots = FALSE,
                                 transform_scaled = FALSE, n_prior_samples = 10000){

  # check input
  if(!inherits(model, "BayesTools_fit"))
    stop("'model' must be a 'BayesTools_fit'")
  check_char(parameters, "parameters", check_length = FALSE)
  check_char(conditional, "conditional", check_length = FALSE, allow_values = c(parameters, "PET", "PEESE", "PETPEESE", "omega", "phacking", "alpha", "pi_null"), allow_NULL = TRUE)
  check_char(conditional_rule, "conditional_rule", allow_values = c("AND", "OR"))
  check_bool(transform_scaled, "transform_scaled")
  check_int(n_prior_samples, "n_prior_samples", lower = 1)

  # extract the list of priors
  priors <- attr(model, "prior_list")
  prior_density_priors <- priors
  formula_scale <- attr(model, "formula_scale")
  condition_event <- .condition_event(
    prior_list        = priors,
    conditional       = conditional,
    conditional_rule  = conditional_rule
  )

  # extract the samples
  model_samples <- .extract_posterior_samples(model, as_list = FALSE)
  if(!is.matrix(model_samples)){
    # deal with automatic coercion into a vector in case of a single predictor
    model_samples <- matrix(model_samples, ncol = 1)
    colnames(model_samples) <- model$monitor
  }
  posterior_density_sources <- .posterior_density_sources(model, model_samples)
  posterior_ordinate_sources <- .posterior_ordinate_sources(model, model_samples)

  # apply conditioning
  if(length(condition_event[["conditional"]]) > 0){

    # subset the posterior distribution
    conditioning_samples <- .condition_event_posterior_mask(
      event        = condition_event,
      prior_list   = priors,
      model_samples = model_samples
    )

    if(sum(conditioning_samples) == 0){
      warning("No samples left after conditioning.", call. = FALSE, immediate. = TRUE)
      return(list())
    }


    model_samples <- model_samples[conditioning_samples,,drop=FALSE]
  }

  # apply scale transformation to posterior samples if requested
  if(transform_scaled && !is.null(formula_scale) && length(formula_scale) > 0){
    model_samples <- transform_scale_samples(model_samples, formula_scale)
    posterior_density_sources <- list()
    posterior_ordinate_sources <- list()
  }

  out    <- list()

  for(p in seq_along(parameters)){

    # prepare parameter specific values
    temp_parameter <- parameters[p]
    temp_prior     <- priors[[temp_parameter]]

    if(is.prior.spike_and_slab(temp_prior)){
      # spike and slab priors
      out[[temp_parameter]] <- .as_mixed_posteriors.spike_and_slab(model_samples, temp_prior, temp_parameter)

    }else if(is.prior.mixture(temp_prior)){
      # mixture priors
      out[[temp_parameter]] <- .as_mixed_posteriors.mixture(model_samples, temp_prior, temp_parameter, condition_event[["conditional"]])

    }else if(is_prior_phacking(temp_prior)){
      # p-hacking priors
      out[[temp_parameter]] <- .as_mixed_posteriors.phacking(model_samples, temp_prior, temp_parameter)

    }else if(is_prior_bias(temp_prior)){
      # composed publication-bias priors
      out[[temp_parameter]] <- .as_mixed_posteriors.bias(model_samples, temp_prior, temp_parameter, condition_event[["conditional"]])

    }else if(is.prior.weightfunction(temp_prior)){
      # weight functions
      out[[temp_parameter]] <- .as_mixed_posteriors.weightfunction(model_samples, temp_prior, temp_parameter)

    }else if(is.prior.factor(temp_prior)){
      # factor priors
      out[[temp_parameter]] <- .as_mixed_posteriors.factor(model_samples, temp_prior, temp_parameter)

    }else if(is.prior.vector(temp_prior)){
      # vector priors
      out[[temp_parameter]] <- .as_mixed_posteriors.vector(model_samples, temp_prior, temp_parameter)

    }else if(is.prior.simple(temp_prior)){
      # simple priors
      out[[temp_parameter]] <- .as_mixed_posteriors.simple(model_samples, temp_prior, temp_parameter)

    }else{
      stop("The posterior samples cannot be mixed: unsupported prior distributions.")
    }

    # add formula relevant information
    if(!is.null(attr(temp_prior, which = "parameter"))){
      class(out[[temp_parameter]]) <- c(class(out[[temp_parameter]]), "mixed_posteriors.formula")
      attr(out[[temp_parameter]], "formula_parameter")  <- attr(temp_prior, which = "parameter")
    }
    if(transform_scaled && !is.null(formula_scale) && length(formula_scale) > 0){
      out[[temp_parameter]] <- .posterior_support_drop(
        out[[temp_parameter]],
        recursive = TRUE
      )
    }

    # add conditioning information
    out[[temp_parameter]] <- .condition_event_set_attributes(
      out[[temp_parameter]],
      condition_event
    )

    out[[temp_parameter]] <- .posterior_density_attach(
      samples            = out[[temp_parameter]],
      sources            = posterior_density_sources,
      parameter          = temp_parameter,
      conditional        = condition_event[["conditional"]],
      conditional_rule   = conditional_rule,
      condition_key      = condition_event[["condition_key"]],
      allow_unlabeled    = length(parameters) == 1L
    )
    out[[temp_parameter]] <- .posterior_ordinate_attach(
      samples            = out[[temp_parameter]],
      sources            = posterior_ordinate_sources,
      parameter          = temp_parameter,
      conditional        = condition_event[["conditional"]],
      conditional_rule   = conditional_rule,
      condition_key      = condition_event[["condition_key"]],
      allow_unlabeled    = length(parameters) == 1L
    )

  }

  attr(out, "prior_list")       <- priors
  out <- .condition_event_set_attributes(out, condition_event)
  if(length(posterior_density_sources) > 0L){
    attr(out, "posterior_density") <- posterior_density_sources[[1]]
    if(length(posterior_density_sources) > 1L){
      attr(out, "posterior_densities") <- posterior_density_sources[-1]
    }
  }
  if(length(posterior_ordinate_sources) > 0L){
    attr(out, "posterior_ordinate") <- posterior_ordinate_sources[[1]]
    if(length(posterior_ordinate_sources) > 1L){
      attr(out, "posterior_ordinates") <- posterior_ordinate_sources[-1]
    }
  }

  # propagate formula_scale attribute for transform_scaled support
  if(!is.null(formula_scale)){
    attr(out, "formula_scale") <- formula_scale
  }

  # generate and store transformed prior densities if requested
  if(transform_scaled && !is.null(formula_scale) && length(formula_scale) > 0){
    prior_densities <- .generate_transformed_prior_densities(
      prior_list       = prior_density_priors,
      column_names     = colnames(model_samples),
      n_grid           = n_prior_samples,
      formula_scale    = formula_scale,
      conditional      = condition_event[["conditional"]],
      conditional_rule = conditional_rule,
      condition_event  = condition_event
    )
    attr(out, "prior_densities")       <- prior_densities
    attr(out, "prior_density_context") <- attr(prior_densities, "context")
    attr(out, "transform_scaled")      <- TRUE
  }else{
    attr(out, "prior_density_context") <- .prior_density_build_context(
      prior_list       = prior_density_priors,
      column_names     = colnames(model_samples),
      n_grid           = n_prior_samples,
      conditional      = condition_event[["conditional"]],
      conditional_rule = conditional_rule,
      condition_event  = condition_event
    )
  }
  if(length(condition_event[["conditional"]]) > 0L ||
     (transform_scaled && !is.null(formula_scale) && length(formula_scale) > 0L)){
    prior_density_context <- attr(out, "prior_density_context", exact = TRUE)
    for(parameter in names(out)){
      out[[parameter]] <- .posterior_support_set_from_prior_context(
        out[[parameter]],
        prior_density_context
      )
    }
  }

  class(out) <- c(class(out), "as_mixed_posteriors", "mixed_posteriors")
  return(out)
}

.as_mixed_posteriors.simple         <- function(model_samples, prior, parameter){

  # check input
  check_char(parameter, "parameter", check_length = FALSE)

  # gather information about the prior distribution
  prior_info <- list(
    "interaction"       = .is_prior_interaction(prior),
    "interaction_terms" = attr(prior, "interaction_terms")
  )

  # prepare output objects
  samples <- model_samples[, parameter]

  # format the output
  samples <- unname(samples)
  attr(samples, "sample_ind") <- FALSE
  attr(samples, "models_ind") <- rep(1, length(samples))
  attr(samples, "parameter")  <- parameter
  attr(samples, "prior_list") <- prior
  attr(samples, "interaction")       <- if(length(prior_info) == 0) FALSE else prior_info[["interaction"]]
  attr(samples, "interaction_terms") <- prior_info[["interaction_terms"]]
  samples <- .posterior_support_set_from_prior_list(samples, prior)
  samples <- .posterior_atoms_set(
    samples,
    .posterior_atoms_from_priors(
      prior,
      1,
      n_columns = 1L,
      column_names = parameter,
      source = "single_model_structure"
    )
  )
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.simple")

  return(samples)
}
.as_mixed_posteriors.vector         <- function(model_samples, prior, parameter){

  # check input
  check_char(parameter, "parameter", check_length = FALSE)

  # gather information about the prior distribution
  K <- prior$parameters[["K"]]
  if(length(K) != 1)
    stop("all vector prior must be of the same length")

  # prepare output objects
  if(K == 1){
    samples <- model_samples[, parameter, drop = FALSE]
  }else{
    samples <- model_samples[, paste0(parameter,"[",1:K,"]"), drop = FALSE]
  }

  rownames(samples) <- NULL
  colnames(samples) <- paste0(parameter,"[",1:K,"]")
  attr(samples, "sample_ind") <- FALSE
  attr(samples, "models_ind") <- rep(1, nrow(samples))
  attr(samples, "parameter")  <- parameter
  attr(samples, "prior_list") <- prior
  samples <- .posterior_support_set_columns_from_prior_list(samples, prior)
  samples <- .posterior_atoms_set(
    samples,
    .posterior_atoms_from_priors(
      prior,
      1,
      n_columns = K,
      column_names = colnames(samples),
      source = "single_model_structure"
    )
  )
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.vector")

  return(samples)
}
.as_mixed_posteriors.factor         <- function(model_samples, prior, parameter){

  # check input
  check_char(parameter, "parameter", check_length = FALSE)
  prior <- .complete_factor_metadata(prior, parameter)

  # gather information about the prior distribution
  prior_info <- list(
    "levels"            = .get_prior_factor_levels(prior),
    "level_names"       = .get_prior_factor_level_names(prior),
    "interaction"       = .is_prior_interaction(prior),
    "interaction_terms" = attr(prior, "interaction_terms"),
    "term_components"   = attr(prior, "term_components"),
    "factor_terms"      = attr(prior, "factor_terms"),
    "factor_contrasts"  = attr(prior, "factor_contrasts"),
    "factor_design"     = attr(prior, "factor_design"),
    "factor_cell_names" = attr(prior, "factor_cell_names"),
    "treatment"         = is.prior.treatment(prior),
    "independent"       = is.prior.independent(prior),
    "orthonormal"       = is.prior.orthonormal(prior),
    "meandif"           = is.prior.meandif(prior),
    "ordered"           = is.prior.ordered(prior)
  )


  if(prior_info[["ordered"]]){

    coefficient_names <- .JAGS_prior_factor_names(parameter, prior)
    samples <- model_samples[, coefficient_names, drop = FALSE]
    ordered_total_indicator <- NULL
    if(is.prior.spike_and_slab(prior$total)){
      indicator_name <- paste0(
        .prior_ordered_total_name(parameter),
        "_indicator"
      )
      if(!indicator_name %in% colnames(model_samples)){
        stop(
          "The fitted samples for ordered factor '", parameter,
          "' do not contain the required total-prior indicator '",
          indicator_name, "'. Refit the model with this package version.",
          call. = FALSE
        )
      }
      ordered_total_indicator <- model_samples[, indicator_name]
    }

    rownames(samples) <- NULL
    colnames(samples) <- coefficient_names
    attr(samples, "sample_ind") <- FALSE
    attr(samples, "models_ind") <- rep(1, nrow(samples))
    attr(samples, "parameter")  <- parameter
    attr(samples, "prior_list") <- prior
    if(!is.null(ordered_total_indicator)){
      attr(samples, "ordered_total_indicator") <-
        as.integer(ordered_total_indicator)
    }
    class(samples) <- c("mixed_posteriors", "mixed_posteriors.factor", "mixed_posteriors.vector")

  }else if(prior_info[["treatment"]]){

    if(prior_info[["levels"]] == 1){

      samples <- .as_mixed_posteriors.simple(model_samples, prior, parameter)
      samples <- matrix(samples, ncol = 1)

    }else{

      samples <- lapply(1:prior_info[["levels"]], function(i) .as_mixed_posteriors.simple(model_samples, prior, paste0(parameter, "[", i, "]")))
      samples <- do.call(cbind, samples)

    }

    level_names <- prior_info[["level_names"]]
    if(is.list(level_names)){
      level_names <- lapply(level_names, function(x) x[-1])
    }else{
      level_names <- level_names[-1]
    }

    rownames(samples) <- NULL
    colnames(samples) <- .format_factor_level_parameter_names(parameter, level_names, ncol(samples))
    attr(samples, "sample_ind") <- FALSE
    attr(samples, "models_ind") <- rep(1, nrow(samples))
    attr(samples, "parameter")  <- parameter
    attr(samples, "prior_list") <- prior
    class(samples) <- c("mixed_posteriors", "mixed_posteriors.factor", "mixed_posteriors.vector")

  }else if(prior_info[["independent"]]){

    if(prior_info[["levels"]] == 1){

      samples <- .as_mixed_posteriors.simple(model_samples, prior, parameter)
      samples <- matrix(samples, ncol = 1)

    }else{

      samples <- lapply(1:prior_info[["levels"]], function(i) .as_mixed_posteriors.simple(model_samples, prior, paste0(parameter, "[", i, "]")))
      samples <- do.call(cbind, samples)

    }

    rownames(samples) <- NULL
    colnames(samples) <- .format_factor_level_parameter_names(parameter, prior_info[["level_names"]], ncol(samples))
    attr(samples, "sample_ind") <- FALSE
    attr(samples, "models_ind") <- rep(1, nrow(samples))
    attr(samples, "parameter")  <- parameter
    attr(samples, "prior_list") <- prior
    class(samples) <- c("mixed_posteriors", "mixed_posteriors.factor", "mixed_posteriors.vector")

  }else if(prior_info[["orthonormal"]] | prior_info[["meandif"]]){

    prior$parameters[["K"]] <- prior_info[["levels"]]
    samples <- .as_mixed_posteriors.vector(model_samples, prior, parameter)
    class(samples) <- c(class(samples), "mixed_posteriors.factor")

  }

  attr(samples, "levels")            <- prior_info[["levels"]]
  attr(samples, "level_names")       <- prior_info[["level_names"]]
  attr(samples, "interaction")       <- if(length(prior_info) == 0) FALSE else prior_info[["interaction"]]
  attr(samples, "interaction_terms") <- prior_info[["interaction_terms"]]
  attr(samples, "term_components")   <- prior_info[["term_components"]]
  attr(samples, "factor_terms")      <- prior_info[["factor_terms"]]
  attr(samples, "factor_contrasts")  <- prior_info[["factor_contrasts"]]
  attr(samples, "factor_design")     <- prior_info[["factor_design"]]
  attr(samples, "factor_cell_names") <- prior_info[["factor_cell_names"]]
  attr(samples, "treatment")         <- prior_info[["treatment"]]
  attr(samples, "independent")       <- prior_info[["independent"]]
  attr(samples, "orthonormal")       <- prior_info[["orthonormal"]]
  attr(samples, "meandif")           <- prior_info[["meandif"]]
  attr(samples, "ordered")           <- prior_info[["ordered"]]
  attr(samples, "ordered_metadata")  <- attr(prior, "ordered_metadata")

  if(isTRUE(prior_info[["treatment"]]) || isTRUE(prior_info[["independent"]])){
    factor_support <- .posterior_support_from_prior_list(prior)
    if(!is.null(factor_support) && !is.null(colnames(samples))){
      attr(samples, "posterior_support") <- stats::setNames(
        rep(list(factor_support), ncol(samples)),
        colnames(samples)
      )
    }
  }

  ordered_atoms <- .posterior_atoms_from_ordered_total(
    prior,
    n_columns = ncol(samples),
    column_names = colnames(samples),
    indicator = if(prior_info[["ordered"]]){
      attr(samples, "ordered_total_indicator", exact = TRUE)
    }else{
      NULL
    },
    source = "ordered_total_posterior_indicator"
  )
  samples <- .posterior_atoms_set(
    samples,
    if(is.null(ordered_atoms)){
      .posterior_atoms_from_priors(
        prior,
        1,
        n_columns = ncol(samples),
        column_names = colnames(samples),
        source = "single_model_structure"
      )
    }else{
      ordered_atoms
    }
  )

  return(samples)
}
.as_mixed_posteriors.weightfunction <- function(model_samples, prior, parameter){

  # check input
  check_char(parameter, "parameter", check_length = FALSE)


  # obtain mapping for the weight coefficients
  omega_info    <- .weightfunction_mapping_info(list(prior))
  omega_mapping <- omega_info$mapping
  omega_names   <- omega_info$names
  omega_par     <- omega_info$pars

  # prepare output objects
  samples <- model_samples[, omega_par, drop = FALSE]

  rownames(samples) <- NULL
  colnames(samples) <- omega_names
  attr(samples, "sample_ind") <- FALSE
  attr(samples, "models_ind") <- rep(1, nrow(samples))
  attr(samples, "parameter")  <- parameter
  attr(samples, "prior_list") <- prior
  samples <- .weightfunction_set_omega_context(samples, omega_info)
  samples <- .posterior_support_set_weightfunction_columns(samples, prior, omega_info)
  samples <- .posterior_atoms_set(
    samples,
    .posterior_atoms_from_priors(
      prior,
      1,
      n_columns = ncol(samples),
      column_names = colnames(samples),
      source = "single_model_structure",
      null_location = 1
    )
  )
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.weightfunction")

  return(samples)
}
.as_mixed_posteriors.phacking      <- function(model_samples, prior, parameter){

  # check input
  check_char(parameter, "parameter", check_length = FALSE)

  par_names <- intersect(.phacking_report_parameter(prior), colnames(model_samples))
  samples   <- model_samples[, par_names, drop = FALSE]

  rownames(samples) <- NULL
  attr(samples, "sample_ind") <- FALSE
  attr(samples, "models_ind") <- rep(1, nrow(samples))
  attr(samples, "parameter")  <- parameter
  attr(samples, "prior_list") <- prior
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.phacking")

  return(samples)
}
.as_mixed_posteriors.bias          <- function(model_samples, prior, parameter, conditional){

  # check input
  check_char(parameter, "parameter", check_length = FALSE)

  branch_info   <- .selection_prior_branch_info(prior)
  has_selection <- vapply(branch_info, function(x) !is.null(x$selection), logical(1))
  has_phacking  <- vapply(branch_info, function(x) !is.null(x$phacking),  logical(1))

  out_names <- NULL
  par_names <- NULL

  if(any(has_selection)){
    selection_priors <- lapply(branch_info[has_selection], function(x) x$selection)
    omega_info       <- .weightfunction_mapping_info(selection_priors, one_sided = TRUE)
    omega_names      <- omega_info$names
    omega_par        <- omega_info$pars
  }
  if(any(has_phacking)){
    phacking_priors <- lapply(branch_info[has_phacking], function(x) x$phacking)
    phacking_par    <- .selection_phacking_report_parameters(phacking_priors)
    phacking_names  <- phacking_par
  }

  if(length(conditional) > 0 && any(c("omega", "phacking", "alpha", "pi_null") %in% conditional)){
    if("omega" %in% conditional && any(has_selection)){
      out_names <- c(out_names, omega_names)
      par_names <- c(par_names, omega_par)
    }
    if("phacking" %in% conditional && any(has_phacking)){
      out_names <- c(out_names, phacking_names)
      par_names <- c(par_names, phacking_par)
    }
    if("alpha" %in% conditional && any(has_phacking)){
      out_names <- c(out_names, "alpha")
      par_names <- c(par_names, "alpha")
    }
    if("pi_null" %in% conditional && any(has_phacking)){
      out_names <- c(out_names, "pi_null")
      par_names <- c(par_names, "pi_null")
    }
  }else{
    if(any(has_selection)){
      out_names <- c(out_names, omega_names)
      par_names <- c(par_names, omega_par)
    }
    if(any(has_phacking)){
      out_names <- c(out_names, phacking_names)
      par_names <- c(par_names, phacking_par)
    }
  }

  if(is.null(par_names)){
    par_names <- character()
    out_names <- character()
  }
  keep_unique <- !duplicated(par_names)
  par_names <- par_names[keep_unique]
  out_names <- out_names[keep_unique]
  keep_par <- par_names %in% colnames(model_samples)
  par_names <- par_names[keep_par]
  out_names <- out_names[keep_par]
  samples   <- model_samples[, par_names, drop = FALSE]

  rownames(samples) <- NULL
  colnames(samples) <- out_names
  attr(samples, "sample_ind") <- FALSE
  attr(samples, "models_ind") <- rep(1, nrow(samples))
  attr(samples, "parameter")  <- parameter
  attr(samples, "prior_list") <- prior
  if(any(has_selection)){
    samples <- .weightfunction_set_omega_context(samples, omega_info)
  }
  class(samples) <- c("mixed_posteriors", "mixed_posteriors.bias")

  return(samples)
}
.as_mixed_posteriors.spike_and_slab <- function(model_samples, prior, parameter){

  # check input
  check_char(parameter, "parameter", check_length = FALSE)

  prior_variable <- .get_spike_and_slab_variable(prior)

  # prepare output objects
  if(is.prior.factor(prior_variable)){

    samples <- .as_mixed_posteriors.factor(model_samples, prior_variable, parameter)
    attr(samples, "models_ind") <- as.vector(model_samples[,paste0(parameter, "_indicator")])

  }else if(is.prior.simple(prior_variable)){

    samples <- .as_mixed_posteriors.simple(model_samples, prior_variable, parameter)
    attr(samples, "models_ind") <- as.vector(model_samples[,paste0(parameter, "_indicator")])

  }

  class(samples) <- c("mixed_posteriors.spike_and_slab", class(samples))
  attr(samples, "prior_list") <- prior
  if(!is.null(dim(samples))){
    samples <- .posterior_support_set_columns_from_prior_list(samples, prior)
  }else{
    samples <- .posterior_support_set_from_prior_list(samples, prior)
  }
  samples <- .posterior_atoms_set(
    samples,
    .posterior_atoms_from_indicator(
      prior = prior,
      indicator = attr(samples, "models_ind"),
      n_columns = if(is.null(dim(samples))) 1L else ncol(samples),
      column_names = if(is.null(dim(samples))) parameter else colnames(samples),
      spike_and_slab = TRUE
    )
  )

  return(samples)
}
.as_mixed_posteriors.mixture        <- function(model_samples, prior, parameter, conditional){

  # check input
  check_char(parameter, "parameter", check_length = FALSE)


  # prepare output objects
  if(inherits(prior, "prior.bias_mixture")){

    is_PET            <- sapply(prior, is.prior.PET)
    is_PEESE          <- sapply(prior, is.prior.PEESE)
    branch_info       <- .selection_prior_branch_info(prior)
    has_selection     <- vapply(branch_info, function(x) !is.null(x$selection), logical(1))
    has_phacking      <- vapply(branch_info, function(x) !is.null(x$phacking), logical(1))

    # prepare weightfunction parameter names
    if(any(has_selection)){
      selection_priors <- lapply(branch_info[has_selection], function(x) x$selection)
      omega_info    <- .weightfunction_mapping_info(selection_priors, one_sided = TRUE)
      omega_mapping <- omega_info$mapping
      omega_names   <- omega_info$names
      omega_par     <- omega_info$pars
    }
    if(any(has_phacking)){
      phacking_priors <- lapply(branch_info[has_phacking], function(x) x$phacking)
      phacking_par    <- .selection_phacking_report_parameters(phacking_priors)
      phacking_names  <- phacking_par
    }

    # deal with conditional parameters
    if(length(conditional) > 0 && any(c("PET", "PEESE", "PETPEESE", "omega", "phacking", "alpha", "pi_null") %in% conditional)){

      out_names <- NULL
      par_names <- NULL

      if("omega" %in% conditional && any(has_selection)){
        out_names <- c(out_names, omega_names)
        par_names <- c(par_names, omega_par)
      }
      if("phacking" %in% conditional && any(has_phacking)){
        out_names <- c(out_names, phacking_names)
        par_names <- c(par_names, phacking_par)
      }
      if("alpha" %in% conditional && any(has_phacking)){
        out_names <- c(out_names, "alpha")
        par_names <- c(par_names, "alpha")
      }
      if("pi_null" %in% conditional && any(has_phacking)){
        out_names <- c(out_names, "pi_null")
        par_names <- c(par_names, "pi_null")
      }
      if("PETPEESE" %in% conditional){
        # subset in case only PET/PEESE is supplied
        out_names <- c(out_names, colnames(model_samples)[colnames(model_samples) %in% c("PET", "PEESE")])
        par_names <- c(par_names, colnames(model_samples)[colnames(model_samples) %in% c("PET", "PEESE")])
      }
      if("PET" %in% conditional && any(is_PET)){
        out_names <- c(out_names, "PET")
        par_names <- c(par_names, "PET")
      }
      if("PEESE" %in% conditional && any(is_PEESE)){
        out_names <- c(out_names, "PEESE")
        par_names <- c(par_names, "PEESE")
      }

    }else{

      out_names <- NULL
      par_names <- NULL

      if(any(has_selection)){
        out_names <- c(out_names, omega_names)
        par_names <- c(par_names, omega_par)
      }
      if(any(has_phacking)){
        out_names <- c(out_names, phacking_names)
        par_names <- c(par_names, phacking_par)
      }
      if(any(is_PET)){
        out_names <- c(out_names, "PET")
        par_names <- c(par_names, "PET")
      }
      if(any(is_PEESE)){
        out_names <- c(out_names, "PEESE")
        par_names <- c(par_names, "PEESE")
      }
    }

    # select samples
    if(is.null(par_names)){
      par_names <- character()
      out_names <- character()
    }
    keep_unique <- !duplicated(par_names)
    par_names <- par_names[keep_unique]
    out_names <- out_names[keep_unique]
    keep_par <- par_names %in% colnames(model_samples)
    par_names <- par_names[keep_par]
    out_names <- out_names[keep_par]
    samples   <- model_samples[, par_names,drop=FALSE]
    indicator <- model_samples[,paste0(parameter, "_indicator")]

    rownames(samples) <- NULL
    colnames(samples) <- out_names
    attr(samples, "sample_ind") <- FALSE
    attr(samples, "models_ind") <- as.vector(indicator)
    attr(samples, "parameter")  <- parameter
    attr(samples, "prior_list") <- prior
    if(any(has_selection)){
      samples <- .weightfunction_set_omega_context(samples, omega_info)
    }
    class(samples) <- c("mixed_posteriors", "mixed_posteriors.bias")

  }else{

    if(inherits(prior, "prior.simple_mixture")){
      samples <- .as_mixed_posteriors.simple(model_samples, prior, parameter)
    }else if(inherits(prior, "prior.factor_mixture")){
      samples <- .as_mixed_posteriors.factor(model_samples, prior, parameter)
    }
    attr(samples, "models_ind") <- as.vector(model_samples[,paste0(parameter, "_indicator")])

  }

  class(samples) <- c("mixed_posteriors.mixture", class(samples))
  attr(samples, "prior_list") <- prior
  samples <- .posterior_atoms_set(
    samples,
    .posterior_atoms_from_indicator(
      prior = prior,
      indicator = attr(samples, "models_ind"),
      n_columns = if(is.null(dim(samples))) 1L else ncol(samples),
      column_names = if(is.null(dim(samples))) parameter else colnames(samples)
    )
  )

  return(samples)
}
