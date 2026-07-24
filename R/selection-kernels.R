#' @title Creates a prior distribution for a p-hacking selection kernel
#'
#' @description \code{prior_phacking()} creates a prior distribution for a
#' mass-preserving p-hacking redistribution kernel. The initial backend supports
#' linear and quadratic power depletion between \code{source} and
#' \code{target}, with depleted probability mass redistributed between
#' \code{target} and \code{destination}.
#'
#' @param side side geometry. Currently only \code{"one-sided"} is supported.
#' @param target target p-value cut point.
#' @param source source p-value cut point. Must be larger than \code{target}.
#' @param destination destination p-value cut point. Must be smaller than
#' \code{target}.
#' @param form power depletion form, either \code{"linear"} or
#' \code{"quadratic"}.
#' @param alpha prior distribution for the p-hacking severity parameter.
#' @param report_scale reporting scale for deterministic summaries.
#' @param prior_weights prior odds associated with a given distribution.
#'
#' @return \code{prior_phacking()} returns an object of class \code{"prior"}.
#'
#' @export
prior_phacking <- function(side = "one-sided",
                           target = .025,
                           source = .25,
                           destination = .005,
                           form = c("linear", "quadratic"),
                           alpha = prior("beta", list(1, 1)),
                           report_scale = "pi_null",
                           prior_weights = 1){

  check_char(side, "side", allow_NA = FALSE)
  check_real(target, "target", lower = 0, upper = 1, allow_bound = FALSE, allow_NA = FALSE)
  check_real(source, "source", lower = 0, upper = 1, allow_bound = FALSE, allow_NA = FALSE)
  check_real(destination, "destination", lower = 0, upper = 1, allow_bound = FALSE, allow_NA = FALSE)
  form <- match.arg(form)
  check_char(report_scale, "report_scale", allow_values = c("pi_null", "alpha"), allow_NA = FALSE)
  .check_prior_weight(prior_weights)

  side <- .weightfunction_normalize_side(side)
  if(side != "one-sided"){
    stop("P-hacking priors currently support only one-sided p-value geometry.", call. = FALSE)
  }

  if(!(destination < target && target < source)){
    stop("P-hacking p-value cut points must satisfy destination < target < source.", call. = FALSE)
  }

  .phack_validate_alpha_prior(alpha)
  constants <- phack_backend_constants(form = form, source = source, destination = destination, target = target)

  output <- list(
    distribution  = "phacking",
    side          = side,
    target        = target,
    source        = source,
    destination   = destination,
    form          = form,
    q             = constants$q,
    alpha         = alpha,
    report_scale  = report_scale,
    parameters    = list(
      target      = target,
      source      = source,
      destination = destination,
      form        = form
    ),
    truncation    = list(lower = 0, upper = 1),
    prior_weights = prior_weights
  )

  class(output) <- c("prior", "prior.phacking")

  return(output)
}

#' @title Creates a composed publication-bias prior
#'
#' @description \code{prior_bias()} composes a step-selection
#' \code{prior_weightfunction()} and/or a \code{prior_phacking()} object. It does
#' not introduce new math; the backend compiler treats the composition as the
#' product of both kernels.
#'
#' @param selection optional \code{prior_weightfunction()} object.
#' @param phacking optional \code{prior_phacking()} object.
#' @inheritParams prior_phacking
#'
#' @return \code{prior_bias()} returns an object of class \code{"prior"}.
#'
#' @export
prior_bias <- function(selection = NULL, phacking = NULL, prior_weights = 1){

  .check_prior_weight(prior_weights)

  if(is.null(selection) && is.null(phacking)){
    stop("At least one of 'selection' or 'phacking' must be specified.", call. = FALSE)
  }
  if(!is.null(selection) && !is.prior.weightfunction(selection)){
    stop("'selection' must be a weightfunction prior created by prior_weightfunction().", call. = FALSE)
  }
  if(!is.null(phacking) && !is_prior_phacking(phacking)){
    stop("'phacking' must be a p-hacking prior created by prior_phacking().", call. = FALSE)
  }
  if(!is.null(selection) && !is.null(phacking) && selection$side != phacking$side){
    stop("'selection' and 'phacking' must use the same side geometry.", call. = FALSE)
  }

  output <- list(
    distribution  = "bias",
    selection     = selection,
    phacking      = phacking,
    parameters    = list(),
    truncation    = list(lower = 0, upper = Inf),
    prior_weights = prior_weights
  )

  class(output) <- c("prior", "prior.bias")

  return(output)
}

#' @title Reports whether x is a p-hacking prior
#'
#' @param x object to test.
#'
#' @return A logical value.
#'
#' @export
is_prior_phacking <- function(x){
  inherits(x, "prior.phacking")
}

#' @title Reports whether x is a composed publication-bias prior
#'
#' @inheritParams is_prior_phacking
#'
#' @return A logical value.
#'
#' @export
is_prior_bias <- function(x){
  inherits(x, "prior.bias")
}

is.prior.phacking <- function(x){
  is_prior_phacking(x)
}

is.prior.bias <- function(x){
  is_prior_bias(x)
}

.phacking_report_parameter <- function(prior){

  if(!is_prior_phacking(prior)){
    stop("'prior' must be a p-hacking prior.", call. = FALSE)
  }

  if(is.null(prior$report_scale)){
    return("pi_null")
  }

  prior$report_scale
}

.phacking_unreported_parameters <- function(prior){

  setdiff(c("alpha", "pi_null"), .phacking_report_parameter(prior))
}

.selection_phacking_report_parameters <- function(phacking_priors){

  if(length(phacking_priors) == 0L){
    return(character())
  }

  unique(vapply(phacking_priors, .phacking_report_parameter, character(1)))
}

.selection_phacking_unreported_parameters <- function(phacking_priors){

  setdiff(c("alpha", "pi_null"), .selection_phacking_report_parameters(phacking_priors))
}

.selection_prior_phacking_report_parameters <- function(prior){

  .selection_phacking_report_parameters(.selection_prior_phacking_priors(prior))
}

#' @title P-hacking calibration helpers
#'
#' @description \code{phack_pi_null()} converts the sampled severity
#' \code{alpha} to the amount of null probability mass depleted from the source
#' interval. \code{phack_alpha_from_pi_null()} applies the inverse calibration.
#' \code{phack_backend_constants()} returns the z-scale constants used by the
#' backend.
#'
#' @param alpha p-hacking severity parameter.
#' @param pi_null null probability mass depleted from the source interval.
#' @param form power depletion form, either \code{"linear"} or
#' \code{"quadratic"}.
#' @param source source p-value cut point.
#' @param destination destination p-value cut point.
#' @param target target p-value cut point.
#'
#' @return Numeric vector for the calibration helpers and a named list for
#' \code{phack_backend_constants()}.
#'
#' @export
phack_pi_null <- function(alpha, form, source, destination, target = .025){

  check_real(alpha, "alpha", lower = 0, upper = 1, check_length = 0, allow_NA = FALSE)
  if(any(alpha >= 1)){
    stop("'alpha' must be lower than 1.", call. = FALSE)
  }
  constants <- phack_backend_constants(form = form, source = source, destination = destination, target = target)

  return(alpha * constants$pi_null_per_alpha)
}

#' @rdname phack_pi_null
#' @export
phack_alpha_from_pi_null <- function(pi_null, form, source, destination, target = .025){

  check_real(pi_null, "pi_null", lower = 0, check_length = 0, allow_NA = FALSE)
  constants <- phack_backend_constants(form = form, source = source, destination = destination, target = target)

  if(any(pi_null >= constants$pi_null_per_alpha)){
    stop("'pi_null' is too large for alpha <= 1 under the specified p-hacking geometry.", call. = FALSE)
  }

  return(pi_null / constants$pi_null_per_alpha)
}

#' @rdname phack_pi_null
#' @export
phack_backend_constants <- function(form, source, destination, target = .025){

  form <- match.arg(form, choices = c("linear", "quadratic"))
  check_real(target, "target", lower = 0, upper = 1, allow_bound = FALSE, allow_NA = FALSE)
  check_real(source, "source", lower = 0, upper = 1, allow_bound = FALSE, allow_NA = FALSE)
  check_real(destination, "destination", lower = 0, upper = 1, allow_bound = FALSE, allow_NA = FALSE)

  if(!(destination < target && target < source)){
    stop("P-hacking p-value cut points must satisfy destination < target < source.", call. = FALSE)
  }

  q <- switch(
    form,
    "linear"    = 1L,
    "quadratic" = 2L
  )

  z_source_lower <- stats::qnorm(1 - source)
  z_target       <- stats::qnorm(1 - target)
  z_dest_upper   <- stats::qnorm(1 - destination)

  source_mass <- .phack_power_null_moment(z_source_lower, z_target, q, anchor = z_source_lower, reverse = FALSE)
  dest_mass   <- .phack_power_null_moment(z_target, z_dest_upper, q, anchor = z_dest_upper, reverse = TRUE)

  return(list(
    form                    = form,
    q                       = q,
    phack_kind              = .phack_kind(form),
    target                  = target,
    source                  = source,
    destination             = destination,
    z_source                = c(z_source_lower, z_target),
    z_destination           = c(z_target, z_dest_upper),
    source_null_mass        = source_mass,
    destination_null_mass   = dest_mass,
    pi_null_per_alpha       = source_mass,
    beta_null_per_alpha     = source_mass / dest_mass
  ))
}

#' @title Compile selection priors for backend consumers
#'
#' @description \code{selection_backend_spec()} compiles step-selection and
#' p-hacking prior objects into active backend parameters. The returned object
#' contains stable p/z geometry, JAGS prior/transform code, monitor names,
#' initial values, and data constants.
#'
#' @param priors a selection prior, p-hacking prior, composed bias prior,
#' \code{prior_none()}, \code{prior_mixture()}, or a list of those priors.
#' @param backend backend target. Currently only \code{"jags"} is supported.
#' @param names list of backend parameter names.
#' @param global_breaks optional global p-value break grid.
#'
#' @return A list describing the compiled backend specification.
#'
#' @export
selection_backend_spec <- function(priors,
                                   backend = "jags",
                                   names = list(omega = "omega", alpha = "alpha"),
                                   global_breaks = NULL){

  check_char(backend, "backend", allow_values = "jags")
  check_list(names, "names", check_names = c("omega", "alpha", "pi_null", "beta_null", "phack_kind", "phack_z_source", "phack_z_dest", "phack_z_destination"), allow_other = FALSE)
  names <- .selection_backend_names(names)

  branches <- .selection_normalize_priors(priors)
  branch_info <- lapply(branches, .selection_branch_info)

  has_selection <- vapply(branch_info, function(x) !is.null(x$selection), logical(1))
  has_phacking  <- vapply(branch_info, function(x) !is.null(x$phacking),  logical(1))
  branch_type   <- vapply(branch_info, function(x) x$type, character(1))

  mode <- .selection_backend_mode(any(has_selection), any(has_phacking))

  step_priors <- lapply(branch_info[has_selection], function(x) x$selection)
  if(is.null(global_breaks)){
    breaks <- if(length(step_priors) > 0L){
      weightfunctions_mapping(step_priors, cuts_only = TRUE, one_sided = TRUE)
    }else{
      c(0, 1)
    }
  }else{
    breaks <- .selection_validate_global_breaks(global_breaks)
    if(length(step_priors) > 0L){
      required_breaks <- weightfunctions_mapping(step_priors, cuts_only = TRUE, one_sided = TRUE)
      if(!all(vapply(required_breaks, function(x) any(abs(x - breaks) < sqrt(.Machine$double.eps)), logical(1)))){
        stop("'global_breaks' must contain all step-selection p-value breaks.", call. = FALSE)
      }
    }
  }
  n_bins <- length(breaks) - 1L

  prior_weights <- vapply(branches, function(x) x$prior_weights, numeric(1))
  uses_indicator <- length(branches) > 1L
  indicator_terms <- if(uses_indicator){
    paste0("equals(bias_indicator, ", seq_along(branches), ")")
  }else{
    rep("1", length(branches))
  }

  prior_code <- character()
  transform_code <- character()

  if(uses_indicator){
    prior_code <- c(prior_code, paste0("bias_indicator ~ dcat(c(", paste0(prior_weights, collapse = ", "), "))"))
  }

  for(i in seq_along(branches)){
    component_id <- if(uses_indicator) i else NULL
    step_code <- if(uses_indicator || !is.null(branch_info[[i]]$selection) || !is.null(branch_info[[i]]$phacking)){
      .selection_jags_step_component_code(branch_info[[i]]$selection, component_id = component_id, n_bins = n_bins, global_cuts = breaks)
    }else{
      character()
    }
    phacking_code <- if(uses_indicator || !is.null(branch_info[[i]]$phacking)){
      .JAGS_phacking_component_syntax(branch_info[[i]]$phacking, component_id = component_id)
    }else{
      character()
    }
    prior_code <- c(
      prior_code,
      step_code,
      phacking_code
    )
  }

  if(uses_indicator){
    for(j in seq_len(n_bins)){
      transform_code <- c(
        transform_code,
        paste0(
          names$omega, "[", j, "] <- ",
          paste0("omega_component_", seq_along(branches), "[", j, "] * ", indicator_terms, collapse = " + ")
        )
      )
    }

    transform_code <- c(
      transform_code,
      .selection_jags_active_scalar("alpha_component_", names$alpha, indicator_terms, seq_along(branches)),
      .selection_jags_active_scalar("phack_kind_component_", names$phack_kind, indicator_terms, seq_along(branches)),
      .selection_jags_active_scalar("pi_null_component_", names$pi_null, indicator_terms, seq_along(branches)),
      .selection_jags_active_scalar("beta_null_component_", names$beta_null, indicator_terms, seq_along(branches))
    )

    if(any(has_phacking)){
      for(k in 1:2){
        transform_code <- c(
          transform_code,
          .selection_jags_active_vector("phack_z_source_component_", names$phack_z_source, k, indicator_terms, seq_along(branches)),
          .selection_jags_active_vector("phack_z_dest_component_", names$phack_z_dest, k, indicator_terms, seq_along(branches))
        )
      }
    }
  }

  monitor <- character()
  if(uses_indicator){
    monitor <- c(monitor, "bias_indicator")
  }
  if(any(has_selection) || any(has_phacking)){
    monitor <- c(monitor, names$omega)
  }
  if(!uses_indicator && any(has_selection)){
    monitor <- c(monitor, .JAGS_monitor_private.weightfunction(branch_info[[which(has_selection)[1L]]]$selection))
  }
  if(any(has_phacking)){
    monitor <- c(
      monitor,
      names$alpha,
      .selection_backend_phacking_auxiliary_monitors(branch_info, has_phacking, names, uses_indicator),
      names$phack_kind,
      names$pi_null
    )
  }

  phacking_priors <- lapply(branch_info[has_phacking], function(x) x$phacking)
  phacking <- .selection_backend_phacking_info(phacking_priors, names)

  init <- .selection_backend_init(branch_info, breaks, prior_weights, names, uses_indicator)

  return(list(
    mode           = mode,
    branch_type    = branch_type,
    prior_weights  = prior_weights,
    jags_omega     = names$omega,
    jags_alpha     = names$alpha,
    jags_pi_null   = names$pi_null,
    jags_beta_null = names$beta_null,
    jags_phack_kind = names$phack_kind,
    jags_phack_z_source = names$phack_z_source,
    jags_phack_z_dest   = names$phack_z_dest,
    step           = list(
      breaks          = breaks,
      n_bins          = as.integer(n_bins),
      coefficient     = names$omega,
      coefficient_ids = paste0(names$omega, "[", seq_len(n_bins), "]"),
      z_lower         = stats::qnorm(1 - breaks[-1]),
      z_upper         = stats::qnorm(1 - breaks[-length(breaks)])
    ),
    phacking      = phacking,
    prior_code    = paste0(prior_code[nzchar(prior_code)], collapse = "\n"),
    transform_code = paste0(transform_code[nzchar(transform_code)], collapse = "\n"),
    monitor       = unique(monitor),
    init          = init,
    data          = list(
      sel_p_cuts       = breaks,
      sel_z_lower      = stats::qnorm(1 - breaks[-1]),
      sel_z_upper      = stats::qnorm(1 - breaks[-length(breaks)]),
      sel_n_bins       = as.integer(n_bins),
      phack_component_z_source = phacking$branch_z_source,
      phack_component_z_dest   = phacking$branch_z_destination,
      phack_component_q        = phacking$branch_q,
      phack_component_beta_null_per_alpha = phacking$branch_beta_null_per_alpha,
      kernel_mode      = .selection_mode_code(mode)
    )
  ))
}
