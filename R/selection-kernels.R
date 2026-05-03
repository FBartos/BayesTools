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

  check_char(side, "side")
  check_real(target, "target", lower = 0, upper = 1, allow_bound = FALSE, allow_NA = FALSE)
  check_real(source, "source", lower = 0, upper = 1, allow_bound = FALSE, allow_NA = FALSE)
  check_real(destination, "destination", lower = 0, upper = 1, allow_bound = FALSE, allow_NA = FALSE)
  form <- match.arg(form)
  check_char(report_scale, "report_scale", allow_values = c("pi_null", "alpha"))
  check_real(prior_weights, "prior_weights", lower = 0, allow_bound = FALSE)

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

  check_real(prior_weights, "prior_weights", lower = 0, allow_bound = FALSE)

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
    monitor <- c(monitor, names$alpha, names$phack_kind, names$pi_null)
  }

  phacking_priors <- lapply(branch_info[has_phacking], function(x) x$phacking)
  phacking <- .selection_backend_phacking_info(phacking_priors, names)

  init <- .selection_backend_init(branch_info, breaks, prior_weights, names, uses_indicator)

  return(list(
    mode           = mode,
    branch_type    = branch_type,
    prior_weights  = prior_weights,
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

.phack_validate_alpha_prior <- function(alpha){

  .check_prior(alpha)

  if(!is.prior.simple(alpha)){
    stop("'alpha' must be a simple prior distribution.", call. = FALSE)
  }
  if(is.prior.discrete(alpha) && !is.prior.point(alpha)){
    stop("'alpha' must be a continuous prior distribution or a point prior.", call. = FALSE)
  }

  if(is.prior.point(alpha)){
    if(alpha$parameters[["location"]] < 0 || alpha$parameters[["location"]] >= 1){
      stop("Point p-hacking alpha priors must be in the interval [0, 1).", call. = FALSE)
    }
  }else{
    if(alpha$truncation[["lower"]] < 0 || alpha$truncation[["upper"]] > 1){
      stop("P-hacking alpha priors must have support within [0, 1].", call. = FALSE)
    }
  }

  return()
}

.phack_kind <- function(form){
  switch(
    form,
    "linear"    = 1L,
    "quadratic" = 2L
  )
}

.phack_power_null_moment <- function(lower, upper, q, anchor, reverse){

  width <- upper - lower
  m0 <- stats::pnorm(upper) - stats::pnorm(lower)
  m1 <- stats::dnorm(lower) - stats::dnorm(upper)

  if(q == 1L){
    if(reverse){
      return((anchor * m0 - m1) / width)
    }else{
      return((m1 - anchor * m0) / width)
    }
  }

  m2 <- m0 + lower * stats::dnorm(lower) - upper * stats::dnorm(upper)
  if(reverse){
    return((m2 - 2 * anchor * m1 + anchor^2 * m0) / width^2)
  }

  return((m2 - 2 * anchor * m1 + anchor^2 * m0) / width^2)
}

.selection_backend_names <- function(names){

  backend_names <- names
  if(!is.null(backend_names$phack_z_destination) && is.null(backend_names$phack_z_dest)){
    backend_names$phack_z_dest <- backend_names$phack_z_destination
  }
  backend_names$phack_z_destination <- NULL

  defaults <- list(
    omega               = "omega",
    alpha               = "alpha",
    pi_null             = "pi_null",
    beta_null           = "beta_null",
    phack_kind          = "phack_kind",
    phack_z_source      = "phack_z_source",
    phack_z_dest        = "phack_z_dest"
  )
  defaults[base::names(backend_names)] <- backend_names

  for(i in seq_along(defaults)){
    check_char(defaults[[i]], paste0("names$", names(defaults)[i]))
  }

  return(defaults)
}

.selection_normalize_priors <- function(priors){

  if(is.null(priors)){
    priors <- list(prior_none())
  }else if(is.prior(priors) && is.prior.mixture(priors)){
    priors <- as.list(priors)
  }else if(is.prior(priors)){
    priors <- list(priors)
  }else if(is.list(priors)){
    if(length(priors) == 0L){
      priors <- list(prior_none())
    }
  }else{
    stop("'priors' must be a prior, prior mixture, list of priors, or NULL.", call. = FALSE)
  }

  if(!all(vapply(priors, is.prior, logical(1)))){
    stop("'priors' must contain only prior objects.", call. = FALSE)
  }

  allowed <- vapply(priors, function(x){
    is.prior.none(x) || is.prior.PET(x) || is.prior.PEESE(x) ||
      is.prior.weightfunction(x) || is_prior_phacking(x) || is_prior_bias(x)
  }, logical(1))

  if(!all(allowed)){
    stop("'priors' contains unsupported selection prior objects.", call. = FALSE)
  }

  return(priors)
}

.selection_branch_info <- function(prior){

  if(is.prior.none(prior)){
    return(list(type = "none", selection = NULL, phacking = NULL))
  }
  if(is.prior.PET(prior)){
    return(list(type = "PET", selection = NULL, phacking = NULL))
  }
  if(is.prior.PEESE(prior)){
    return(list(type = "PEESE", selection = NULL, phacking = NULL))
  }
  if(is.prior.weightfunction(prior)){
    return(list(type = "weightfunction", selection = prior, phacking = NULL))
  }
  if(is_prior_phacking(prior)){
    return(list(type = "phack", selection = NULL, phacking = prior))
  }
  if(is_prior_bias(prior)){
    has_selection <- !is.null(prior$selection)
    has_phacking  <- !is.null(prior$phacking)
    type <- if(has_selection && has_phacking){
      "combined"
    }else if(has_selection){
      "weightfunction"
    }else{
      "phack"
    }
    return(list(type = type, selection = prior$selection, phacking = prior$phacking))
  }

  stop("Unsupported selection prior object.", call. = FALSE)
}

.selection_backend_mode <- function(has_selection, has_phacking){

  if(has_selection && has_phacking){
    return("step_phack_power")
  }
  if(has_selection){
    return("step")
  }
  if(has_phacking){
    return("phack_power")
  }

  return("none")
}

.selection_mode_code <- function(mode){
  switch(
    mode,
    "none"             = 0L,
    "step"             = 1L,
    "phack_power"      = 2L,
    "step_phack_power" = 3L
  )
}

.selection_validate_global_breaks <- function(global_breaks){

  check_real(global_breaks, "global_breaks", check_length = 0, allow_NA = FALSE)
  if(length(global_breaks) < 2L){
    stop("'global_breaks' must contain at least 0 and 1.", call. = FALSE)
  }
  if(any(global_breaks < 0 | global_breaks > 1)){
    stop("'global_breaks' must be within [0, 1].", call. = FALSE)
  }
  if(anyDuplicated(global_breaks)){
    stop("'global_breaks' must not contain duplicate values.", call. = FALSE)
  }
  if(!all(global_breaks == cummax(global_breaks))){
    stop("'global_breaks' must be monotonically increasing.", call. = FALSE)
  }
  if(!isTRUE(all.equal(global_breaks[1], 0)) || !isTRUE(all.equal(global_breaks[length(global_breaks)], 1))){
    stop("'global_breaks' must start at 0 and end at 1.", call. = FALSE)
  }

  return(global_breaks)
}

.selection_jags_step_component_code <- function(selection, component_id, n_bins, global_cuts){

  if(is.null(selection)){
    return(.JAGS_weightfunction_none_component_syntax(component_id = component_id, n_bins = n_bins))
  }

  return(.JAGS_weightfunction_component_syntax(
    prior           = selection,
    component_id    = component_id,
    global_cuts     = global_cuts,
    force_one_sided = TRUE
  ))
}

.selection_jags_active_scalar <- function(prefix, target, indicator_terms, component_ids){

  paste0(
    target, " <- ",
    paste0(prefix, component_ids, " * ", indicator_terms, collapse = " + ")
  )
}

.selection_jags_active_vector <- function(prefix, target, index, indicator_terms, component_ids){

  paste0(
    target, "[", index, "] <- ",
    paste0(prefix, component_ids, "[", index, "] * ", indicator_terms, collapse = " + ")
  )
}

.selection_backend_phacking_info <- function(phacking_priors, names){

  if(length(phacking_priors) == 0L){
    return(list(
      form                 = "none",
      q                    = 0L,
      z_source             = c(0, 0),
      z_destination        = c(0, 0),
      coefficient          = names$alpha,
      coefficient_ids      = names$alpha,
      branch_form          = character(),
      branch_q             = integer(),
      branch_phack_kind    = integer(),
      branch_beta_null_per_alpha = numeric(),
      branch_z_source      = matrix(nrow = 0L, ncol = 2L),
      branch_z_destination = matrix(nrow = 0L, ncol = 2L)
    ))
  }

  constants <- lapply(phacking_priors, function(x){
    phack_backend_constants(x$form, x$source, x$destination, target = x$target)
  })

  forms <- vapply(constants, function(x) x$form, character(1))
  q <- vapply(constants, function(x) x$q, integer(1))
  z_source <- do.call(rbind, lapply(constants, function(x) x$z_source))
  z_destination <- do.call(rbind, lapply(constants, function(x) x$z_destination))

  return(list(
    form                 = if(length(unique(forms)) == 1L) forms[1] else unique(forms),
    q                    = if(length(unique(q)) == 1L) q[1] else unique(q),
    z_source             = if(nrow(z_source) == 1L) as.numeric(z_source[1,]) else z_source,
    z_destination        = if(nrow(z_destination) == 1L) as.numeric(z_destination[1,]) else z_destination,
    coefficient          = names$alpha,
    coefficient_ids      = names$alpha,
    branch_form          = forms,
    branch_q             = q,
    branch_phack_kind    = vapply(constants, function(x) x$phack_kind, integer(1)),
    branch_beta_null_per_alpha = vapply(constants, function(x) x$beta_null_per_alpha, numeric(1)),
    branch_z_source      = z_source,
    branch_z_destination = z_destination
  ))
}

.selection_backend_init <- function(branch_info, breaks, prior_weights, names, uses_indicator){

  active_branch <- which.max(prior_weights)

  init <- list()
  for(i in seq_along(branch_info)){
    component_id <- if(uses_indicator) i else NULL
    if(!is.null(branch_info[[i]]$selection)){
      init <- c(init, .selection_JAGS_init_weightfunction_component(branch_info[[i]]$selection, component_id = component_id))
    }
    if(!is.null(branch_info[[i]]$phacking)){
      init <- c(init, .selection_JAGS_init_phacking_component(branch_info[[i]]$phacking, component_id = component_id))
    }
  }
  if(uses_indicator){
    init[["bias_indicator"]] <- active_branch
  }

  return(init)
}

.selection_JAGS_init_weightfunction_component <- function(prior, component_id = NULL){

  init <- list()
  if(prior$weights$type == "fixed"){
    return()
  }else if(prior$weights$type == "cumulative"){
    eta_name <- if(is.null(component_id)) "eta" else paste0("eta_component_", component_id)
    init[[eta_name]] <- stats::rgamma(length(prior$weights[["alpha"]]), shape = prior$weights[["alpha"]], rate = 1)
  }

  return(init)
}

.selection_JAGS_init_phacking_component <- function(prior, component_id = NULL){

  alpha_name <- if(is.null(component_id)) "alpha" else paste0("alpha_component_", component_id)
  .JAGS_init.simple(prior$alpha, alpha_name)
}

.selection_prior_branch_info <- function(prior){

  if(is.prior.mixture(prior)){
    return(lapply(prior, .selection_branch_info))
  }

  list(.selection_branch_info(prior))
}

.selection_prior_has_selection <- function(prior){

  branch_info <- .selection_prior_branch_info(prior)
  any(vapply(branch_info, function(x) !is.null(x$selection), logical(1)))
}

.selection_prior_has_phacking <- function(prior){

  branch_info <- .selection_prior_branch_info(prior)
  any(vapply(branch_info, function(x) !is.null(x$phacking), logical(1)))
}

.selection_prior_selection_priors <- function(prior){

  branch_info <- .selection_prior_branch_info(prior)
  selection_priors <- lapply(branch_info, function(x) x$selection)
  selection_priors[!vapply(selection_priors, is.null, logical(1))]
}

.selection_prior_phacking_priors <- function(prior){

  branch_info <- .selection_prior_branch_info(prior)
  phacking_priors <- lapply(branch_info, function(x) x$phacking)
  phacking_priors[!vapply(phacking_priors, is.null, logical(1))]
}

.selection_prior_stop_unsupported_generic <- function(generic, prior){

  if(is_prior_bias(prior)){
    selection_backend_spec(prior)
    stop(
      sprintf(
        "No %s is implemented for composed bias priors; use the selection or p-hacking component explicitly.",
        generic
      ),
      call. = FALSE
    )
  }

  if(is_prior_phacking(prior)){
    selection_backend_spec(prior)
    stop(
      sprintf(
        "No %s is implemented for p-hacking priors; use the alpha prior component explicitly.",
        generic
      ),
      call. = FALSE
    )
  }

  stop(sprintf("No %s is implemented for this prior.", generic), call. = FALSE)
}

.selection_prior_has_PET <- function(prior){

  if(is.prior.mixture(prior)){
    return(any(sapply(prior, is.prior.PET)))
  }

  is.prior.PET(prior)
}

.selection_prior_has_PEESE <- function(prior){

  if(is.prior.mixture(prior)){
    return(any(sapply(prior, is.prior.PEESE)))
  }

  is.prior.PEESE(prior)
}

.selection_bias_parameter_names <- function(prior, include_kind = TRUE){

  c(
    if(.selection_prior_has_PET(prior)) "PET",
    if(.selection_prior_has_PEESE(prior)) "PEESE",
    if(.selection_prior_has_selection(prior)) "omega",
    if(.selection_prior_has_phacking(prior)) c(
      .selection_prior_phacking_report_parameters(prior),
      if(include_kind) c("beta_null", "phack_kind")
    )
  )
}

.selection_initial_omega <- function(selection, breaks){

  local <- .selection_weightfunction_mean(selection)
  expansion <- .weightfunction_mapping_expansion(selection, force_one_sided = TRUE)
  global_bin_indices <- .weightfunction_global_bin_indices(breaks, expansion)

  vapply(seq_len(length(breaks) - 1L), function(i){
    ind <- global_bin_indices[i]
    local[expansion$index[ind]]
  }, numeric(1))
}

.selection_weightfunction_mean <- function(selection){

  components <- .weightfunction_marginal_components(selection)
  vapply(components, .prior_weightfunction_component_mean, numeric(1))
}

.selection_prior_mean <- function(prior){

  if(is.prior.point(prior)){
    return(prior$parameters[["location"]])
  }

  return(mean(prior))
}

.selection_format_number <- function(x){
  format(x, scientific = FALSE, digits = 16, trim = TRUE)
}

.JAGS_phacking_component_syntax <- function(prior, component_id = NULL){

  alpha_name <- if(is.null(component_id)) "alpha" else paste0("alpha_component_", component_id)
  kind_name <- if(is.null(component_id)) "phack_kind" else paste0("phack_kind_component_", component_id)
  pi_null_name <- if(is.null(component_id)) "pi_null" else paste0("pi_null_component_", component_id)
  beta_null_name <- if(is.null(component_id)) "beta_null" else paste0("beta_null_component_", component_id)
  z_source_name <- if(is.null(component_id)) "phack_z_source" else paste0("phack_z_source_component_", component_id)
  z_destination_name <- if(is.null(component_id)) "phack_z_dest" else paste0("phack_z_dest_component_", component_id)

  if(is.null(prior)){
    return(paste0(
      alpha_name, " <- 0\n",
      kind_name, " <- 0\n",
      pi_null_name, " <- 0\n",
      beta_null_name, " <- 0\n",
      z_source_name, "[1] <- 0\n",
      z_source_name, "[2] <- 0\n",
      z_destination_name, "[1] <- 0\n",
      z_destination_name, "[2] <- 0\n"
    ))
  }

  if(!is_prior_phacking(prior)){
    stop("'prior' must be a p-hacking prior.", call. = FALSE)
  }

  constants <- phack_backend_constants(prior$form, prior$source, prior$destination, target = prior$target)

  paste0(
    .JAGS_prior.simple(prior$alpha, alpha_name),
    kind_name, " <- ", constants$phack_kind, "\n",
    pi_null_name, " <- ", alpha_name, " * ", .selection_format_number(constants$pi_null_per_alpha), "\n",
    beta_null_name, " <- ", alpha_name, " * ", .selection_format_number(constants$beta_null_per_alpha), "\n",
    z_source_name, "[1] <- ", .selection_format_number(constants$z_source[1]), "\n",
    z_source_name, "[2] <- ", .selection_format_number(constants$z_source[2]), "\n",
    z_destination_name, "[1] <- ", .selection_format_number(constants$z_destination[1]), "\n",
    z_destination_name, "[2] <- ", .selection_format_number(constants$z_destination[2]), "\n"
  )
}
