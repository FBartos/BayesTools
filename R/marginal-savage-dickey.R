#' @title Compute Savage-Dickey inclusion Bayes factors
#'
#' @description Computes Savage-Dickey (density ratio) inclusion Bayes factors
#' based the change of height from prior to posterior distribution at the test value.
#'
#' @param posterior marginal posterior distribution generated via the
#' \code{marginal_posterior} function
#' @param null_hypothesis point null hypothesis to test. Defaults to \code{0}
#' @param normal_approximation whether the height of prior and posterior density should be
#' approximated via a normal distribution (rather than kernel density). Defaults to \code{FALSE}.
#' @param silent whether warnings should be returned silently. Defaults to \code{FALSE}
#' @param density_method density source for the posterior ordinate. \code{"KDE"}
#' computes the Gaussian kernel density estimate at the null exactly (the
#' kernel sum over the continuous draws, with bandwidth
#' \code{stats::bw.nrd0()} of those draws; no evaluation grid or
#' interpolation), using boundary reflection when exact posterior-support
#' metadata is available. Finite sample ranges are not treated as exact
#' support; Gaussian kernel tails are evaluated at finite null values.
#' \code{"precomputed"} requires valid \code{posterior_ordinate} metadata
#' when present, or otherwise valid \code{posterior_density} metadata (see
#' [posterior_metadata()]).
#' Posterior atom status must be declared through package-generated
#' \code{atoms} metadata or \code{posterior_atom_attribute()}. An
#' ordinary Savage-Dickey density ratio is rejected when atom status is unknown
#' or when any positive atom is located exactly at the null. For a scalar
#' marginal posterior, a positive atom at the null is an error; for a list of
#' marginal posteriors (levels), such a level returns \code{NA} with the reason
#' in its \code{"warnings"} attribute ("fixed at the null hypothesis value" for
#' an atom of mass one, e.g., the reference level of a treatment-coded factor)
#' and the other levels are computed. Declared atoms at
#' other locations are removed from the continuous posterior ordinate and the
#' resulting density is scaled by the remaining continuous mass.
#'
#' @details Marginal posterior vectors may carry \code{posterior_ordinate}
#' metadata with exact \code{value} and \code{ordinate} entries. When
#' \code{density_method = "precomputed"} and
#' \code{normal_approximation = FALSE}, a matching ordinate is used for the
#' Savage-Dickey ratio. If no matching ordinate is available, a valid
#' \code{posterior_density} grid with \code{x} and \code{y} coordinates is
#' used. If neither valid precomputed source is available, an error is thrown.
#' Exact support metadata is also checked before accepting precomputed posterior
#' ordinates or densities; support is used to exclude the null only after it
#' is validated against the posterior samples. A stale precomputed value is
#' ignored when compatible exact support excludes the null; returned Bayes
#' factors label this path as exact support exclusion rather than as a KDE
#' estimate. Support is not inferred from prior-density grids,
#' posterior-density grids, or plotting ranges, which may be finite numerical
#' integration ranges rather than true support boundaries. Exact support
#' metadata is created with [posterior_support_attribute()] (\code{exact =
#' TRUE}); \code{type = "points"} support is not treated as a continuous
#' interval for KDE boundary reflection. The same object may be supplied as
#' \code{support} of \code{posterior_density} metadata. A prior density that is zero or
#' infinite at the null makes the density ratio a 0/0 or singular limit that a
#' kernel estimate of the posterior ordinate cannot estimate (e.g., with a
#' log-singular prior the posterior carries the same singularity and the true
#' ratio is a finite limit): a scalar call stops (see Value), and for a list of
#' marginal posteriors (levels) such a level returns \code{NA} with the reason
#' in its \code{"warnings"} attribute ("The prior density at the null
#' hypothesis value is zero" or "infinite") and the other levels are computed.
#'
#' Marginal posteriors created by \code{marginal_posterior()} record the
#' mixture component of each draw and each component's exact support: the
#' model for \code{mix_posteriors()} ensembles, and, for single fits with
#' mixture or spike-and-slab priors (\code{as_mixed_posteriors()}), the
#' combination of the component indicators of the terms entering the parameter
#' or level. When continuous components have different exact supports (for
#' example, a truncated prior in only some models or mixture components), the
#' KDE posterior ordinate is estimated per component and mixed by the
#' components' shares of the continuous draws: each component's ordinate uses
#' its own support (zero when the support excludes the null, the one-sided
#' limit on a support bound, as for the prior ordinate). The
#' \code{"posterior_density_components"} attribute of the Bayes factor lists the
#' components. When all continuous components share their support, the pooled
#' KDE is used.
#'
#' When the null hypothesis lies outside the continuous posterior draws (and is
#' not an exact support bound), the KDE posterior density at the null is an
#' extrapolation from Gaussian kernel tails. The finite Bayes factor is then
#' returned with a warning that it is not reliable evidence; for a model
#' mixture, only the draws of models whose support contains the null count.
#' Diagnostic messages are stored in the \code{"warnings"} attribute of each
#' Bayes factor and, unless \code{silent = TRUE}, emitted once per parameter
#' or level, prefixed with its label (for example \code{mu[A]}).
#'
#' @return \code{Savage_Dickey_BF} returns a Bayes factor. The prior ordinate
#' at the null follows the exactness rule of [hypothesis_BF()] point
#' hypotheses, with the same condition classes: a prior point mass at the null
#' value stops with an error of class \code{BayesTools_point_mass_at_null}, a
#' zero, infinite or undefined ordinate with \code{BayesTools_zero_ordinate},
#' \code{BayesTools_infinite_ordinate} or \code{BayesTools_undefined_ordinate},
#' and an ordinate without an exact structural value (a prior-density
#' combination evaluated only on a numerical grid, a quadrature rejected by
#' its diagnostics, or a density grid without recorded provenance) with
#' \code{BayesTools_inexact_ordinate}; each also has class
#' \code{BayesTools_hypothesis_ordinate}. A declared posterior point mass at the
#' null value of a scalar marginal posterior stops with an error of class
#' \code{BayesTools_posterior_point_mass_at_null} (also
#' \code{BayesTools_hypothesis_ordinate}). For a list of marginal posteriors
#' (and in marginal inference), a zero or infinite prior ordinate or a
#' posterior point mass at the null gives the level an \code{NA} Bayes factor
#' with its reason instead (see Details).
#'
#' @export
Savage_Dickey_BF <- function(posterior, null_hypothesis = 0, normal_approximation = FALSE, silent = FALSE,
                             density_method = c("KDE", "precomputed")){

  .Savage_Dickey_BF.checked(
    posterior            = posterior,
    null_hypothesis      = null_hypothesis,
    normal_approximation = normal_approximation,
    silent               = silent,
    density_method       = density_method,
    null_mass_NA         = is.list(posterior)
  )
}

# 'null_mass_NA': a level with a declared posterior point mass at the null, or
# with a zero or infinite prior ordinate there, gets an NA Bayes factor with
# its reason instead of stopping the whole evaluation (list posteriors and
# marginal inference; a scalar call keeps the error).
.Savage_Dickey_BF.checked <- function(posterior, null_hypothesis, normal_approximation,
                                      silent, density_method, null_mass_NA){

  if(!inherits(posterior, "marginal_posterior")){
    if(is.numeric(posterior)){
      .bt_draws_stop_plain(
        "'Savage_Dickey_BF' requires an object of class 'marginal_posterior', not plain numeric draws"
      )
    }
    stop("'Savage_Dickey_BF' requires an object of class 'marginal_posterior'.")
  }
  check_real(null_hypothesis, "null_hypothesis", allow_NA = FALSE)
  if(!is.finite(null_hypothesis)){
    stop("The 'null_hypothesis' argument must be finite.", call. = FALSE)
  }
  check_bool(normal_approximation, "normal_approximation", allow_NA = FALSE)
  check_bool(silent, "silent", allow_NA = FALSE)
  density_method <- .posterior_density_method(density_method)

  .Savage_Dickey_BF.marginal(
    posterior            = posterior,
    null_hypothesis      = null_hypothesis,
    normal_approximation = normal_approximation,
    silent               = silent,
    density_method       = density_method,
    null_mass_NA         = null_mass_NA
  )
}

# Savage-Dickey Bayes factors of a scalar marginal posterior or of each level of
# a list of marginal posteriors.
.Savage_Dickey_BF.marginal <- function(posterior, null_hypothesis, normal_approximation,
                                       silent, density_method, null_mass_NA = FALSE){

  if(is.list(posterior)){
    bf <- list()
    for(i in seq_along(posterior)){
      posterior_i <- .posterior_precomputed_child(
        parent          = posterior,
        child           = posterior[[i]],
        index           = i,
        null_hypothesis = null_hypothesis,
        density_method  = density_method
      )
      bf[[i]] <- .Savage_Dickey_BF.fun(
        posterior_i, null_hypothesis, normal_approximation, silent, density_method,
        label        = .Savage_Dickey_BF.level_label(posterior, i),
        null_mass_NA = null_mass_NA
      )
    }
    names(bf) <- names(posterior)
  }else{
    bf <- .Savage_Dickey_BF.fun(
      posterior, null_hypothesis, normal_approximation, silent, density_method,
      label        = .Savage_Dickey_BF.parameter_label(posterior),
      null_mass_NA = null_mass_NA
    )
  }

  return(bf)
}

# Labels of emitted Savage-Dickey diagnostics follow the marginal summary
# tables: 'parameter[level]', or 'parameter' for a single intercept level.
.Savage_Dickey_BF.parameter_label <- function(posterior){

  parameter <- attr(posterior, "parameter", exact = TRUE)
  if(!is.character(parameter) || length(parameter) != 1L ||
     is.na(parameter) || !nzchar(parameter)){
    return(NULL)
  }

  parameter
}

.Savage_Dickey_BF.level_label <- function(posterior, index){

  parameter <- .Savage_Dickey_BF.parameter_label(posterior)
  if(is.null(parameter)){
    parameter <- .Savage_Dickey_BF.parameter_label(posterior[[index]])
  }
  level <- names(posterior)[index]
  if(is.null(level) || is.na(level) || !nzchar(level)){
    return(parameter)
  }
  if(is.null(parameter)){
    return(level)
  }
  if(length(posterior) == 1L && identical(level, "intercept")){
    return(parameter)
  }

  paste0(parameter, "[", level, "]")
}

.Savage_Dickey_BF_extrapolation_warning <- paste0(
  "Posterior samples do not span both sides of the null hypothesis. ",
  "The posterior density at the null hypothesis is an extrapolation from ",
  "Gaussian kernel tails; the Bayes factor is not reliable evidence."
)

.Savage_Dickey_BF.fun    <- function(posterior, null_hypothesis, normal_approximation, silent, density_method,
                                     label = NULL, null_mass_NA = FALSE){

  prior <- .bt_meta_get(posterior, "prior_density")
  if(is.null(prior))
    stop("there are no prior densities for the posterior distribution", call. = FALSE)

  warnings <- NULL
  stored_posterior_density <- NULL
  stored_posterior_ordinate <- NULL
  posterior_density_source <- if(isTRUE(normal_approximation)) "normal" else "KDE"
  posterior_density_boundary_reflection <- FALSE
  posterior_density_support_bounds <- NULL
  posterior_density_fallback_warnings <- NULL
  BF_error_percent <- NA_real_
  posterior_ordinate_status <- NULL
  posterior_density_status <- NULL
  if(!normal_approximation){
    posterior_ordinate_status <- .posterior_ordinate_direct_status(
      posterior,
      null_hypothesis = null_hypothesis
    )
    posterior_density_status <- .posterior_density_direct_status(posterior)
    if(identical(density_method, "precomputed") &&
       isTRUE(posterior_ordinate_status[["valid"]])){
      stored_posterior_ordinate <- posterior_ordinate_status[["value"]]
    }
    if(identical(density_method, "precomputed") &&
       isTRUE(posterior_density_status[["valid"]])){
      stored_posterior_density <- posterior_density_status[["value"]]
    }
  }
  stored_posterior_density_support <- if(!is.null(stored_posterior_density)){
    stored_posterior_density[["support"]]
  }else{
    NULL
  }
  if(!isTRUE(normal_approximation) &&
     identical(density_method, "precomputed") &&
     is.null(stored_posterior_ordinate) &&
     is.null(stored_posterior_density)){
    if(!is.null(posterior_ordinate_status) &&
       isTRUE(posterior_ordinate_status[["present"]]) &&
       isTRUE(posterior_ordinate_status[["relevant"]]) &&
       !isTRUE(posterior_ordinate_status[["valid"]])){
      stop(
        "Precomputed posterior ordinate metadata is present but invalid ",
        "for the requested null hypothesis.",
        call. = FALSE
      )
    }
    if(!is.null(posterior_density_status) &&
       isTRUE(posterior_density_status[["present"]]) &&
       isTRUE(posterior_density_status[["relevant"]]) &&
       !isTRUE(posterior_density_status[["valid"]])){
      stop(
        "Precomputed posterior density metadata is present but invalid.",
        call. = FALSE
      )
    }
    stop(
      "'density_method = \"precomputed\"' requires valid posterior ordinate ",
      "or posterior density metadata for the requested null hypothesis.",
      call. = FALSE
    )
  }

  posterior_atoms <- .posterior_atoms_get(posterior)
  if(is.null(posterior_atoms)){
    stop(
      "Posterior atom status is unknown. Savage-Dickey evaluation requires an ",
      "explicit atom/no-atom declaration; attach posterior_atom_attribute() ",
      "metadata or use a BayesTools posterior producer that records it. ",
      "Marginal posteriors created by BayesTools 0.3.0 do not record it: ",
      "recompute them with marginal_posterior() from mixed posteriors created ",
      "by the current version (refitting models fitted with 0.3.0).",
      call. = FALSE
    )
  }
  if(ncol(posterior_atoms$locations) != 1L){
    stop("Posterior atom metadata do not describe a scalar marginal posterior.",
         call. = FALSE)
  }
  null_point_mass <- sum(
    posterior_atoms$mass[
      posterior_atoms$locations[, 1L] == null_hypothesis
    ]
  )
  if(null_point_mass > 0){
    point_mass_reason <- paste0(
      "The posterior contains a declared point mass at the exact null ",
      "hypothesis value. The ordinary Savage-Dickey density ratio is invalid."
    )
    if(!isTRUE(null_mass_NA)){
      .hypothesis_stop_ordinate(
        "BayesTools_posterior_point_mass_at_null",
        point_mass_reason
      )
    }
    reason <- if(null_point_mass >= 1 - sqrt(.Machine$double.eps)){
      paste0(
        "The posterior is fixed at the null hypothesis value. The ",
        "Savage-Dickey Bayes factor is undefined."
      )
    }else{
      point_mass_reason
    }
    BF <- NA_real_
    attr(BF, "warnings") <- reason
    attr(BF, "posterior_density_source") <- "null_point_mass"
    if(!silent){
      .Savage_Dickey_BF.emit_warnings(reason, label)
    }
    return(BF)
  }
  # the one exactness rule of point hypotheses (hypothesis_BF()): a prior
  # point mass at the null, a zero, infinite or undefined prior ordinate, and
  # a prior ordinate without an exact structural value (an unknown route, a
  # quadrature rejected by its diagnostics, a density without provenance)
  # stop with their classed conditions. A zero or infinite ordinate makes the
  # density ratio a 0/0 or singular limit that a posterior kernel estimate
  # cannot estimate; a level of a list posterior (marginal inference) then
  # gets an NA Bayes factor with the reason, as a level fixed at the null.
  prior_ordinate <- tryCatch(
    .hypothesis_check_prior_ordinate(
      prior, null_hypothesis,
      label = paste0(if(is.null(label)) "parameter" else label, " = ", null_hypothesis)
    ),
    BayesTools_zero_ordinate     = function(condition) condition,
    BayesTools_infinite_ordinate = function(condition) condition
  )
  if(inherits(prior_ordinate, "condition")){
    if(!isTRUE(null_mass_NA)){
      stop(prior_ordinate)
    }
    behavior <- if(inherits(prior_ordinate, "BayesTools_zero_ordinate")) "zero" else "infinite"
    reason <- paste0(
      "The prior density at the null hypothesis value is ", behavior, ". The ",
      "Savage-Dickey Bayes factor is undefined."
    )
    BF <- NA_real_
    attr(BF, "warnings") <- reason
    attr(BF, "posterior_density_source") <- paste0(behavior, "_prior_ordinate")
    if(!silent){
      .Savage_Dickey_BF.emit_warnings(reason, label)
    }
    return(BF)
  }
  continuous_posterior <- .Savage_Dickey_BF.continuous_posterior(
    posterior = posterior,
    posterior_atoms = posterior_atoms
  )
  continuous_samples <- continuous_posterior$samples
  continuous_mass <- continuous_posterior$continuous_mass

  # A model mixture whose continuous components have different exact supports
  # mixes per-component KDE ordinates (KDE path only).
  component_plan <- NULL
  posterior_density_components <- NULL
  if(!isTRUE(normal_approximation) && is.null(stored_posterior_ordinate) &&
     is.null(stored_posterior_density) && continuous_mass > 0){
    component_plan <- .Savage_Dickey_BF.component_plan(
      samples    = continuous_samples,
      index      = continuous_posterior$component_index,
      components = continuous_posterior$components
    )
    if(!is.null(component_plan[["warning"]])){
      warnings <- c(warnings, component_plan[["warning"]])
      component_plan <- NULL
    }
  }

  # (a null outside the exact support of the prior has a zero prior ordinate
  # and has stopped above)
  posterior_range <- if(length(continuous_samples) > 0L){
    range(continuous_samples)
  }else{
    range(posterior)
  }
  posterior_support_bounds <- .posterior_support_bounds(
    posterior,
    interval_only = TRUE
  )
  if(is.null(posterior_support_bounds)){
    stored_support <- .posterior_support_from_attribute(stored_posterior_density_support)
    if(!is.null(stored_support) && isTRUE(stored_support$exact) &&
       .posterior_support_has_interval(stored_support)){
      posterior_support_bounds <- stored_support$bounds
    }
  }
  null_at_support_boundary <- FALSE
  if(!is.null(posterior_support_bounds)){
    null_at_support_boundary <- any(
      is.finite(posterior_support_bounds) &
        null_hypothesis == posterior_support_bounds
    )
  }
  if(!is.null(component_plan)){
    # only components whose support contains the null contribute an ordinate
    contains <- .Savage_Dickey_BF.component_contains(component_plan, null_hypothesis)
    if(any(contains)){
      posterior_range <- range(unlist(component_plan$draws[contains], use.names = FALSE))
      null_at_support_boundary <- any(vapply(
        component_plan$bounds[contains],
        function(bounds) any(is.finite(bounds) & null_hypothesis == bounds),
        logical(1)
      ))
    }
  }
  if(!is.null(stored_posterior_density) && is.null(stored_posterior_ordinate)){
    posterior_range <- range(stored_posterior_density[["x"]], finite = TRUE)
    if(null_hypothesis < posterior_range[1] || null_hypothesis > posterior_range[2]){
      stop(
        "Stored posterior density does not span both sides of the null hypothesis.",
        call. = FALSE
      )
    }
  }
  if(is.null(stored_posterior_ordinate) &&
     (null_hypothesis < posterior_range[1] || null_hypothesis > posterior_range[2]) &&
     !isTRUE(null_at_support_boundary)){
    warnings <- c(warnings, .Savage_Dickey_BF_extrapolation_warning)
  }

  kde_height <- function(support = NULL){
    if(continuous_mass <= 0 || length(continuous_samples) == 0L){
      height <- 0
      attr(height, "continuous_mass") <- continuous_mass
      return(height)
    }
    height <- .Savage_Dickey_BF.kd(
      samples         = continuous_samples,
      null_hypothesis = null_hypothesis,
      support         = support,
      warn_extrapolation = FALSE
    )
    height <- height * continuous_mass
    attr(height, "continuous_mass") <- continuous_mass
    if(!is.null(attr(height, "kde_extrapolation", exact = TRUE))){
      warnings <<- c(warnings, .Savage_Dickey_BF_extrapolation_warning)
    }
    support_warning <- attr(height, "posterior_support_warning", exact = TRUE)
    if(!is.null(support_warning)){
      warnings <<- c(warnings, support_warning)
    }
    if(isTRUE(attr(height, "boundary_reflection", exact = TRUE))){
      posterior_density_boundary_reflection <<- TRUE
    }
    if(isTRUE(attr(height, "posterior_support_exclusion", exact = TRUE))){
      posterior_density_source <<- "exact_support_exclusion"
    }
    support_bounds <- attr(height, "posterior_support_bounds", exact = TRUE)
    if(!is.null(support_bounds)){
      posterior_density_support_bounds <<- support_bounds
    }
    height
  }

  if(normal_approximation){
    if(continuous_mass <= 0 || length(continuous_samples) == 0L){
      posterior_height <- 0
    }else{
      posterior_height <- .Savage_Dickey_BF.normal(
        continuous_samples,
        null_hypothesis
      ) * continuous_mass
    }
  }else if(!is.null(stored_posterior_ordinate)){
    support_exclusion <- .Savage_Dickey_BF.support_exclusion(
      posterior,
      null_hypothesis = null_hypothesis,
      source_support  = stored_posterior_density_support
    )
    warnings <- c(warnings, support_exclusion[["warnings"]])
    if(isTRUE(support_exclusion[["excluded"]])){
      fallback_warning <- paste0(
        "Exact ", support_exclusion[["source_label"]],
        " excludes the null hypothesis. Ignoring the precomputed posterior ordinate."
      )
      warnings <- c(warnings, fallback_warning)
      posterior_density_fallback_warnings <- c(posterior_density_fallback_warnings, fallback_warning)
      posterior_height <- 0
      posterior_density_source <- "exact_support_exclusion"
      posterior_density_support_bounds <- support_exclusion[["bounds"]]
    }else{
      posterior_height <- stored_posterior_ordinate[["y"]]
      posterior_density_source <- "precomputed"
      BF_error_percent <- .posterior_ordinate_bf_error_percent(stored_posterior_ordinate)
    }
  }else if(!is.null(stored_posterior_density)){
    support_exclusion <- .Savage_Dickey_BF.support_exclusion(
      posterior,
      null_hypothesis = null_hypothesis,
      source_support  = stored_posterior_density[["support"]]
    )
    warnings <- c(warnings, support_exclusion[["warnings"]])
    if(isTRUE(support_exclusion[["excluded"]])){
      fallback_warning <- paste0(
        "Exact ", support_exclusion[["source_label"]],
        " excludes the null hypothesis. Ignoring the precomputed posterior density."
      )
      warnings <- c(warnings, fallback_warning)
      posterior_density_fallback_warnings <- c(posterior_density_fallback_warnings, fallback_warning)
      posterior_height <- 0
      posterior_density_source <- "exact_support_exclusion"
      posterior_density_support_bounds <- support_exclusion[["bounds"]]
    }else{
      posterior_height <- .posterior_density_height(stored_posterior_density, null_hypothesis)
    }
    if(!isTRUE(support_exclusion[["excluded"]]) &&
       (!is.finite(posterior_height) || posterior_height <= 0)){
      stop(
        "Stored posterior density has zero or non-finite height at the null hypothesis.",
        call. = FALSE
      )
    }else if(!isTRUE(support_exclusion[["excluded"]])){
      posterior_density_source <- "precomputed"
      BF_error_percent <- .posterior_density_bf_error_percent(stored_posterior_density, null_hypothesis)
    }
  }else if(!is.null(component_plan)){
    posterior_height <- .Savage_Dickey_BF.component_height(
      plan            = component_plan,
      null_hypothesis = null_hypothesis,
      pooled_samples  = continuous_samples
    ) * continuous_mass
    posterior_density_boundary_reflection <- isTRUE(attr(posterior_height, "boundary_reflection", exact = TRUE))
    if(isTRUE(attr(posterior_height, "posterior_support_exclusion", exact = TRUE))){
      posterior_density_source <- "exact_support_exclusion"
    }
    posterior_density_components <- attr(posterior_height, "components", exact = TRUE)
  }else{
    posterior_height <- kde_height(stored_posterior_density_support)
  }
  if(identical(posterior_density_source, "exact_support_exclusion")){
    # an exact zero posterior density is not a kernel-tail extrapolation
    warnings <- warnings[warnings != .Savage_Dickey_BF_extrapolation_warning]
  }
  warnings <- unique(warnings)
  if(!silent){
    .Savage_Dickey_BF.emit_warnings(warnings, label)
  }

  # the checked prior ordinate is regular with a finite log density
  prior_height <- .prior_linear_density_exact_height(prior_ordinate)
  BF <- exp(log(prior_height) - log(posterior_height))

  if(!is.null(warnings)){
    attr(BF, "warnings") <- warnings
  }
  if(is.finite(BF_error_percent)){
    attr(BF, "BF_error_percent") <- BF_error_percent
  }
  attr(BF, "posterior_density_source") <- posterior_density_source
  if(isTRUE(posterior_density_boundary_reflection)){
    attr(BF, "posterior_density_boundary_reflection") <- TRUE
  }
  if(!is.null(posterior_density_support_bounds)){
    attr(BF, "posterior_density_support") <- posterior_density_support_bounds
  }
  if(!is.null(posterior_density_components)){
    attr(BF, "posterior_density_components") <- posterior_density_components
  }
  if(length(posterior_density_fallback_warnings) > 0L){
    attr(BF, "posterior_density_fallback") <- TRUE
    attr(BF, "posterior_density_fallback_warnings") <- posterior_density_fallback_warnings
  }

  return(BF)
}

.Savage_Dickey_BF.emit_warnings <- function(warnings, label = NULL){

  for(message in warnings){
    if(!is.null(label)){
      message <- paste0(label, ": ", message)
    }
    warning(message, call. = FALSE)
  }

  invisible(NULL)
}

.Savage_Dickey_BF.continuous_posterior <- function(posterior, posterior_atoms){

  sample_values <- as.numeric(posterior)
  finite <- is.finite(sample_values)
  components <- .posterior_components_get(posterior)
  component_index <- NULL
  if(!is.null(components) && length(components$index) == length(sample_values)){
    component_index <- components$index[finite]
  }
  sample_values <- sample_values[finite]
  atom_mass <- sum(posterior_atoms$mass)
  if(!is.finite(atom_mass) || atom_mass < 0 || atom_mass > 1){
    stop("Posterior atom masses must form a finite probability in [0, 1].",
         call. = FALSE)
  }
  continuous_mass <- 1 - atom_mass
  atom_locations <- posterior_atoms$locations[, 1L]
  keep <- rep(TRUE, length(sample_values))
  if(length(atom_locations) > 0L && atom_mass > 0){
    for(location in atom_locations){
      keep <- keep & sample_values != location
    }
  }
  continuous_samples <- sample_values[keep]
  # Preserve support and density metadata used by KDE / boundary reflection.
  kept_metadata <- .bt_meta_get_fields(
    posterior,
    c("support", "posterior_density", "posterior_ordinate", "prior_density")
  )
  kept_metadata <- kept_metadata[!vapply(kept_metadata, is.null, logical(1))]
  if(length(kept_metadata) > 0L){
    continuous_samples <- .bt_meta_assign(continuous_samples, kept_metadata)
  }
  prior_list <- attr(posterior, "prior_list", exact = TRUE)
  if(!is.null(prior_list)){
    attr(continuous_samples, "prior_list") <- prior_list
  }

  list(
    samples         = continuous_samples,
    continuous_mass = continuous_mass,
    atom_mass       = atom_mass,
    component_index = if(!is.null(component_index)) component_index[keep],
    components      = components
  )
}

# Per-component posterior ordinates of a model mixture whose continuous
# components have different exact supports. Returns NULL for the pooled
# estimate: without component metadata, with a single continuous component, or
# when all continuous components share their support. Components whose support
# or draws cannot be used also fall back to the pooled estimate, with a warning.
.Savage_Dickey_BF.component_plan <- function(samples, index, components){

  if(is.null(index) || is.null(components) || length(samples) < 2L){
    return(NULL)
  }

  ids <- sort(unique(index))
  if(length(ids) < 2L){
    return(NULL)
  }

  draws  <- vector("list", length(ids))
  bounds <- vector("list", length(ids))
  for(i in seq_along(ids)){
    support <- components$supports[[ids[i]]]
    if(is.null(support) || !isTRUE(support$exact)){
      return(NULL)
    }
    draws[[i]]   <- as.numeric(samples)[index == ids[i]]
    support_info <- .posterior_support_for_kde(draws[[i]], support = support)
    if(is.null(support_info[["bounds"]])){
      return(list(warning = paste0(
        "Exact posterior support metadata of the mixed model components is ",
        "unusable for their posterior samples. Falling back to the pooled ",
        "kernel density estimate."
      )))
    }
    bounds[[i]] <- support_info[["bounds"]]
  }

  if(all(vapply(bounds, function(x) all(x == bounds[[1L]]), logical(1)))){
    return(NULL)
  }

  list(
    ids    = ids,
    counts = vapply(draws, length, integer(1)),
    draws  = draws,
    bounds = bounds
  )
}

.Savage_Dickey_BF.component_contains <- function(plan, null_hypothesis){

  vapply(plan$bounds, function(bounds){
    null_hypothesis >= bounds[1] && null_hypothesis <= bounds[2]
  }, logical(1))
}

# height = sum_m (n_m / n_c) f_m(null) over the continuous components, where
# f_m is the boundary-reflected KDE of component m's draws with its exact
# support (zero outside the support, the one-sided limit on its bound). A
# component with a single draw uses the pooled bandwidth.
.Savage_Dickey_BF.component_height <- function(plan, null_hypothesis, pooled_samples){

  contains <- .Savage_Dickey_BF.component_contains(plan, null_hypothesis)
  heights <- numeric(length(plan$ids))
  boundary_reflection <- FALSE
  pooled_bw <- NULL

  for(i in which(contains)){
    if(plan$counts[i] >= 2L){
      height <- .Savage_Dickey_BF.kd(
        samples            = plan$draws[[i]],
        null_hypothesis    = null_hypothesis,
        support            = plan$bounds[[i]],
        warn_extrapolation = FALSE
      )
      boundary_reflection <- boundary_reflection ||
        isTRUE(attr(height, "boundary_reflection", exact = TRUE))
    }else{
      if(is.null(pooled_bw)){
        pooled_bw <- stats::bw.nrd0(as.numeric(pooled_samples))
      }
      height <- .density_kde_gaussian_height(
        x      = plan$draws[[i]],
        value  = null_hypothesis,
        bw     = pooled_bw,
        bounds = plan$bounds[[i]]
      )
      boundary_reflection <- boundary_reflection || any(is.finite(plan$bounds[[i]]))
    }
    heights[i] <- as.numeric(height)
  }

  height <- sum(plan$counts / sum(plan$counts) * heights)
  attr(height, "boundary_reflection") <- boundary_reflection
  attr(height, "posterior_support_exclusion") <- !any(contains)
  attr(height, "components") <- data.frame(
    component = plan$ids,
    n         = plan$counts,
    lower     = vapply(plan$bounds, `[`, numeric(1), 1L),
    upper     = vapply(plan$bounds, `[`, numeric(1), 2L),
    ordinate  = heights
  )

  height
}

.Savage_Dickey_BF.normal <- function(samples, null_hypothesis){

  height <- stats::dnorm(null_hypothesis, mean = mean(samples), sd = stats::sd(samples))

  return(height)
}

.Savage_Dickey_BF.support_exclusion <- function(samples, null_hypothesis,
                                                 source_support = NULL){

  supports <- list(source_support, .posterior_support_get(samples))
  support_labels <- c("stored posterior density support", "posterior support")
  valid_bounds <- NULL
  support_warnings <- NULL

  for(i in seq_along(supports)){
    support <- supports[[i]]
    support <- .posterior_support_from_attribute(support)
    if(is.null(support) || !isTRUE(support$exact)){
      next
    }

    support_info <- .posterior_support_for_kde(samples, support = support)
    if(!is.null(support_info[["warning"]])){
      support_warnings <- c(support_warnings, support_info[["warning"]])
    }
    if(is.null(support_info[["bounds"]])){
      next
    }
    if(is.null(valid_bounds)){
      valid_bounds <- support_info[["bounds"]]
    }

    if(.posterior_support_excludes_value(
      .posterior_support_new(support_info[["bounds"]]),
      null_hypothesis
    )){
      return(list(
        excluded    = TRUE,
        bounds      = support_info[["bounds"]],
        source_label = support_labels[[i]],
        warnings    = support_warnings
      ))
    }
  }

  list(excluded = FALSE, bounds = valid_bounds, source_label = NULL,
       warnings = support_warnings)
}

# The posterior ordinate at the null: the exact Gaussian kernel sum of the
# continuous draws at the null, with bandwidth bw.nrd0() of those draws,
# reflected at finite exact support bounds (the one-sided limit on a bound);
# no evaluation grid, interpolation or binning. A null outside the draws (and
# not on a support bound) is a kernel-tail extrapolation.
.Savage_Dickey_BF.kd     <- function(samples, null_hypothesis, support = NULL,
                                     warn_extrapolation = TRUE){

  sample_values <- as.numeric(samples)
  sample_values <- sample_values[is.finite(sample_values)]
  support_info <- .posterior_support_for_kde(samples, support = support)
  support_bounds <- support_info[["bounds"]]
  bounds <- c(-Inf, Inf)

  if(!is.null(support_bounds)){
    if(null_hypothesis < support_bounds[1] || null_hypothesis > support_bounds[2]){
      height <- 0
      attr(height, "posterior_support_bounds") <- support_bounds
      attr(height, "posterior_support_exclusion") <- TRUE
      return(height)
    }
    bounds <- support_bounds
  }

  height <- .density_kde_gaussian_height(
    x      = sample_values,
    value  = null_hypothesis,
    bw     = stats::bw.nrd0(sample_values),
    bounds = bounds
  )
  if(any(is.finite(bounds) & null_hypothesis == bounds)){
    attr(height, "kde_extrapolation") <- NULL
  }
  if(!is.null(attr(height, "kde_extrapolation", exact = TRUE)) &&
     isTRUE(warn_extrapolation)){
    warning(.Savage_Dickey_BF_extrapolation_warning, call. = FALSE)
  }

  if(any(is.finite(bounds))){
    attr(height, "boundary_reflection") <- TRUE
    attr(height, "posterior_support_bounds") <- support_bounds
  }else if(!is.null(support_info[["warning"]])){
    attr(height, "posterior_support_warning") <- support_info[["warning"]]
  }

  return(height)
}
