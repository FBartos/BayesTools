#' @title Transform standardized posterior samples back to original scale
#'
#' @description Transforms posterior samples from standardized continuous
#' predictors back to the original scale. This function is used when predictors
#' were standardized during model fitting via the \code{formula_scale} parameter.
#'
#' @param fit a model fitted with [JAGS_fit()]. Matrices of posterior samples
#' and other fit objects are not supported.
#' @param formula_scale optional nested list containing standardization
#' information keyed by parameter name, replacing the \code{formula_scale}
#' attribute of \code{fit}. Each parameter entry contains scaling info (mean and
#' sd) for each standardized predictor, e.g.,
#' \code{list(mu = list(mu_x1 = list(mean = 0, sd = 1)))}. The transformation is
#' derived from the fitted structure stored in \code{fit}: the fixed-effect
#' design, a log intercept, and the random-effect SD structure.
#'
#' @details The function transforms regression coefficients and intercepts
#' to account for predictor standardization using a combinatorial approach that
#' correctly handles interactions of any order.
#'
#' Internal standardized random-effect coordinates (latent \code{xRE_Zx}
#' columns and realized \code{xRE_COEFx} columns) are omitted from the returned
#' matrix. They do not have the same transformation as fixed coefficients or
#' random-effect covariance summaries and must not be presented as
#' original-scale coefficients. This omission also applies when no scaling
#' information is supplied.
#'
#' For a k-way interaction between standardized predictors, the expansion of
#' \eqn{\prod_{i} (x_i - \mu_i)/\sigma_i} contributes to all lower-order terms.
#' The contribution to a target term T from a source term S (where T is a subset
#' of S's scaled components) is:
#' \deqn{(-1)^{|extra|} \cdot \prod_{i \in extra} \mu_i / \prod_{i \in S_{scaled}} \sigma_i}
#' where \eqn{extra = S_{scaled} \setminus T_{scaled}}.
#' The posterior must contain every lower-order coefficient with a nonzero
#' contribution from this expansion. If it does not, the function stops rather
#' than returning an incomplete transformation.
#'
#' For a fitted object, or formula-scale metadata created by [JAGS_formula()] or
#' [JAGS_fit()], the fixed-coefficient transformation is derived from the
#' stored fixed-effect design (persisted factor levels and contrasts) and
#' verified to reproduce the fitted linear predictor exactly, which also covers
#' nested terms such as \code{~ f/x}. Formulas whose centered terms induce
#' lower-order effects that the fitted formula does not contain (e.g.,
#' \code{~ x + x:f} with \code{x} standardized) have no original-scale
#' representation, and the function stops.
#'
#' @return A numeric matrix of posterior samples transformed back to the
#' original predictor scale, with chains merged in their existing order.
#' When no scaling information is supplied, the function returns the sample
#' matrix without changing its values, rather than returning the fitted object.
#'
#' @seealso [JAGS_formula()] [JAGS_fit()]
#'
#' @export
transform_scale_samples <- function(fit, formula_scale = NULL){

  if(!inherits(fit, "BayesTools_fit")){
    stop(
      "'fit' must be a model fitted with JAGS_fit(); matrices of posterior ",
      "samples and other fit objects are not supported.",
      call. = FALSE
    )
  }
  coordinates <- parameter_coordinates(fit)

  # extract formula_scale from fit if available
  if(is.null(formula_scale) && !is.null(attr(fit, "formula_scale"))){
    formula_scale <- attr(fit, "formula_scale")
  }

  if(is.null(formula_scale) || length(formula_scale) == 0){
    # extract unscaled samples under the same internal-coordinate policy
    return(.bt_remove_internal_random_coordinates(
      posterior = .fit_to_posterior(fit),
      coordinates = coordinates
    ))
  }
  # the fitted structure is the model's, whatever standardization is passed
  formula_scale <- .bt_formula_scale_list_complete(formula_scale, fit)

  .bt_transform_scale_posterior(
    posterior      = as.matrix(.fit_to_posterior(fit)),
    formula_scale  = formula_scale,
    formula_design = attr(fit, "formula_design", exact = TRUE),
    coordinates    = coordinates
  )
}

# Combinatorial unscaling of a posterior sample matrix with fitted
# formula-scale information (which carries the fitted unscale designs).
.bt_transform_scale_posterior <- function(posterior, formula_scale,
                                          formula_design = NULL,
                                          coordinates = NULL, targets = NULL){

  .check_formula_scale_info(formula_scale)
  formula_scale <- .bt_formula_scale_list_with_unscale_designs(
    formula_scale,
    formula_design
  )

  posterior <- .apply_unscale_transform(as.matrix(posterior), formula_scale, targets = targets)
  .bt_remove_internal_random_coordinates(
    posterior = posterior,
    coordinates = coordinates
  )
}

.bt_remove_internal_random_coordinates <- function(posterior,
                                                   coordinates = NULL){

  posterior <- as.matrix(posterior)
  column_names <- colnames(posterior)
  if(is.null(column_names) || length(column_names) == 0L){
    return(posterior)
  }

  remove <- grepl(
    "_xRE_(?:GROUP_Z|UNIT_COEF|COEF|Z)x(?:\\[|$)",
    column_names,
    perl = TRUE
  )

  if(!is.null(coordinates)){
    .bt_validate_parameter_coordinates(coordinates)
    coordinate_rows <- match(
      column_names,
      coordinates$coordinate_name
    )
    registered <- !is.na(coordinate_rows)
    remove[registered] <- coordinates$internal[coordinate_rows[registered]] &
      coordinates$role[coordinate_rows[registered]] %in% c(
        "random_latent",
        "random_group_coefficient"
      )
  }

  posterior[, !remove, drop = FALSE]
}


#' @title Transform prior samples to original scale
#'
#' @description Generate prior samples and transform them using the same
#' matrix transformation as posterior samples. This is the correct approach for
#' visualizing priors on the original (unscaled) scale, especially for the intercept
#' which depends on contributions from multiple coefficient priors.
#'
#' @param fit a fitted model object with \code{prior_list} and optionally
#' \code{formula_scale} attributes
#' @param n_samples number of samples to generate (default: 10000)
#' @param seed random seed for reproducibility (optional). The caller's
#' random-number state (\code{.Random.seed} and \code{RNGkind()}) is restored
#' afterwards. With \code{NULL}, the caller's random-number stream is used.
#' @param formula_scale optional nested list containing standardization information.
#' If not provided, extracted from \code{fit} attribute. The transformation is
#' derived from the fitted structure stored in \code{fit}: the fixed-effect
#' design, a log intercept, and the random-effect SD structure.
#'
#' @details When models use auto-scaling (standardizing predictors), the posterior
#' samples are on the standardized scale. To correctly visualize priors on the
#' original scale, we cannot simply apply a linear transformation to individual
#' priors because the intercept on the original scale is a weighted sum of
#' multiple priors:
#'
#' \deqn{\beta_0^{orig} = \beta_0^* - \sum_i \frac{\mu_i}{\sigma_i} \beta_i^*}
#'

#' This function generates samples from ALL priors simultaneously and applies
#' the same matrix transformation used for posterior samples, which correctly
#' handles the intercept and all other parameters.
#'
#' Monitored random-effect nodes that the model defines deterministically from
#' other nodes are computed from the prior draws of those nodes with the
#' model's own definitions: standard deviations derived from a variance
#' allocation ([random_variance_allocation()]), scalar correlations
#' sampled on the Fisher-z or logit scale, and LKJ Cholesky factors,
#' correlation matrices, and partial correlations. Every parameter-catalog
#' quantity that depends only on such nodes and on the priors can therefore be
#' evaluated on the prior draws (see [parameter_draws()]). The auxiliary
#' nodes of spike-and-slab and mixture priors are included from the components
#' of their [rng()] draws: the component indicator (\code{<parameter>_indicator})
#' of both, and the inclusion probability (\code{<parameter>_inclusion}) and
#' slab draws (\code{<parameter>_variable}) of spike-and-slab priors; the
#' random-number stream of the other columns is unchanged. The totals of
#' ordered priors (\code{<parameter>_ordered_total}) are included with the
#' same auxiliary nodes of a spike-and-slab or mixture total; the slices of a
#' spike-and-slab total of an interaction share its inclusion probability and
#' indicator, as in the fitted model. Standardized latent random effects,
#' nodes derived from them, such as realized group coefficients, the
#' component nodes of mixture priors, and the auxiliary nodes of Dirichlet
#' priors are not included, so catalog quantities that depend on them cannot
#' be evaluated on the prior draws.
#' Variance-allocation inclusion indicators are included, drawn from their
#' inclusion probabilities.
#' Valid LKJ concentrations do not guarantee representable primitive Beta draws.
#' Nonfinite draws or draws outside the open interval `(0, 1)` stop before
#' deterministic replay with `BayesTools_lkj_rng_unavailable` (parent
#' `BayesTools_prior_rng_unavailable`). The condition records `block`,
#' `primitive`, `K`, `eta`, original-row `failed_draws` and `n_failed`.
#' Choose a concentration whose primitive draws remain representable; no retry,
#' clamping, endpoint substitution or prior change is performed. Explicit scalar
#' correlation transforms that numerically saturate a boundary retain their
#' compiler refusal rather than changing the requested prior or scale.
#'
#' @return A matrix of prior samples on the original (unscaled) scale, with
#' columns matching the structure of posterior samples.
#'
#' @seealso [transform_scale_samples()] [plot_posterior()]
#'
#' @examples
#' # With a fitted model that used formula_scale:
#' # prior_samples <- transform_prior_samples(fit, n_samples = 10000)
#' # This can then be used with density() or for custom plotting
#'
#' @export
transform_prior_samples <- function(fit, n_samples = 10000, seed = NULL, formula_scale = NULL){

  # A seeded call leaves the caller's random-number state as it found it.
  if(!is.null(seed)){
    rng_state <- .bt_rng_state()
    on.exit(.bt_rng_restore(rng_state), add = TRUE)
  }
  check_int(n_samples, "n_samples", lower = 1, allow_NA = FALSE)
  check_int(
    seed,
    "seed",
    lower = 0,
    upper = .Machine$integer.max,
    allow_NULL = TRUE,
    allow_NA = FALSE
  )

  # Extract prior_list from fit

  prior_list <- attr(fit, "prior_list")

  if(is.null(prior_list)){
    stop("'fit' must have 'prior_list' attribute.")
  }

  # Extract formula_scale from fit if not provided; the fitted structure is
  # the model's, whatever standardization is passed
  if(is.null(formula_scale)){
    formula_scale <- attr(fit, "formula_scale")
  }
  formula_scale <- .bt_formula_scale_list_complete(formula_scale, fit)

  # Get posterior column names for structure matching
  if(inherits(fit, "runjags") || inherits(fit, "BayesTools_fit")){
    posterior <- as.matrix(.fit_to_posterior(fit))
  }else{
    stop("'fit' must be a fitted model object.")
  }

  prior_samples <- .generate_transformed_prior_samples(
    prior_list   = prior_list,
    column_names = colnames(posterior),
    n_samples    = n_samples,
    seed         = seed,
    formula_scale = formula_scale,
    formula_design = attr(fit, "formula_design", exact = TRUE)
  )

  return(prior_samples)
}


# Helper: Generate prior samples and apply the same unscaling transform used
# for posterior samples.
#
# @param prior_list Named list of prior objects
# @param column_names Column names to match from the posterior structure
# @param n_samples Number of samples to generate
# @param seed Optional random seed
# @param formula_scale Optional nested scaling information
# @return Matrix with transformed prior samples
.bt_expression_point_prior_names <- function(prior_list){

  names(prior_list)[vapply(prior_list, function(prior){
    is.prior.point(prior) && is.expression(prior$parameters$location)
  }, logical(1))]
}

.generate_transformed_prior_samples <- function(
    prior_list, column_names, n_samples, seed = NULL, formula_scale = NULL,
    formula_design = NULL, retain_state = FALSE){

  if(!is.null(formula_scale) && length(formula_scale) > 0){
    .check_formula_scale_info(formula_scale)
  }

  expression_points <- .bt_expression_point_prior_names(prior_list)
  needed <- column_names
  for(prefix in names(formula_scale)){
    scale <- formula_scale[[prefix]]
    spec <- attr(scale, "unscale_design", exact = TRUE)
    if(is.null(spec)) next
    transform <- .bt_formula_coefficient_transform(names(spec$multipliers), scale, prefix)
    requested <- intersect(column_names, transform$target_names)
    needed <- union(needed, transform$dependencies$source[transform$dependencies$target %in% requested])
    needed <- union(needed, transform$state_dependencies$source[transform$state_dependencies$target %in% requested])
  }
  unsupported <- expression_points[vapply(expression_points, function(owner){
    any(.prior_linear_prior_columns(owner, prior_list[[owner]]) %in% needed)
  }, logical(1))]
  if(length(unsupported)) .bt_formula_transform_stop(
    "Fresh expression-point replay is unavailable without a certified persisted recipe.",
    reason = "uncertified_point_replay", missing = unsupported)

  generation_priors <- prior_list[!names(prior_list) %in% expression_points]
  prior_samples <- if(length(generation_priors)) .generate_prior_sample_matrix(
    prior_list     = generation_priors,
    n_samples      = n_samples,
    column_names   = NULL,
    seed           = seed
  ) else matrix(numeric(), n_samples, 0L, dimnames = list(NULL, character()))
  retained_primitives <- prior_samples
  prior_samples <- .bt_add_lkj_prior_samples(
    samples        = prior_samples,
    formula_design = formula_design,
    column_names   = column_names,
    n_samples      = n_samples
  )
  retained_columns <- setdiff(colnames(retained_primitives), colnames(prior_samples))
  if(length(retained_columns)) prior_samples <- cbind(prior_samples,
    retained_primitives[, retained_columns, drop = FALSE])
  retained_primitives <- prior_samples
  prior_samples <- .bt_add_random_allocation_indicator_prior_samples(
    samples      = prior_samples,
    prior_list   = prior_list,
    column_names = column_names,
    n_samples    = n_samples
  )
  retained_columns <- setdiff(colnames(retained_primitives), colnames(prior_samples))
  if(length(retained_columns)) prior_samples <- cbind(prior_samples,
    retained_primitives[, retained_columns, drop = FALSE])
  retained_primitives <- prior_samples
  prior_samples <- .bt_add_random_deterministic_prior_samples(
    samples        = prior_samples,
    prior_list     = prior_list,
    formula_design = formula_design,
    column_names   = column_names
  )
  retained_columns <- setdiff(colnames(retained_primitives), colnames(prior_samples))
  if(length(retained_columns)) prior_samples <- cbind(prior_samples,
    retained_primitives[, retained_columns, drop = FALSE])

  if(!is.null(formula_scale) && length(formula_scale) > 0){
    prior_samples <- .apply_unscale_transform(prior_samples, formula_scale, targets = column_names)
  }

  if(retain_state) return(prior_samples)
  return(prior_samples[, intersect(column_names, colnames(prior_samples)), drop = FALSE])
}

.bt_formula_sample_contributions <- function(context, samples){

  if(!identical(context$linear_weight_space, "formula_contribution")) return(samples)
  fitted <- samples
  for(owner in names(context$prior_list)){
    prior <- context$prior_list[[owner]]
    if(!is.prior(prior)) next
    columns <- intersect(.prior_linear_prior_columns(owner, prior), context$column_names)
    if(!length(columns)) next
    multiplier <- attr(prior, "multiply_by", exact = TRUE)
    if(is.null(multiplier)) next
    if(is.prior.point(prior) && is.numeric(prior$parameters$location) && isTRUE(prior$parameters$location == 0)){
      samples[, columns] <- 0
      next
    }
    if(is.character(multiplier)){
      if(!multiplier %in% colnames(fitted)) .bt_formula_transform_stop(
        "Fresh formula contribution state is unavailable.", reason = "missing_multiplier_state", state = multiplier)
      multiplier <- fitted[, multiplier]
    }
    if(any(!is.finite(multiplier))) .bt_formula_transform_stop(
      "Fresh formula contribution state is nonfinite.", reason = "nonfinite_transform")
    if(isTRUE(all(multiplier == 0))){
      samples[, columns] <- 0
    }else samples[, columns] <- fitted[, columns, drop = FALSE] * multiplier
    if(any(!is.finite(samples[, columns, drop = FALSE]))) .bt_formula_transform_stop(
      "Fresh formula contribution arithmetic is nonfinite.", reason = "nonfinite_transform")
  }
  samples
}


.bt_add_random_allocation_indicator_prior_samples <- function(
    samples, prior_list, column_names, n_samples){

  for(parameter in names(prior_list)){
    prior <- prior_list[[parameter]]
    indicator <- attr(prior, "random_allocation_indicator", exact = TRUE)
    if(!is.character(indicator) || length(indicator) != 1L ||
       is.na(indicator) || !nzchar(indicator) ||
       !indicator %in% column_names || indicator %in% colnames(samples)){
      next
    }
    if(!parameter %in% colnames(samples)){
      stop(
        "Random-effect allocation gate prior samples are missing probability coordinate '",
        parameter, "'.",
        call. = FALSE
      )
    }

    probability <- samples[, parameter]
    if(length(probability) != n_samples || any(!is.finite(probability)) ||
       any(probability < 0 | probability > 1)){
      stop(
        "Random-effect allocation gate prior probabilities are invalid for '",
        parameter, "'.",
        call. = FALSE
      )
    }
    samples <- cbind(
      samples,
      stats::rbinom(n_samples, size = 1L, prob = probability)
    )
    colnames(samples)[ncol(samples)] <- indicator
  }

  ordered <- intersect(column_names, colnames(samples))
  samples[, ordered, drop = FALSE]
}


# Generated deterministic monitors (allocation-derived SDs, Fisher-z and logit
# scalar correlations, LKJ factors, matrices, and partial correlations) are
# computed from the prior draws of their dependencies with the registered node
# evaluators, which the posterior draws use as well, so that every quantity
# available from the posterior draws is also available from the prior draws.
.bt_add_random_deterministic_prior_samples <- function(samples, prior_list,
                                                       formula_design,
                                                       column_names){

  if(!is.list(formula_design)){
    return(samples)
  }
  samples <- .bt_add_deterministic_prior_samples(
    samples        = samples,
    prior_list     = prior_list,
    formula_design = formula_design,
    column_names   = column_names
  )

  ordered <- intersect(column_names, colnames(samples))
  samples[, ordered, drop = FALSE]
}

# Monitored generated deterministic nodes missing from the prior draws are
# computed from the prior draws of their dependencies with the registered node
# evaluators. A node whose dependencies have no prior draws (for example an
# external SD source defined in the model syntax) stays unavailable, as it is
# in fitted draws without it; like the monitored JAGS node, the values are not
# validated here (the catalog evaluators check supports).
.bt_add_deterministic_prior_samples <- function(samples, prior_list,
                                                formula_design, column_names){

  nodes <- .bt_deterministic_nodes(
    prior_list = prior_list,
    formula_design = formula_design
  )
  for(node in nodes){
    targets <- node$coordinates[
      node$coordinates %in% column_names
    ]
    if(length(targets) == 0L){
      next
    }
    values <- .bt_deterministic_node_evaluate(
      node,
      .bt_deterministic_lookup(samples, prior_list)
    )
    if(is.null(values)){
      next
    }
    present <- intersect(targets, colnames(samples))
    absent <- setdiff(targets, colnames(samples))
    if(length(present)) samples[, present] <- values[, present, drop = FALSE]
    if(length(absent)) samples <- cbind(samples, values[, absent, drop = FALSE])
  }

  samples
}


.bt_add_lkj_prior_samples <- function(samples, formula_design, column_names,
                                      n_samples){

  if(!is.list(formula_design)){
    return(samples)
  }
  terms <- unlist(lapply(formula_design, `[[`, "random_effects"), recursive = FALSE)
  for(random_term in terms){
    correlation <- random_term$correlation
    if(!is.list(correlation) || !identical(correlation$type, "lkj")){
      next
    }
    K <- random_term$n_columns
    eta <- correlation$eta
    primitive_names <- correlation$primitive_names
    alpha <- .bt_lkj_cholesky_alpha(K = K, eta = eta)
    if(length(primitive_names) != length(alpha)){
      stop(
        "Stored LKJ primitive metadata do not match the random-effect dimension.",
        call. = FALSE
      )
    }
    keep <- primitive_names %in% column_names &
      !primitive_names %in% colnames(samples)
    for(i in which(keep)){
      primitive_draws <- stats::rbeta(n_samples, shape1 = alpha[[i]], shape2 = alpha[[i]])
      failed_draws <- which(!is.finite(primitive_draws) | primitive_draws <= 0 | primitive_draws >= 1)
      if(length(failed_draws) > 0L){
        stop(errorCondition(
          paste0("LKJ prior sampling for block '", random_term$block_name,
            "' is numerically unavailable: primitive '", primitive_names[[i]],
            "' produced ", length(failed_draws),
            " nonfinite or out-of-support draws outside (0, 1) at K = ", as.integer(K),
            " and eta = ", format(eta, digits = 17, trim = TRUE),
            ". Choose a concentration whose primitive draws remain representable inside this interval."),
          call = NULL,
          class = c("BayesTools_lkj_rng_unavailable", "BayesTools_prior_rng_unavailable"),
          block = random_term$block_name, primitive = primitive_names[[i]],
          K = K, eta = eta, failed_draws = failed_draws, n_failed = length(failed_draws)
        ))
      }
      samples <- cbind(samples, primitive_draws)
      colnames(samples)[ncol(samples)] <- primitive_names[[i]]
    }
    # The monitored Cholesky factors, correlation matrices, and partial
    # correlations of the primitives are added with the other generated
    # deterministic nodes (.bt_add_deterministic_prior_samples()).
  }

  ordered <- intersect(column_names, colnames(samples))
  samples[, ordered, drop = FALSE]
}

# 'auxiliary': also return the columns of the fitted auxiliary nodes of a
# spike-and-slab or mixture prior (its indicator, and the inclusion
# probability and slab draws of a spike-and-slab prior).
.generate_factor_prior_sample_matrix <- function(prior, parameter, n_samples,
                                                 auxiliary = FALSE, allocation_registry = NULL){

  K <- .get_prior_factor_levels(prior)
  if(is.null(K) || is.na(K)){
    stop("The number of factor coefficients for prior '", parameter, "' is unknown.", call. = FALSE)
  }

  auxiliary_samples <- NULL
  if(is.prior.spike_and_slab(prior)){
    prior_variable  <- .get_spike_and_slab_variable(prior)
    prior_inclusion <- .get_spike_and_slab_inclusion(prior)
    variable <- .generate_factor_prior_sample_matrix(prior_variable, parameter, n_samples)
    inclusion_probability <- rng(prior_inclusion, n_samples)
    inclusion <- stats::rbinom(n_samples, size = 1, prob = inclusion_probability)
    samples <- variable * inclusion
    colnames(variable) <- .JAGS_prior_factor_names(paste0(parameter, "_variable"), prior_variable)
    auxiliary_samples <- cbind(
      matrix(inclusion, ncol = 1L, dimnames = list(NULL, paste0(parameter, "_indicator"))),
      matrix(inclusion_probability, ncol = 1L, dimnames = list(NULL, paste0(parameter, "_inclusion"))),
      variable
    )
  }else if(is.prior.mixture(prior)){
    prior_weights <- attr(prior, "prior_weights")
    prior_weights <- prior_weights / sum(prior_weights)
    components <- sample(seq_along(prior_weights), size = n_samples, replace = TRUE, prob = prior_weights)
    samples <- matrix(NA_real_, nrow = n_samples, ncol = K)

    for(component in unique(components)){
      samples[components == component, ] <- .generate_factor_prior_sample_matrix(
        prior      = prior[[component]],
        parameter  = parameter,
        n_samples  = sum(components == component)
      )
    }
    auxiliary_samples <- matrix(
      components, ncol = 1L, dimnames = list(NULL, paste0(parameter, "_indicator"))
    )
  }else if(is.prior.point(prior)){
    location <- prior$parameters[["location"]]
    samples <- matrix(rep(location, length.out = K), nrow = n_samples, ncol = K, byrow = TRUE)
  }else if(is.prior.ordered(prior)){
    # the coefficients in the random-number stream of rng(), with the fitted
    # nodes of the total
    draws <- .prior_ordered_draws(prior, n_samples, allocation_registry = allocation_registry)
    samples <- draws$coefficients
    auxiliary_samples <- .prior_ordered_total_samples(draws, parameter)
    if(!is.null(allocation_registry)){
      for(record in .prior_ordered_dirichlet_records(prior)){
        if(isTRUE(allocation_registry[[record$key]]$exported)) next
        allocation <- draws$allocation_samples[[record$key]]
        colnames(allocation) <- paste0(record$node, "[", seq_len(record$dim), "]")
        auxiliary_samples <- cbind(auxiliary_samples, allocation)
        allocation_registry[[record$key]]$exported <- TRUE
      }
    }
  }else if(is.prior.orthonormal(prior) || is.prior.meandif(prior)){
    prior$parameters[["K"]] <- K
    samples <- rng(prior, n_samples, transform_factor_samples = FALSE)
    samples <- matrix(samples, nrow = n_samples, ncol = K)
  }else if(is.prior.treatment(prior) || is.prior.independent(prior)){
    samples <- replicate(K, rng(prior, n_samples, transform_factor_samples = FALSE))
    samples <- matrix(samples, nrow = n_samples, ncol = K)
  }else{
    samples <- rng(prior, n_samples, transform_factor_samples = FALSE)
    samples <- matrix(samples, nrow = n_samples, ncol = K)
  }

  colnames(samples) <- .JAGS_prior_factor_names(parameter, prior)
  .prior_numerical_finite(samples, prior$distribution)
  if(isTRUE(auxiliary) && !is.null(auxiliary_samples)){
    samples <- cbind(samples, auxiliary_samples)
  }
  return(samples)
}


# Helper: Generate a matrix of prior samples matching posterior structure
#
# @param prior_list Named list of prior objects
# @param n_samples Number of samples to generate
# @param column_names Optional vector of column names to match (filters output)
# @param seed Optional random seed
# @return Matrix with prior samples (rows = samples, columns = parameters)
.generate_prior_sample_matrix <- function(prior_list, n_samples, column_names = NULL, seed = NULL){

  expression_points <- .bt_expression_point_prior_names(prior_list)
  if(length(expression_points)) .bt_formula_transform_stop(
    "Fresh expression-point replay is unavailable without a certified persisted recipe.",
    reason = "uncertified_point_replay", missing = expression_points)
  if(!is.null(seed)){
    set.seed(seed)
  }

  # Determine which parameters to sample
  param_names <- names(prior_list)

  if(is.null(param_names) || length(param_names) == 0){
    stop("'prior_list' must be a named list of priors.")
  }

  # Initialize list to collect samples (handles varying column counts per prior)
  samples_list <- list()
  ordered_priors <- prior_list[vapply(prior_list, function(prior){
    is.prior.ordered(prior) && !is.prior.mixture(prior) &&
      !is.null(attr(prior, "ordered_metadata", exact = TRUE))
  }, logical(1))]
  .bt_validate_ordered_shared_allocations(ordered_priors)
  allocation_registry <- new.env(parent = emptyenv())

  for(param_name in param_names){
    prior <- prior_list[[param_name]]

    if(is.null(prior)){
      # No prior for this parameter - use zeros
      samples_list[[param_name]] <- matrix(0, nrow = n_samples, ncol = 1)
      colnames(samples_list[[param_name]]) <- param_name

    }else if(is.prior.none(prior)){
      # No effect prior - use zeros
      samples_list[[param_name]] <- matrix(0, nrow = n_samples, ncol = 1)
      colnames(samples_list[[param_name]]) <- param_name

    }else if(is.prior.factor(prior) || inherits(prior, "prior.factor_mixture") || inherits(prior, "prior.factor_spike_and_slab")){
      samples_list[[param_name]] <- .generate_factor_prior_sample_matrix(
        prior, param_name, n_samples, auxiliary = TRUE,
        allocation_registry = if(param_name %in% names(ordered_priors)) allocation_registry else NULL
      )

    }else if(is.prior.spike_and_slab(prior)){
      # the value and the fitted auxiliary nodes, from rng()'s components
      parts <- .rng_spike_and_slab_parts(prior, n_samples)
      samples_list[[param_name]] <- cbind(
        as.numeric(parts$value),
        parts$inclusion,
        parts$inclusion_probability,
        as.numeric(parts$variable)
      )
      colnames(samples_list[[param_name]]) <- paste0(
        param_name, c("", "_indicator", "_inclusion", "_variable")
      )

    }else if(is.prior.point(prior)){
      # Point prior - constant values
      samples_list[[param_name]] <- matrix(
        prior$parameters[["location"]],
        nrow = n_samples,
        ncol = 1
      )
      colnames(samples_list[[param_name]]) <- param_name

    }else if(is.prior.simple(prior)){
      # Simple priors - single column
      samples_list[[param_name]] <- matrix(
        rng(prior, n_samples),
        nrow = n_samples,
        ncol = 1
      )
      colnames(samples_list[[param_name]]) <- param_name

    }else if(is.prior.vector(prior)){
      # Vector priors return matrix from rng
      temp_samples <- rng(prior, n_samples)
      if(is.matrix(temp_samples)){
        n_cols <- ncol(temp_samples)
        col_names <- paste0(param_name, "[", 1:n_cols, "]")
        colnames(temp_samples) <- col_names
        samples_list[[param_name]] <- temp_samples
      }else{
        samples_list[[param_name]] <- matrix(temp_samples, nrow = n_samples, ncol = 1)
        colnames(samples_list[[param_name]]) <- param_name
      }

    }else{
      temp_samples <- tryCatch(
        rng(prior, n_samples),
        error = function(e){
          if(inherits(e, "BayesTools_prior_rng_unavailable")){
            stop(e)
          }
          stop(
            "Could not generate samples for prior '", param_name, "': ",
            e$message,
            call. = FALSE
          )
        }
      )
      if(!is.numeric(temp_samples)){
        stop(
          "Could not generate samples for prior '", param_name,
          "': rng() did not return numeric samples.",
          call. = FALSE
        )
      }
      .prior_numerical_finite(temp_samples, prior$distribution)

      components <- attr(temp_samples, "components", exact = TRUE)
      if(is.matrix(temp_samples)){
        n_cols <- ncol(temp_samples)
        col_names <- paste0(param_name, "[", 1:n_cols, "]")
        colnames(temp_samples) <- col_names
        samples_list[[param_name]] <- temp_samples
      }else{
        samples_list[[param_name]] <- matrix(temp_samples, nrow = n_samples, ncol = 1)
        colnames(samples_list[[param_name]]) <- param_name
      }
      # the component indicator of a mixture prior (a fitted auxiliary node)
      if(is.prior.mixture(prior) && !is.null(components)){
        samples_list[[param_name]] <- cbind(
          samples_list[[param_name]],
          matrix(components, ncol = 1L,
                 dimnames = list(NULL, paste0(param_name, "_indicator")))
        )
      }
    }
  }

  # Combine all samples into one matrix
  samples <- do.call(cbind, samples_list)
  .prior_numerical_finite(samples, "joint prior")

  # Filter to match column_names if provided
  if(!is.null(column_names)){
    available_cols <- intersect(column_names, colnames(samples))
    if(length(available_cols) > 0){
      samples <- samples[, available_cols, drop = FALSE]
    }
  }

  return(samples)
}
