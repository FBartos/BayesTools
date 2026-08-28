#' @title Transform standardized posterior samples back to original scale
#'
#' @description Transforms posterior samples from standardized continuous
#' predictors back to the original scale. This function is used when predictors
#' were standardized during model fitting via the \code{formula_scale} parameter.
#'
#' @param fit a fitted model object with \code{formula_scale} attribute, or
#' a matrix of posterior samples
#' @param formula_scale nested list containing standardization information keyed by
#' parameter name. Each parameter entry contains scaling info (mean and sd) for
#' each standardized predictor, e.g., \code{list(mu = list(mu_x1 = list(mean = 0, sd = 1)))}.
#' If \code{fit} is provided and has a \code{formula_scale} attribute, this will be used automatically.
#'
#' @details The function transforms regression coefficients and intercepts
#' to account for predictor standardization using a combinatorial approach that
#' correctly handles interactions of any order.
#'
#' Internal standardized random-effect coordinates (latent \code{xRE_Zx}
#' columns and realized \code{xRE_COEFx} columns) are omitted from the returned
#' matrix. They do not have the same transformation as fixed coefficients or
#' random-effect covariance summaries and must not be presented as
#' original-scale coefficients.
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
#' @return \code{transform_scale_samples} returns posterior samples transformed
#' back to the original predictor scale.
#'
#' @seealso [JAGS_formula()] [JAGS_fit()]
#'
#' @export
transform_scale_samples <- function(fit, formula_scale = NULL){

  coordinates <- NULL
  if(inherits(fit, "BayesTools_fit")){
    coordinates <- parameter_coordinates(fit)
  }

  # extract formula_scale from fit if available
  if(is.null(formula_scale) && !is.null(attr(fit, "formula_scale"))){
    formula_scale <- attr(fit, "formula_scale")
  }

  if(is.null(formula_scale) || length(formula_scale) == 0){
    # no scaling information, return as is
    return(fit)
  }

  .check_formula_scale_info(formula_scale)

  # extract posterior samples
  if(inherits(fit, "runjags") || inherits(fit, "BayesTools_fit")){
    posterior <- as.matrix(.fit_to_posterior(fit))
  }else if(is.matrix(fit)){
    posterior <- fit
  }else{
    stop("'fit' must be a fitted model object or a matrix of posterior samples.")
  }

  # Apply the combinatorial unscaling transformation
  posterior <- .apply_unscale_transform(posterior, formula_scale)
  posterior <- .bt_remove_internal_random_coordinates(
    posterior = posterior,
    coordinates = coordinates
  )

  return(posterior)
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
#' @param seed random seed for reproducibility (optional)
#' @param formula_scale optional nested list containing standardization information.
#' If not provided, extracted from \code{fit} attribute.
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

  # Extract formula_scale from fit if not provided
  if(is.null(formula_scale)){
    formula_scale <- attr(fit, "formula_scale")
  }

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
.generate_transformed_prior_samples <- function(
    prior_list, column_names, n_samples, seed = NULL, formula_scale = NULL,
    formula_design = NULL){

  if(!is.null(formula_scale) && length(formula_scale) > 0){
    .check_formula_scale_info(formula_scale)
  }

  prior_samples <- .generate_prior_sample_matrix(
    prior_list     = prior_list,
    n_samples      = n_samples,
    column_names   = column_names,
    seed           = seed
  )
  prior_samples <- .bt_add_lkj_prior_samples(
    samples        = prior_samples,
    formula_design = formula_design,
    column_names   = column_names,
    n_samples      = n_samples
  )
  prior_samples <- .bt_add_random_allocation_indicator_prior_samples(
    samples      = prior_samples,
    prior_list   = prior_list,
    column_names = column_names,
    n_samples    = n_samples
  )

  if(!is.null(formula_scale) && length(formula_scale) > 0){
    prior_samples <- .apply_unscale_transform(prior_samples, formula_scale)
  }

  return(prior_samples)
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
      samples <- cbind(
        samples,
        stats::rbeta(n_samples, shape1 = alpha[[i]], shape2 = alpha[[i]])
      )
      colnames(samples)[ncol(samples)] <- primitive_names[[i]]
    }
  }

  ordered <- intersect(column_names, colnames(samples))
  samples[, ordered, drop = FALSE]
}

.generate_factor_prior_sample_matrix <- function(prior, parameter, n_samples){

  K <- .get_prior_factor_levels(prior)
  if(is.null(K) || is.na(K)){
    stop("The number of factor coefficients for prior '", parameter, "' is unknown.", call. = FALSE)
  }

  if(is.prior.spike_and_slab(prior)){
    prior_variable  <- .get_spike_and_slab_variable(prior)
    prior_inclusion <- .get_spike_and_slab_inclusion(prior)
    samples <- .generate_factor_prior_sample_matrix(prior_variable, parameter, n_samples)
    inclusion <- stats::rbinom(n_samples, size = 1, prob = rng(prior_inclusion, n_samples))
    samples <- samples * inclusion
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
  }else if(is.prior.point(prior)){
    location <- prior$parameters[["location"]]
    samples <- matrix(rep(location, length.out = K), nrow = n_samples, ncol = K, byrow = TRUE)
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
      samples_list[[param_name]] <- .generate_factor_prior_sample_matrix(prior, param_name, n_samples)

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

      if(is.matrix(temp_samples)){
        n_cols <- ncol(temp_samples)
        col_names <- paste0(param_name, "[", 1:n_cols, "]")
        colnames(temp_samples) <- col_names
        samples_list[[param_name]] <- temp_samples
      }else{
        samples_list[[param_name]] <- matrix(temp_samples, nrow = n_samples, ncol = 1)
        colnames(samples_list[[param_name]]) <- param_name
      }
    }
  }

  # Combine all samples into one matrix
  samples <- do.call(cbind, samples_list)

  # Filter to match column_names if provided
  if(!is.null(column_names)){
    available_cols <- intersect(column_names, colnames(samples))
    if(length(available_cols) > 0){
      samples <- samples[, available_cols, drop = FALSE]
    }
  }

  return(samples)
}
