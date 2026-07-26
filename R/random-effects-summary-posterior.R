#' Extract random-effect summary posterior distributions
#'
#' @description Extracts posterior draws for derived random-effect summary
#' quantities and attaches analytic marginal prior densities where available.
#' For allocation summaries, the helper keeps raw Dirichlet allocation weights
#' internal and exposes interpretable scalar summaries such as mean-variance
#' SD-component variance ratios.
#'
#' @param fit model fit created by [JAGS_fit].
#' @param summary summary quantity to extract. `"variance_ratio"` returns
#'   `K * w` for mean-variance SD-component allocations, `"variance_fraction"`
#'   returns true total-variance fractions `w`, and `"sd_multiplier"` returns
#'   SD multipliers.
#' @param allocation optional allocation label filter.
#' @param component optional allocation component label filter.
#' @param formula_parameter optional formula parameter filter.
#' @param n_prior_points number of grid points used for the attached analytic
#'   prior density.
#'
#' @return A named list of posterior vectors with class `mixed_posteriors` and
#'   `marginal_posterior`. The list can be passed to [plot_posterior], and the
#'   individual vectors carry a `prior_density` attribute.
#'
#' @export
random_effects_summary_posterior <- function(
    fit,
    summary = c("variance_ratio", "variance_fraction", "sd_multiplier"),
    allocation = NULL,
    component = NULL,
    formula_parameter = NULL,
    n_prior_points = 4096){

  if(!inherits(fit, "runjags") || !inherits(fit, "BayesTools_fit")){
    stop("'fit' must be a BayesTools JAGS fit.", call. = FALSE)
  }
  summary <- .bt_random_effect_summary_posterior_type(summary)
  check_char(allocation, "allocation", check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)
  check_char(component, "component", check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)
  check_char(formula_parameter, "formula_parameter", check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)
  check_int(n_prior_points, "n_prior_points", lower = 16)

  prior_list <- attr(fit, "prior_list", exact = TRUE)
  check_list(prior_list, "prior_list")
  if(!all(vapply(prior_list, is.prior, logical(1)))){
    stop("'prior_list' must be a list of priors.", call. = FALSE)
  }

  model_samples <- .extract_posterior_samples(fit, as_list = FALSE)
  summary_mode <- if(identical(summary$summary, "sd_multiplier")) "full" else "standard"
  random_summary <- .bt_random_effect_summary_samples(
    model_samples = model_samples,
    prior_list = prior_list,
    formula_design = attr(fit, "formula_design", exact = TRUE),
    mode = summary_mode
  )

  summary_samples <- random_summary$model_samples
  summary_priors <- random_summary$prior_list
  if(ncol(summary_samples) == 0L || length(summary_priors) == 0L){
    stop("No random-effect summaries are available in 'fit'.", call. = FALSE)
  }

  keep <- vapply(summary_priors, function(prior){
    identical(attr(prior, "random_summary", exact = TRUE), summary$summary)
  }, logical(1))
  if(!is.null(allocation)){
    keep <- keep & vapply(summary_priors, function(prior){
      .bt_random_effect_summary_posterior_attr_matches(
        prior,
        "random_allocation",
        allocation
      )
    }, logical(1))
  }
  if(!is.null(component)){
    keep <- keep & vapply(summary_priors, function(prior){
      .bt_random_effect_summary_posterior_attr_matches(
        prior,
        "random_component",
        component
      )
    }, logical(1))
  }
  if(!is.null(formula_parameter)){
    keep <- keep & vapply(summary_priors, function(prior){
      .bt_random_effect_summary_posterior_attr_matches(
        prior,
        "parameter",
        formula_parameter
      )
    }, logical(1))
  }

  if(!any(keep)){
    .bt_random_effect_summary_posterior_no_match_stop(
      summary = summary,
      allocation = allocation,
      component = component,
      formula_parameter = formula_parameter
    )
  }

  raw_names <- names(summary_priors)[keep]
  display_names <- .bt_random_effect_summary_display_names(
    names = raw_names,
    raw_names = raw_names,
    prior_list = summary_priors,
    formula_prefix = TRUE
  )
  display_names <- make.unique(display_names)

  out <- vector("list", length(raw_names))
  names(out) <- display_names
  out_prior_list <- vector("list", length(raw_names))
  names(out_prior_list) <- display_names

  for(i in seq_along(raw_names)){
    raw_name <- raw_names[i]
    summary_prior <- summary_priors[[raw_name]]
    values <- unname(as.numeric(summary_samples[, raw_name]))
    attr(values, "sample_ind") <- FALSE
    attr(values, "models_ind") <- rep(1, length(values))
    attr(values, "parameter") <- display_names[i]
    attr(values, "summary_name") <- raw_name
    attr(values, "random_summary") <- summary$summary
    attr(values, "random_summary_label") <- attr(summary_prior, "random_summary_label", exact = TRUE)
    attr(values, "random_allocation") <- attr(summary_prior, "random_allocation", exact = TRUE)
    attr(values, "random_component") <- attr(summary_prior, "random_component", exact = TRUE)
    attr(values, "formula_parameter") <- attr(summary_prior, "parameter", exact = TRUE)
    attr(values, "prior_list") <- prior_none()

    prior_density <- .bt_random_effect_summary_posterior_prior_density(
      summary_prior = summary_prior,
      prior_list = prior_list,
      summary_type = summary$summary,
      n_grid = n_prior_points
    )
    if(!is.null(prior_density)){
      attr(values, "prior_density") <- prior_density
      values <- .posterior_support_set(
        values,
        .posterior_support_new(
          bounds = attr(prior_density, "support", exact = TRUE),
          exact = TRUE,
          source = "prior",
          type = "interval"
        )
      )
    }

    class(values) <- c(
      "mixed_posteriors",
      "mixed_posteriors.simple",
      "marginal_posterior.simple",
      "marginal_posterior"
    )
    out[[i]] <- values
    out_prior_list[[i]] <- attr(values, "prior_list", exact = TRUE)
  }

  attr(out, "prior_list") <- out_prior_list
  attr(out, "summary_names") <- raw_names
  class(out) <- c("as_mixed_posteriors", "mixed_posteriors", "list")
  out
}

.bt_random_effect_summary_posterior_attr_matches <- function(prior, attribute,
                                                             values){

  value <- attr(prior, attribute, exact = TRUE)
  if(is.null(value)){
    return(FALSE)
  }

  any(value %in% values)
}

.bt_random_effect_summary_posterior_type <- function(summary){

  summary <- match.arg(
    summary,
    c("variance_ratio", "variance_fraction", "sd_multiplier")
  )

  switch(
    summary,
    "variance_ratio" = list(input = summary, summary = "var_ratio"),
    "variance_fraction" = list(input = summary, summary = "var_frac"),
    "sd_multiplier" = list(input = summary, summary = "sd_multiplier")
  )
}

.bt_random_effect_summary_posterior_no_match_stop <- function(summary,
                                                              allocation,
                                                              component,
                                                              formula_parameter){

  has_filters <- !is.null(allocation) || !is.null(component) ||
    !is.null(formula_parameter)
  filter_text <- if(has_filters) " matched the requested filters" else " are available"

  detail <- switch(
    summary$summary,
    "var_ratio" = paste0(
      "Variance-ratio summaries are created only for ",
      "random_variance_allocation(..., target = \"sd_component\", ",
      "scale = \"mean_variance\"). Total-variance allocations are returned ",
      "by summary = \"variance_fraction\"."
    ),
    "var_frac" = paste0(
      "Variance-fraction summaries are created for true total-variance ",
      "allocations. Mean-variance SD-component allocations are returned by ",
      "summary = \"variance_ratio\"."
    ),
    "sd_multiplier" = paste0(
      "SD-multiplier summaries are available only for SD-component ",
      "variance allocations."
    ),
    "Requested random-effect summaries are not available."
  )

  stop(
    "No random-effect ",
    summary$input,
    " summaries",
    filter_text,
    ". ",
    detail,
    call. = FALSE
  )
}

.bt_random_effect_summary_posterior_prior_density <- function(
    summary_prior,
    prior_list,
    summary_type,
    n_grid){

  allocation <- attr(summary_prior, "random_allocation_metadata", exact = TRUE)
  index <- attr(summary_prior, "random_allocation_index", exact = TRUE)
  if(is.null(allocation) || is.null(index) || length(index) != 1L || is.na(index)){
    return(NULL)
  }

  weight_prior <- prior_list[[allocation$weight_name]]
  if(is.null(weight_prior) ||
     !is.prior.simplex(weight_prior) ||
     !identical(weight_prior$distribution, "dirichlet")){
    return(NULL)
  }

  alpha <- weight_prior$parameters[["alpha"]]
  if(index < 1L || index > length(alpha)){
    return(NULL)
  }

  alpha_i <- alpha[index]
  beta_i <- sum(alpha) - alpha_i
  if(!is.finite(alpha_i) || !is.finite(beta_i) || alpha_i <= 0 || beta_i <= 0){
    return(NULL)
  }

  scale <- .bt_random_effect_summary_posterior_prior_scale(
    allocation = allocation,
    summary_type = summary_type,
    K = length(alpha)
  )

  .bt_random_effect_summary_posterior_scaled_beta_density(
    alpha = alpha_i,
    beta = beta_i,
    scale = scale,
    transform = if(identical(summary_type, "sd_multiplier")) "sqrt" else "linear",
    n_grid = n_grid
  )
}

.bt_random_effect_summary_posterior_prior_scale <- function(allocation,
                                                            summary_type,
                                                            K){

  if(identical(summary_type, "var_ratio")){
    return(.bt_random_effect_summary_allocation_n_targets(allocation, K))
  }

  if(identical(summary_type, "sd_multiplier")){
    allocation_scale <- .bt_random_effect_allocation_scale_metadata(
      allocation,
      context = "Random-effect allocation posterior metadata"
    )
    if(identical(allocation_scale, "mean_variance")){
      return(.bt_random_effect_summary_allocation_n_targets(allocation, K))
    }
  }

  1
}

.bt_random_effect_summary_posterior_scaled_beta_density <- function(alpha, beta,
                                                                    scale,
                                                                    transform,
                                                                    n_grid){

  density_evaluator <- function(x){
    out <- numeric(length(x))
    outside <- !is.finite(x) | x < 0

    if(identical(transform, "sqrt")){
      upper_support <- sqrt(scale)
      outside <- outside | x > upper_support
      interior <- !outside & x > 0 & x < upper_support
      out[interior] <- stats::dbeta(
        x[interior]^2 / scale,
        alpha,
        beta
      ) * (2 * x[interior] / scale)

      at_lower <- !outside & x == 0
      if(any(at_lower)){
        lower_power <- 2 * alpha - 1
        out[at_lower] <- if(lower_power < 0){
          Inf
        }else if(lower_power == 0){
          2 * scale^(-alpha) / beta(alpha, beta)
        }else{
          0
        }
      }

      at_upper <- !outside & x == upper_support
      if(any(at_upper)){
        out[at_upper] <- if(beta < 1){
          Inf
        }else if(beta == 1){
          2 * alpha / upper_support
        }else{
          0
        }
      }
    }else{
      upper_support <- scale
      outside <- outside | x > upper_support
      interior <- !outside & x > 0 & x < upper_support
      out[interior] <- stats::dbeta(
        x[interior] / scale,
        alpha,
        beta
      ) / scale

      at_lower <- !outside & x == 0
      if(any(at_lower)){
        out[at_lower] <- if(alpha < 1){
          Inf
        }else if(alpha == 1){
          beta / scale
        }else{
          0
        }
      }

      at_upper <- !outside & x == upper_support
      if(any(at_upper)){
        out[at_upper] <- if(beta < 1){
          Inf
        }else if(beta == 1){
          alpha / scale
        }else{
          0
        }
      }
    }

    out[outside] <- 0
    out
  }

  cdf_evaluator <- function(x){
    if(identical(transform, "sqrt")){
      out <- stats::pbeta(x^2 / scale, alpha, beta)
      out[x <= 0] <- 0
      out[x >= sqrt(scale)] <- 1
      out
    }else{
      stats::pbeta(x / scale, alpha, beta)
    }
  }

  tail_prob <- .prior_linear_density_tail_prob()
  if(identical(transform, "sqrt")){
    support <- c(0, sqrt(scale))
    lower <- if(alpha < 0.5) sqrt(scale * stats::qbeta(tail_prob / 2, alpha, beta)) else support[1]
    upper <- if(beta < 1) sqrt(scale * stats::qbeta(1 - tail_prob / 2, alpha, beta)) else support[2]
  }else{
    support <- c(0, scale)
    lower <- if(alpha < 1) scale * stats::qbeta(tail_prob / 2, alpha, beta) else support[1]
    upper <- if(beta < 1) scale * stats::qbeta(1 - tail_prob / 2, alpha, beta) else support[2]
  }

  if(!is.finite(lower) || !is.finite(upper) || lower >= upper){
    stop(
      "Could not construct a finite plotting grid for the scaled-Beta prior.",
      call. = FALSE
    )
  }

  x <- seq(lower, upper, length.out = n_grid)
  if(support[1] <= 1 && support[2] >= 1){
    x <- sort(unique(c(x, 1)))
  }

  y <- density_evaluator(x)
  if(anyNA(y) || any(!is.finite(y))){
    stop(
      "The scaled-Beta plotting grid includes a singular support boundary.",
      call. = FALSE
    )
  }

  out <- list(
    density = list(
      x = x,
      y = y,
      mass = 1
    ),
    points = .prior_linear_density_empty_points(),
    n_grid = length(x)
  )
  class(out) <- c("prior_linear_density", "prior_density")
  attr(out, "support") <- support
  attr(out, "density_evaluator") <- density_evaluator
  attr(out, "cdf_evaluator") <- cdf_evaluator
  attr(out, "singular_boundaries") <- support[
    is.infinite(density_evaluator(support))
  ]
  out
}
