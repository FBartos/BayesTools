#' Extract random-effect summary posterior distributions
#'
#' @description Extracts posterior draws for derived random-effect summary
#' quantities and attaches analytic marginal prior densities where available.
#' For allocation summaries, the helper keeps raw Dirichlet allocation weights
#' internal and exposes interpretable scalar summaries such as mean-variance
#' SD-component variance multipliers. Gated total-variance proportions are
#' normalized over active components and omit draws on which every component
#' is excluded because a variance share is undefined there.
#'
#' @param fit model fit created by [JAGS_fit].
#' @param summary semantic quantity to extract: `"var_mult"`, `"var_prop"`,
#'   `"sd_mult"`, `"sd_total"`, `"var_total"`, `"sd_common"`, or
#'   `"var_common"`.
#' @param allocation optional allocation label filter.
#' @param component optional allocation component label filter.
#' @param formula_parameter optional formula parameter filter.
#' @param simplify_names whether returned quantities should use centrally
#'   generated simplified display labels. Defaults to `FALSE`.
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
    summary = c(
      "var_mult", "var_prop", "sd_mult",
      "sd_total", "var_total", "sd_common", "var_common"
    ),
    allocation = NULL,
    component = NULL,
    formula_parameter = NULL,
    simplify_names = FALSE,
    n_prior_points = 4096){

  if(!inherits(fit, "runjags") || !inherits(fit, "BayesTools_fit")){
    stop("'fit' must be a BayesTools JAGS fit.", call. = FALSE)
  }
  summary <- .bt_random_effect_summary_posterior_type(summary)
  check_char(allocation, "allocation", check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)
  check_char(component, "component", check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)
  check_char(formula_parameter, "formula_parameter", check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)
  check_bool(simplify_names, "simplify_names", allow_NA = FALSE)
  check_int(n_prior_points, "n_prior_points", lower = 16)

  prior_list <- attr(fit, "prior_list", exact = TRUE)
  check_list(prior_list, "prior_list")
  if(!all(vapply(prior_list, is.prior, logical(1)))){
    stop("'prior_list' must be a list of priors.", call. = FALSE)
  }

  catalog <- parameter_catalog(fit)
  quantities <- catalog$quantities
  keep <- startsWith(quantities$role, "random_") &
    !quantities$internal & quantities$status != "unavailable" &
    quantities$quantity == summary$summary
  keys <- quantities$extraction_key
  if(!is.null(allocation)){
    keep <- keep & vapply(keys, function(key){
      !is.null(key$allocation_label) && key$allocation_label %in% allocation
    }, logical(1))
  }
  if(!is.null(component)){
    keep <- keep & quantities$component %in% component
  }
  if(!is.null(formula_parameter)){
    keep <- keep & quantities$formula_parameter %in% formula_parameter
  }
  quantities <- quantities[keep, , drop = FALSE]
  if(nrow(quantities) == 0L && !any(startsWith(
    catalog$quantities$role,
    "random_"
  ))){
    stop("No random-effect summaries are available in 'fit'.", call. = FALSE)
  }
  if(nrow(quantities) == 0L){
    .bt_random_effect_summary_posterior_no_match_stop(
      summary = summary,
      allocation = allocation,
      component = component,
      formula_parameter = formula_parameter
    )
  }

  display_names <- if(simplify_names){
    quantities$display_label
  }else{
    quantities$canonical_name
  }

  out <- vector("list", nrow(quantities))
  names(out) <- display_names
  out_prior_list <- vector("list", nrow(quantities))
  names(out_prior_list) <- display_names

  for(i in seq_len(nrow(quantities))){
    quantity <- quantities[i, , drop = FALSE]
    key      <- quantity$extraction_key[[1L]]
    draws    <- .bt_parameter_draws_from_quantities(fit, quantity)
    values   <- unname(as.numeric(as.matrix(draws)[, 1L]))
    if(identical(quantity$quantity, "var_prop")){
      values <- values[!is.na(values)]
      if(length(values) == 0L){
        stop(
          "The selected variance proportion is undefined because no posterior draw has positive realized allocation variance.",
          call. = FALSE
        )
      }
    }
    attr(values, "sample_ind") <- FALSE
    attr(values, "models_ind") <- rep(1, length(values))
    attr(values, "parameter") <- display_names[i]
    attr(values, "summary_name") <- key$summary_name
    attr(values, "random_summary") <- summary$summary
    attr(values, "random_summary_label") <- display_names[i]
    attr(values, "random_allocation") <- key$allocation_label
    attr(values, "random_component") <- quantity$component
    attr(values, "formula_parameter") <- quantity$formula_parameter
    attr(values, "prior_list") <- prior_none()

    prior_density <- .bt_random_effect_summary_posterior_prior_density(
      fit      = fit,
      quantity = quantity,
      n_grid   = n_prior_points
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
  attr(out, "summary_names") <- vapply(
    quantities$extraction_key,
    `[[`,
    character(1),
    "summary_name"
  )
  class(out) <- c("as_mixed_posteriors", "mixed_posteriors", "list")
  out
}

.bt_random_effect_summary_posterior_type <- function(summary){

  summary <- match.arg(
    summary,
    c(
      "var_mult", "var_prop", "sd_mult",
      "sd_total", "var_total", "sd_common", "var_common"
    )
  )
  list(input = summary, summary = summary)
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
    "var_mult" = paste0(
      "Variance-multiplier summaries are created only for ",
      "random_variance_allocation(..., target = \"sd_component\", ",
      "scale = \"mean_variance\"). Total-variance allocations are returned ",
      "by summary = \"var_prop\"."
    ),
    "var_prop" = paste0(
      "Variance-proportion summaries are created for true total-variance ",
      "allocations. Mean-variance SD-component allocations are returned by ",
      "summary = \"var_mult\"."
    ),
    "sd_mult" = paste0(
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
    fit,
    quantity,
    n_grid){

  key <- quantity$extraction_key[[1L]]
  if(is.null(key$allocation_label) || is.null(key$index)){
    return(NULL)
  }
  random_term <- if(nzchar(key$random_block)){
    .bt_parameter_catalog_find_random_term(fit, key)
  }else{
    NULL
  }
  allocation <- .bt_parameter_catalog_find_allocation(fit, key, random_term)
  index <- key$index
  if(identical(quantity$quantity, "sd_mult") &&
     index > allocation$n_targets){
    index <- index - allocation$n_targets
  }
  if(is.null(allocation) || is.null(index) || length(index) != 1L || is.na(index)){
    return(NULL)
  }

  prior_list <- attr(fit, "prior_list", exact = TRUE)
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

  semantic_transform <- .bt_parameter_transform_from_quantity(fit, quantity)
  density_transform <- if(identical(
    semantic_transform,
    list(type = "identity")
  )){
    list(scale = 1, type = "linear")
  }else if(is.list(semantic_transform) &&
           identical(semantic_transform$type, "affine") &&
           identical(semantic_transform$offset, 0) &&
           semantic_transform$scale > 0){
    list(scale = semantic_transform$scale, type = "linear")
  }else if(is.list(semantic_transform) &&
           identical(semantic_transform$type, "sqrt_scale")){
    list(scale = semantic_transform$scale, type = "sqrt")
  }else{
    return(NULL)
  }

  .bt_random_effect_summary_posterior_scaled_beta_density(
    alpha = alpha_i,
    beta = beta_i,
    scale = density_transform$scale,
    transform = density_transform$type,
    n_grid = n_grid
  )
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
