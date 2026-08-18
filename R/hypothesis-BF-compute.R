

.hypothesis_BF_compute <- function(quantity, statement, density_method) {

  left  <- statement[["left"]]
  right <- statement[["right"]]
  explicit <- isTRUE(statement[["explicit"]])

  if(.hypothesis_sides_point_complement(left, right)){
    if(explicit){
      return(.hypothesis_BF_result_labels(.hypothesis_point_BF(
        quantity       = quantity,
        side           = left,
        density_method = density_method,
        inverse        = FALSE
      ), alternative = left, null = right))
    }
    return(.hypothesis_BF_result_labels(.hypothesis_point_BF(
      quantity       = quantity,
      side           = left,
      density_method = density_method,
      inverse        = TRUE
    ), alternative = right, null = left))
  }
  if(.hypothesis_sides_point_complement(right, left)){
    return(.hypothesis_BF_result_labels(.hypothesis_point_BF(
      quantity       = quantity,
      side           = right,
      density_method = density_method,
      inverse        = TRUE
    ), alternative = left, null = right))
  }

  if(identical(left[["type"]], "region") &&
     identical(right[["type"]], "region")){
    return(.hypothesis_BF_result_labels(
      .hypothesis_region_odds_BF(quantity, left, right),
      alternative = left,
      null        = right
    ))
  }

  if(identical(left[["type"]], "point") &&
     identical(right[["type"]], "region")){
    return(.hypothesis_BF_result_labels(.hypothesis_transitive_BF(
      quantity       = quantity,
      point_side     = left,
      region_side    = right,
      density_method = density_method,
      inverse        = FALSE
    ), alternative = left, null = right))
  }
  if(identical(left[["type"]], "region") &&
     identical(right[["type"]], "point")){
    return(.hypothesis_BF_result_labels(.hypothesis_transitive_BF(
      quantity       = quantity,
      point_side     = right,
      region_side    = left,
      density_method = density_method,
      inverse        = TRUE
    ), alternative = left, null = right))
  }

  stop("Unsupported hypothesis comparison.", call. = FALSE)
}


.hypothesis_BF_result_labels <- function(result, alternative, null) {

  result[["alternative"]] <- alternative[["label"]]
  result[["null"]]        <- null[["label"]]

  return(result)
}


.hypothesis_sides_point_complement <- function(point_side, other_side) {

  if(!identical(point_side[["type"]], "point") ||
     !identical(other_side[["type"]], "not_point")){
    return(FALSE)
  }

  identical(
    .hypothesis_expression_key(.hypothesis_side_expression(point_side)),
    .hypothesis_expression_key(.hypothesis_side_expression(other_side))
  ) &&
    isTRUE(all.equal(point_side[["value"]], other_side[["value"]]))
}


.hypothesis_point_BF <- function(quantity, side, density_method,
                                 inverse = FALSE) {

  marginal <- .hypothesis_point_marginal(quantity, side)
  if(!is.null(marginal)){
    normal_approximation <- identical(density_method, "normal")
    posterior <- .posterior_precomputed_child(
      parent          = marginal[["posterior_parent"]],
      child           = marginal[["posterior"]],
      index           = marginal[["posterior_index"]],
      null_hypothesis = side[["value"]],
      density_method  = if(normal_approximation) "KDE" else density_method
    )
    inclusion_BF <- Savage_Dickey_BF(
      posterior            = posterior,
      null_hypothesis      = side[["value"]],
      normal_approximation = normal_approximation,
      silent               = TRUE,
      density_method       = if(normal_approximation) "KDE" else density_method
    )
    prior_value <- .hypothesis_prior_density_height(
      marginal[["prior_density"]],
      side[["value"]]
    )
    .hypothesis_check_prior_density(prior_value, side[["label"]])
    posterior_value <- prior_value / as.numeric(inclusion_BF)
    BF <- posterior_value / prior_value
    bf_warnings <- attr(inclusion_BF, "warnings", exact = TRUE)
    BF_error <- attr(inclusion_BF, "BF_error_percent", exact = TRUE)
    density_source <- attr(inclusion_BF, "posterior_density_source", exact = TRUE)
    fallback_warnings <- attr(inclusion_BF, "posterior_density_fallback_warnings", exact = TRUE)
    if(length(fallback_warnings) > 0L){
      warning(
        paste(unique(fallback_warnings), collapse = " "),
        call. = FALSE
      )
    }
    method <- if(identical(density_source, "normal")){
      "Savage-Dickey (normal)"
    }else if(identical(density_method, "precomputed") &&
             identical(density_source, "precomputed")){
      "Savage-Dickey (precomputed)"
    }else{
      "Savage-Dickey"
    }
  }else{
    posterior <- .hypothesis_eval_expression(
      .hypothesis_side_expression(side),
      quantity[["posterior_draws"]]
    )
    prior_value <- .hypothesis_prior_object_density_height(quantity, side)
    if(is.null(prior_value)){
      if(!is.null(quantity[["prior_density"]]) &&
         .hypothesis_expression_is_parameter(.hypothesis_side_expression(side),
                                             quantity[["parameter"]])){
        prior_value <- .hypothesis_prior_density_height(
          quantity[["prior_density"]],
          side[["value"]]
        )
      }else{
        prior <- .hypothesis_eval_expression(
          .hypothesis_side_expression(side),
          .hypothesis_prior_draws(quantity)
        )
        prior_value <- .hypothesis_draw_density_height(
          prior,
          side[["value"]],
          "prior",
          density_method
        )
      }
    }
    .hypothesis_check_prior_density(prior_value, side[["label"]])
    posterior_value <- .hypothesis_draw_density_height(
      posterior,
      side[["value"]],
      "posterior",
      density_method
    )
    BF <- posterior_value / prior_value
    bf_warnings <- NULL
    BF_error <- NA_real_
    method <- if(identical(density_method, "normal")){
      "Savage-Dickey (normal)"
    }else{
      "kernel Savage-Dickey"
    }
  }

  if(inverse){
    BF <- 1 / as.numeric(BF)
  }

  return(list(
    BF        = as.numeric(BF),
    prior     = prior_value,
    posterior = posterior_value,
    method    = method,
    BF_error  = if(is.null(BF_error)) NA_real_ else as.numeric(BF_error),
    warning   = .hypothesis_collapse_warning(bf_warnings)
  ))
}


.hypothesis_point_marginal <- function(quantity, side) {

  if(!is.null(quantity[["posterior_marginal"]]) &&
     .hypothesis_expression_is_parameter(.hypothesis_side_expression(side),
                                         quantity[["parameter"]])){
    return(list(
      posterior        = quantity[["posterior_marginal"]],
      posterior_parent = quantity[["posterior_marginal_parent"]],
      posterior_index  = quantity[["posterior_marginal_index"]],
      prior_density    = quantity[["prior_density"]]
    ))
  }

  symbol <- .hypothesis_direct_symbol(.hypothesis_side_expression(side))
  if(is.null(symbol) || is.null(quantity[["posterior_marginals"]]) ||
     !symbol %in% names(quantity[["posterior_marginals"]])){
    return(NULL)
  }

  return(list(
    posterior        = quantity[["posterior_marginals"]][[symbol]],
    posterior_parent = quantity[["posterior_marginal_parent"]],
    posterior_index  = quantity[["posterior_marginal_indices"]][[symbol]],
    prior_density    = quantity[["prior_densities"]][[symbol]]
  ))
}


.hypothesis_region_odds_BF <- function(quantity, left, right) {

  prior_left      <- .hypothesis_region_mass(quantity, left, prior = TRUE)
  prior_right     <- .hypothesis_region_mass(quantity, right, prior = TRUE)
  posterior_left  <- .hypothesis_region_mass(quantity, left, prior = FALSE)
  posterior_right <- .hypothesis_region_mass(quantity, right, prior = FALSE)

  .hypothesis_check_prior_mass(prior_left, left[["label"]])
  .hypothesis_check_prior_mass(prior_right, right[["label"]])

  if(posterior_left == 0 && posterior_right == 0){
    return(list(
      BF        = NA_real_,
      prior     = prior_left / prior_right,
      posterior = NA_real_,
      method    = "prior-posterior odds",
      BF_error  = NA_real_,
      warning   = "Both posterior region masses are zero; region-odds BF is undefined."
    ))
  }

  BF       <- (posterior_left / posterior_right) / (prior_left / prior_right)
  BF_error <- .hypothesis_region_odds_BF_error_percent(quantity, left, right)
  warning  <- NULL
  if(posterior_left == 0 || posterior_right == 0){
    warning <- "Posterior region mass is zero; reported BF is boundary-valued."
  }

  return(list(
    BF        = BF,
    prior     = prior_left / prior_right,
    posterior = posterior_left / posterior_right,
    method    = "prior-posterior odds",
    BF_error  = BF_error,
    warning   = warning
  ))
}


.hypothesis_transitive_BF <- function(quantity, point_side, region_side,
                                      density_method, inverse) {

  if(!.hypothesis_point_region_compatible(point_side, region_side)){
    stop("Point-vs-region hypotheses must use the same scalar expression.",
         call. = FALSE)
  }

  point_BF <- .hypothesis_point_BF(
    quantity       = quantity,
    side           = point_side,
    density_method = density_method,
    inverse        = FALSE
  )
  region_prior     <- .hypothesis_region_mass(quantity, region_side, prior = TRUE)
  region_posterior <- .hypothesis_region_mass(quantity, region_side, prior = FALSE)
  .hypothesis_check_prior_mass(region_prior, region_side[["label"]])

  region_BF       <- region_posterior / region_prior
  region_BF_error <- .hypothesis_region_BF_error_percent(quantity, region_side)
  BF              <- point_BF[["BF"]] / region_BF
  if(inverse){
    BF <- 1 / BF
  }

  warning <- point_BF[["warning"]]
  if(region_posterior == 0){
    warning <- .hypothesis_collapse_warning(c(
      warning,
      "Posterior region mass is zero; reported BF is boundary-valued."
    ))
  }

  return(list(
    BF        = BF,
    prior     = NA_real_,
    posterior = NA_real_,
    method    = "transitive Savage-Dickey",
    BF_error  = .hypothesis_combine_BF_error_percent(
      point_BF[["BF_error"]],
      region_BF_error
    ),
    warning   = warning
  ))
}


.hypothesis_region_odds_BF_error_percent <- function(quantity, left, right) {

  posterior_var <- .hypothesis_region_log_odds_mc_var(
    quantity = quantity,
    left     = left,
    right    = right,
    prior    = FALSE
  )
  prior_var     <- .hypothesis_region_log_odds_mc_var(
    quantity = quantity,
    left     = left,
    right    = right,
    prior    = TRUE
  )

  .hypothesis_log_BF_error_percent(c(posterior_var, prior_var))
}


.hypothesis_region_BF_error_percent <- function(quantity, side) {

  posterior_var <- .hypothesis_region_log_mass_mc_var(
    quantity = quantity,
    side     = side,
    prior    = FALSE
  )
  prior_var     <- .hypothesis_region_log_mass_mc_var(
    quantity = quantity,
    side     = side,
    prior    = TRUE
  )

  .hypothesis_log_BF_error_percent(c(posterior_var, prior_var))
}


.hypothesis_region_log_odds_mc_var <- function(quantity, left, right, prior) {

  if(prior && is.null(quantity[["prior_draws"]])){
    return(0)
  }

  draws <- if(prior){
    quantity[["prior_draws"]]
  }else{
    quantity[["posterior_draws"]]
  }
  left_values  <- .hypothesis_draw_region_indicator(left, draws)
  right_values <- .hypothesis_draw_region_indicator(right, draws)

  .hypothesis_log_odds_indicator_mc_var(left_values, right_values)
}


.hypothesis_region_log_mass_mc_var <- function(quantity, side, prior) {

  if(prior && is.null(quantity[["prior_draws"]])){
    return(0)
  }

  draws <- if(prior){
    quantity[["prior_draws"]]
  }else{
    quantity[["posterior_draws"]]
  }
  values <- .hypothesis_draw_region_indicator(side, draws)

  .hypothesis_log_prob_indicator_mc_var(values)
}


# NOTE: These MC variances intentionally use the iid delta-method
# approximation because hypothesis_BF() currently receives flattened draws.
# Revisit after the draw backend moves to posterior objects so chain/iteration
# metadata can support autocorrelation-aware SEs for the derived indicators.
.hypothesis_log_odds_indicator_mc_var <- function(left, right) {

  n <- length(left)
  if(n <= 1L || length(right) != n){
    return(NA_real_)
  }

  left    <- as.numeric(left)
  right   <- as.numeric(right)
  p_left  <- mean(left)
  p_right <- mean(right)
  if(!is.finite(p_left) || !is.finite(p_right) ||
     p_left <= 0 || p_right <= 0){
    return(NA_real_)
  }

  log_var <- stats::var(left) / (n * p_left^2) +
    stats::var(right) / (n * p_right^2) -
    2 * stats::cov(left, right) / (n * p_left * p_right)

  .hypothesis_normalize_log_BF_var(log_var)
}


.hypothesis_log_prob_indicator_mc_var <- function(values) {

  n <- length(values)
  if(n <= 1L){
    return(NA_real_)
  }

  values <- as.numeric(values)
  p      <- mean(values)
  if(!is.finite(p) || p <= 0){
    return(NA_real_)
  }

  log_var <- stats::var(values) / (n * p^2)

  .hypothesis_normalize_log_BF_var(log_var)
}


.hypothesis_normalize_log_BF_var <- function(log_var) {

  if(!is.finite(log_var)){
    return(NA_real_)
  }
  if(log_var < 0 && log_var > -sqrt(.Machine$double.eps)){
    log_var <- 0
  }
  if(log_var < 0){
    return(NA_real_)
  }

  return(log_var)
}


.hypothesis_log_BF_error_percent <- function(log_var) {

  log_var <- as.numeric(log_var)
  if(any(!is.finite(log_var) | log_var < 0)){
    return(NA_real_)
  }

  total <- sum(log_var)
  if(!is.finite(total) || total < 0){
    return(NA_real_)
  }

  return(100 * sqrt(total))
}


.hypothesis_combine_BF_error_percent <- function(...) {

  BF_error <- as.numeric(unlist(list(...), use.names = FALSE))
  if(length(BF_error) == 0L ||
     any(!is.finite(BF_error) | BF_error < 0)){
    return(NA_real_)
  }

  return(100 * sqrt(sum((BF_error / 100)^2)))
}


.hypothesis_point_region_compatible <- function(point_side, region_side) {

  point_key <- .hypothesis_expression_key(
    .hypothesis_side_expression(point_side)
  )
  region_exprs <- .hypothesis_region_scalar_expressions(
    .hypothesis_side_expression(region_side)
  )
  if(length(region_exprs) == 0L || anyNA(region_exprs)){
    return(FALSE)
  }

  identical(unique(region_exprs), point_key)
}
