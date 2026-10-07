

.hypothesis_BF_compute <- function(quantity, statement, density_method,
                                    allow_unavailable = FALSE) {

  left  <- statement[["left"]]
  right <- statement[["right"]]
  explicit <- isTRUE(statement[["explicit"]])

  if(.hypothesis_sides_point_complement(left, right)){
    if(explicit){
      return(.hypothesis_BF_result_labels(.hypothesis_BF_available(.hypothesis_point_BF(
        quantity       = quantity,
        side           = left,
        density_method = density_method,
        inverse        = FALSE
      ), allow_unavailable), alternative = left, null = right))
    }
    return(.hypothesis_BF_result_labels(.hypothesis_BF_available(.hypothesis_point_BF(
      quantity       = quantity,
      side           = left,
      density_method = density_method,
      inverse        = TRUE
    ), allow_unavailable), alternative = right, null = left))
  }
  if(.hypothesis_sides_point_complement(right, left)){
    return(.hypothesis_BF_result_labels(.hypothesis_BF_available(.hypothesis_point_BF(
      quantity       = quantity,
      side           = right,
      density_method = density_method,
      inverse        = TRUE
    ), allow_unavailable), alternative = left, null = right))
  }

  if(identical(left[["type"]], "region") &&
     identical(right[["type"]], "region")){
    return(.hypothesis_BF_result_labels(
      .hypothesis_BF_available(.hypothesis_region_odds_BF(
        quantity, left, right, explicit = explicit), allow_unavailable),
      alternative = left,
      null        = right
    ))
  }

  if(identical(left[["type"]], "point") &&
     identical(right[["type"]], "region")){
    return(.hypothesis_BF_result_labels(.hypothesis_BF_available(.hypothesis_transitive_BF(
      quantity       = quantity,
      point_side     = left,
      region_side    = right,
      density_method = density_method,
      inverse        = FALSE,
      explicit       = explicit
    ), allow_unavailable), alternative = left, null = right))
  }
  if(identical(left[["type"]], "region") &&
     identical(right[["type"]], "point")){
    return(.hypothesis_BF_result_labels(.hypothesis_BF_available(.hypothesis_transitive_BF(
      quantity       = quantity,
      point_side     = right,
      region_side    = left,
      density_method = density_method,
      inverse        = TRUE,
      explicit       = explicit
    ), allow_unavailable), alternative = left, null = right))
  }

  stop("Unsupported hypothesis comparison.", call. = FALSE)
}

.hypothesis_BF_available <- function(result, allow_unavailable){

  if(!isTRUE(allow_unavailable)){
    return(result)
  }
  tryCatch(result, error = function(condition){
    if(!inherits(condition, c("BayesTools_hypothesis_ordinate",
                              "BayesTools_hypothesis_region",
                              "BayesTools_transformation_image_unavailable"))){
      stop(condition)
    }
    list(BF = NA_real_, BF_error = NA_real_, prior = NA_real_,
         posterior = NA_real_, method = "unavailable",
         warning = conditionMessage(condition), inexact_prior = NULL)
  })
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
    point_side[["value"]] == other_side[["value"]]
}


.hypothesis_point_BF <- function(quantity, side, density_method,
                                 inverse = FALSE) {

  inexact_prior <- NULL
  evaluation_value <- side$value
  log_jacobian <- 0
  marginal <- .hypothesis_point_marginal(quantity, side, density_method)
  linear <- NULL
  if(is.null(marginal)){
    expression <- .hypothesis_side_expression(side)
    symbols <- .hypothesis_expression_symbols(expression)
    linear <- .hypothesis_linear_coefficients(expression, symbols, quantity$posterior_draws)
    if(!is.null(linear)){
      evaluation_value <- .hypothesis_affine_null(linear, side$value)
      log_jacobian <- log(abs(linear$divisor))
    }
  }
  if(is.null(marginal) && density_method %in% c("KDE", "normal")){
    marginal <- .hypothesis_linear_point_marginal(quantity, side)
  }
  if(is.null(marginal) && identical(density_method, "precomputed")){
    if(!is.null(linear) && length(linear$coefficients) == 0L){
      .hypothesis_check_prior_ordinate(.prior_linear_combination_density(
        list(value = prior("point", list(0))), c(value = 1), n_grid = 64L),
        evaluation_value, side$label)
    }
    .hypothesis_draw_density_height(numeric(), evaluation_value, "posterior", density_method)
  }
  if(!is.null(marginal)){
    if(!is.null(marginal$evaluation_value)) evaluation_value <- marginal$evaluation_value
    if(!is.null(marginal$log_jacobian)) log_jacobian <- marginal$log_jacobian
    prior_ordinate <- .hypothesis_check_prior_ordinate(marginal$prior_density,
      evaluation_value, side$label)
    normal_approximation <- identical(density_method, "normal")
    posterior <- .posterior_precomputed_child(
      parent = marginal$posterior_parent, child = marginal$posterior,
      index = marginal$posterior_index, null_hypothesis = evaluation_value,
      density_method = if(normal_approximation) "KDE" else density_method)
    inclusion_BF <- Savage_Dickey_BF(posterior, null_hypothesis = evaluation_value,
      normal_approximation = normal_approximation, silent = TRUE,
      density_method = if(normal_approximation) "KDE" else density_method)
    log_prior <- attr(inclusion_BF, "log_prior_height", exact = TRUE)
    log_posterior <- attr(inclusion_BF, "log_posterior_height", exact = TRUE)
    log_BF <- -attr(inclusion_BF, "log_BF", exact = TRUE)
    numerical_diagnostics <- attr(inclusion_BF, "numerical_diagnostics", exact = TRUE)
    bf_warnings <- attr(inclusion_BF, "warnings", exact = TRUE)
    BF_error <- attr(inclusion_BF, "BF_error_percent", exact = TRUE)
    density_source <- attr(inclusion_BF, "posterior_density_source", exact = TRUE)
    fallback_warnings <- attr(inclusion_BF, "posterior_density_fallback_warnings", exact = TRUE)
    if(length(fallback_warnings) > 0L) warning(paste(unique(fallback_warnings), collapse = " "), call. = FALSE)
    method <- if(identical(density_source, "normal")) "Savage-Dickey (normal)" else
      if(identical(density_method, "precomputed") && identical(density_source, "precomputed"))
        "Savage-Dickey (precomputed)" else "Savage-Dickey"
  }else{
    posterior <- if(is.null(linear)){
      .hypothesis_eval_expression(.hypothesis_side_expression(side), quantity$posterior_draws)
    }else .hypothesis_affine_value(linear, quantity$posterior_draws)
    prior_density <- .hypothesis_expression_prior_density(quantity, side)
    if(!is.null(prior_density)){
      prior_ordinate <- .hypothesis_check_prior_ordinate(prior_density, evaluation_value, side$label)
      log_prior <- prior_ordinate$log_density
    }else if(.hypothesis_quantity_has_prior_structure(quantity)){
      .hypothesis_stop_inexact_ordinate(side$label,
        "the expression is not a linear combination with an exact prior density")
    }else{
      draws <- .hypothesis_prior_draws(quantity)
      prior <- if(is.null(linear)) .hypothesis_eval_expression(.hypothesis_side_expression(side), draws) else
        .hypothesis_affine_value(linear, draws)
      log_prior <- .hypothesis_draw_density_log_height(prior, evaluation_value, "prior", density_method)
      if(!identical(density_method, "normal")) .hypothesis_check_prior_density(exp(log_prior), side$label)
      inexact_prior <- side$label
    }
    log_posterior <- .hypothesis_draw_density_log_height(posterior, evaluation_value, "posterior", density_method)
    ratio <- .hypothesis_log_ratio(log_prior, log_posterior, normal = identical(density_method, "normal"))
    log_BF <- -ratio$log_BF
    numerical_diagnostics <- ratio$diagnostics
    bf_warnings <- NULL
    BF_error <- NA_real_
    method <- if(identical(density_method, "normal")) "Savage-Dickey (normal)" else "kernel Savage-Dickey"
  }
  if(inverse) log_BF <- -log_BF
  log_prior <- as.numeric(log_prior)
  log_posterior <- as.numeric(log_posterior)
  log_prior <- log_prior + log_jacobian
  log_posterior <- log_posterior + log_jacobian
  numerical_diagnostics$log_prior_height <- log_prior
  numerical_diagnostics$log_posterior_height <- log_posterior
  numerical_diagnostics$evaluation_value <- evaluation_value
  numerical_diagnostics$requested_value <- side$value
  numerical_diagnostics$log_jacobian <- log_jacobian
  list(BF = exp(log_BF), log_BF = log_BF, prior = exp(log_prior),
    posterior = exp(log_posterior), log_prior_height = log_prior,
    log_posterior_height = log_posterior, numerical_diagnostics = numerical_diagnostics,
    method = method, BF_error = if(is.null(BF_error)) NA_real_ else as.numeric(BF_error),
    warning = .hypothesis_collapse_warning(bf_warnings), inexact_prior = inexact_prior)
}


.hypothesis_point_marginal <- function(quantity, side, density_method = "KDE") {

  if(!is.null(quantity[["posterior_marginal"]]) &&
     .hypothesis_expression_is_parameter(.hypothesis_side_expression(side),
                                         quantity[["parameter"]])){
    posterior <- quantity$posterior_marginal
    evaluation <- .bt_meta_get(posterior, "hypothesis_evaluation")
    if(!is.null(evaluation) && !identical(density_method, "precomputed")){
      if(!identical(as.numeric(posterior), evaluation$numerator / evaluation$divisor)){
        stop("Affine hypothesis evaluation coordinates do not match the reported draws.", call. = FALSE)
      }
      context <- .bt_meta_get(posterior, "prior_context")
      density <- .prior_density_from_context(context, evaluation$weights,
        output_transformation = "lin",
        output_transformation_arguments = list(a = evaluation$offset, b = 1))
      unshifted <- evaluation$numerator
      class(unshifted) <- c("marginal_posterior.simple", "marginal_posterior")
      unshifted <- .bt_meta_assign(unshifted, list(prior_density = density,
        atoms = .posterior_atoms_new(column_names = quantity$parameter, source = "linear_target"),
        support = .posterior_support_from_prior_context_weights(context, evaluation$weights,
          "lin", list(a = evaluation$offset, b = 1)),
        condition = .bt_meta_get(posterior, "condition")))
      components <- .posterior_components_get(posterior)
      if(!is.null(components)){
        unshifted <- .posterior_components_set(unshifted, .posterior_components_new(
          components$index, lapply(components$supports, .posterior_support_transform,
            transformation = "lin", transformation_arguments = list(a = 0, b = evaluation$divisor)),
          components$keys))
      }
      return(list(posterior = unshifted, prior_density = density,
        evaluation_value = .hypothesis_affine_product(side$value, evaluation$divisor),
        log_jacobian = log(evaluation$divisor)))
    }
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


# A linear expression of marginal posteriors (e.g., mu[b] - mu[a] or
# 2 * mu[a] + mu[b]) is evaluated as the marginal posterior of the linear
# combination: the prior density and exact support of the combination come
# from the joint prior context (structurally fixed levels add their value) and
# the posterior ordinate is the (per-component) boundary-reflected KDE, or the
# normal approximation, of Savage_Dickey_BF(). NULL for nonlinear
# expressions, levels with other point masses, row-varying level weights,
# nonlinear transformed levels, unusable exact supports, or missing joint
# metadata; the point hypothesis then has no exact prior ordinate.
.hypothesis_linear_point_marginal <- function(quantity, side) {

  expr <- .hypothesis_side_expression(side)
  symbols <- unique(.hypothesis_expression_symbols(expr))
  marginals <- .hypothesis_symbol_marginals(quantity, symbols)
  if(is.null(marginals)){
    return(NULL)
  }
  linear <- .hypothesis_linear_coefficients(expr, symbols, quantity[["posterior_draws"]])
  if(is.null(linear)){
    return(NULL)
  }
  active <- names(linear[["coefficients"]])[linear[["coefficients"]] != 0]
  if(length(active) == 0L){
    constant_density <- .prior_linear_combination_density(list(value = prior("point", list(0))),
      c(value = 1), n_grid = 64L)
    .hypothesis_check_prior_ordinate(constant_density,
      .hypothesis_affine_null(linear, side$value), side$label)
  }

  context <- NULL
  weights <- NULL
  space <- NULL
  offset  <- 0
  for(symbol in active){
    level <- marginals[[symbol]]
    .bt_formula_measure_check(level, "prior_density")
    level_space <- .bt_linear_weight_space(level)
    if(!is.null(space) && !identical(space, level_space)){
      .hypothesis_linear_target_stop("Linear target levels use incompatible weight spaces.", "weight_space")
    }
    space <- level_space
    atoms <- .posterior_atoms_get(level)
    if(is.null(atoms)){
      return(NULL)
    }
    if(any(atoms$mass > 0)){
      # a structurally fixed level (a reference level: its prior and
      # posterior are one point) adds a constant to the combination
      fixed <- .hypothesis_fixed_level_value(level, atoms)
      if(is.null(fixed)){
        return(NULL)
      }
      offset <- .hypothesis_affine_sum(offset,
        .hypothesis_affine_product(linear[["coefficients"]][[symbol]], fixed))
      next
    }
    level_context <- .bt_meta_get(level, "prior_context")
    if(!.hypothesis_is_prior_density_context(level_context) ||
       (!is.null(context) && !identical(level_context, context))){
      return(NULL)
    }
    context <- level_context
    level_weights <- .bt_meta_get(level, "linear_weights")
    if(is.null(level_weights) ||
       !is.null(.bt_meta_get(level, "joint_prior_transformation"))){
      return(NULL)
    }
    if(!is.null(dim(level_weights))){
      if(nrow(level_weights) != 1L){
        return(NULL)
      }
      level_weights <- stats::setNames(as.numeric(level_weights[1L, ]), colnames(level_weights))
    }
    if(is.null(names(level_weights))){
      return(NULL)
    }
    coefficient <- linear[["coefficients"]][[symbol]]
    weights <- .hypothesis_add_linear_weights(weights,
      .hypothesis_affine_product(coefficient, level_weights))
    offset <- .hypothesis_affine_sum(offset,
      .hypothesis_affine_product(coefficient, .hypothesis_level_linear_offset(level)))
  }

  shift <- list(a = offset, b = 1)
  if(length(weights) == 0L || all(weights == 0)){
    .hypothesis_check_prior_ordinate(.prior_density_from_context(context, weights,
      output_transformation = "lin", output_transformation_arguments = shift),
      .hypothesis_affine_null(linear, side$value), side$label)
  }
  support <- .posterior_support_from_prior_context_weights(
    context,
    weights,
    output_transformation           = "lin",
    output_transformation_arguments = shift
  )
  if(!is.null(support) && !.hypothesis_linear_support_usable(support)){
    if(identical(support$type, "interval") && length(support$points) == 0L){
      support <- NULL
    }else return(NULL)
  }

  n_draws <- nrow(quantity[["posterior_draws"]])
  components <- .hypothesis_linear_components(marginals[active], context, weights, shift, n_draws)
  component_supports <- if(is.null(components)) list() else
    components$supports[!vapply(components$supports, is.null, logical(1))]
  if(!all(vapply(component_supports, .hypothesis_linear_support_usable, logical(1)))){
    return(NULL)
  }

  prior_density <- .prior_density_from_context(
    context,
    weights,
    output_transformation           = "lin",
    output_transformation_arguments = shift
  )
  posterior <- .hypothesis_affine_value(linear, quantity[["posterior_draws"]])
  class(posterior) <- c("marginal_posterior.simple", "marginal_posterior")
  posterior <- .bt_meta_assign(posterior, list(prior_density = prior_density,
    prior_context = context, linear_weights = weights, linear_weight_space = space,
    linear_offset = offset))
  posterior <- .posterior_support_set(posterior, support)
  posterior <- .posterior_atoms_set(posterior, .posterior_atoms_new(
    column_names = "value",
    source       = "linear_combination",
    declared     = TRUE
  ))
  posterior <- .posterior_components_set(posterior, components)

  return(list(
    posterior        = posterior,
    posterior_parent = NULL,
    posterior_index  = NULL,
    prior_density    = prior_density,
    evaluation_value = .hypothesis_affine_null(linear, side$value),
    log_jacobian     = log(abs(linear$divisor))
  ))
}


# The value of a structurally fixed level: its posterior atom holds all of
# the mass and its prior density is the point mass at the same location.
# NULL otherwise.
.hypothesis_fixed_level_value <- function(level, atoms) {

  if(length(atoms$mass) != 1L || ncol(atoms$locations) != 1L ||
     abs(atoms$mass - 1) > 8 * .Machine$double.eps){
    return(NULL)
  }
  location <- atoms$locations[1L, 1L]
  prior_density <- .bt_meta_get(level, "prior_density")
  if(!inherits(prior_density, "prior_linear_density") ||
     abs(.prior_linear_density_point_mass(prior_density, location) - 1) >
       8 * .Machine$double.eps){
    return(NULL)
  }
  location
}


.hypothesis_symbol_marginals <- function(quantity, symbols) {

  if(length(symbols) == 0L){
    return(NULL)
  }

  out <- list()
  for(symbol in symbols){
    if(!is.null(quantity[["posterior_marginals"]]) &&
       symbol %in% names(quantity[["posterior_marginals"]])){
      out[[symbol]] <- quantity[["posterior_marginals"]][[symbol]]
    }else if(!is.null(quantity[["posterior_marginal"]]) &&
             identical(symbol, quantity[["parameter"]])){
      out[[symbol]] <- quantity[["posterior_marginal"]]
    }else{
      return(NULL)
    }
  }

  out
}


# Coefficients of an expression that is linear in its symbols, from its values
# at the origin, the unit vectors, and two checking points with mixed signs;
# NULL for a nonlinear expression. Probing alone cannot see a kink outside the
# probe points (abs(x - 2) is x - 2 at all of them), so the expression must
# also have a linear form.
.hypothesis_linear_coefficients <- function(expr, symbols, draws) {

  missing <- setdiff(symbols, names(draws))
  if(length(missing)) .hypothesis_stop_unknown_quantity(missing, names(draws))
  .hypothesis_affine_read(expr, symbols)
}

.hypothesis_expression_linear_form <- function(expr) {

  !is.null(.hypothesis_affine_read(expr, .hypothesis_expression_symbols(expr)))
}


.hypothesis_add_linear_weights <- function(total, weights) {

  if(is.null(total)){
    return(weights)
  }

  columns <- union(names(total), names(weights))
  out <- stats::setNames(numeric(length(columns)), columns)
  out[names(total)] <- out[names(total)] + total
  out[names(weights)] <- out[names(weights)] + weights

  out
}


# Exact interval support without point masses.
.hypothesis_linear_support_usable <- function(support) {

  support <- .posterior_support_from_attribute(support)
  !is.null(support) && isTRUE(support$exact) &&
    identical(support$type, "interval") && length(support$points) == 0L
}


# Components of a linear combination of marginal posteriors: the component of a
# draw combines the component keys of its terms (the model of a model-mixture
# ensemble, or the mixture indicators of a single fit), and each component's
# exact support is that of the combination within the component.
.hypothesis_linear_components <- function(levels, context, weights, shift, n_draws) {

  level_components <- lapply(levels, .posterior_components_get)
  level_components <- level_components[!vapply(level_components, is.null, logical(1))]
  if(length(level_components) == 0L){
    return(NULL)
  }

  draw_keys <- list()
  for(components in level_components){
    if(is.null(components$keys) || length(components$index) != n_draws){
      stop("The mixture components of the marginal posterior levels do not match their draws.",
           call. = FALSE)
    }
    keys <- components$keys[components$index, , drop = FALSE]
    for(column in colnames(keys)){
      draw_keys[[column]] <- keys[, column]
    }
  }
  draw_keys <- do.call(cbind, draw_keys)

  if(!".model" %in% colnames(draw_keys)){
    # every mixture prior entering the combination needs its component;
    # without one (the mixture terms cancel), the pooled ordinate applies
    parameters <- .posterior_components_mixture_parameters(context, weights)
    if(length(parameters) == 0L){
      return(NULL)
    }
    if(!all(parameters %in% colnames(draw_keys))){
      stop("The mixture components of the marginal posterior levels do not cover the mixture priors of the linear combination.",
           call. = FALSE)
    }
    draw_keys <- draw_keys[, parameters, drop = FALSE]
  }

  key_text <- do.call(paste, c(as.data.frame(draw_keys), sep = "\r"))
  index <- match(key_text, unique(key_text))
  keys <- draw_keys[!duplicated(key_text), , drop = FALSE]
  rownames(keys) <- NULL

  .posterior_components_new(
    index    = index,
    supports = .posterior_components_supports(
      context                         = context,
      keys                            = keys,
      weights                         = weights,
      output_transformation           = "lin",
      output_transformation_arguments = shift
    ),
    keys     = keys
  )
}


.hypothesis_region_odds_BF <- function(quantity, left, right,
                                       explicit = FALSE) {

  prior_left      <- .hypothesis_region_mass(quantity, left, prior = TRUE)
  prior_right     <- .hypothesis_region_mass(quantity, right, prior = TRUE)
  posterior_left  <- .hypothesis_region_mass(quantity, left, prior = FALSE)
  posterior_right <- .hypothesis_region_mass(quantity, right, prior = FALSE)

  .hypothesis_check_prior_mass(prior_left, left[["label"]],
                               allow_one = explicit)
  .hypothesis_check_prior_mass(prior_right, right[["label"]],
                               allow_one = explicit)

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
                                      density_method, inverse,
                                      explicit = TRUE) {

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
  .hypothesis_check_prior_mass(region_prior, region_side[["label"]],
                               allow_one = explicit)

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
    warning   = warning,
    inexact_prior = point_BF[["inexact_prior"]]
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
  if(prior){
    # Prior masses computed from the prior object's distribution function are
    # exact and contribute no Monte Carlo variance.
    left_exact  <- .hypothesis_region_prior_mass_exact(quantity, left)
    right_exact <- .hypothesis_region_prior_mass_exact(quantity, right)
    if(left_exact && right_exact){
      return(0)
    }
    if(left_exact || right_exact){
      estimated <- if(left_exact) right else left
      return(.hypothesis_log_prob_indicator_mc_var(
        .hypothesis_draw_region_indicator(estimated, draws)
      ))
    }
  }
  left_values  <- .hypothesis_draw_region_indicator(left, draws)
  right_values <- .hypothesis_draw_region_indicator(right, draws)

  .hypothesis_log_odds_indicator_mc_var(left_values, right_values)
}


.hypothesis_region_prior_mass_exact <- function(quantity, side) {

  # Mirrors .hypothesis_region_mass(): the prior object's distribution
  # function is used whenever it defines the region mass.
  !is.null(.hypothesis_prior_object_region_mass(quantity, side))
}


.hypothesis_region_log_mass_mc_var <- function(quantity, side, prior) {

  if(prior && is.null(quantity[["prior_draws"]])){
    return(0)
  }
  if(prior && .hypothesis_region_prior_mass_exact(quantity, side)){
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
