

.hypothesis_expression_key <- function(text) {

  .hypothesis_expression_ast_key(.hypothesis_parse_expression(text))
}

.hypothesis_expression_ast_key <- function(expr){

  canonicalize <- function(node){
    if(identical(.hypothesis_call_name(node), "(")){
      return(canonicalize(node[[2L]]))
    }
    if(!is.call(node)){
      return(node)
    }
    as.call(c(
      list(node[[1L]]),
      lapply(as.list(node[-1L]), canonicalize)
    ))
  }
  .hypothesis_expression_text(canonicalize(expr))
}


.hypothesis_side_expression <- function(side) {

  .bt_hypothesis_node_language(side[["expression"]])
}


.hypothesis_region_scalar_expressions <- function(expr) {

  if(!is.call(expr)){
    return(NA_character_)
  }

  fun <- .hypothesis_call_name(expr)
  if(is.null(fun)){
    return(NA_character_)
  }
  if(fun == "("){
    return(.hypothesis_region_scalar_expressions(expr[[2L]]))
  }
  if(fun %in% c("&", "|")){
    return(unlist(lapply(as.list(expr[-1L]),
                         .hypothesis_region_scalar_expressions),
                  use.names = FALSE))
  }
  if(fun == "!"){
    return(.hypothesis_region_scalar_expressions(expr[[2L]]))
  }
  if(fun %in% c("<", "<=", ">", ">=")){
    lhs <- expr[[2L]]
    rhs <- expr[[3L]]
    lhs_symbols <- .hypothesis_expression_symbols(lhs)
    rhs_symbols <- .hypothesis_expression_symbols(rhs)

    if(length(lhs_symbols) > 0L && length(rhs_symbols) == 0L){
      return(.hypothesis_expression_ast_key(lhs))
    }
    if(length(lhs_symbols) == 0L && length(rhs_symbols) > 0L){
      return(.hypothesis_expression_ast_key(rhs))
    }
  }

  return(NA_character_)
}


.hypothesis_direct_symbol <- function(expr_text) {

  expr <- .hypothesis_parse_expression(expr_text)
  if(is.name(expr)){
    return(.hypothesis_decode_escaped_constant(as.character(expr)))
  }

  return(NULL)
}


.hypothesis_region_mass <- function(quantity, side, prior) {

  if(prior){
    prior_object_mass <- .hypothesis_prior_object_region_mass(quantity, side)
    if(!is.null(prior_object_mass)){
      return(prior_object_mass)
    }
  }

  if(prior && is.null(quantity[["prior_draws"]])){
    if(is.null(quantity[["prior_density"]])){
      stop("Prior information is required for region hypotheses.",
           call. = FALSE)
    }
    mass <- .hypothesis_prior_density_prob(
      prior_density = quantity[["prior_density"]],
      side          = side,
      parameter     = quantity[["parameter"]]
    )
  }else{
    draws <- if(prior){
      quantity[["prior_draws"]]
    }else{
      quantity[["posterior_draws"]]
    }
    mass <- .hypothesis_draw_region_mass(side, draws)
  }

  return(mass)
}


.hypothesis_prior_object_region_mass <- function(quantity, side) {

  prior_object <- quantity[["prior_object"]]
  if(is.null(prior_object) || !is.prior.simple(prior_object) ||
     is.prior.point(prior_object) || is.prior.discrete(prior_object)){
    return(NULL)
  }

  comparison <- .hypothesis_simple_parameter_comparison(
    side      = side,
    parameter = quantity[["parameter"]]
  )
  if(is.null(comparison)){
    return(NULL)
  }

  mass <- tryCatch(
    switch(
      comparison[["operator"]],
      "<"  = cdf(prior_object, comparison[["value"]]),
      "<=" = cdf(prior_object, comparison[["value"]]),
      ">"  = ccdf(prior_object, comparison[["value"]]),
      ">=" = ccdf(prior_object, comparison[["value"]])
    ),
    error = function(e) NULL
  )
  if(is.null(mass) || length(mass) != 1L || !is.finite(mass)){
    return(NULL)
  }

  mass <- max(0, min(1, as.numeric(mass)))
  mass
}


.hypothesis_prior_object_density_height <- function(quantity, side) {

  prior_object <- quantity[["prior_object"]]
  if(is.null(prior_object) || !is.prior.simple(prior_object) ||
     is.prior.point(prior_object) || is.prior.discrete(prior_object) ||
     !.hypothesis_expression_is_parameter(.hypothesis_side_expression(side),
                                          quantity[["parameter"]])){
    return(NULL)
  }

  height <- tryCatch(
    pdf(prior_object, side[["value"]]),
    error = function(e) NULL
  )
  if(is.null(height) || !is.numeric(height) || length(height) != 1L){
    return(NULL)
  }

  as.numeric(height)
}


.hypothesis_simple_parameter_comparison <- function(side, parameter) {

  if(is.null(parameter)){
    return(NULL)
  }

  expr <- .hypothesis_side_expression(side)
  expr <- .hypothesis_unwrap_parentheses(expr)
  if(!is.call(expr)){
    return(NULL)
  }
  op <- as.character(expr[[1L]])
  if(!op %in% c("<", "<=", ">", ">=")){
    return(NULL)
  }

  lhs <- expr[[2L]]
  rhs <- expr[[3L]]
  lhs_symbols <- .hypothesis_expression_symbols(lhs)
  rhs_symbols <- .hypothesis_expression_symbols(rhs)

  if(length(lhs_symbols) == 1L && identical(lhs_symbols, parameter) &&
     length(rhs_symbols) == 0L && is.name(lhs)){
    return(list(
      operator = op,
      value    = .hypothesis_parse_number(paste(deparse(rhs, width.cutoff = 500L),
                                                collapse = ""))
    ))
  }

  if(length(rhs_symbols) == 1L && identical(rhs_symbols, parameter) &&
     length(lhs_symbols) == 0L && is.name(rhs)){
    return(list(
      operator = switch(op, "<" = ">", "<=" = ">=", ">" = "<", ">=" = "<="),
      value    = .hypothesis_parse_number(paste(deparse(lhs, width.cutoff = 500L),
                                                collapse = ""))
    ))
  }

  NULL
}


.hypothesis_draw_region_mass <- function(side, draws) {

  values <- .hypothesis_draw_region_indicator(side, draws)

  return(mean(values))
}


.hypothesis_draw_region_indicator <- function(side, draws) {

  values <- .hypothesis_eval_condition(
    .hypothesis_side_expression(side),
    draws
  )

  return(values)
}


.hypothesis_prior_density_prob <- function(prior_density, side, parameter) {

  if(is.null(parameter)){
    stop("A quantity name is required for prior-density region tests.",
         call. = FALSE)
  }
  if(!inherits(prior_density, "prior_linear_density")){
    stop("Prior density is not a deterministic scalar prior density.",
         call. = FALSE)
  }

  evaluate_probability <- function(density_object){
    prob <- 0
    if(!is.null(density_object[["density"]])){
      density <- density_object[["density"]]
      x <- density[["x"]]
      y <- density[["y"]] * density[["mass"]]
      draws <- data.frame(x, check.names = FALSE)
      names(draws) <- parameter
      inside <- .hypothesis_eval_condition(
        .hypothesis_side_expression(side),
        draws
      )
      if(length(x) > 1L){
        prob <- prob + .hypothesis_trapz(x, y * as.numeric(inside))
      }
    }

    points <- density_object[["points"]]
    if(!is.null(points) && nrow(points) > 0L){
      draws <- data.frame(points[["x"]], check.names = FALSE)
      names(draws) <- parameter
      inside <- .hypothesis_eval_condition(
        .hypothesis_side_expression(side),
        draws
      )
      prob <- prob + sum(points[["p"]][inside])
    }
    prob
  }

  prob <- evaluate_probability(prior_density)
  refinements <- .prior_linear_density_refinements(prior_density)
  if(length(refinements) > 0L){
    tolerance <- .prior_linear_density_refinement_tolerance()
    previous <- prob
    converged <- FALSE
    for(i in seq_along(refinements)){
      current <- evaluate_probability(refinements[[i]])
      change <- abs(current - previous)
      bound <- tolerance$absolute +
        tolerance$relative * max(abs(current), abs(previous))
      if(is.finite(current) && change <= bound){
        prob <- current
        converged <- TRUE
        break
      }
      previous <- current
    }
    if(!converged){
      stop(
        "Adaptive prior-probability evaluation did not converge within the ",
        "documented grid-refinement error criterion.",
        call. = FALSE
      )
    }
  }

  probability_bound <- 16 * .Machine$double.eps * max(1, abs(prob))
  if(!is.finite(prob) || prob < -probability_bound ||
     prob > 1 + probability_bound){
    stop("Computed prior probability lies materially outside [0, 1].",
         call. = FALSE)
  }
  prob <- max(0, min(1, prob))

  return(prob)
}

.hypothesis_prior_draws <- function(quantity) {

  if(is.null(quantity[["prior_draws"]])){
    stop("Prior draws are required for this hypothesis expression.",
         call. = FALSE)
  }

  return(quantity[["prior_draws"]])
}


.hypothesis_eval_expression <- function(text, draws) {

  expr <- .hypothesis_parse_expression(text)
  .hypothesis_validate_expression(expr, condition = FALSE)

  draws <- as.data.frame(draws, check.names = FALSE)
  missing <- setdiff(.hypothesis_expression_symbols(expr), names(draws))
  if(length(missing) > 0L){
    stop("Hypothesis expression references unknown quantity '",
         paste(missing, collapse = "', '"), "'.", call. = FALSE)
  }

  env <- .hypothesis_draw_environment(draws)
  values <- eval(expr, envir = env)
  check_real(values, "hypothesis expression", check_length = 0,
             allow_NA = FALSE)
  if(any(!is.finite(values))){
    stop("Hypothesis expression produced non-finite values.", call. = FALSE)
  }

  return(as.numeric(values))
}


.hypothesis_eval_condition <- function(text, draws) {

  expr <- .hypothesis_parse_expression(text)
  .hypothesis_validate_expression(expr, condition = TRUE)

  draws <- as.data.frame(draws, check.names = FALSE)
  missing <- setdiff(.hypothesis_expression_symbols(expr), names(draws))
  if(length(missing) > 0L){
    stop("Hypothesis expression references unknown quantity '",
         paste(missing, collapse = "', '"), "'.", call. = FALSE)
  }

  env <- .hypothesis_draw_environment(draws)
  values <- eval(expr, envir = env)
  if(!is.logical(values)){
    stop("Region hypothesis must evaluate to logical values.", call. = FALSE)
  }
  if(anyNA(values)){
    stop("Region hypothesis produced missing values.", call. = FALSE)
  }

  return(values)
}

.hypothesis_draw_environment <- function(draws){

  values <- as.list(draws)
  for(name in intersect(
    names(values),
    .hypothesis_escaped_constant_names()
  )){
    values[[.hypothesis_escaped_constant_symbol(name)]] <- values[[name]]
  }
  list2env(values, parent = .hypothesis_eval_parent())
}


.hypothesis_expression_is_parameter <- function(expr_text, parameter) {

  identical(.hypothesis_direct_symbol(expr_text), parameter)
}


.hypothesis_sample_density_height <- function(samples, value, label) {

  sample_sd <- stats::sd(samples)
  if(length(samples) < 2L || !is.finite(sample_sd) || sample_sd <= 0){
    stop("Cannot estimate ", label, " density from degenerate samples.",
         call. = FALSE)
  }
  .hypothesis_warn_point_draw_cluster(samples, value, label)

  if(value < min(samples) || value > max(samples)){
    warning(
      "The ", label, " samples do not span the point hypothesis. The Gaussian ",
      "KDE height is estimated from kernel tails.",
      call. = FALSE,
      immediate. = TRUE
    )
  }

  density <- stats::density(samples)
  height <- stats::approx(
    density[["x"]],
    density[["y"]],
    xout = value
  )[["y"]]
  if(is.na(height)){
    height <- .density_kde_gaussian_height(
      x     = samples,
      value = value,
      bw    = density[["bw"]]
    )
  }

  if(!is.finite(height) || height < 0){
    stop("Could not estimate ", label, " density at the point hypothesis.",
         call. = FALSE)
  }

  return(height)
}


.hypothesis_warn_point_draw_cluster <- function(samples, value, label,
                                               threshold = 0.01) {

  point_matches <- samples == value
  point_matches[is.na(point_matches)] <- FALSE
  point_share <- mean(point_matches)
  if(is.finite(point_share) && point_share > threshold){
    warning(
      "More than 1% of ", label,
      " draws exactly match the point hypothesis value. Raw-draw ",
      "point-null Bayes factors use continuous-density approximations and ",
      "may be unreliable for spike-and-slab or other point-mass draws.",
      call. = FALSE
    )
  }

  invisible(point_share)
}


.hypothesis_draw_density_height <- function(samples, value, label,
                                            density_method) {

  if(identical(density_method, "precomputed")){
    stop(
      "'density_method = \"precomputed\"' requires valid posterior density ",
      "metadata and cannot be used with raw draws or compound expressions.",
      call. = FALSE
    )
  }

  if(identical(density_method, "normal")){
    return(.hypothesis_normal_density_height(samples, value, label))
  }

  return(.hypothesis_sample_density_height(samples, value, label))
}


.hypothesis_normal_density_height <- function(samples, value, label) {

  sample_sd <- stats::sd(samples)
  if(length(samples) < 2L || !is.finite(sample_sd) || sample_sd <= 0){
    stop("Cannot estimate ", label, " normal density from degenerate samples.",
         call. = FALSE)
  }
  .hypothesis_warn_point_draw_cluster(samples, value, label)

  height <- stats::dnorm(value, mean = mean(samples), sd = sample_sd)
  if(!is.finite(height) || height < 0){
    stop("Could not estimate ", label,
         " normal density at the point hypothesis.", call. = FALSE)
  }

  return(height)
}


.hypothesis_prior_density_height <- function(prior_density, value) {

  if(is.null(prior_density)){
    stop("Prior density is required for point hypotheses.", call. = FALSE)
  }

  if(.prior_linear_density_point_mass(prior_density, value) > 0){
    stop(
      "There is a point mass in the prior at the exact null hypothesis value. The Savage-Dickey density ratio is invalid.",
      call. = FALSE
    )
  }

  .prior_linear_density_height(prior_density, value)
}


.hypothesis_check_prior_mass <- function(mass, label) {

  if(!is.finite(mass) || mass <= 0 || mass >= 1){
    stop("Prior region mass for hypothesis '", label,
         "' is zero, one, or non-finite.", call. = FALSE)
  }

  return(invisible(TRUE))
}

.hypothesis_check_prior_density <- function(density, label) {

  if(!is.finite(density) || density <= 0){
    stop("Prior density at point hypothesis '", label,
         "' is zero or non-finite.", call. = FALSE)
  }

  return(invisible(TRUE))
}
