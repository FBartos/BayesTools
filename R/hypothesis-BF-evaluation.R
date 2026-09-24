

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
      value    = .hypothesis_parse_number(rhs)
    ))
  }

  if(length(rhs_symbols) == 1L && identical(rhs_symbols, parameter) &&
     length(lhs_symbols) == 0L && is.name(rhs)){
    return(list(
      operator = switch(op, "<" = ">", "<=" = ">=", ">" = "<", ">=" = "<="),
      value    = .hypothesis_parse_number(lhs)
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

  comparison <- .hypothesis_simple_parameter_comparison(side, parameter)
  condition  <- .hypothesis_side_expression(side)

  # Densities with a structural representation (point masses, scalar and
  # normal distribution functions, conditional-normal quadratures, and finite
  # mixtures of them) give the region probability directly; the grid below
  # is used only for the other combinations and for regions that are not
  # unions of intervals.
  region <- .hypothesis_prior_region(condition, parameter)
  exact <- if(is.null(region)){
    NULL
  }else{
    .prior_linear_density_region_probability(prior_density, region)
  }
  prob <- if(is.null(exact)){
    .hypothesis_prior_density_grid_prob(prior_density, comparison, condition,
                                        parameter)
  }else{
    as.numeric(exact)
  }

  # a quadrature total may exceed [0, 1] by at most its absolute error
  numerical_error <- if(is.null(exact)){
    0
  }else{
    attr(exact, "numerical_diagnostics", exact = TRUE)$absolute_error
  }
  probability_bound <- max(16 * .Machine$double.eps * max(1, abs(prob)),
                           numerical_error)
  if(!is.finite(prob) || prob < -probability_bound ||
     prob > 1 + probability_bound){
    stop("Computed prior probability lies materially outside [0, 1].",
         call. = FALSE)
  }
  prob <- max(0, min(1, prob))

  return(prob)
}


# Region probability of one density grid: the continuous part is the linear
# interpolant of the grid ordinates, and point masses use the exact condition.
.hypothesis_prior_grid_probability <- function(density_object, comparison,
                                               condition, parameter) {

  prob <- 0
  if(!is.null(density_object[["density"]])){
    density <- density_object[["density"]]
    x <- density[["x"]]
    y <- density[["y"]]
    # Region and total masses both integrate the interpolant exactly
    # (trapezoid rule), so their ratio is independent of the grid's
    # Riemann-sum normalization.
    total <- .hypothesis_trapz(x, y)
    if(length(x) < 2L || !is.finite(total) || total <= 0){
      stop("The continuous prior density grid has no positive integrable ",
           "mass; the prior region probability is unavailable.",
           call. = FALSE)
    }
    inside <- if(!is.null(comparison)){
      .hypothesis_comparison_grid_integral(x, y, comparison)
    }else{
      .hypothesis_condition_grid_integral(x, y, condition, parameter)
    }
    prob <- prob + density[["mass"]] * inside / total
  }

  points <- density_object[["points"]]
  if(!is.null(points) && nrow(points) > 0L){
    inside <- .hypothesis_condition_indicator(
      condition, parameter, points[["x"]]
    )
    prob <- prob + sum(points[["p"]][inside])
  }
  prob
}

# Region probability from the prior-density grid of a density with
# provenance, refined until the documented criterion is met.
.hypothesis_prior_density_grid_prob <- function(prior_density, comparison,
                                                condition, parameter) {

  evaluate_probability <- function(density_object){
    .hypothesis_prior_grid_probability(
      density_object, comparison, condition, parameter
    )
  }

  # Point masses alone are exact; a continuous grid needs the provenance
  # that refines it.
  if(is.null(attr(prior_density, "adaptive_evaluation", exact = TRUE))){
    if(is.null(prior_density[["density"]]) ||
       !isTRUE(prior_density[["density"]][["mass"]] > 0)){
      return(evaluate_probability(prior_density))
    }
    .prior_linear_density_stop_no_provenance("region probability")
  }
  .prior_linear_density_check_grid(.prior_density_route_from_adaptive(
    attr(prior_density, "adaptive_evaluation", exact = TRUE)
  ))
  prob <- evaluate_probability(prior_density)
  refined <- .prior_linear_density_refinement(prior_density)
  if(is.null(refined)){
    stop(
      "Adaptive prior-probability evaluation did not converge within the ",
      "documented grid-refinement error criterion.",
      call. = FALSE
    )
  }
  if(!is.null(refined)){
    tolerance <- .prior_linear_density_refinement_tolerance()
    previous <- prob
    converged <- FALSE
    for(i in seq_len(4L)){
      current <- evaluate_probability(refined)
      change <- abs(current - previous)
      bound <- tolerance$absolute +
        tolerance$relative * max(abs(current), abs(previous))
      if(is.finite(current) && change <= bound){
        prob <- current
        converged <- TRUE
        break
      }
      previous <- current
      if(i < 4L){
        next_refined <- .prior_linear_density_refinement(refined)
        if(is.null(next_refined)){
          break
        }
        refined <- next_refined
      }
    }
    if(!converged){
      stop(
        "Adaptive prior-probability evaluation did not converge within the ",
        "documented grid-refinement error criterion.",
        call. = FALSE
      )
    }
  }

  prob
}


# The region of a condition on 'parameter' whose relations are linear in it
# (e.g. 'theta > 0', '2 * theta < 1', combined with '&', '|', '!'): its
# continuous part as disjoint intervals, and the exact condition as the
# indicator of point masses. NULL for other conditions.
.hypothesis_prior_region <- function(condition, parameter) {

  expr <- .hypothesis_parse_expression(condition)
  template <- stats::setNames(data.frame(0), parameter)
  intervals <- .hypothesis_region_intervals(expr, parameter, template)
  if(is.null(intervals)){
    return(NULL)
  }

  list(
    intervals = intervals,
    indicator = function(values){
      .hypothesis_condition_indicator(condition, parameter, values)
    }
  )
}


.hypothesis_region_intervals <- function(expr, parameter, template) {

  fun <- .hypothesis_call_name(expr)
  if(is.null(fun)){
    return(NULL)
  }
  if(fun == "("){
    return(.hypothesis_region_intervals(expr[[2L]], parameter, template))
  }
  if(fun == "!"){
    inner <- .hypothesis_region_intervals(expr[[2L]], parameter, template)
    return(if(is.null(inner)) NULL else .prior_region_intervals_complement(inner))
  }
  if(fun %in% c("&", "|")){
    left  <- .hypothesis_region_intervals(expr[[2L]], parameter, template)
    right <- .hypothesis_region_intervals(expr[[3L]], parameter, template)
    if(is.null(left) || is.null(right)){
      return(NULL)
    }
    return(if(fun == "&"){
      .prior_region_intervals_intersect(left, right)
    }else{
      .prior_region_intervals_union(left, right)
    })
  }
  if(!fun %in% c("<", "<=", ">", ">=")){
    return(NULL)
  }

  lhs <- .hypothesis_unwrap_parentheses(expr[[2L]])
  rhs <- .hypothesis_unwrap_parentheses(expr[[3L]])
  if(!all(.hypothesis_expression_symbols(expr) %in% parameter)){
    return(NULL)
  }
  # lhs - rhs = constant + coefficient * parameter; the relation holds below
  # or above its root. A bare parameter compared with a constant keeps the
  # constant as the exact boundary.
  if(is.name(lhs) && length(.hypothesis_expression_symbols(rhs)) == 0L){
    constant <- -tryCatch(.hypothesis_parse_number(rhs), error = function(e) NA_real_)
    coefficient <- 1
  }else if(is.name(rhs) && length(.hypothesis_expression_symbols(lhs)) == 0L){
    constant <- tryCatch(.hypothesis_parse_number(lhs), error = function(e) NA_real_)
    coefficient <- -1
  }else{
    linear <- tryCatch(
      .hypothesis_linear_coefficients(call("-", lhs, rhs), parameter, template),
      error = function(e) NULL
    )
    if(is.null(linear)){
      return(NULL)
    }
    constant <- linear[["constant"]]
    coefficient <- unname(linear[["coefficients"]][[1L]])
  }
  if(!is.finite(constant) || !is.finite(coefficient)){
    return(NULL)
  }
  below <- fun %in% c("<", "<=")
  if(coefficient == 0){
    holds <- if(below) constant < 0 || (fun == "<=" && constant == 0) else
      constant > 0 || (fun == ">=" && constant == 0)
    return(.prior_region_intervals(if(holds) -Inf else numeric(),
                                   if(holds) Inf else numeric()))
  }
  boundary <- -constant / coefficient
  if(below == (coefficient > 0)){
    .prior_region_intervals(-Inf, boundary)
  }else{
    .prior_region_intervals(boundary, Inf)
  }
}


.hypothesis_condition_indicator <- function(condition, parameter, values) {

  draws <- data.frame(values, check.names = FALSE)
  names(draws) <- parameter

  .hypothesis_eval_condition(condition, draws)
}


.hypothesis_comparison_grid_integral <- function(x, y, comparison) {

  # Integrate the continuous interpolant up to the actual boundary.
  # Multiplying grid ordinates by a step indicator moves that boundary
  # to neighbouring knots and introduces first-order grid error.
  value <- comparison[["value"]]
  lower <- comparison[["operator"]] %in% c("<", "<=")
  keep <- if(lower) x <= value else x >= value
  boundary <- value > min(x) && value < max(x) && !any(x == value)
  boundary_y <- if(boundary) stats::approx(x, y, xout = value)$y else NULL
  x <- x[keep]
  y <- y[keep]
  if(boundary){
    x <- if(lower) c(x, value) else c(value, x)
    y <- if(lower) c(y, boundary_y) else c(boundary_y, y)
  }

  .hypothesis_trapz(x, y)
}


.hypothesis_condition_grid_integral <- function(x, y, condition, parameter) {

  # Integrate the linear interpolant of the grid density over the region.
  # Grid cells whose endpoints disagree on the condition contain a region
  # boundary; it is located by bisection on the condition and inserted as
  # a knot with a linearly interpolated ordinate, so the integral has no
  # first-order boundary error.
  n <- length(x)
  inside <- .hypothesis_condition_indicator(condition, parameter, x)
  left_inside  <- inside[-n]
  right_inside <- inside[-1L]
  width        <- diff(x)
  cell_area    <- width * (y[-n] + y[-1L]) / 2

  unchanged <- left_inside == right_inside
  integral  <- sum(cell_area[unchanged & left_inside])

  changed <- which(!unchanged)
  if(length(changed) > 0L){
    boundary <- .hypothesis_condition_boundary(
      condition    = condition,
      parameter    = parameter,
      lower        = x[changed],
      upper        = x[changed + 1L],
      lower_inside = left_inside[changed]
    )
    boundary_y <- y[changed] + (y[changed + 1L] - y[changed]) *
      (boundary - x[changed]) / width[changed]
    left_area  <- (boundary - x[changed]) * (y[changed] + boundary_y) / 2
    right_area <- (x[changed + 1L] - boundary) *
      (boundary_y + y[changed + 1L]) / 2
    integral <- integral +
      sum(left_area[left_inside[changed]]) +
      sum(right_area[right_inside[changed]])
  }

  integral
}


.hypothesis_condition_boundary <- function(condition, parameter, lower, upper,
                                           lower_inside) {

  # Vectorized bisection: each interval keeps one endpoint on each side of
  # the condition boundary until it is resolved to about 1e-12 relative
  # precision (or to 2^-40 of the grid cell width near zero).
  cell_width <- upper - lower
  for(i in seq_len(64L)){
    width  <- upper - lower
    active <- width > 1e-12 * pmax(abs(lower), abs(upper)) &
      width > 2^-40 * cell_width
    if(!any(active)){
      break
    }
    middle <- (lower[active] + upper[active]) / 2
    stalled <- middle <= lower[active] | middle >= upper[active]
    if(all(stalled)){
      break
    }
    middle_inside <- .hypothesis_condition_indicator(
      condition, parameter, middle
    )
    move_lower <- middle_inside == lower_inside[active] & !stalled
    move_upper <- middle_inside != lower_inside[active] & !stalled
    active_i <- which(active)
    lower[active_i[move_lower]] <- middle[move_lower]
    upper[active_i[move_upper]] <- middle[move_upper]
  }

  (lower + upper) / 2
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

  # the exact Gaussian kernel sum at the value (bandwidth bw.nrd0 of the
  # draws); no evaluation grid, interpolation or binning
  height <- .density_kde_gaussian_height(
    x     = samples,
    value = value,
    bw    = stats::bw.nrd0(samples)
  )

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


.hypothesis_check_prior_mass <- function(mass, label, allow_one = FALSE) {

  # A region with prior mass one is a valid encompassing hypothesis in an
  # explicit comparison. An implicit statement compares a region with its
  # complement, which then has zero prior mass.
  if(!is.finite(mass) || mass <= 0){
    stop("Prior region mass for hypothesis '", label,
         "' is zero or non-finite.", call. = FALSE)
  }
  if(!isTRUE(allow_one) && mass >= 1){
    stop("Prior region mass for hypothesis '", label,
         "' is one, so its complement has zero prior mass.", call. = FALSE)
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
