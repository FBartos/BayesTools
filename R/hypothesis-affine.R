# A bounded binary64 affine numerator and constant divisor. Keeping the
# divisor avoids rounding a reciprocal before certifying a comparison.
.hypothesis_numerical_stop <- function(reason, operation, inputs = NULL,
                                        diagnostics = NULL, region = FALSE){

  leaf <- if(region) "BayesTools_hypothesis_region_numerical_unavailable" else
    "BayesTools_hypothesis_numerical_unavailable"
  parent <- if(region) "BayesTools_hypothesis_region" else "BayesTools_hypothesis_ordinate"
  observed <- if(identical(reason, "log_resolution")){
    paste0("the log-ordinate resolution estimate is ", format(diagnostics$resolution, digits = 6))
  }else if(identical(reason, "normal_scale_unavailable")){
    paste0("the Normal ", diagnostics$label, " standard deviation is ",
      format(diagnostics$standard_deviation, digits = 6))
  }else{
    paste0("", gsub("_", " ", reason), " in ", operation)
  }
  stop(errorCondition(
    paste0("Hypothesis inference was rejected by diagnostics: ", observed,
           ". Inspect the numerical diagnostics or use a supported region hypothesis."),
    class = c(leaf, parent, "BayesTools_numerical_unavailable",
              "BayesTools_numerical_condition"), call = NULL,
    operation = operation, family = "hypothesis", requested_scale = if(region) "probability" else "log",
    indices = 1L, reason = reason, inputs = inputs, diagnostics = diagnostics
  ))
}

.hypothesis_affine_product <- function(left, right, operation = "affine multiplication",
                                       region = FALSE){

  out <- left * right
  bad <- !is.finite(out) | (left != 0 & right != 0 & out == 0)
  if(any(bad)) .hypothesis_numerical_stop("coefficient_range", operation,
    list(left = left, right = right), region = region)
  out
}

.hypothesis_affine_sum <- function(left, right, operation = "affine addition",
                                   region = FALSE){

  out <- left + right
  bad <- !is.finite(out) | (left != 0 & right != 0 &
    ((out == left & right != -left) | (out == right & left != -right)))
  if(any(bad)) .hypothesis_numerical_stop("absorbed_nonzero_term", operation,
    list(left = left, right = right), region = region)
  out
}

.hypothesis_affine_new <- function(constant = 0, coefficients = numeric(), divisor = 1){

  if(!is.finite(divisor) || divisor == 0 || !is.finite(constant) ||
     any(!is.finite(coefficients))){
    .hypothesis_numerical_stop("coefficient_range", "affine representation")
  }
  if(divisor < 0){
    constant <- -constant
    coefficients <- -coefficients
    divisor <- -divisor
  }
  coefficients <- coefficients[coefficients != 0]
  if(length(coefficients)) coefficients <- coefficients[order(names(coefficients))]
  list(constant = constant, coefficients = coefficients, divisor = divisor)
}

.hypothesis_affine_add <- function(left, right, sign = 1, region = FALSE){

  if(identical(left$divisor, right$divisor)){
    left_scale <- right_scale <- 1
    divisor <- left$divisor
  }else{
    left_scale <- right$divisor
    right_scale <- left$divisor
    divisor <- .hypothesis_affine_product(left$divisor, right$divisor, region = region)
  }
  columns <- sort(union(names(left$coefficients), names(right$coefficients)))
  a <- b <- stats::setNames(numeric(length(columns)), columns)
  a[names(left$coefficients)] <- .hypothesis_affine_product(left$coefficients, left_scale, region = region)
  b[names(right$coefficients)] <- .hypothesis_affine_product(right$coefficients, sign * right_scale, region = region)
  .hypothesis_affine_new(
    .hypothesis_affine_sum(.hypothesis_affine_product(left$constant, left_scale, region = region),
      .hypothesis_affine_product(right$constant, sign * right_scale, region = region), region = region),
    .hypothesis_affine_sum(a, b, region = region), divisor
  )
}

.hypothesis_affine_scale <- function(form, constant, divide = FALSE, region = FALSE){

  if(length(constant$coefficients) != 0L) return(NULL)
  if(divide && constant$constant == 0) return(NULL)
  factor <- if(divide) constant$divisor else constant$constant
  denominator <- if(divide) constant$constant else constant$divisor
  .hypothesis_affine_new(
    .hypothesis_affine_product(form$constant, factor, region = region),
    .hypothesis_affine_product(form$coefficients, factor, region = region),
    .hypothesis_affine_product(form$divisor, denominator, region = region)
  )
}

.hypothesis_affine_read <- function(expr, symbols, region = FALSE){

  if(is.character(expr)) expr <- .hypothesis_parse_expression(expr)
  if(is.numeric(expr) && length(expr) == 1L) return(.hypothesis_affine_new(expr))
  if(is.name(expr)){
    name <- .hypothesis_decode_escaped_constant(as.character(expr))
    if(!name %in% symbols) return(NULL)
    return(.hypothesis_affine_new(coefficients = stats::setNames(1, name)))
  }
  if(!is.call(expr)) return(NULL)
  fun <- .hypothesis_call_name(expr)
  if(is.null(fun)) return(NULL)
  if(!fun %in% c("(", "+", "-", "*", "/")){
    if(!fun %in% c("^", "abs", "exp", "log", "sqrt", "plogis", "qlogis")) return(NULL)
    if(length(.hypothesis_expression_symbols(expr)) != 0L) return(NULL)
    arguments <- lapply(as.list(expr[-1L]), .hypothesis_affine_read,
      symbols = symbols, region = region)
    if(any(vapply(arguments, is.null, logical(1)))) return(NULL)
    values <- vapply(arguments, function(form){
      value <- form$constant / form$divisor
      if(!is.finite(value) || (form$constant != 0 && value == 0)){
        .hypothesis_numerical_stop("coefficient_range", "constant function argument", region = region)
      }
      value
    }, numeric(1))
    value <- do.call(fun, as.list(values), envir = .hypothesis_eval_parent())
    collapsed <- (fun == "exp" && (value == 0 || (values[[1L]] != 0 && value == 1))) ||
      (fun == "plogis" && (value %in% c(0, 1) || (values[[1L]] != 0 && value == .5))) ||
      (fun == "^" && values[[1L]] != 0 &&
       (value == 0 || (abs(values[[1L]]) != 1 && values[[2L]] != 0 && value == 1)))
    if(!is.finite(value) || collapsed){
      .hypothesis_numerical_stop("coefficient_range", "constant function", inputs = values, region = region)
    }
    return(.hypothesis_affine_new(value))
  }
  args <- lapply(as.list(expr[-1L]), .hypothesis_affine_read, symbols = symbols, region = region)
  if(any(vapply(args, is.null, logical(1)))) return(NULL)
  if(length(args) == 1L){
    if(fun == "(") return(args[[1L]])
    if(fun %in% c("+", "-")) return(.hypothesis_affine_scale(args[[1L]],
      .hypothesis_affine_new(if(fun == "-") -1 else 1), region = region))
    return(NULL)
  }
  if(length(args) != 2L) return(NULL)
  left <- args[[1L]]
  right <- args[[2L]]
  switch(fun,
    "+" = .hypothesis_affine_add(left, right, region = region),
    "-" = .hypothesis_affine_add(left, right, -1, region = region),
    "*" = if(length(left$coefficients) == 0L){
      .hypothesis_affine_scale(right, left, region = region)
    }else .hypothesis_affine_scale(left, right, region = region),
    "/" = .hypothesis_affine_scale(left, right, divide = TRUE, region = region)
  )
}

.hypothesis_affine_value <- function(form, draws, offset = FALSE, region = FALSE){

  out <- rep(if(offset) form$constant else 0, nrow(draws))
  for(symbol in names(form$coefficients)){
    term <- .hypothesis_affine_product(form$coefficients[[symbol]], draws[[symbol]], region = region)
    out <- .hypothesis_affine_sum(out, term, region = region)
  }
  out
}

.hypothesis_affine_null <- function(form, value, region = FALSE){

  .hypothesis_affine_sum(.hypothesis_affine_product(value, form$divisor, region = region),
                         -form$constant, "translated null", region = region)
}

.hypothesis_affine_comparison <- function(expr, symbols){

  if(!is.call(expr)) return(NULL)
  op <- .hypothesis_call_name(expr)
  if(is.null(op) || !op %in% c("<", "<=", ">", ">=", "==", "!=")) return(NULL)
  left <- .hypothesis_affine_read(expr[[2L]], symbols, region = TRUE)
  right <- .hypothesis_affine_read(expr[[3L]], symbols, region = TRUE)
  if(is.null(left) || is.null(right)) return(NULL)
  form <- .hypothesis_affine_add(left, right, -1, region = TRUE)
  list(form = form, operator = op, value = -form$constant)
}

.hypothesis_structural_condition <- function(expr, draws){

  if(is.character(expr)) expr <- .hypothesis_parse_expression(expr)
  fun <- .hypothesis_call_name(expr)
  if(is.null(fun)) return(NULL)
  if(fun == "(") return(.hypothesis_structural_condition(expr[[2L]], draws))
  if(fun == "!"){
    inner <- .hypothesis_structural_condition(expr[[2L]], draws)
    return(if(is.null(inner)) NULL else !inner)
  }
  if(fun %in% c("&", "|")){
    left <- .hypothesis_structural_condition(expr[[2L]], draws)
    right <- .hypothesis_structural_condition(expr[[3L]], draws)
    if(is.null(left) || is.null(right)) return(NULL)
    return(if(fun == "&") left & right else left | right)
  }
  comparison <- .hypothesis_affine_comparison(expr, names(draws))
  if(is.null(comparison)) return(NULL)
  value <- .hypothesis_affine_value(comparison$form, draws, region = TRUE)
  do.call(comparison$operator, list(value, comparison$value))
}

# The scalar tolerance governs relative density/BF accuracy. This screen is
# an operand-resolution heuristic, not an estimator or backend error bound.
.hypothesis_log_ratio <- function(log_prior, log_posterior, normal = FALSE){

  log_prior <- as.numeric(log_prior)
  log_posterior <- as.numeric(log_posterior)
  diagnostics <- list(log_prior = log_prior, log_posterior = log_posterior,
    resolution = 8 * .Machine$double.eps * max(1, abs(log_prior), abs(log_posterior)),
    budget = log1p(.prior_linear_density_refinement_tolerance()$relative))
  if(!is.finite(log_prior) || (normal && !is.finite(log_posterior)) || is.na(log_posterior)){
    .hypothesis_numerical_stop("nonfinite_log_ordinate", "log density ratio",
      list(log_prior = log_prior, log_posterior = log_posterior), diagnostics)
  }
  if(is.finite(log_posterior) && diagnostics$resolution > diagnostics$budget){
    .hypothesis_numerical_stop("log_resolution", "log density ratio",
      list(log_prior = log_prior, log_posterior = log_posterior), diagnostics)
  }
  list(log_BF = log_prior - log_posterior, diagnostics = diagnostics)
}
