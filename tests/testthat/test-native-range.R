skip_if_not_test_profile("unit")

test_that("numerical range advice matches the requested operation", {
  reason <- function(operation, scale){
    condition <- NULL
    value <- withCallingHandlers(.prior_numerical_result(Inf, 1, TRUE, operation, "invgamma", scale),
      warning = function(w){condition <<- w; invokeRestart("muffleWarning")})
    expect_identical(value, Inf)
    expect_s3_class(condition, "BayesTools_numerical_range_limit")
    expect_null(conditionCall(condition))
    expect_identical(condition$operation, operation)
    expect_identical(condition$requested_scale, scale)
    condition$reason
  }
  for(operation in c("density", "distribution")){
    expect_identical(reason(operation, "natural"), "Use an available logarithmic result or inspect the declared numerical limit")
  }
  for(operation in c("quantile", "sampling")){
    expect_identical(reason(operation, "natural"), "Inspect the declared numerical limit and the numerical condition")
  }
  for(scale in c("log", "finite")){
    expect_identical(reason("density", scale), "Inspect the declared numerical limit and the numerical condition")
  }
})

test_that("initialization retains parent range meaning and consumes one draw", {
  for(range in c(FALSE, TRUE)){
    parent <- .prior_numerical_condition("density", "moment", "log", 1L, "Controlled parent", range = range, error = TRUE)
    parent$log_density <- if(range) Inf else NaN
    count <- 0L
    condition <- testthat::with_mocked_bindings(
      expect_error(.JAGS_init.simple(prior("moment", list(location = 0, tau = 1, order = 1)), "theta"), class = "BayesTools_numerical_condition"),
      rng = function(...){count <<- count + 1L; .5}, .prior_simple_lpdf_evaluator = function(...) function(value) stop(parent), .package = "BayesTools")
    expect_identical(count, 1L)
    expect_identical(condition$parent, parent)
    expect_identical(condition$values, .5)
    expect_identical(condition$log_density, parent$log_density)
    expect_identical(condition$indices, 1L)
    expect_null(conditionCall(condition))
    expect_identical(inherits(condition, "BayesTools_numerical_range_limit"), range)
    expect_identical(condition$reason, if(range){
      "This consumer requires finite usable normalized log density; inspect the declared numerical limit and the numerical condition"
    }else "The declared normalized log density could not be resolved at supported precision")
  }
})

native_range_capture <- function(expression){

  conditions <- list()
  value <- withCallingHandlers(expression,
    BayesTools_numerical_condition = function(condition){
      conditions[[length(conditions) + 1L]] <<- condition
      if(inherits(condition, "warning")) invokeRestart("muffleWarning")
    })
  list(value = value, conditions = conditions)
}

native_range_available <- function(){

  .BayesTools_require_native_nonlocal()
  .Call("BayesTools_native_range_environment", PACKAGE = "BayesTools")$available
}

test_that("native range eligibility exposes exact internal runtime facts", {

  .BayesTools_require_native_nonlocal()
  facts <- .Call("BayesTools_native_range_environment", PACKAGE = "BayesTools")
  expect_identical(names(facts), c("compiled_profile", "documented_cpu", "hypervisor", "control_word", "mxcsr",
    "ln2_raw", "ln2_lower_raw", "ln2_upper_raw", "log2_raw", "threshold_raw", "available", "representation"))
  expect_type(facts$available, "logical")
  if(facts$compiled_profile){
    expect_identical(facts$representation, "x87-extended-little-endian-10")
    expect_identical(lengths(facts[c("ln2_raw", "ln2_lower_raw", "ln2_upper_raw", "log2_raw", "threshold_raw")]),
      c(ln2_raw = 10L, ln2_lower_raw = 10L, ln2_upper_raw = 10L, log2_raw = 8L, threshold_raw = 8L))
    if(facts$available){
      expect_identical(bitwAnd(facts$control_word, 3072L), 0L)
      expect_true(bitwAnd(facts$control_word, 768L) %in% c(512L, 768L))
      expect_identical(bitwAnd(facts$mxcsr, 57408L), 0L)
    }
  }else{
    expect_false(facts$available)
    expect_identical(facts$representation, "unsupported")
    expect_identical(facts$ln2_raw, raw())
    expect_true(is.na(facts$control_word))
  }
})

test_that("native FAR densities distinguish range from unavailable backends and support zeros", {

  priors <- list(prior("invgamma", list(shape = 2, scale = 1)),
    prior("moment", list(tau = 1)), prior("invmoment", list(tau = 1, df = 3)))
  values <- c(1e-310, 1e200, 1e-200)
  available <- native_range_available()
  for(i in seq_along(priors)){
    result <- native_range_capture(lpdf(priors[[i]], c(values[i], 0, NA_real_, NaN)))
    expect_identical(result$value[2:4], c(-Inf, NA_real_, NaN))
    expect_length(result$conditions, 1L)
    condition <- result$conditions[[1L]]
    expect_s3_class(condition, if(available) "BayesTools_numerical_range_limit" else "BayesTools_numerical_unavailable")
    expect_identical(condition$indices, 1L)
    expect_identical(condition$family, priors[[i]]$distribution)
    expect_identical(condition$operation, "density")
    expect_identical(condition$requested_scale, "log")
    expect_null(condition$call)
    if(available){
      expect_identical(conditionMessage(condition), paste0("The ", priors[[i]]$distribution,
        " prior density is numerically outside the requested representable range. ",
        "Inspect the declared numerical limit and the numerical condition."))
    }
    if(available) expect_identical(result$value[1L], -Inf) else expect_true(is.nan(result$value[1L]))
    natural <- native_range_capture(pdf(priors[[i]], values[i]))
    if(available) expect_identical(natural$value, 0) else expect_true(is.nan(natural$value))
    if(available){
      expect_identical(conditionMessage(natural$conditions[[1L]]), paste0("The ", priors[[i]]$distribution,
        " prior density is numerically outside the requested representable range. ",
        "Use an available logarithmic result or inspect the declared numerical limit."))
      expect_identical(natural$conditions[[1L]]$indices, 1L)
      expect_identical(natural$conditions[[1L]]$requested_scale, "natural")
      expect_null(conditionCall(natural$conditions[[1L]]))
    }
  }
})

test_that("finite range requests preserve the condition while recommending inspection", {
  result <- native_range_capture(.prior_numerical_result(Inf, 1, TRUE, "density", "moment", "finite"))
  expect_identical(result$value, Inf)
  expect_identical(class(result$conditions[[1L]]), c("BayesTools_numerical_range_limit",
    "BayesTools_numerical_condition", "warning", "condition"))
  expect_identical(conditionMessage(result$conditions[[1L]]), paste0(
    "The moment prior density is numerically outside the requested representable range. ",
    "Inspect the declared numerical limit and the numerical condition."))
  expect_identical(result$conditions[[1L]]$indices, 1L)
  expect_identical(result$conditions[[1L]]$requested_scale, "finite")
  expect_null(conditionCall(result$conditions[[1L]]))
  condition <- expect_error(.prior_numerical_result(Inf, 1, TRUE, "density", "moment", "finite", bounded = TRUE),
    class = "BayesTools_numerical_range_limit")
  expect_identical(class(condition), c("BayesTools_numerical_range_limit", "BayesTools_numerical_condition", "error", "condition"))
  expect_identical(condition$indices, 1L)
  expect_null(conditionCall(condition))
})

test_that("truncated native density warnings retain original eligible vector indices", {

  priors <- list(prior("invgamma", list(shape = 2, scale = 1), list(lower = 1e-310, upper = 1)),
    prior("moment", list(tau = 1), list(lower = 0, upper = 1e200)),
    prior("invmoment", list(tau = 1, df = 3), list(lower = 1e-200, upper = Inf)))
  values <- list(c(.5, 1e-310, 0, 1e-311, NA_real_, NaN),
    c(1, 1e200, 0, -1e200, NA_real_, NaN), c(1, 1e-200, 0, 1e-210, NA_real_, NaN))
  if(!native_range_available()){
    for(p in priors) expect_error(.prior_simple_lpdf_evaluator(p), class = "BayesTools_numerical_unavailable")
    return(invisible(NULL))
  }
  for(i in seq_along(priors)){
    for(log in c(TRUE, FALSE)){
      result <- native_range_capture(if(log) lpdf(priors[[i]], values[[i]]) else pdf(priors[[i]], values[[i]]))
      expect_true(is.finite(result$value[1L]))
      expect_identical(result$value[2:4], rep(if(log) -Inf else 0, 3L))
      expect_identical(result$value[5:6], c(NA_real_, NaN))
      expect_length(result$conditions, 1L)
      expect_s3_class(result$conditions[[1L]], "BayesTools_numerical_range_limit")
      expect_identical(result$conditions[[1L]]$indices, 2L)
      expect_identical(result$conditions[[1L]]$requested_scale, if(log) "log" else "natural")
    }
  }
})

test_that("mixed native vectors report both unavailable and certified range subsets", {

  skip_if_not(native_range_available(), "The native FAR range certificate is unavailable on this runtime.")
  # One finite original distance certifies FAR; the other original subtraction
  # overflows and has no exp/logaddexp certificate.
  result <- native_range_capture(.dmoment_prior(c(-1e308, .Machine$double.xmax, NA_real_, NaN),
    -.Machine$double.xmax, 1, 1, log = TRUE))
  expect_identical(result$value[1L], -Inf)
  expect_true(is.nan(result$value[2L]))
  expect_identical(result$value[3:4], c(NA_real_, NaN))
  expect_identical(vapply(result$conditions, function(e) class(e)[1L], character(1)),
    c("BayesTools_numerical_unavailable", "BayesTools_numerical_range_limit"))
  expect_identical(lapply(result$conditions, `[[`, "indices"), list(2L, 1L))
})

test_that("finite consumers refuse certified FAR logs while natural integration retains only tail zeros", {

  p <- prior("moment", list(tau = 1), list(lower = 0, upper = 1e200))
  if(!native_range_available()){
    expect_error(.prior_simple_lpdf_evaluator(p), class = "BayesTools_numerical_unavailable")
    return(invisible(NULL))
  }
  strict <- .prior_simple_lpdf_evaluator(p)
  condition <- tryCatch(strict(c(0, 1, 1e200)), error = identity)
  expect_s3_class(condition, "BayesTools_numerical_range_limit")
  expect_identical(condition$indices, 3L)
  expect_identical(condition$log_density, -Inf)
  expect_s3_class(condition$parent, "BayesTools_numerical_range_limit")
  expect_warning(expect_identical(strict(0), -Inf), NA)
  natural <- .prior_simple_lpdf_evaluator(p, purpose = "natural_integral")
  expect_warning(value <- natural(c(0, 1, 1e200)), NA)
  expect_identical(attr(value, "BayesTools_certified_log_tail"), 3L)
  expect_identical(unname(exp(value[c(1L, 3L)])), c(0, 0))
  expect_warning(.prior_quadrature_log_density(p, c(1, 1e200), natural), NA)
  conditional <- .prior_conditional_normal_integrand_parts(list(multiplier = p,
    additive_mean = 0, additive_sd = 1, product_mean = 1e-200, product_sd = 0))
  expect_warning(shared <- conditional$shared(1e200), NA)
  expect_identical(conditional$log_value(shared, 1), -Inf)
  expect_identical(conditional$value(shared, 1), 0)
  product <- .prior_scale_product_integrand_parts(list(factor = p,
    multiplier = prior("uniform", list(1, 2)), map = NULL, scale = 1))
  expect_warning(product_shared <- product$shared(1), NA)
  expect_identical(product$log_value(product_shared, 1e200), -Inf)
  expect_identical(product$value(product_shared, 1e200), 0)
  expect_error(.bt_JAGS_bridge_compile_simple_log_prior(p, "theta")(list(theta = 1e200)),
    class = "BayesTools_numerical_range_limit")
  expect_error(.bt_JAGS_marglik_compile_prior_rows_component(p, "theta")(
    matrix(c(1, 1e200), ncol = 1, dimnames = list(NULL, "theta"))), class = "BayesTools_numerical_range_limit")
  factor <- prior_factor("moment", list(tau = 1), contrast = "independent")
  factor_error <- tryCatch(.bt_JAGS_bridge_compile_simple_log_prior(factor, c("b[1]", "b[2]"))(
    list("b[1]" = 1, "b[2]" = 1e200)), error = identity)
  expect_s3_class(factor_error, "BayesTools_numerical_range_limit")
  expect_identical(factor_error$indices, 2L)
  invgamma <- prior("invgamma", list(shape = 2, scale = 1))
  expect_error(.bt_JAGS_marglik_compile_invgamma_prior_rows(invgamma, "theta")(
    matrix(c(1, 1e-310), ncol = 1, dimnames = list(NULL, "theta"))), class = "BayesTools_numerical_range_limit")
  ordinate <- prior_density_ordinate(p, 1e200)
  expect_identical(ordinate$behavior, "regular")
  expect_false(ordinate$exact)
  expect_true(is.na(ordinate$log_density))
  unresolved <- prior("invmoment", list(tau = 1, order = 1000,
    df = .Machine$double.xmin * .Machine$double.eps))
  expect_error(.prior_simple_lpdf_evaluator(unresolved, purpose = "natural_integral")(.5),
    class = "BayesTools_numerical_unavailable")
})

test_that("natural integration preserves unrelated failures and checked kernel images", {

  p <- prior("normal", list(0, 1))
  natural <- .prior_simple_lpdf_evaluator(p, purpose = "natural_integral")
  testthat::with_mocked_bindings({
    expect_warning(natural(1), "unrelated warning", fixed = TRUE)
  }, .prior_simple_lpdf = function(...) {warning("unrelated warning"); 0})
  for(value in c(NaN, Inf)){
    testthat::with_mocked_bindings({
      expect_error(natural(1), class = "BayesTools_numerical_unavailable")
    }, .prior_simple_lpdf = function(...) value)
  }
  spec <- list(factor = prior("invgamma", list(shape = 2, scale = 1)),
    multiplier = prior("uniform", list(1, 2)), bounds = c(1, 2), scale = 1,
    offset = 0, sources = list(), map = NULL)
  parts <- .prior_scale_product_integrand_parts(spec)
  expect_error(parts$value(parts$shared(1), Inf), class = "BayesTools_numerical_condition")
})

test_that("one FAR intervals normalize and invert while positive both FAR intervals refuse", {

  skip_if_not(native_range_available(), "The native FAR range certificate is unavailable on this runtime.")
  priors <- list(prior("invgamma", list(shape = 2, scale = 1), list(lower = 1e-310, upper = 1)),
    prior("moment", list(tau = 1), list(lower = 0, upper = 1e200)),
    prior("invmoment", list(tau = 1, df = 3), list(lower = 1e-200, upper = 1)))
  for(p in priors){
    expect_true(is.finite(.prior_simple_log_C(p)))
    q <- quant(p, c(.1, .5, .9))
    expect_true(all(is.finite(q)))
    expect_equal(cdf(p, q), c(.1, .5, .9), tolerance = 3e-12)
    set.seed(278)
    u <- stats::runif(5)
    seed <- .Random.seed
    expected <- quant(p, u)
    set.seed(278)
    expect_identical(rng(p, 5), expected)
    expect_identical(.Random.seed, seed)
  }
  expect_equal(.prior_simple_log_C(priors[[1L]]), stats::pgamma(1, 2, lower.tail = FALSE, log.p = TRUE), tolerance = 2e-13)
  expect_equal(.prior_simple_log_C(priors[[2L]]), -log(2), tolerance = 2e-13)
  both <- list(prior("invgamma", list(shape = 2, scale = 1), list(lower = 1e-310, upper = 2e-310)),
    prior("moment", list(tau = 1), list(lower = 1e200, upper = 2e200)),
    prior("invmoment", list(tau = 1, df = 3), list(lower = 1e-200, upper = 2e-200)))
  for(p in both) expect_error(.prior_simple_log_C(p), class = "BayesTools_numerical_unavailable")
})

test_that("nonlocal FAR starts retain the range parent without replacement draws", {

  skip_if_not(native_range_available(), "The native FAR range certificate is unavailable on this runtime.")
  p <- prior("moment", list(tau = 1))
  calls <- 0L
  testthat::with_mocked_bindings({
    condition <- tryCatch(.JAGS_init.simple(p, "theta"), error = identity)
    expect_s3_class(condition, "BayesTools_numerical_range_limit")
    expect_identical(condition$operation, "initialization")
    expect_identical(condition$indices, 1L)
    expect_identical(condition$values, 1e200)
    expect_identical(condition$log_density, -Inf)
    expect_s3_class(condition$parent, "BayesTools_numerical_range_limit")
  }, rng = function(...){calls <<- calls + 1L; 1e200})
  expect_identical(calls, 1L)
})
