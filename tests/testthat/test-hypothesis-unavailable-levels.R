skip_if_not_test_profile("unit")

.unavailable_factor_posterior <- function(){

  formula <- JAGS_formula(~ fac, "mu", data.frame(fac = factor(c("A", "B", "C"))),
    prior_list = list(intercept = prior("normal", list(0, 1)),
                      fac = prior_factor("normal", list(0, 1), contrast = "treatment")))
  fit <- coda::mcmc(cbind(mu_intercept = seq(-1, 1, length.out = 201),
    "mu_fac[1]" = seq(-1, 2, length.out = 201),
    "mu_fac[2]" = seq(-2, 1, length.out = 201)))
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- formula$prior_list
  fit <- attach_test_parameter_map(fit)
  marginal_posterior(as_mixed_posteriors(fit, "mu_fac"), "mu_fac",
                    use_formula = FALSE, prior_samples = TRUE)
}

test_that("D3 whole-factor unavailable rows retain labels and valid siblings", {

  posterior <- .unavailable_factor_posterior()
  for(statement in c("mu_fac = 0", "mu_fac > 0", "mu_fac = 0 vs mu_fac > 0")){
    result <- tryCatch(hypothesis_BF(posterior, hypothesis = statement,
      parameter = "mu_fac", columns = "all", seed = 8), error = identity)
    expect_s3_class(result, "BayesTools_hypothesis_BF")
    if(inherits(result, "error")) next
    expect_identical(result$method[1L], "unavailable")
    expect_true(all(is.na(result[1L, c("BF", "BF_error", "prior", "posterior")])))
    expect_true(nzchar(attr(result, "warnings")[[rownames(result)[1L]]]))
    expect_identical(rownames(result), paste0("mu_fac[", c("A", "B", "C"), "]"))
    for(level in c("B", "C")){
      scalar <- hypothesis_BF(.hypothesis_marginal_child(posterior[[level]]),
                              hypothesis = statement, parameter = "mu_fac", columns = "all", seed = 8)
      expect_equal(as.numeric(result[level == c("A", "B", "C"), "BF"]), as.numeric(scalar$BF), tolerance = 1e-14)
    }
    expect_identical(result$Alternative, rep(if(statement == "mu_fac = 0") "mu_fac != 0" else sub(" vs.*", "", statement), 3))
  }
  expect_error(hypothesis_BF(posterior, hypothesis = "mu_fac[A] = 0"),
               class = "BayesTools_hypothesis_ordinate")
  expect_error(hypothesis_BF(posterior, hypothesis = "mu_fac[A] > 0"),
               class = "BayesTools_hypothesis_region")
  expect_error(hypothesis_BF(posterior, hypothesis = "mu_fac[missing] = 0"),
               class = "BayesTools_parameter_not_found")
})

test_that("D3 opted-in routes catch errors while warnings and generic errors propagate", {

  posterior <- .unavailable_factor_posterior()
  check_ordinate <- .hypothesis_check_prior_ordinate
  testthat::local_mocked_bindings(.package = "BayesTools", .hypothesis_check_prior_ordinate = function(...){
    warning(structure(list(message = "Known ordinate warning.", call = NULL),
      class = c("BayesTools_inexact_ordinate", "BayesTools_hypothesis_ordinate", "warning", "condition")))
    check_ordinate(...)
  })
  continuous <- posterior
  continuous$A <- NULL
  warnings <- list()
  result <- withCallingHandlers(hypothesis_BF(continuous, hypothesis = "mu_fac = 0", parameter = "mu_fac"),
    warning = function(w){ warnings[[length(warnings) + 1L]] <<- w; invokeRestart("muffleWarning") })
  expect_true(length(warnings) > 0L)
  expect_true(all(vapply(warnings, inherits, logical(1), "BayesTools_hypothesis_ordinate")))
  expect_true(all(is.finite(as.numeric(result$BF))))
  testthat::local_mocked_bindings(.package = "BayesTools", .hypothesis_check_prior_ordinate = function(...) stop("Metadata is malformed.", call. = FALSE))
  expect_error(hypothesis_BF(posterior, hypothesis = "mu_fac = 0", parameter = "mu_fac"),
               "Metadata is malformed.", fixed = TRUE)
})

test_that("D3 known region diagnostic refusals form complete unavailable tables", {

  posterior <- .unavailable_factor_posterior()
  testthat::local_mocked_bindings(.package = "BayesTools",
    .prior_region_from_adaptive = function(...) list(available = TRUE,
      converged = FALSE, probability = .5, messages = "known integration failure",
      absolute_error = .001))
  result <- hypothesis_BF(posterior, hypothesis = "mu_fac > 0", parameter = "mu_fac", columns = "all")
  expect_true(all(is.na(as.numeric(result$BF))))
  expect_identical(result$method, rep("unavailable", 3L))
  expect_match(attr(result, "warnings")[["mu_fac[B]"]],
               "Conditional-normal prior probability was rejected by diagnostics:", fixed = TRUE)
  exported <- as.data.frame(result)
  expect_true(all(nzchar(exported$warning)))
  expect_true(any(grepl("known integration failure", capture.output(print(result)), fixed = TRUE)))
  condition <- tryCatch(hypothesis_BF(.hypothesis_marginal_child(posterior$B),
    hypothesis = "mu_fac > 0", parameter = "mu_fac"), error = identity)
  expect_s3_class(condition, "BayesTools_prior_region_probability_rejected")
  expect_identical(condition$diagnostics$absolute_error, .001)
  expect_error(.prior_linear_density_stop_refinement(NULL, quantity = "probability"),
    "Adaptive prior-probability evaluation did not converge within the documented grid-refinement error criterion.",
    class = "BayesTools_prior_region_probability_rejected", fixed = TRUE)
})

test_that("D3 region missing provenance and invalid probability metadata remain strict", {

  posterior <- .unavailable_factor_posterior()
  density <- .bt_meta_get(posterior$B, "prior_density")
  attr(density, "adaptive_evaluation") <- NULL
  posterior$B <- .bt_meta_set(posterior$B, "prior_density", density)
  expect_error(hypothesis_BF(posterior, hypothesis = "mu_fac > 0", parameter = "mu_fac"),
               "provenance", fixed = TRUE)
  testthat::local_mocked_bindings(.package = "BayesTools",
    .prior_linear_density_region_probability = function(...){
      probability <- 2
      attr(probability, "numerical_diagnostics") <- list(absolute_error = .001)
      probability
    })
  condition <- tryCatch(hypothesis_BF(.hypothesis_marginal_child(posterior$C),
    hypothesis = "mu_fac > 0", parameter = "mu_fac"), error = identity)
  expect_s3_class(condition, "BayesTools_prior_region_probability_rejected")
  expect_identical(conditionMessage(condition), "Computed prior probability lies materially outside [0, 1].")
  expect_identical(condition$probability, 2)
  expect_identical(condition$absolute_error, .001)
})

test_that("D3 marginal table known prior refusals preserve reasons and scalar strictness", {

  fit <- coda::mcmc(cbind(beta = seq(.1, 2, length.out = 20),
                          beta_indicator = rep(1L, 20)))
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- list(beta = prior_spike_and_slab(
    prior("normal", list(0, 1)), prior_inclusion = prior("point", list(.5))))
  fit <- attach_test_parameter_map(fit)
  scalar <- marginal_posterior(as_mixed_posteriors(fit, "beta"), "beta",
                               NULL, prior_samples = TRUE)
  expect_error(Savage_Dickey_BF(scalar, silent = TRUE), class = "BayesTools_point_mass_at_null")
  result <- tryCatch(.Savage_Dickey_BF.checked(scalar, null_hypothesis = 0,
    normal_approximation = FALSE, density_method = "KDE", silent = TRUE,
    null_mass_NA = TRUE), error = identity)
  expect_false(inherits(result, "error"))
  if(inherits(result, "error")) return()
  expect_true(is.na(result))
  expect_identical(attr(result, "posterior_density_source"), "point_mass_prior_ordinate")
  expect_match(attr(result, "warnings"), "point mass", fixed = TRUE)
})

test_that("R133 D3 transports the specific image leaf while computing valid siblings", {
  # Synthetic compiled factor metadata and mocked prior transport, not fitted numerical evidence.
  posterior <- .unavailable_factor_posterior()
  posterior$A <- NULL
  region_mass <- .hypothesis_region_mass
  testthat::local_mocked_bindings(.package = "BayesTools", .hypothesis_region_mass = function(quantity, side, prior){
    if(prior && identical(quantity$label, "mu_fac[B]")){
      region <- list(intervals = .prior_region_intervals(-Inf, 0), indicator = function(x) x <= 0)
      result <- .prior_region_transformed("exp_lin", list(b = 2), region,
        function(r) .prior_region_atoms(1e-300, 1, r), function() c(0, Inf))
      return(result$probability)
    }
    region_mass(quantity, side, prior)
  })
  result <- hypothesis_BF(posterior, hypothesis = "mu_fac > 0", parameter = "mu_fac", columns = "all")
  expect_identical(result$method[1L], "unavailable")
  expect_true(is.na(as.numeric(result$BF[1L])))
  expect_true(is.finite(as.numeric(result$BF[2L])))
  expect_match(attr(result, "warnings")[["mu_fac[B]"]], "finite source values", fixed = TRUE)
  expect_error(hypothesis_BF(.hypothesis_marginal_child(posterior$B), hypothesis = "mu_fac > 0", parameter = "mu_fac"),
    class = "BayesTools_transformation_image_unavailable")
})
