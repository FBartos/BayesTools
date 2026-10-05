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
    expect_true(nzchar(result$warning[1L]))
    expect_identical(rownames(result), paste0("mu_fac[", c("A", "B", "C"), "]"))
    for(level in c("B", "C")){
      scalar <- hypothesis_BF(posterior, hypothesis = gsub("mu_fac", paste0("mu_fac[", level, "]"), statement, fixed = TRUE),
                              columns = "all", seed = 8)
      expect_equal(as.numeric(result[level == c("A", "B", "C"), "BF"]), as.numeric(scalar$BF), tolerance = 1e-14)
    }
    expect_identical(result$Alternative, rep(if(statement == "mu_fac = 0") "mu_fac != 0" else sub(" vs.*", "", statement), 3))
  }
  expect_error(hypothesis_BF(posterior, hypothesis = "mu_fac[A] = 0"),
               class = "BayesTools_hypothesis_ordinate")
  expect_error(hypothesis_BF(posterior, hypothesis = "mu_fac[A] > 0"),
               class = "BayesTools_hypothesis_region")
  expect_error(hypothesis_BF(posterior, hypothesis = "unknown = 0"))
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
  expect_warning(result <- hypothesis_BF(continuous, hypothesis = "mu_fac = 0", parameter = "mu_fac"),
                  class = "BayesTools_hypothesis_ordinate")
  expect_true(all(is.finite(as.numeric(result$BF))))
  testthat::local_mocked_bindings(.package = "BayesTools", .hypothesis_check_prior_ordinate = function(...) stop("Metadata is malformed.", call. = FALSE))
  expect_error(hypothesis_BF(posterior, hypothesis = "mu_fac = 0", parameter = "mu_fac"),
               "Metadata is malformed.", fixed = TRUE)
})

test_that("D3 marginal table known prior refusals preserve reasons and scalar strictness", {

  posterior <- .unavailable_factor_posterior()
  density <- .prior_linear_combination_density(list(theta = prior("point", list(0))), c(theta = 1))
  scalar <- .hypothesis_marginal_child(posterior$B)
  scalar <- .bt_meta_set(scalar, "prior_density", density)
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
