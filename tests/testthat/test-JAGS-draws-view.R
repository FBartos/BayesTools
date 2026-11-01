skip_if_not_test_profile("unit")

.draws_view_fit_for_test <- function(){

  chains <- coda::mcmc.list(lapply(c(0, 3), function(offset){
    coda::mcmc(cbind(theta = offset + stats::qnorm((seq_len(60) - .5) / 60)), start = 11, thin = 2)
  }))
  fit <- structure(list(mcmc = chains, sample = 120L, model = "original model",
                         data = "original data", end.state = list("state"),
                         monitor = "theta",
                         summary.pars = list(mutate = NULL)),
                    class = c("runjags", "BayesTools_fit"))
  attr(fit, "prior_list") <- list(theta = prior("normal", list(0, 1)))
  fit <- attach_test_parameter_map(fit)
  attr(fit, "runtime_setup") <- function(...) stop("Runtime setup must not run.")
  attr(fit, "runtime_cache") <- function(...) stop("Runtime cache must not run.")
  attr(fit, "runtime_state") <- new.env(parent = emptyenv())
  fit
}

.draws_view_replacement <- function(){

  coda::mcmc.list(lapply(c(-2, 2), function(offset){
    coda::mcmc(cbind(theta = offset + stats::qnorm((seq_len(30) - .5) / 30)), start = 51, thin = 3)
  }))
}

test_that("a derived-draw view prints without a leading space on its second line (F221)", {

  view <- JAGS_with_draws(.draws_view_fit_for_test(), .draws_view_replacement())
  output <- capture.output(print(view))
  expect_identical(output, c(
    "Derived-draw view: 2 chains, 60 draws.",
    "The complete original sampling fit is available as 'fit$original_fit'."
  ))
})

test_that("D9 derived views preserve original state and local metadata without sampler slots", {

  fit <- .draws_view_fit_for_test()
  replacement <- .draws_view_replacement()
  view <- JAGS_with_draws(fit, replacement)
  expect_identical(names(view), c("original_fit", "mcmc"))
  expect_identical(class(view), c("BayesTools_draws_view", "BayesTools_fit"))
  expect_identical(view$original_fit, fit)
  expect_identical(attr(view$original_fit, "runtime_state"), attr(fit, "runtime_state"))
  expect_identical(attr(view, "parameter_map"), attr(fit, "parameter_map"))
  expect_null(attr(view, "runtime_setup"))
  expect_null(attr(view, "runtime_cache"))
  expect_null(attr(view, "runtime_state"))
  attr(view, "prior_list") <- list(theta = prior("normal", list(1, 2)))
  renewed <- JAGS_with_draws(view, replacement)
  expect_identical(renewed$original_fit, fit)
  expect_identical(attr(renewed, "prior_list"), attr(view, "prior_list"))
  expect_identical(coda::as.mcmc.list(view), replacement)
  expect_identical(as.matrix(coda::as.mcmc(view)), do.call(rbind, lapply(replacement, as.matrix)))
  expect_identical(attr(coda::as.mcmc(view), "mcpar"), c(1, 60, 1))
  expect_identical(JAGS_draw_geometry(view)$chains$iterations, c(30L, 30L))
  selection <- parameter_catalog_resolve(parameter_catalog(view), "theta")
  expect_equal(as.numeric(as.matrix(parameter_draws(view, selection))), as.numeric(coda::as.mcmc(view)))
  expect_equal(as.numeric(as_mixed_posteriors(view, "theta")$theta), as.numeric(coda::as.mcmc(view)))
})

test_that("D9 view boundaries refuse before evaluating runtime or inference arguments", {

  view <- JAGS_with_draws(.draws_view_fit_for_test(), .draws_view_replacement())
  expected <- c(
    extension = "Sampling extension is unavailable for a derived-draw view. Extend 'fit$original_fit' with 'JAGS_extend()' and regenerate the view with 'JAGS_with_draws()'.",
    convergence = "Model convergence is unavailable for a derived-draw view. Use 'JAGS_check_convergence()' on 'fit$original_fit'.",
    bridge = "Bridge sampling is unavailable for a derived-draw view. Use 'JAGS_bridgesampling()' on 'fit$original_fit'.")
  calls <- list(extension = function() JAGS_extend(view, runtime_setup = stop("forced")),
                convergence = function() JAGS_check_convergence(view, max_Rhat = stop("forced")),
                bridge = function() JAGS_bridgesampling(view, log_posterior = stop("forced"), seed = stop("forced")))
  for(operation in names(calls)){
    condition <- tryCatch(calls[[operation]](), error = identity)
    expect_s3_class(condition, "BayesTools_draws_view_unavailable")
    expect_identical(conditionMessage(condition), unname(expected[operation]))
    expect_identical(condition$operation, operation)
    expect_null(condition$call)
  }
  partial <- tryCatch(JAGS_with_draws(view, coda::mcmc.list(lapply(.draws_view_replacement(), function(chain){
    coda::mcmc(matrix(numeric(), nrow = nrow(chain), ncol = 0), start = 51, thin = 3)
  }))), error = identity)
  expect_s3_class(partial, "BayesTools_draws_view")
  if(inherits(partial, "error")) return()
  expect_identical(dim(coda::as.mcmc(partial)), c(60L, 0L))
  expect_identical(dim(.extract_posterior_samples(partial)), c(60L, 0L))
  expect_identical(dim(as.matrix(.fit_to_posterior(partial))), c(60L, 0L))
  expect_identical(coda::as.mcmc.list(partial), partial$mcmc)
  condition <- tryCatch(JAGS_materialize_draws(partial, "theta"), error = identity)
  expect_s3_class(condition, "BayesTools_draws_view_coordinate_unavailable")
  expect_identical(condition$missing, "theta")
  expect_identical(conditionMessage(condition),
    "The derived-draw view does not contain requested sampled coordinates: 'theta'. Regenerate the view with those coordinates or use 'fit$original_fit'.")
})

test_that("D9 estimates use supplied two-chain geometry and preserve diagnostic early returns", {

  fit <- .draws_view_fit_for_test()
  replacement <- .draws_view_replacement()
  view <- JAGS_with_draws(fit, replacement)
  table <- tryCatch(JAGS_estimates_table(view), error = identity)
  expect_s3_class(table, "BayesTools_table")
  if(inherits(table, "error")) return()
  summary <- summary(replacement, quantiles = NULL)$statistics
  expected <- c(MCMC_error = summary[["Time-series SE"]],
    MCMC_SD_error = summary[["Time-series SE"]] / summary[["SD"]],
    ESS = unname(coda::effectiveSize(replacement)),
    R_hat = coda::gelman.diag(replacement, multivariate = FALSE, autoburnin = FALSE)$psrf[1, 1])
  expect_equal(unlist(table[1, names(expected)]), expected, tolerance = 1e-14, ignore_attr = TRUE)
  expect_equal(table$Mean, mean(as.numeric(coda::as.mcmc(view))), tolerance = 1e-14)
  expect_false("ESS" %in% names(JAGS_estimates_table(view, conditional = TRUE)))
  expect_false("ESS" %in% names(JAGS_estimates_table(view, remove_diagnostics = TRUE)))
  expect_identical(view$original_fit, fit)
})

test_that("estimates diagnostics match known-answer chains parameter by parameter", {

  # four stationary AR(1) chains per parameter: theta with phi = 0.5 and equal
  # chains, eta with phi = 0.9 and chain means 0, 0.5, 1, 1.5. With unit
  # marginal variance, a chain's spectral density at zero is
  # (1 + phi) / (1 - phi), so the effective sample size is n (1 - phi) / (1 + phi)
  # per chain and the MCMC error sqrt((1 + phi) / (1 - phi) / (n chains)); chain
  # shifts change neither. Over 100 replicates of this design coda's ESS
  # scatters by 3-4.5% and its MCMC error by 2-3.5% (SDs), so 0.2 and 0.15 are
  # about five SDs, while the IID standard error is 0.58 and 0.23 of the target.
  # Equal chains give R-hat <= 1.002 in those replicates; the shifted means give
  # a population R-hat of at least sqrt(1 + 0.3125). A wrong chain split,
  # column mapping, or error definition fails
  set.seed(1)
  n   <- 5000L
  ar1 <- function(phi) {
    x <- numeric(n)
    x[1L] <- stats::rnorm(1L)
    innovation <- stats::rnorm(n, sd = sqrt(1 - phi^2))
    for (i in 2:n) x[i] <- phi * x[i - 1L] + innovation[i]
    x
  }
  draws <- coda::mcmc.list(lapply(0:3, function(chain) {
    coda::mcmc(cbind(theta = ar1(0.5), eta = ar1(0.9) + 0.5 * chain))
  }))
  fit <- .draws_view_fit_for_test()
  fit$mcmc    <- draws
  fit$monitor <- c("theta", "eta")
  attr(fit, "prior_list") <- list(theta = prior("normal", list(0, 1)), eta = prior("normal", list(0, 1)))
  fit   <- attach_test_parameter_map(fit)
  table <- JAGS_estimates_table(JAGS_with_draws(fit, draws))

  phi          <- c(theta = 0.5, eta = 0.9)
  expected_ess <- 4 * n * (1 - phi) / (1 + phi)
  expected_se  <- sqrt((1 + phi) / (1 - phi) / (4 * n))
  expect_equal_each(unname(table[c("theta", "eta"), "ESS"] / expected_ess), c(1, 1), tolerance = 0.2)
  expect_equal_each(unname(table[c("theta", "eta"), "MCMC_error"] / expected_se), c(1, 1), tolerance = 0.15)
  expect_lt(table["theta", "R_hat"], 1.01)
  expect_gt(table["eta", "R_hat"], sqrt(1 + 0.3125))
  expect_equal(table[["MCMC_SD_error"]], table[["MCMC_error"]] / table[["SD"]], tolerance = 1e-12)
})

test_that("D9 fresh views clear fit-level posterior estimates and repeated-view metadata", {

  fit <- .draws_view_fit_for_test()
  density <- posterior_density_attribute(c(-4, 0, 4), c(.01, .2, .01), method = "iwmde", density_method = "precomputed")
  ordinate <- posterior_ordinate_attribute(0, .2, method = "qCMDE", density_method = "precomputed")
  for(field in c("posterior_density", "posterior_densities", "posterior_ordinate", "posterior_ordinates")){
    value <- if(grepl("density|densities", field)) density else ordinate
    if(field %in% c("posterior_densities", "posterior_ordinates")) value <- list(theta = value)
    fit <- .bt_meta_set(fit, field, value)
  }
  view <- JAGS_with_draws(fit, .draws_view_replacement())
  expect_null(attr(view, "bayestools_meta"))
  expect_identical(view$original_fit, fit)
  mixed <- as_mixed_posteriors(view, "theta")
  marginal <- marginal_posterior(mixed, "theta", NULL, prior_samples = TRUE)
  expect_null(.bt_meta_get(mixed$theta, "posterior_density"))
  expect_null(.bt_meta_get(marginal, "posterior_ordinate"))
  expect_error(Savage_Dickey_BF(marginal, density_method = "precomputed", silent = TRUE),
               "requires valid posterior ordinate", fixed = TRUE)
  view <- .bt_meta_set(view, "posterior_density", density)
  renewed <- JAGS_with_draws(view, .draws_view_replacement())
  expect_null(attr(renewed, "bayestools_meta"))
  expect_identical(renewed$original_fit, fit)
})
